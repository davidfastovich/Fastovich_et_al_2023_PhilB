"""
SCRIPT TAKES ALL MODEL OUTPUT FROM ~/MODEL_TEMPERATURE_ANOMALY_50KM_ENA AND 
CALCULATES MODEL SKILL (HARGREAVES ET AL., 2013, JC) BY BILINEARLY
INTERPOLATING MODEL OUTPUT TO THE PROXY DATA. SINCE THE PROXY DATA IS MORE
ABUNDANT IN THE NORTH, WE CORRECT FOR THIS BY RESAMPLING AT A GRID CELL
LEVEL.

MODEL DATA IS FROM KAGEYAMA ET AL. (2013), IVANOVIC ET AL. (2016), BROWN AND
GALBRAITH (2016), AND HE ET AL. (2013). PROXY DATA IS FROM FASTOVICH ET AL.
(2020) AND FASTOVICH ET AL. (2022).

INPUT:
    MODEL
        * MODEL_TEMPERATURE_ANOMALY_50KM_ENA/CCSM-MARUM_TAS_50KM_ENA.NC
        * MODEL_TEMPERATURE_ANOMALY_50KM_ENA/CM2MC_TAS_50KM_ENA.NC
        * MODEL_TEMPERATURE_ANOMALY_50KM_ENA/COSMOS-W_TAS_50KM_ENA.NC
        * MODEL_TEMPERATURE_ANOMALY_50KM_ENA/IPSL_TAS_50KM_ENA.NC
        * MODEL_TEMPERATURE_ANOMALY_50KM_ENA/MIROC-W_TAS_50KM_ENA.NC
        * MODEL_TEMPERATURE_ANOMALY_50KM_ENA/TRACE_MWF_TAS_50KM_ENA.NC
        * MODEL_TEMPERATURE_ANOMALY_50KM_ENA/CCSM-NCAR_TAS_50KM_ENA.NC
        * MODEL_TEMPERATURE_ANOMALY_50KM_ENA/COSMOS-S_TAS_50KM_ENA.NC
        * MODEL_TEMPERATURE_ANOMALY_50KM_ENA/HADCM3_TAS_50KM_ENA.NC
        * MODEL_TEMPERATURE_ANOMALY_50KM_ENA/MIROC-S_TAS_50KM_ENA.NC
        * MODEL_TEMPERATURE_ANOMALY_50KM_ENA/TRACE_LORENZ_TAS_50KM_ENA.NC
        * MODEL_TEMPERATURE_ANOMALY_50KM_ENA/TRACE_TAS_50KM_ENA.NC
        
    PROXY
        * PROXY_DATA/TAVE_ANOMALIES.CSV
        * PROXY_DATA/GDGT_MAT.CSV
        
OUTPUT:
    * SKILL_SCORES/TAS_SKILL.CSV

MINICONDA ENVIRONMENT:
    PYTHON: 3.8.10
    SYS: BUILT INTO PYTHON
    PARALLEL_SKILL: HOMEBREWED FUNCTION
    PANDAS: 1.2.4
    NUMPY: 1.20.3
    JOBLIB: 1.1.0
    XARRAY: 0.18.2
    CARTOPY: 0.18.0
    MATPLOTLIB: 3.4.2
    GLOB: BUILT INTO PYTHON
    GEOPANDAS: 0.10.2
    
@author: David Fastovich
"""

import sys
sys.path.append('modules')
from parallel_skill import parallel_skill
import pandas as pd
import numpy as np
from joblib import Parallel, delayed
import xarray as xr
import glob
import geopandas

################
# READ IN MODELS
################

# Models to read in
model_filenames = glob.glob('model_temperature_anomaly_50km_ena/*.nc')

# Create an empty xarray that will serve as the extra dimension for indexing
# the models into a single xarray file
model_concat_xr = xr.DataArray(model_filenames, dims=['model'], name='model')

# Empty list to hold the models
models = list()

# Loop through and append each model onto the empty list
for ii, file in enumerate(model_filenames):
    if ii == 10:
        # Fill Lorenz with NAs
        models.append(xr.open_dataset(file, decode_times=False)['tas'].fillna(-99999))
    else:
        models.append(xr.open_dataset(file, decode_times=False)['tas'])

# Use the empty xarray to reshape the list into a single xarray where the 
# 'model' dimension is associated with each unique model
models = xr.concat(models, dim=model_concat_xr)

###########################################
# READ IN PROXY TEMPERATURE RECONSTRUCTIONS
###########################################

proxy = pd.read_csv('proxy_data/tave_anomalies.csv')
gdgt = pd.read_csv('proxy_data/gdgt_mat.csv')
proxy = pd.concat([proxy, gdgt])

# Drop NAs
proxy = proxy.dropna(axis=0)

##############################################################
# REPROJECT PROXY DATA TO THE 50KM GRID SHARED WITH THE MODELS
##############################################################

# Using geopandas to make quick work of reprojection
# Create geodataframe
proxy_gdf = geopandas.GeoDataFrame(proxy,
                                   geometry=geopandas.points_from_xy(proxy.long, proxy.lat))

# Set CRS of data frame and I know that its 4326
proxy_gdf = proxy_gdf.set_crs(epsg=4326)

# Reproject to North America Albers Equal Area Conic
# EPSG: 102008
# https://epsg.io/102008
# Origin: 40 N, 96 W
# Standard parallels = 20 N, 60 N
# Datum: NAD83
proxy_gdf = proxy_gdf.to_crs('+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs')

# Adding new column of reprojected lat and lon
proxy_gdf['lat_reproject'] = proxy_gdf.geometry.y
proxy_gdf['long_reproject'] = proxy_gdf.geometry.x

# Drop geometry column because its not necessary
proxy_gdf.drop('geometry', axis=1, inplace=True)

######################################################################
# CALCULATE MODEL SKILL BY LOOPING THROUGH VALUES - NO CHANGE BASELINE
######################################################################

# Skill calculation in parallel on 12 cores, with 10000 replicates - 10 minutes
# List needs to be flattened because joblib feeds 1 list item to each core.
flat_clim = list()
for ii, model in enumerate(models.model.values):
    flat_clim.append({'values': models.sel(model=model).values, 'lat':models.lat.values, 'lon':models.lon.values, 'model_name':model})

skill_ens_list = Parallel(n_jobs=12)(delayed(parallel_skill)(i['values'], model_lat=i['lat'], model_lon=i['lon'], reps=10000, ni=0, proxy=proxy_gdf, model_name=i['model_name']) for i in flat_clim)

# Dataframe to hold skill
skill = pd.DataFrame(columns=['model', 'median_skill', 'lower', 'upper'])

for result in skill_ens_list:
    skill = skill.append({'model': result['model_name'],
                          'median_skill': np.median(result['skill_ens']),
                          'lower': np.quantile(result['skill_ens'], q=0.025),
                          'upper': np.quantile(result['skill_ens'], q=0.975)}, ignore_index=True)

#################################################################
# CALCULATE MODEL SKILL BY LOOPING THROUGH VALUES - MEAN BASELINE
#################################################################

# Skill calculation in parallel on 10 cores
# Takes ~7 minutes
skill_ens_list_mean_baseline = Parallel(n_jobs=12)(delayed(parallel_skill)(i['values'], model_lat=i['lat'], model_lon=i['lon'], reps=10000, ni=proxy_gdf['mean'].mean(), proxy=proxy_gdf, model_name=i['model_name']) for i in flat_clim)

# Dataframe to hold skill
skill_mean_baseline = pd.DataFrame(columns=['model', 'median_skill', 'lower', 'upper'])

for result in skill_ens_list_mean_baseline:
    skill_mean_baseline = skill_mean_baseline.append({'model': result['model_name'],
                                                      'median_skill': np.median(result['skill_ens']),
                                                      'lower': np.quantile(result['skill_ens'], q=0.025),
                                                      'upper': np.quantile(result['skill_ens'], q=0.975)}, ignore_index=True)
#############    
# SAVE AS CSV
#############

# Assign column to identify the two alternative hypotheses
skill['null'] = 'Mean Change Skill'
skill_mean_baseline['null'] = 'Spatial Configuration Skill'

# Combine the two data frames and save as csv
pd.concat([skill, skill_mean_baseline]).to_csv('skill_scores/tas_skill.csv', index=False)