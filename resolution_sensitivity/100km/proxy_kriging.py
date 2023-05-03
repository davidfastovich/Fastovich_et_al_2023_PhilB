"""
**********THESE ARE SPATIAL GRAIN SENSITIVITY TESTS ON A 100KM GRID***********

OUR PROXY RECONSTRUCTIONS ARE POINT ESTIMATES FROM FASTOVICH ET AL., (2022, 
QUATERNARY SCIENCE REVIEWS) AND OUR ANALYSES REQUIRE FULL CLIMATE FIELDS.
THEREFORE, WE USE UNIVERSAL KRIGING TO SPATIALLY INTERPOLATE OUR PROXY
RECONSTRUCTIONS TO A FULL SPATIAL FIELD. SILL, RANGE, AND NUGGET PARAMETERS
ARE TAKEN FROM VARIOGRAM ESTIMATES OF THE GPCC AND NSEP REANALYSIS DATASETS.

INPUT:
    PRECIPITATION
        * PROXY_DATA/ANNP_ANOMALIES.CSV
    
    TEMPERATURE
        * PROXY_DATA/TAVE_ANOMALIES.CSV
        * PROXY_DATA/GDGT_MAT.CSV

OUTPUT:
    PRECIPITATION
        * RESOLUTION_SENSITIVITY/100KM/MODEL_PRECIPITATION_ANOMALY_100KM_ENA/PROXY_KRIGING_PR_ANOM_100KM_ENA.NC
    
    TEMPERATURE
        * RESOLUTION_SENSITIVITY/100KM/MODEL_TEMPERATURE_ANOMALY_100KM_ENA/PROXY_KRIGING_TAS_ANOM_100KM_ENA.NC

MINICONDA ENVIRONMENT:
    PYTHON: 3.8.10
    PANDAS: 1.2.4
    XARRAY: 0.18.2
    RIOXARRAY: 0.5.0
    GEOPANDAS: 0.10.2
    PYKRIGE: 1.6.1

@author: David Fastovich
"""

import pandas as pd
import numpy as np
import xarray as xr
import rioxarray
import geopandas
from pykrige.uk import UniversalKriging

###########################################
# READ IN PROXY TEMPERATURE RECONSTRUCTIONS
###########################################

proxy_tas = pd.read_csv('proxy_data/tave_anomalies.csv')
gdgt = pd.read_csv('proxy_data/gdgt_mat.csv')
proxy_tas = pd.concat([proxy_tas, gdgt])

# Drop NAs
proxy_tas = proxy_tas.dropna(axis=0)

##############################################################
# REPROJECT PROXY DATA TO THE 50KM GRID SHARED WITH THE MODELS
##############################################################

# Using geopandas to make quick work of reprojection
# Create geodataframe
proxy_tas_gdf = geopandas.GeoDataFrame(proxy_tas,
                                   geometry=geopandas.points_from_xy(proxy_tas.long, proxy_tas.lat))

# Set CRS of data frame and I know that its 4326
proxy_tas_gdf = proxy_tas_gdf.set_crs(epsg=4326)

# Reproject to North America Albers Equal Area Conic
# EPSG: 102008
# https://epsg.io/102008
# Origin: 40 N, 96 W
# Standard parallels = 20 N, 60 N
# Datum: NAD83
proxy_tas_gdf = proxy_tas_gdf.to_crs('+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs')

# Adding new column of reprojected lat and lon
proxy_tas_gdf['lat_reproject'] = proxy_tas_gdf.geometry.y
proxy_tas_gdf['long_reproject'] = proxy_tas_gdf.geometry.x

# Drop geometry column because its not necessary
proxy_tas_gdf.drop('geometry', axis=1, inplace=True)

###########################################
# READ IN PROXY TEMPERATURE RECONSTRUCTIONS
###########################################

proxy_pr = pd.read_csv('proxy_data/annp_anomalies.csv')

# Changing units to mm/day from mm/year by dividing by 365.25 days
proxy_pr['mean'] = proxy_pr['mean']/365.25
proxy_pr['sd'] = proxy_pr['sd']/365.25

# Drop NAs
proxy_pr = proxy_pr.dropna(axis=0)

##############################################################
# REPROJECT PROXY DATA TO THE 50KM GRID SHARED WITH THE MODELS
##############################################################

# Using geopandas to make quick work of reprojection
# Create geodataframe
proxy_pr_gdf = geopandas.GeoDataFrame(proxy_pr,
                                   geometry=geopandas.points_from_xy(proxy_pr.long, proxy_pr.lat))

# Set CRS of data frame and I know that its 4326
proxy_pr_gdf = proxy_pr_gdf.set_crs(epsg=4326)

# Reproject to North America Albers Equal Area Conic
# EPSG: 102008
# https://epsg.io/102008
# Origin: 40 N, 96 W
# Standard parallels = 20 N, 60 N
# Datum: NAD83
proxy_pr_gdf = proxy_pr_gdf.to_crs('+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs')

# Adding new column of reprojected lat and lon
proxy_pr_gdf['lat_reproject'] = proxy_pr_gdf.geometry.y
proxy_pr_gdf['long_reproject'] = proxy_pr_gdf.geometry.x

# Drop geometry column because its not necessary
proxy_pr_gdf.drop('geometry', axis=1, inplace=True)

#########
# KRIGING
#########

# The same grid that everything has been interpolated to
cell_width = 100000
Y = np.arange(-3000000, 3000000, cell_width).astype('float32')
X = np.arange(-4000000, 4000000, cell_width).astype('float32')

# Model fit using automap in R
UK_pr = UniversalKriging(
    proxy_pr_gdf['long_reproject'],
    proxy_pr_gdf['lat_reproject'],
    proxy_pr_gdf['mean'],
    variogram_model="exponential",
    variogram_parameters={'sill': 301813, 'range': 892257, 'nugget': 0}, # Based on GPCC reanalysis data over ENA from 1951 to 1980
    drift_terms=["regional_linear"],
)

UK_pr_pred, UK_pr_ss = UK_pr.execute("grid", X, Y)

# Model fit using automap in R
UK_tas = UniversalKriging(
    proxy_tas_gdf['long_reproject'],
    proxy_tas_gdf['lat_reproject'],
    proxy_tas_gdf['mean'],
    variogram_model="gaussian",    
    variogram_parameters={'sill': 3.5, 'range': 791996, 'nugget': 1.5}, # Based on NCEP reanalysis data over ENA from 1951 to 1980
    drift_terms=["regional_linear"],
)
UK_tas_pred, UK_tas_ss = UK_tas.execute("grid", X, Y)

##################################
# CREATE XARRAY AND SAVE AS NETCDF
##################################

tas_ds = xr.Dataset(
    data_vars=dict(
        tas=(["lat", "lon"], UK_tas_pred),
    ),
    coords=dict(
        lon=(["lon"], X),
        lat=(["lat"], Y),
    ),
)
tas_ds = tas_ds.rio.write_crs('+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs')
tas_ds.to_netcdf("resolution_sensitivity/100km/model_temperature_anomaly_100km_ena/PROXY_KRIGING_tas_anom_100km_ena.nc")

pr_ds = xr.Dataset(
    data_vars=dict(
        pr=(["lat", "lon"], UK_pr_pred),
    ),
    coords=dict(
        lon=(["lon"], X),
        lat=(["lat"], Y),
    ),
)
pr_ds = pr_ds.rio.write_crs('+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs')
pr_ds.to_netcdf("resolution_sensitivity/100km/model_precipitation_anomaly_100km_ena/PROXY_KRIGING_pr_anom_100km_ena.nc")