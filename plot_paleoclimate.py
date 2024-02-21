"""
SCRIPT PLOTS ALL MODEL OUTPUT FROM  ~/MODEL_PRECIPITATION_ANOMALY_50KM_ENA AND 
~/MODEL_PRECIPITATION_ANOMALY_50KM_ENA.

INPUT:
    PRECIPITATION
      * MODEL_PRECIPITATION_ANOMALY_50KM_ENA/CCSM-MARUM_PR_50KM_ENA.NC
      * MODEL_PRECIPITATION_ANOMALY_50KM_ENA/CM2MC_PR_50KM_ENA.NC
      * MODEL_PRECIPITATION_ANOMALY_50KM_ENA/COSMOS-W_PR_50KM_ENA.NC
      * MODEL_PRECIPITATION_ANOMALY_50KM_ENA/IPSL_PR_50KM_ENA.NC
      * MODEL_PRECIPITATION_ANOMALY_50KM_ENA/MIROC-W_PR_50KM_ENA.NC
      * MODEL_PRECIPITATION_ANOMALY_50KM_ENA/TRACE_MWF_PR_50KM_ENA.NC
      * MODEL_PRECIPITATION_ANOMALY_50KM_ENA/CCSM-NCAR_PR_50KM_ENA.NC
      * MODEL_PRECIPITATION_ANOMALY_50KM_ENA/COSMOS-S_PR_50KM_ENA.NC
      * MODEL_PRECIPITATION_ANOMALY_50KM_ENA/HADCM3_PR_50KM_ENA.NC
      * MODEL_PRECIPITATION_ANOMALY_50KM_ENA/MIROC-S_PR_50KM_ENA.NC
      * MODEL_PRECIPITATION_ANOMALY_50KM_ENA/TRACE_LORENZ_PR_50KM_ENA.NC
      * MODEL_PRECIPITATION_ANOMALY_50KM_ENA/TRACE_PR_50KM_ENA.NC
    
    TEMPERATURE
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
      
OUTPUT:
  FIGURES/ALL_MODEL_TAS.PDF
  FIGURES/ALL_MODEL_PR.PDF

MINICONDA ENVIRONMENT:
    PYTHON: 3.8.10
    XARRAY: 0.18.2
    NUMPY: 1.20.3
    CARTOPY: 0.18.0
    MATPLOTLIB: 3.4.2
    GLOB: BUILT INTO PYTHON

@author: David Fastovich
"""

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import xarray as xr
import numpy as np
import pandas as pd
import cartopy.crs as ccrs
import glob
import cartopy.feature as cfeature

##################################
# READ IN MODELS FOR PRECIPITATION
##################################

# Models to read in
model_pr_filenames = glob.glob('model_precipitation_anomaly_50km_ena/*.nc')

# Create an empty xarray that will serve as the extra dimension for indexing
# the models into a single xarray file
model_pr_concat_xr = xr.DataArray(model_pr_filenames, dims=['model'], name='model')

# Empty list to hold the models
models_pr = list()
    
# Loop through and append each model onto the empty list
for ii, file in enumerate(model_pr_filenames):
    # if 'LORENZ' in file:
        # Fill Lorenz with NAs
        # models_pr.append(xr.open_dataset(file, decode_times=False)['pr'].fillna(-99999))
    # else:
        models_pr.append(xr.open_dataset(file, decode_times=False)['pr'])

# Use the empty xarray to reshape the list into a single xarray where the 
# 'model' dimension is associated with each unique model
models_pr = xr.concat(models_pr, dim=model_pr_concat_xr)

# Stripping name for plotting
models_pr = models_pr.assign_coords({'model':np.array(models_pr.model.to_dataframe()['model'].str.replace('model_precipitation_anomaly_50km_ena/','').str.replace("_pr_anom_50km_ena.nc","").to_list())})

# Alphebetize
models_pr = models_pr.sortby('model')

################################
# READ IN MODELS FOR TEMPERATURE
################################

# Models to read in
model_tas_filenames = glob.glob('model_temperature_anomaly_50km_ena/*.nc')

# Create an empty xarray that will serve as the extra dimension for indexing
# the models into a single xarray file
model_tas_concat_xr = xr.DataArray(model_tas_filenames, dims=['model'], name='model')

# Empty list to hold the models
models_tas = list()
    
# Loop through and append each model onto the empty list
for ii, file in enumerate(model_tas_filenames):
    # if 'LORENZ' in file:
        # Fill Lorenz with NAs
        # models_tas.append(xr.open_dataset(file, decode_times=False)['tas'].fillna(-99999))
    # else:
        models_tas.append(xr.open_dataset(file, decode_times=False)['tas'])

# Use the empty xarray to reshape the list into a single xarray where the 
# 'model' dimension is associated with each unique model
models_tas = xr.concat(models_tas, dim=model_tas_concat_xr)

# Stripping name for plotting
models_tas = models_tas.assign_coords({'model':np.array(models_tas.model.to_dataframe()['model'].str.replace('model_temperature_anomaly_50km_ena/','').str.replace("_tas_anom_50km_ena.nc","").to_list())})

# Alphebetize
models_tas = models_tas.sortby('model')

###############
# SET FONT SIZE
###############

SMALL_SIZE = 16
MEDIUM_SIZE = 18
BIGGER_SIZE = 20
plt.rc('font', size=SMALL_SIZE)  # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)  # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)  # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
plt.rc('legend', fontsize=BIGGER_SIZE)  # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

######
# PLOT
######

#%%

# Temperature plot
proj = ccrs.AlbersEqualArea(central_longitude=-96, central_latitude=40, standard_parallels=(20.0, 60.0))

# GridSpec configuration
fig = plt.figure(figsize=(13, 10), constrained_layout=True)
gs = gridspec.GridSpec(nrows=5, ncols=4, figure=fig, height_ratios=(1, 1, 1, 1, 0.1))

for ii, model in enumerate(models_tas.model.values):
    ax = fig.add_subplot(gs[ii], projection=proj)
    cb = ax.contourf(models_tas.lon,
                        models_tas.lat,
                        models_tas.sel(model=model).values,
                        cmap=plt.get_cmap('seismic'),
                        levels=np.arange(-10, 2, 0.5),
                        vmin=-10,
                        vmax=10,
                        extend='both'
                        )
    ax.set_title(models_tas.sel(model=model).model.values)
    ax.set_extent([-120, -70, 26, 50], crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.COASTLINE, zorder=2, edgecolor='#000000')

cb_ax = fig.add_subplot(gs[4,1:3])
cbar = plt.colorbar(cb, orientation='horizontal', cax=cb_ax)
cbar.set_label('Temperature Anomaly (K)')

plt.savefig('figures/all_tas_models.pdf', bbox_inches='tight', dpi=300)
plt.savefig('figures/all_tas_models.png', bbox_inches='tight', dpi=300)


#%%

# Precipitation plot
# GridSpec configuration
fig = plt.figure(figsize=(13, 10), constrained_layout=True)
gs = gridspec.GridSpec(nrows=5, ncols=4, figure=fig, height_ratios=(1, 1, 1, 1, 0.1))

for ii, model in enumerate(models_pr.model.values):
    ax = fig.add_subplot(gs[ii], projection=proj)
    cb = ax.contourf(models_pr.lon,
                        models_pr.lat,
                        models_pr.sel(model=model).values,
                        cmap=plt.cm.get_cmap('BrBG'),
                        transform=proj,
                        levels=np.arange(-1.2, 1.3, 0.1),
                        extend='both'
                        )
    ax.set_title(models_pr.sel(model=model).model.values)
    ax.set_extent([-120, -70, 26, 50], crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.COASTLINE, zorder=2, edgecolor='#000000')

cb_ax = fig.add_subplot(gs[4,1:3])
cbar = plt.colorbar(cb, orientation='horizontal', cax=cb_ax)
cbar.set_label('Precipitation Anomaly (mm/day)')

plt.savefig('figures/all_pr_models.pdf', bbox_inches='tight', dpi=300)
plt.savefig('figures/all_pr_models.png', bbox_inches='tight', dpi=300)
