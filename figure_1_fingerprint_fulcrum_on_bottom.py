"""
SCRIPT TAKES ALL MODEL OUTPUT FROM ~/MODEL_TEMPERATURE_ANOMALY AND 
~/MODEL_PRECIPITATION_ANOMALY AND PLOTS THE TEMPERATURE AND PRECIPITATION FIELDS
FROM THE CLIMATE SIMULATIONS. ADDITIONALLY, FIGURE 1C IS ADDED BY TAKING 1000
YEAR BINS AND REGRESSING TEMPERATURE AS A FUNCTION OF LATITUDE.

IMPORTANT: NO ANALYSES ARE PERFORMED WITH THIS DATA. THIS IS STRICTLY CREATED
FOR PLOTTING PURPOSES FOR FIGURE 1.

INPUT:
    PRECIPITATION
      * MODEL_PRECIPITATION_ANOMALY_05DEGREE/CCSM-MARUM_PR_05DEGREE.NC
      * MODEL_PRECIPITATION_ANOMALY_05DEGREE/CM2MC_PR_05DEGREE.NC
      * MODEL_PRECIPITATION_ANOMALY_05DEGREE/COSMOS-W_PR_05DEGREE.NC
      * MODEL_PRECIPITATION_ANOMALY_05DEGREE/IPSL_PR_05DEGREE.NC
      * MODEL_PRECIPITATION_ANOMALY_05DEGREE/MIROC-W_PR_05DEGREE.NC
      * MODEL_PRECIPITATION_ANOMALY_05DEGREE/TRACE_MWF_PR_05DEGREE.NC
      * MODEL_PRECIPITATION_ANOMALY_05DEGREE/CCSM-NCAR_PR_05DEGREE.NC
      * MODEL_PRECIPITATION_ANOMALY_05DEGREE/COSMOS-S_PR_05DEGREE.NC
      * MODEL_PRECIPITATION_ANOMALY_05DEGREE/HADCM3_PR_05DEGREE.NC
      * MODEL_PRECIPITATION_ANOMALY_05DEGREE/MIROC-S_PR_05DEGREE.NC
      * MODEL_PRECIPITATION_ANOMALY_05DEGREE/TRACE_LORENZ_PR_05DEGREE.NC
      * MODEL_PRECIPITATION_ANOMALY_05DEGREE/TRACE_PR_05DEGREE.NC
    
    TEMPERATURE
      * MODEL_TEMPERATURE_ANOMALY_05DEGREE/CCSM-MARUM_TAS_05DEGREE.NC
      * MODEL_TEMPERATURE_ANOMALY_05DEGREE/CM2MC_TAS_05DEGREE.NC
      * MODEL_TEMPERATURE_ANOMALY_05DEGREE/COSMOS-W_TAS_05DEGREE.NC
      * MODEL_TEMPERATURE_ANOMALY_05DEGREE/IPSL_TAS_05DEGREE.NC
      * MODEL_TEMPERATURE_ANOMALY_05DEGREE/MIROC-W_TAS_05DEGREE.NC
      * MODEL_TEMPERATURE_ANOMALY_05DEGREE/TRACE_MWF_TAS_05DEGREE.NC
      * MODEL_TEMPERATURE_ANOMALY_05DEGREE/CCSM-NCAR_TAS_05DEGREE.NC
      * MODEL_TEMPERATURE_ANOMALY_05DEGREE/COSMOS-S_TAS_05DEGREE.NC
      * MODEL_TEMPERATURE_ANOMALY_05DEGREE/HADCM3_TAS_05DEGREE.NC
      * MODEL_TEMPERATURE_ANOMALY_05DEGREE/MIROC-S_TAS_05DEGREE.NC
      * MODEL_TEMPERATURE_ANOMALY_05DEGREE/TRACE_LORENZ_TAS_05DEGREE.NC
      * MODEL_TEMPERATURE_ANOMALY_05DEGREE/TRACE_TAS_05DEGREE.NC
        
    PROXY TEMPERATURE
      * PROXY_DATA/TAVE_ANOMALIES.CSV
      * PROXY_DATA/GDGT_MAT.CSV
      * PROXY_DATA/ENS_TAVE_BINNED1000.CSV
    
    PROXY PRECIPITATION
      * PROXY_DATA/ANNP_ANOMALIES.CSV
      
OUTPUT:
  FIGURES/FINGERPRINT_FULCRUM_ON_BOTTOM.PDF

MINICONDA ENVIRONMENT:
  PYTHON: 3.8.10
  MATPLOTLIB: 3.5.1
  XARRAY: 0.18.2
  RIOXARRAY: 0.5.0
  NUMPY: 1.20.3
  PANDAS: 1.2.4
  CARTOPY: 0.20.2
  GLOB: INCLUDED WITH PYTHON

@author: David Fastovich
"""

import glob
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker as mticker
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from sklearn.linear_model import LinearRegression
import xarray as xr
import rioxarray
import numpy as np
import pandas as pd
import cartopy.crs as ccrs
import cartopy
import cartopy.mpl.geoaxes
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

####################################################
# READ IN MODELS FOR TEMPERATURE - PLATE CARREE GRID
####################################################

# Models to read in
model_tas_filenames = glob.glob('model_temperature_anomaly_05degree/*.nc')

# Excluding TraCE-21KA Lorenz since that has data only over North America, but
# the full TraCE simulation is included
model_tas_filenames = [item for item in model_tas_filenames if 'LORENZ' not in item]

# Create an empty xarray that will serve as the extra dimension for indexing
# the models into a single xarray file
model_tas_concat_xr = xr.DataArray(model_tas_filenames, dims=['model'], name='model')

# Empty list to hold the models
models_tas = list()

# Loop through and append each model onto the empty list
for file in model_tas_filenames:
    models_tas.append(xr.open_dataset(file)['tas'])

# Use the empty xarray to reshape the list into a single xarray where the 
# 'model' dimension is associated with each unique model
models_tas = xr.concat(models_tas, dim=model_tas_concat_xr)

# Take model mean and model standard deviation
model_tas_mean = models_tas.mean(dim='model')
model_tas_sd = models_tas.std(dim='model')

######################################################
# READ IN MODELS FOR PRECIPITATION - PLATE CARREE GRID
######################################################

# Models to read in
model_pr_filenames = glob.glob('model_precipitation_anomaly_05degree/*.nc')

# Excluding TraCE-21KA Lorenz since that has data only over North America, but
# the full TraCE simulation is included
model_pr_filenames = [item for item in model_pr_filenames if 'LORENZ' not in item]

# Create an empty xarray that will serve as the extra dimension for indexing
# the models into a single xarray file
model_pr_concat_xr = xr.DataArray(model_pr_filenames, dims=['model'], name='model')

# Empty list to hold the models
models_pr = list()

# Loop through and append each model onto the empty list
for file in model_pr_filenames:
    models_pr.append(xr.open_dataset(file)['pr'])

# Use the empty xarray to reshape the list into a single xarray where the 
# 'model' dimension is associated with each unique model
models_pr = xr.concat(models_pr, dim=model_pr_concat_xr)

# Take model mean and model standard deviation
model_pr_mean = models_pr.mean(dim='model')
model_pr_sd = models_pr.std(dim='model')

##################################
# READ IN PROXY DATA - TEMPERATURE
##################################

proxy_tave = pd.read_csv('proxy_data/tave_anomalies.csv')
proxy_tave['proxy'] = 'Fossil-Pollen'
gdgt = pd.read_csv('proxy_data/gdgt_mat.csv')
gdgt['proxy'] = 'brGDGT'
proxy_tave = pd.concat([proxy_tave, gdgt])

# Drop NAs
proxy_tave = proxy_tave.dropna(axis=0)

####################################
# READ IN PROXY DATA - PRECIPITATION
####################################

proxy_pr = pd.read_csv('proxy_data/annp_anomalies.csv')
proxy_pr['proxy'] = 'Fossil-Pollen'

# Changing units to mm/day from mm/year by dividing by 365.25 days
proxy_pr['mean'] = proxy_pr['mean']/365.25
proxy_pr['sd'] = proxy_pr['sd']/365.25

# Drop NAs
proxy_pr = proxy_pr.dropna(axis=0)

#####################################
# READ IN PROXY DATA - 1000 YEAR BINS
#####################################

# From results of Fastovich et al., (2020, Geophysical Research Letters)

kilo_bins = pd.read_csv('proxy_data/ens_tave_binned1000.csv')

###############
# SET FONT SIZE
###############

SMALL_SIZE = 20
MEDIUM_SIZE = 22
BIGGER_SIZE = 24
plt.rc('font', size=SMALL_SIZE)  # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)  # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)  # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
plt.rc('legend', fontsize=BIGGER_SIZE)  # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
plt.rcParams["axes.prop_cycle"] = plt.cycler("color", plt.cm.tab20.colors) # Configure the default colors used for the spaghetti plot

############################
# PLOT DATA - SCATTER ON TOP
############################

(13.2/12)*12.75
#%%

# Marker dictionary, plot title
proxy_marker_dict = {'Fossil-Pollen': 's', 'brGDGT': '^'}

# GridSpec configuration
fig = plt.figure(figsize=(12.75, 14.025), constrained_layout=True)
gs = gridspec.GridSpec(nrows=3, ncols=2, figure=fig, height_ratios=(1, 0.05, 1), width_ratios=(1, 1))

# AX1 - Mean Temperature -----------------------------------------------------

ax1 = fig.add_subplot(gs[0, 0], projection=ccrs.PlateCarree())
cb = ax1.contourf(model_tas_mean.lon,
                  model_tas_mean.lat,
                  model_tas_mean.values,
                  cmap=plt.get_cmap('seismic'),
                  levels=np.arange(-10, 2, 0.5),
                  vmin=-10,
                  vmax=10,
                  extend='min',
                  transform=ccrs.PlateCarree())
ax1.set_title('Mean Annual\nTemperature Anomaly')
ax1.set_extent([-95, -51, 16.5, 61.5], crs=ccrs.PlateCarree())
ax1.coastlines()
plt.text(-95, 62, 'A', weight='bold')

# Scatter points on top of map
for kind in proxy_marker_dict:
    d = proxy_tave.loc[proxy_tave['proxy'] == kind,:]
    ax1.scatter(x=d['long'],
               y=d['lat'],
               c=d['mean'],
               s=100,
               marker=proxy_marker_dict[kind],
               edgecolors='black',
               cmap=plt.get_cmap('seismic'),
               vmin=-10,
               vmax=10,
               label=kind,
               transform=ccrs.PlateCarree())

ax1_ins = inset_axes(ax1, width=3.2, height=3.2, loc='center', bbox_to_anchor=(534, 1038.5), axes_class=cartopy.mpl.geoaxes.GeoAxes, axes_kwargs=dict(map_projection=cartopy.crs.PlateCarree()))
cb = ax1_ins.contourf(model_tas_mean.lon,
                  model_tas_mean.lat,
                  model_tas_mean.values,
                  cmap=plt.get_cmap('seismic'),
                  levels=np.arange(-10, 2, 0.5),
                  vmin=-10,
                  vmax=10,
                  extend='min',
                  transform=ccrs.PlateCarree())
ax1_ins.coastlines()

# Add ticks for latitude and longitude
ax1.tick_params(top=True, labeltop=True, bottom=True, labelbottom=True, right=True, labelright=True)
ax1.set_xticks(np.arange(-90, -50, 5), crs=ccrs.PlateCarree())
ax1.set_yticks(np.arange(20, 65, 5), crs=ccrs.PlateCarree())

# Configure and add colorbar
cb_ax = fig.add_subplot(gs[1, 0])
colorbar = plt.colorbar(cb, orientation='horizontal', cax=cb_ax)
colorbar.set_label('Multi-Model Mean (K)', fontsize=16)

# AX2 - Mean Precipitation ---------------------------------------------------

ax2 = fig.add_subplot(gs[0, 1], projection=ccrs.PlateCarree())
cb = ax2.contourf(model_pr_mean.lon,
                  model_pr_mean.lat,
                  model_pr_mean.values,
                  cmap=plt.get_cmap('BrBG'),
                  levels=np.arange(-1.2, 1.3, 0.1),
                  extend='both',
                  transform=ccrs.PlateCarree())
ax2.set_title('Mean Annual\nPrecipitation Anomaly')
ax2.set_extent([-95, -51, 16.5, 61.5], crs=ccrs.PlateCarree())
ax2.coastlines()
plt.text(-95, 62, 'B', weight='bold')

# Scatter points on top of map
ax2.scatter(x=proxy_pr['long'],
           y=proxy_pr['lat'],
           c=proxy_pr['mean'],
           s=100,
           marker='s',
           edgecolors='black',
           cmap=plt.get_cmap('BrBG'),
           vmin=-1.2,
           vmax=1.2,
           label="Fossil-Pollen",
           transform=ccrs.PlateCarree())

ax2_ins = inset_axes(ax2, width=3.2, height=3.2, loc='center', bbox_to_anchor=(1272.5, 1038.5), axes_class=cartopy.mpl.geoaxes.GeoAxes, axes_kwargs=dict(map_projection=cartopy.crs.PlateCarree()))
cb = ax2_ins.contourf(model_pr_mean.lon,
                  model_pr_mean.lat,
                  model_pr_mean.values,
                  cmap=plt.get_cmap('BrBG'),
                  levels=np.arange(-1.2, 1.2, 0.1),
                  extend='both',
                  transform=ccrs.PlateCarree())
ax2_ins.coastlines()

# Add ticks for latitude and longitude
ax2.tick_params(top=True, labeltop=True, bottom=True, labelbottom=True, right=True, labelright=True)
ax2.set_xticks(np.arange(-90, -50, 5), crs=ccrs.PlateCarree())
ax2.set_yticks(np.arange(20, 65, 5), crs=ccrs.PlateCarree())

# Configure and add colorbar
cb_ax = fig.add_subplot(gs[1, 1])
colorbar = plt.colorbar(cb, orientation='horizontal', cax=cb_ax)
colorbar.set_label('Multi-Model Mean (mm/day)', fontsize=16)

# AX3 - Latitudinal Gradient -------------------------------------------------
col = ['#A50026', '#D73027', '#F46D43', '#FDAE61', '#FEE090', '#FFFFBF', '#E0F3F8', '#ABD9E9', '#74ADD1', '#4575B4', '#313695']
ax3 = fig.add_subplot(gs[2, 0:2])

# Perform linear regressions of temperature ~ latitude for each 1000 year bin
# and add to AX3
for ii, age_bin in enumerate(kilo_bins['age'].unique()[(kilo_bins['age'].unique() <= 18500) & (kilo_bins['age'].unique() >= 8000)]):
    # Initialize linear regression instance
    lin = LinearRegression()
    
    # Fit the instance
    lin.fit(kilo_bins['lat'][kilo_bins['age'] == age_bin].to_numpy().reshape(-1, 1), kilo_bins['tave'][kilo_bins['age'] == age_bin])
    
    # Plot
    ax3.scatter(x=kilo_bins['lat'][kilo_bins['age'] == age_bin], y=kilo_bins['tave'][kilo_bins['age'] == age_bin], color=col[ii])
    ax3.plot(np.arange(25, 50, 1), lin.predict(np.arange(25, 50, 1).reshape(-1, 1)), c=col[ii], label=np.floor(age_bin/1000).astype('int'), lw=5)
ax3.set_ylabel('Temperature ($^{\circ}$C)')
ax3.set_xlabel('Latitude ($^{\circ}$)')
ax3.set_xlim(26, 49)
plt.text(26.25, 22.25, 'C', weight='bold')

# First legend for lines
legend1 = plt.legend(bbox_to_anchor=(1.1, 0.5), loc='center', frameon=False, title='Age (ka)', title_fontproperties={'weight':'bold'})

leg = plt.gca().add_artist(legend1)
leg._legend_box.align = "left"
plt.savefig('figures/fingerprint_fulcrum_on_bottom.pdf', bbox_inches='tight', dpi=300)
