"""
SCRIPT TAKES ALL MODEL OUTPUT FROM ~/MODEL_TEMPERATURE_ANOMALY AND 
~/MODEL_PRECIPITATION_ANOMALY AND BILINEARLY INTERPOLATES IT TO A COMMON 0.5 
DEGREE GRID DEFINED BY THE OUTPUT FROM TRACE-21KA (LORENZ ET AL., 2016, 
SCIENTIFIC DATA).

THE PRODUCTS OF THIS SCRIPT ARE NOT USED IN ANY ANALYSES AND ARE ONLY MEANT
TO ILLUSTRATE THE CLIMATIC FIELDS IN FIGURE 1. IN FACT, WE EXCLUDE THIS 
FROM SPATIAL GRAIN SENSTIVITY ANALYSES BECAUSE IT HAS NO BEARING ON THE FINAL
RESULTS, THOUGH THE CLIMATIC FIELDS ARE IMPORTANT FOR THE FINAL CONCLUSIONS.

INPUT:
    PRECIPITATION
        * MODEL_PRECIPITATION_ANOMALY/CCSM-MARUM_PR_ANOM.NC
        * MODEL_PRECIPITATION_ANOMALY/CM2MC_PR_ANOM.NC
        * MODEL_PRECIPITATION_ANOMALY/COSMOS-W_PR_ANOM.NC
        * MODEL_PRECIPITATION_ANOMALY/IPSL_PR_ANOM.NC
        * MODEL_PRECIPITATION_ANOMALY/MIROC-W_PR_ANOM.NC
        * MODEL_PRECIPITATION_ANOMALY/TRACE_MWF_PR_ANOM.NC
        * MODEL_PRECIPITATION_ANOMALY/CCSM-NCAR_PR_ANOM.NC
        * MODEL_PRECIPITATION_ANOMALY/COSMOS-S_PR_ANOM.NC
        * MODEL_PRECIPITATION_ANOMALY/HADCM3_PR_ANOM.NC
        * MODEL_PRECIPITATION_ANOMALY/MIROC-S_PR_ANOM.NC
        * MODEL_PRECIPITATION_ANOMALY/TRACE_LORENZ_PR_ANOM.NC
        * MODEL_PRECIPITATION_ANOMALY/TRACE_PR_ANOM.NC
    
    TEMPERATURE
        * MODEL_TEMPERATURE_ANOMALY/CCSM-MARUM_TAS_ANOM.NC
        * MODEL_TEMPERATURE_ANOMALY/CM2MC_TAS_ANOM.NC
        * MODEL_TEMPERATURE_ANOMALY/COSMOS-W_TAS_ANOM.NC
        * MODEL_TEMPERATURE_ANOMALY/IPSL_TAS_ANOM.NC
        * MODEL_TEMPERATURE_ANOMALY/MIROC-W_TAS_ANOM.NC
        * MODEL_TEMPERATURE_ANOMALY/TRACE_MWF_TAS_ANOM.NC
        * MODEL_TEMPERATURE_ANOMALY/CCSM-NCAR_TAS_ANOM.NC
        * MODEL_TEMPERATURE_ANOMALY/COSMOS-S_TAS_ANOM.NC
        * MODEL_TEMPERATURE_ANOMALY/HADCM3_TAS_ANOM.NC
        * MODEL_TEMPERATURE_ANOMALY/MIROC-S_TAS_ANOM.NC
        * MODEL_TEMPERATURE_ANOMALY/TRACE_LORENZ_TAS_ANOM.NC
        * MODEL_TEMPERATURE_ANOMALY/TRACE_TAS_ANOM.NC

OUTPUT:
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

MINICONDA ENVIRONMENT:
    PYTHON: 3.8.10
    XARRAY: 0.18.2
    RIOXARRAY: 0.5.0
    RASTERIO: 1.2.1
    NUMPY: 1.20.3

@author: David Fastovich
"""

import sys
sys.path.append('modules')
from bilinear_interpolation import gridded_bilinear_interpolation
import xarray as xr
import numpy as np
import glob

# Prepare the target grid from TraCE-21ka 0.5 degree simulation: 
#
# Latitude:
# [10.25 10.75 11.25 11.75 ... 78.25 78.75 79.25 79.75]
#
# Longitude:
# [-172.75, -172.25, -171.75, ...,  -49.25,  -48.75,  -48.25]

target_lat = np.arange(start=-89.75, stop=90, step=0.5)
target_lon = np.arange(start=-179.75, stop=180, step=0.5)

# Loop through temperature simulations and bilinearly interpolat2
for file in glob.glob('model_temperature_anomaly/*.nc'):
    # Read in model
    model = xr.open_dataset(file, decode_times=False).squeeze()
    
    # Renaming to common coordinates
    if 'lat' not in model.tas.dims:
        model = model.rename({'latitude':'lat', 'longitude':'lon'})
    
    # Bilinear interpolation function requires increasing y-values
    regridded_model = gridded_bilinear_interpolation(new_lat=target_lat, 
                                   new_lon=target_lon, 
                                   model_lat=model.lat.values, 
                                   model_lon=model.lon.values, 
                                   climate_var=model.tas.values.transpose(),
                                   varname='tas')
    
    # Save as netCDF4
    regridded_model.to_netcdf(('model_temperature_anomaly_05degree/' 
                               + file[26:len(file)-3] 
                               + '_05degree.nc'))

# Loop through precipitation simulations and bilinearly interpolat2
for file in glob.glob('model_precipitation_anomaly/*.nc'):
    # Read in model
    model = xr.open_dataset(file, decode_times=False).squeeze()
    
    # Renaming to common coordinates
    if 'lat' not in model.pr.dims:
        model = model.rename({'latitude':'lat', 'longitude':'lon'})
    
    # Bilinear interpolation function requires increasing y-values
    regridded_model = gridded_bilinear_interpolation(new_lat=target_lat, 
                                   new_lon=target_lon, 
                                   model_lat=model.lat.values, 
                                   model_lon=model.lon.values, 
                                   climate_var=model.pr.values.transpose(),
                                   varname='pr')
    
    # Save as netCDF4
    regridded_model.to_netcdf(('model_precipitation_anomaly_05degree/' 
                               + file[28:len(file)-3] 
                               + '_05degree.nc'))