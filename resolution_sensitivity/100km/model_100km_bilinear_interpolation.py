"""
**********THESE ARE SPATIAL GRAIN SENSITIVITY TESTS ON A 100KM GRID***********

SCRIPT TAKES ALL MODEL OUTPUT FROM ~/MODEL_TEMPERATURE_ANOMALY AND 
~/MODEL_PRECIPITATION_ANOMALY AND BILINEARLY INTERPOLATES TO A 50KM EQUAL AREA
GRID.

THE PRODUCTS OF THIS SCRIPT ARE THE PALEOCLIMATIC GRIDS USED WITHIN THE SPATIAL
ERROR MODELS THAT FORM THE FOUNDATION OF THE MANUSCRIPT.

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
        * RESOLUTION_SENSITIVITY/100KM/MODEL_PRECIPITATION_ANOMALY_100KM_ENA/CCSM-MARUM_PR_100KM_ENA.NC
        * RESOLUTION_SENSITIVITY/100KM/MODEL_PRECIPITATION_ANOMALY_100KM_ENA/CM2MC_PR_100KM_ENA.NC
        * RESOLUTION_SENSITIVITY/100KM/MODEL_PRECIPITATION_ANOMALY_100KM_ENA/COSMOS-W_PR_100KM_ENA.NC
        * RESOLUTION_SENSITIVITY/100KM/MODEL_PRECIPITATION_ANOMALY_100KM_ENA/IPSL_PR_100KM_ENA.NC
        * RESOLUTION_SENSITIVITY/100KM/MODEL_PRECIPITATION_ANOMALY_100KM_ENA/MIROC-W_PR_100KM_ENA.NC
        * RESOLUTION_SENSITIVITY/100KM/MODEL_PRECIPITATION_ANOMALY_100KM_ENA/TRACE_MWF_PR_100KM_ENA.NC
        * RESOLUTION_SENSITIVITY/100KM/MODEL_PRECIPITATION_ANOMALY_100KM_ENA/CCSM-NCAR_PR_100KM_ENA.NC
        * RESOLUTION_SENSITIVITY/100KM/MODEL_PRECIPITATION_ANOMALY_100KM_ENA/COSMOS-S_PR_100KM_ENA.NC
        * RESOLUTION_SENSITIVITY/100KM/MODEL_PRECIPITATION_ANOMALY_100KM_ENA/HADCM3_PR_100KM_ENA.NC
        * RESOLUTION_SENSITIVITY/100KM/MODEL_PRECIPITATION_ANOMALY_100KM_ENA/MIROC-S_PR_100KM_ENA.NC
        * RESOLUTION_SENSITIVITY/100KM/MODEL_PRECIPITATION_ANOMALY_100KM_ENA/TRACE_LORENZ_PR_100KM_ENA.NC
        * RESOLUTION_SENSITIVITY/100KM/MODEL_PRECIPITATION_ANOMALY_100KM_ENA/TRACE_PR_100KM_ENA.NC
    
    TEMPERATURE
        * RESOLUTION_SENSITIVITY/100KM/MODEL_TEMPERATURE_ANOMALY_100KM_ENA/CCSM-MARUM_TAS_100KM_ENA.NC
        * RESOLUTION_SENSITIVITY/100KM/MODEL_TEMPERATURE_ANOMALY_100KM_ENA/CM2MC_TAS_100KM_ENA.NC
        * RESOLUTION_SENSITIVITY/100KM/MODEL_TEMPERATURE_ANOMALY_100KM_ENA/COSMOS-W_TAS_100KM_ENA.NC
        * RESOLUTION_SENSITIVITY/100KM/MODEL_TEMPERATURE_ANOMALY_100KM_ENA/IPSL_TAS_100KM_ENA.NC
        * RESOLUTION_SENSITIVITY/100KM/MODEL_TEMPERATURE_ANOMALY_100KM_ENA/MIROC-W_TAS_100KM_ENA.NC
        * RESOLUTION_SENSITIVITY/100KM/MODEL_TEMPERATURE_ANOMALY_100KM_ENA/TRACE_MWF_TAS_100KM_ENA.NC
        * RESOLUTION_SENSITIVITY/100KM/MODEL_TEMPERATURE_ANOMALY_100KM_ENA/CCSM-NCAR_TAS_100KM_ENA.NC
        * RESOLUTION_SENSITIVITY/100KM/MODEL_TEMPERATURE_ANOMALY_100KM_ENA/COSMOS-S_TAS_100KM_ENA.NC
        * RESOLUTION_SENSITIVITY/100KM/MODEL_TEMPERATURE_ANOMALY_100KM_ENA/HADCM3_TAS_100KM_ENA.NC
        * RESOLUTION_SENSITIVITY/100KM/MODEL_TEMPERATURE_ANOMALY_100KM_ENA/MIROC-S_TAS_100KM_ENA.NC
        * RESOLUTION_SENSITIVITY/100KM/MODEL_TEMPERATURE_ANOMALY_100KM_ENA/TRACE_LORENZ_TAS_100KM_ENA.NC
        * RESOLUTION_SENSITIVITY/100KM/MODEL_TEMPERATURE_ANOMALY_100KM_ENA/TRACE_TAS_100KM_ENA.NC

MINICONDA ENVIRONMENT:
    PYTHON: 3.8.10
    XARRAY: 0.18.2
    RIOXARRAY: 0.5.0
    RASTERIO: 1.2.1
    NUMPY: 1.20.3

@author: David Fastovich
"""

import xarray as xr
import rioxarray
from rasterio.enums import Resampling
import numpy as np
import glob

# 100km grid from the origin of 40 N, 96 W
cell_width = 100000
target_lat = np.arange(-3000000, 3000000, cell_width)
target_lon = np.arange(-4000000, 4000000, cell_width)

# Use target_lat and target_lon to build target xarray for rioxarray to match
# rioxarray requires 'x' and 'y', not 'lon' and 'lat'
target_array = xr.DataArray(
    dims=['x', 'y'],
    coords=dict(
        x=(['x'], target_lon),
        y=(['y'], target_lat),
    ),
)

# Adding CRS to target array
# Projection: North America Albers Equal Area Conic
# EPSG: 102008
# https://epsg.io/102008
# Origin: 40 N, 96 W
# Standard parallels = 20 N, 60 N
# Datum: NAD83
target_array = target_array.rio.write_crs('+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs')

# Loop through temperature simulations and bilinearly interpolate onto 50km
# grid
for file in glob.glob('model_temperature_anomaly/*.nc'):
    # Read in model
    model = xr.open_dataset(file, decode_times=False).squeeze()
    
    # Add spherical earth CRS to model output since all models assume a spherical
    # earth.
    # This is the Proj4 string for cartopy.crs.Geodetic()
    model = model.tas.rio.write_crs('+proj=longlat +datum=WGS84 +no_defs')
    
    # Renaming to coordinate names rioxarray expects
    if 'lat' not in model.dims:
        model = model.rename({'latitude':'y', 'longitude':'x'})
    else:
        model = model.rename({'lat':'y', 'lon':'x'})
    
    # Resample onto target 50km grid using bilinear interpolation
    # This also crops data to the target region
    ena_regridded_model = model.rio.reproject_match(target_array, resampling=Resampling.bilinear)
        
    # Rename back to cf compliant coordinate names
    ena_regridded_model = ena_regridded_model.rename({'x':'lon',
                                                      'y':'lat'})
    
    # Save as netCDF4
    ena_regridded_model.to_netcdf(('resolution_sensitivity/100km/model_temperature_anomaly_100km_ena/' 
                                   + file[26:len(file)-3] 
                                   + '_100km_ena.nc'))

# Loop through precipitation simulations and bilinearly interpolate onto 50km
# grid
for file in glob.glob('model_precipitation_anomaly/*.nc'):
    # Read in model
    model = xr.open_dataset(file, decode_times=False).squeeze()
    
    # Add spherical earth CRS to model output since all models assume a spherical
    # earth.
    # This is the Proj4 string for cartopy.crs.Geodetic()
    model = model.pr.rio.write_crs('+proj=longlat +datum=WGS84 +no_defs')
    
    # Renaming to coordinate names rioxarray expects
    if 'lat' not in model.dims:
        model = model.rename({'latitude':'y', 'longitude':'x'})
    else:
        model = model.rename({'lat':'y', 'lon':'x'})
    
    # Resample onto target 50km grid using bilinear interpolation
    # This also crops data to the target region
    ena_regridded_model = model.rio.reproject_match(target_array, resampling=Resampling.bilinear)
        
    # Rename back to cf compliant coordinate names
    ena_regridded_model = ena_regridded_model.rename({'x':'lon',
                                                      'y':'lat'})
    
    # Save as netCDF4
    ena_regridded_model.to_netcdf(('resolution_sensitivity/100km/model_precipitation_anomaly_100km_ena/' 
                                   + file[28:len(file)-3] 
                                   + '_100km_ena.nc'))
