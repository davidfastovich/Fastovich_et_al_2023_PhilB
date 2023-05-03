"""
SCRIPT TAKES THE MODERN OBERSVATIONAL DATA OF LAND/SEA MASK, TEMPEATURE, AND
PRECIPITATION AND REPROJECTS IT TO A SHARED 50KM GRID.

INPUT:
    * MODERN_DATA/GLDASP5_LANDMASK_025D.NC4
    * MODERN_DATA/MODERN_DATA/WC2.1_10M_PREC/*.TIF
    * MODERN_DATA/MODERN_DATA/WC2.1_10M_TAVG/*.TIF

OUTPUT:
    * RESOLUTION_SENSITIVITY/100KM/MODERN_DATA/POSTPROCESSED/GLDAS_50KM_ALBERS.NC
    * RESOLUTION_SENSITIVITY/100KM/MODERN_DATA/POSTPROCESSED/WC_PR_ANNUAL_50KM_ENA.TIF
    * RESOLUTION_SENSITIVITY/100KM/MODERN_DATA/POSTPROCESSED/WC_TAVE_ANNUAL_50KM_ENA.TIF

MINICONDA ENVIRONMENT:
    PYTHON: 3.8.10
    PANDAS: 1.2.4
    XARRAY: 0.18.2
    RIOXARRAY: 0.5.0
    RASTERIO: 1.2.1
    NUMPY: 1.20.3
    GLOB: BUILT INTO PYTHON

@author: David Fastovich
"""

import xarray as xr
import rioxarray
from rasterio.enums import Resampling
import numpy as np
import glob


# 50km grid from the origin of 40 N, 96 W
target_lat = np.arange(-3000000, 3000000, 50000)
target_lon = np.arange(-4000000, 4000000, 50000)

# Use target_lat and target_lon to build target xarray for rioxarray to match
# rioxarray requires 'x' and 'y', not 'lon' and 'lat'
target_array = xr.DataArray(
    data=np.random.rand(160,120), # filling it with random data because why not
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

#######
# GLDAS
#######

# Loop through temperature simulations and bilinearly interpolate onto 50km
# grid
mask = xr.open_dataset('modern_data/GLDASp5_landmask_025d.nc4')['GLDAS_mask'].squeeze()

# Add spherical earth CRS to model output as specified by GLDAS
# This is the Proj4 string for cartopy.crs.Geodetic()
mask = mask.rio.write_crs('+proj=longlat +datum=WGS84 +no_defs')

# Renaming to coordinate names rioxarray expects
if 'lat' not in mask.dims:
    mask = mask.rename({'latitude':'y', 'longitude':'x'})
else:
    mask = mask.rename({'lat':'y', 'lon':'x'})

# Resample onto target 50km grid using the nearest grid cell. Do not want to 
# bilinearly interpolate here because that would lead to values that range 
# between 0 and 1.
ena_regridded_mask = mask.rio.reproject_match(target_array, resampling=Resampling.nearest)
    
# Rename back to cf compliant coordinate names
ena_regridded_mask = ena_regridded_mask.rename({'x':'lon',
                                                  'y':'lat'})
# Save as netCDF4
ena_regridded_mask.to_netcdf(('modern_data/postprocessed/GLDAS_50km_albers.nc'))

########################
# WORLD CLIM TEMPERATURE
########################

# Create an empty xarray that will serve as the extra dimension for indexing
# the models into a single xarray file
wc_concat_xr = xr.DataArray(np.arange(1, 13, 1), dims=['month'], name='month')

# Empty list to hold the models
tave_month_list = list()

# Loop through and append each model onto the empty list
tave_files = glob.glob('modern_data/wc2.1_10m_tavg/*.tif')
tave_files.sort() # Want the files in month order

for file in tave_files:
    tave_month_list.append(rioxarray.open_rasterio(file))

# Use the empty xarray to reshape the list into a single xarray where the 
# 'month' dimension is associated with each monthly GeoTIFF read in
wc_tave = xr.concat(tave_month_list, dim=wc_concat_xr)

# Convert to Kelvin to maintain scale consistency with the climate models
wc_tave = wc_tave + 273.15

# Take annual average and get ride of the unnecessary dimension and convert 
wc_tave_annual = wc_tave.mean('month').squeeze()

# Resample onto target 50km grid using the nearest sample since we're
# effectively coarsening the data here
ena_regridded_wc_tave_annual = wc_tave_annual.rio.reproject_match(target_array, resampling=Resampling.nearest)

# For some reason the y-axis flips when saving as a GeoTIFF - fixing here
ena_regridded_wc_tave_annual.values = np.flip(ena_regridded_wc_tave_annual.values, 0)

# Save as geoTIFF
ena_regridded_wc_tave_annual.rio.to_raster(('modern_data/postprocessed/wc_tave_annual_50km_ena.tif'))

##########################
# WORLD CLIM PRECIPITATION
##########################

# We'll be converting to mm/day later on so if we set the xarray to have a
# months as a datetime64 object our life becomes much easier
time_dim = np.arange(np.datetime64('1970-01'),
                     np.datetime64('1971-01'),
                     np.timedelta64(1, 'M')) # the year doesn't matter we're only after the month

# Create an empty xarray that will serve as the extra dimension for indexing
# the models into a single xarray file
wc_concat_xr = xr.DataArray(time_dim, dims=['month'], name='month')

# Empty list to hold the models
pr_month_list = list()

# Loop through and append each model onto the empty list
pr_files = glob.glob('modern_data/wc2.1_10m_prec/*.tif')
pr_files.sort() # Want the files in month order

for file in pr_files:
    pr_month_list.append(rioxarray.open_rasterio(file))

# Use the empty xarray to reshape the list into a single xarray where the 
# 'month' dimension is associated with each monthly GeoTIFF read in
wc_pr = xr.concat(pr_month_list, dim=wc_concat_xr)

# Convert to mm/day to maintain consistency with the climate model output
wc_pr = wc_pr / wc_pr.month.dt.days_in_month

# Take annual average and get ride of the unnecessary dimension and convert 
wc_pr_annual = wc_pr.mean('month').squeeze()

# Resample onto target 50km grid using the nearest sample since we're
# effectively coarsening the data here
ena_regridded_wc_pr_annual = wc_pr_annual.rio.reproject_match(target_array, resampling=Resampling.nearest)

# For some reason the y-axis flips when saving as a GeoTIFF - fixing here
ena_regridded_wc_pr_annual.values = np.flip(ena_regridded_wc_pr_annual.values, 0)

# Save as geoTIFF
ena_regridded_wc_pr_annual.rio.to_raster(('modern_data/postprocessed/wc_pr_annual_50km_ena.tif'))
