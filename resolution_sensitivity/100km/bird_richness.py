"""
**********THESE ARE SPATIAL GRAIN SENSITIVITY TESTS ON A 100KM GRID***********

SCRIPT TAKES IN THE APMHIBIAN RANGE MAPS FROM BIRDLIFE, RASTERIZES THEM, SORTS
THEM INTO FAMILY BASED ON THE SHEARD ET AL., (2020, NATURE COMMUNICATIONS)
TRAIT DATABSED AND THEN SUMS RASTERS TO CREATE SPECIES RICHNESS ESTIMATES. ALL
OF THESE PROCEDURES ARE PERFORMED ON THE NORTH AMERICAN ALBERS EQUAL AREA 
PROJECTION AT A 50KM SPATIAL GRAIN.

THE BIRD RANGE MAPS ARE IN A GDB WHICH PYTHON HAS QUITE AN ISSUE WITH - IT EATS
MEMORY LIKE NO ONTHER. TO GET AROUND THIS, I READ IN 1000 LINES OF THE 17397
LINES OF THE GDB AT A TIME AND STORE THAT IN A LIST. ANALYSIS PROCEEDS AS
NORMAL AFTER THAT.

ALL STATISTICAL ANALYSES THAT FOLLOW ARE BASED ON THE ORDER PASSERIFORMES BUT
THE OTHER ORDERS ARE STILL RETAINED.

INPUT:

    * SPECIES_DATA/BOTW.GDB
    * SPECIES_DATA/SHEARD_DATASET_HWI_20200410_ORIG.XLSX

OUTPUT:

    * RESOLUTION_SENSITIVITY/100KM/RICHNESS/RICHNESS/BIRDS_PASSERIFORMES_RICHNESS_EPSG_102008.TIF

MINICONDA ENVIRONMENT:
    PYTHON: 3.8.10
    GEOPANDAS: 0.10.2
    NUMPY: 1.20.3
    XARRAY: 0.18.2
    RIOXARRAY: 0.5.0
    GEOCUBE: 0.1.2
    FUNCTOOLS: BUILT INTO PYTHON
    MATH: BUILT INTO PYTHON

@author: David Fastovich
"""

import geopandas as gpd
import pandas as pd
import numpy as np
import xarray as xr
import rioxarray
from geocube.api.core import make_geocube
from geocube.rasterize import rasterize_image
from functools import partial
import math

####################
# READ IN SHAPEFILES
####################

# List to shapefiles to be concatenoated
bird_pd_list = list()

# Read in trait database
traits = pd.read_excel('species_data/Sheard_Dataset_HWI_20200410_orig.xlsx')

for ii in np.arange(1, 17397, 1000):
    print(ii)
    one_bird = gpd.read_file('species_data/BOTW.gdb', ignore_fields=('yrcompiled', 'source', 'origin', 'data_sens', 'sens_comm', 'compiler', 'tax_comm', 'dist_comm', 'generalised', 'reviewers', 'citation', 'yrmodified', 'version'), rows=slice(ii,ii+1000))

    one_bird = one_bird.merge(traits, left_on='sci_name', right_on='IUCN name', how='left')
    
    # Only keeping extant ranges (not probably - see BirdLife metadata for description)
    one_bird = one_bird[one_bird.presence == 1]

    # Rename the id column since I'm reusing some code
    one_bird.rename(columns={'OBJECTID':'id_no'})
    
    # Transforming to albers equal area projection
    one_bird = one_bird.to_crs('+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs')
        
    bird_pd_list.append(one_bird)
    
# Create a single GeoDataFrame
birds = pd.concat(bird_pd_list)

# Drop nan life histories
birds = birds.dropna(subset=['Order'])

#############
# TARGET GRID
#############

# Projection: North America Albers Equal Area Conic
# EPSG: 102008
# https://epsg.io/102008
# Origin: 40 N, 96 W
# Standard parallels = 20 N, 60 N
# Datum: NAD83
target_crs = '+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs'

# Make grid
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

target_array = target_array.rio.write_crs(target_crs)

#####################################################
# RASTERIZE POLYGONS BUT IN A FOR LOOP FOR EACH ORDER
#####################################################

for order in birds['Order'].unique():
    
    # Subset the order we way
    one_order = birds[birds.Order == order]
    
    # Using geocube.api to rasterize which has a rasterio and GDAL backend
    # Performing and array split becausethe GDAL backend eats computer
    # memory which will crash python. Breaking it into chunks makes it manageable.
    # We want to feed in about a thousand rows at a time
    birds_array_split = np.array_split(one_order, math.ceil(one_order.shape[0]/1000))
    birds_grid_list = list()
    
    for chunk in birds_array_split:
        geo_grid = make_geocube(
            vector_data=chunk.loc[:, ['id_no','presence', 'geometry']],
            measurements=['presence'],
            like=target_array,
            group_by='id_no', # Tells geocube to keep the vector layers (individual species) separate
            rasterize_function=partial(rasterize_image, all_touched=True), # Passing custom argument through to assign all cells that touch the polygon a value of 1
        )
        birds_grid_list.append(geo_grid)
        
    #####################################################
    # SPECIES RICHNESS AND WEIGHTED ENDEMISM CALCULATIONS
    #####################################################
    
    # Quick functino to calcuate endemism
    def calc_weighted_endemism(x):
        weight = 1/np.sum(~np.isnan(x)) # Weight is 1/the number of cells the species exist in
        x[~np.isnan(x)] = weight # Replacing cell values where the species exist to weight vales
        return(x)
    
    # Lists to hold calculated values
    richness_list = list()
    weighted_endemism_list = list()
    
    # Each metric is calculated on one chunk and the chunks are added together
    # This is allowed because the only operation is addition and order does not
    # matter
    for chunk in birds_grid_list:
        
        # Sum richness for one chunk
        richness_list.append(chunk.presence.sum('id_no'))
        
        # Apply weighted_endemism function - this is applied onto each individual 
        # 'id_no' lat-lon array
        presence_array = chunk.presence.values
        for x in presence_array:
          weighted_endemism_list.append(calc_weighted_endemism(x))
    
    #####################################
    # RESHAPE ARRAYS AFTER CALCULATIONS #
    #####################################
    
    ############# RICHNESS ##############
    
    # xarray to concatenate the species richness lists - lat-lon array for each chunk
    concat_dim = xr.DataArray(np.arange(0,len(richness_list)), dims=['id'], name = 'id')
    
    # Concatenate species richness
    richness_concat = xr.concat(richness_list, concat_dim)
    richness = richness_concat.sum('id')
    
    # Save as GeoTIFF
    # For some reason the y-axis flips when saving as a GeoTIFF - fixing here
    richness.values = np.flip(richness.values, 0)
    richness.rio.to_raster('resolution_sensitivity/100km/richness/birds_' + order + '_richness_epsg_102008.tif')
    
