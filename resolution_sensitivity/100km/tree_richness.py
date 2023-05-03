"""
**********THESE ARE SPATIAL GRAIN SENSITIVITY TESTS ON A 100KM GRID***********

SCRIPT TAKES IN THE TREE RANGE MAPS FROM THE DIGITIZED VERSION OF THE ED LITTLE 
REANGE MAPS (ORIGINALLY PUBLISHED IN 1976 AND DIGITIZED BY THE USGS IN 199),
RASTERIZES THEM, AND THEN SUMS THEM TO CREATE SPECIES RICHNESS ESTIMATES. ALL
OF THESE PROCEDURES ARE PERFORMED ON THE NORTH AMERICAN ALBERS EQUAL AREA
PROJECTION AT A 50KM SPATIAL GRAIN.

INPUT:

    * SPECIES_DATA/E_Little_Jr_Atlas_of_United_States_Trees_Digital_Representation_1999/*
    * SPECIES_DATA/E_Little_Jr_Atlas_of_United_States_Trees_Digital_Representation_1999/Little_datatable.csv

OUTPUT:

    * RESOLUTION_SENSITIVITY/100KM/RICHNESS/TREE_RICHNESS_EPSG_102008.TIF

MINICONDA ENVIRONMENT:
    PYTHON: 3.8.10
    GEOPANDAS: 0.10.2
    NUMPY: 1.20.3
    XARRAY: 0.18.2
    RIOXARRAY: 0.5.0
    GEOCUBE: 0.1.2
    FUNCTOOLS: BUILT INTO PYTHON

@author: David Fastovich
"""

import geopandas as gpd
import pandas as pd
import numpy as np
import xarray as xr
from geocube.api.core import make_geocube
from functools import partial
from geocube.rasterize import rasterize_image

####################
# READ IN SHAPEFILES
####################

# Read in csv that contain shapefile names and subspecies to be combined
little_trees = pd.read_csv('species_data/E_Little_Jr_Atlas_of_United_States_Trees_Digital_Representation_1999/Little_datatable.csv')

# List to shapefiles to be concatenoated
tree_pd_list = list()

# Loop through and read in each shapefile
for index, row in little_trees.iterrows():
    one_tree = gpd.read_file('species_data/E_Little_Jr_Atlas_of_United_States_Trees_Digital_Representation_1999/' + row['SHP/*'] + '.shp')
    one_tree = one_tree[['AREA', 'PERIMETER', 'geometry']]
    one_tree = one_tree.set_crs('EPSG:4267')
    one_tree = one_tree.to_crs('+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs')
    one_tree['tree'] = row['SHP/*']
    one_tree['combine'] = row['combine']
    tree_pd_list.append(one_tree)
    
# Create a single GeoDataFrame
trees_gdf = pd.concat(tree_pd_list)

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

####################
# RASTERIZE POLYGONS
####################

# Add a 'presence' column to the individual and subspecies GeoDataFrames
trees_gdf['presence'] = 1

# Remove subspecies to be dealt with separately
subspecies_to_combine = trees_gdf[trees_gdf['combine'].notna()]
trees_gdf = trees_gdf[trees_gdf['combine'].isna()]

# Before rasterizing we need to add a unique id for each species in trees_gdf
# Not doing this for subspecies because those already have a 'combine' id that
# groups together subspecies
trees_gdf['id'] = trees_gdf['tree'].factorize()[0]+1

trees_raster = make_geocube(
    vector_data=trees_gdf,
    measurements=['presence'],
    like=target_array,
    group_by='id',
    rasterize_function=partial(rasterize_image, all_touched=True),
)

# Rasterize subspecies with the combine column
subspecies_raster =  make_geocube(
    vector_data=subspecies_to_combine,
    measurements=['presence'],
    like = target_array, # maintain identical grids
    group_by='combine',
    rasterize_function=partial(rasterize_image, all_touched=True),
)

# Changing unique id dimension name and continuing numbering from where
# tree_raster leaves off
subspecies_raster = subspecies_raster.rename({'combine':'id'})
subspecies_raster = subspecies_raster.assign_coords({'id':(subspecies_raster.id.values + trees_raster.id.values.max()).astype('int')})

# Concatonate the two xarrays into a single
trees_raster = xr.concat([trees_raster, subspecies_raster], dim = 'id')

##################
# SPECIES RICHNESS
##################

# Sum across id's to get richness
trees_richness = trees_raster.sum('id')

# Save as GeoTIFF
# For some reason the y-axis flips when saving as a GeoTIFF - fixing here
trees_richness.presence.values = np.flip(trees_richness.presence.values, 0)
trees_richness.rio.to_raster('resolution_sensitivity/100km/richness/tree_richness_epsg_102008.tif')