"""
PLOT RICHNESS ESTIMATES FOR N0RTH AMERICA WITH A BOUNDING BOX OVER EASTERN NORTH
AMERICA (STUDY REGION). CORRESPONDS TO FIGURE S1 IN THE MANUSCRIPT.

INPUT:

  * RICHNESS/AMPH_RICHNESS_EPSG_102008.TIF
  * RICHNESS/BIRDS_PASSERIFORMES_RICHNESS_EPSG_102008.TIF
  * RICHNESS/MAMMS_RICHNESS_EPSG_102008.TIF
  * RICHNESS/REPS_RICHNESS_EPSG_102008.TIF
  * RICHNESS/TREE_RICHNESS_EPSG_102008.TIF
      
OUTPUT:
  * FIGURES/RICHNESS_MAPS.PDF

MINICONDA ENVIRONMENT:
  PYTHON: 3.8.10
  XARRAY: 0.18.2
  RIOXARRAY: 0.5.0
  NUMPY: 1.20.3
  MATPLOTLIB: 3.5.1
  CARTOPY: 0.20.2
  SHAPELY: 1.8.0

@author: David Fastovich
"""

import platform
import xarray as xr
import rioxarray
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker as mticker
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from shapely import geometry

###############################
# READ IN SPECIES RICHNESS DATA
###############################

# Holds the xarray instances of species richness
xarray_list = list()

# Xarray to concatenate xarray_list along species
biodiv_concat = xr.DataArray(['Amphibians', 'Birds', 'Mammals', 'Reptiles', 'Trees'], dims=['organism'], name='organism')

# Richness files and sort
biodiv_files = [
    'richness/amph_richness_epsg_102008.tif',
    'richness/birds_Passeriformes_richness_epsg_102008.tif',
    'richness/mamms_richness_epsg_102008.tif',
    'richness/reps_richness_epsg_102008.tif',
    'richness/tree_richness_epsg_102008.tif'
 ]
biodiv_files.sort()

# Loop through richness files and combine them into a single xarray
for file in biodiv_files:
    # Read in richness data
    rio_xr = rioxarray.open_rasterio(file).squeeze()
    rio_xr = rio_xr.assign_coords({'x':np.flip(rio_xr.x.values)})
    rio_xr.values = np.flip(rio_xr.values)
    rio_xr_crop = rio_xr.sel(x=slice(2700000, -2500000)).sel(y=slice(-1650000, 1600000))
    xarray_list.append(rio_xr_crop)
    
# Concatenate xarray_list to get a single multidimensional xarray
species_richness = xr.concat(xarray_list, biodiv_concat)

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

# Projection - North American Albers Equal Area projection
# https://epsg.io/102008
# Proj.4 - +proj=aea +lat_0=40 +lon_0=-96 +lat_1=20 +lat_2=60 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +type=crs
proj = ccrs.AlbersEqualArea(central_longitude=-96, central_latitude=40, standard_parallels=(20.0, 60.0))

# Box of study region
geom = geometry.box(minx=-90.6849314, maxx=-55.3347072, miny=24.4006828, maxy=49.6764456)

# GridSpec configuration
fig = plt.figure(figsize=(6.75, 24), constrained_layout=True)
gs = gridspec.GridSpec(nrows=10, ncols=1, figure=fig, height_ratios=(1, 0.1, 1, 0.1, 1, 0.1, 1, 0.1, 1, 0.1))

# AX0 - Amphibians ------------------------------------------------------------

ax0 = fig.add_subplot(gs[0,0], projection=proj)
cb = ax0.pcolormesh(species_richness.x,
                    species_richness.y,
                    species_richness.sel(organism='Amphibians').values,
                    cmap=plt.cm.get_cmap('plasma'),
                    transform=proj,
                    rasterized=True)
ax0.set_title('Amphibians')
ax0.set_extent([-120, -70, 26, 50], crs=ccrs.PlateCarree())
ax0.add_feature(cfeature.COASTLINE, zorder=2, edgecolor='#000000')
ax0.add_feature(cfeature.LAKES, facecolor='#000000')
ax0.add_feature(cfeature.OCEAN, facecolor='#000000', zorder=1)
ax0.add_geometries([geom], crs=ccrs.PlateCarree(), facecolor='none', edgecolor='red', lw=3)
plt.text(-2320000, 1500000, 'A', weight='bold')

gl = ax0.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False)
gl.top_labels = False
gl.right_labels = False
gl.xlabel_style = {'rotation': 0, 'verticalalignment':'bottom'}
gl.ylabel_style = {'rotation': 0}
gl.xlocator = mticker.FixedLocator([-120, -110, -100, -90, -80, -70])
gl.ylocator = mticker.FixedLocator([30, 40, 50])

ticks = np.arange(species_richness.sel(organism='Amphibians').values.min(), \
                  species_richness.sel(organism='Amphibians').values.max(), \
                  5)
              
cb_ax = fig.add_subplot(gs[1,0])
plt.colorbar(cb, orientation='horizontal', cax=cb_ax, ticks=ticks)

# Getting rid of unneeded labels:
# https://stackoverflow.com/questions/64019387/specify-the-lat-lon-label-location-in-cartopy-remove-at-some-sides
plt.draw()  #enable the use of ._lables()
for ea in gl._labels:
    if platform.platform() in ['macOS-10.16-x86_64-i386-64bit', 'Linux-3.10.0-1160.15.2.el7.x86_64-x86_64-with-glibc2.17']:
        if '120°W'==ea.artist.get_text():
            ea.artist.set_text("")
    else:
        if '120°W'==ea[2].get_text():
            ea[2].set_text("")
        if '70°W'==ea[2].get_text():
            ea[2].set_text("")
# AX1 - Birds ----------------------------------------------------------------

ax1 = fig.add_subplot(gs[2,0], projection=proj)
cb = ax1.pcolormesh(species_richness.x,
                    species_richness.y,
                    species_richness.sel(organism='Birds').values,
                    cmap=plt.cm.get_cmap('plasma'),
                    transform=proj,
                    rasterized=True)

ax1.set_title('Birds')
ax1.set_extent([-120, -70, 26, 50], crs=ccrs.PlateCarree())
ax1.add_feature(cfeature.COASTLINE, zorder=2, edgecolor='#000000')
ax1.add_feature(cfeature.LAKES, facecolor='#000000')
ax1.add_feature(cfeature.OCEAN, facecolor='#000000', zorder=1)
ax1.add_geometries([geom], crs=ccrs.PlateCarree(), facecolor='none', edgecolor='red', lw=3)
plt.text(-2320000, 1500000, 'B', weight='bold')

gl = ax1.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False)
gl.top_labels = False
gl.xlabel_style = {'rotation': 0}
gl.ylabel_style = {'rotation': 0}
gl.xlocator = mticker.FixedLocator([-120, -110, -100, -90, -80, -70])
gl.ylocator = mticker.FixedLocator([30, 40, 50])

ticks = np.arange(0,                                                      \
                  species_richness.sel(organism='Birds').values.max(),    \
                  50)

cb_ax = fig.add_subplot(gs[3,0])
cbar = plt.colorbar(cb, orientation='horizontal', cax=cb_ax, ticks=ticks)

plt.draw()  #enable the use of ._lables()
for ea in gl._labels:
    if platform.platform() in ['macOS-10.16-x86_64-i386-64bit', 'Linux-3.10.0-1160.15.2.el7.x86_64-x86_64-with-glibc2.17']:
        if '120°W'==ea.artist.get_text():
            ea.artist.set_text("")
    else:
        if '120°W'==ea[2].get_text():
            ea[2].set_text("")
        if '70°W'==ea[2].get_text():
            ea[2].set_text("")

# AX2 - Mammals --------------------------------------------------------------

ax2 = fig.add_subplot(gs[4,0], projection=proj)
cb = ax2.pcolormesh(species_richness.x,
                    species_richness.y,
                    species_richness.sel(organism='Mammals').values,
                    cmap=plt.cm.get_cmap('plasma'),
                    transform=proj,
                    rasterized=True)
ax2.set_title('Mammals')
ax2.set_extent([-120, -70, 26, 50], crs=ccrs.PlateCarree())
ax2.add_feature(cfeature.COASTLINE, zorder=2, edgecolor='#000000')
ax2.add_feature(cfeature.LAKES, facecolor='#000000')
ax2.add_feature(cfeature.OCEAN, facecolor='#000000', zorder=1)
ax2.add_geometries([geom], crs=ccrs.PlateCarree(), facecolor='none', edgecolor='red', lw=3)
plt.text(-2320000, 1500000, 'C', weight='bold')

gl = ax2.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False)
gl.top_labels = False
gl.xlabel_style = {'rotation': 0}
gl.ylabel_style = {'rotation': 0}
gl.xlocator = mticker.FixedLocator([-120, -110, -100, -90, -80, -70])
gl.ylocator = mticker.FixedLocator([30, 40, 50])

ticks = np.arange(0,                                                        \
                  species_richness.sel(organism='Mammals').values.max()+10, \
                  10)

cb_ax = fig.add_subplot(gs[5,0])
plt.colorbar(cb, orientation='horizontal', cax=cb_ax, ticks=ticks)

plt.draw()  #enable the use of ._lables()
for ea in gl._labels:
    if platform.platform() in ['macOS-10.16-x86_64-i386-64bit', 'Linux-3.10.0-1160.15.2.el7.x86_64-x86_64-with-glibc2.17']:
        if '120°W'==ea.artist.get_text():
            ea.artist.set_text("")
    else:
        if '120°W'==ea[2].get_text():
            ea[2].set_text("")
        if '70°W'==ea[2].get_text():
            ea[2].set_text("")
        
# AX3 - Reptiles --------------------------------------------------------------

ax3 = fig.add_subplot(gs[6,0], projection=proj)
cb = ax3.pcolormesh(species_richness.x,
                    species_richness.y,
                    species_richness.sel(organism='Reptiles').values,
                    cmap=plt.cm.get_cmap('plasma'),
                    transform=proj,
                    rasterized=True)
ax3.set_title('Reptiles')
ax3.set_extent([-120, -70, 26, 50], crs=ccrs.PlateCarree())
ax3.add_feature(cfeature.COASTLINE, zorder=2, edgecolor='#000000')
ax3.add_feature(cfeature.LAKES, facecolor='#000000')
ax3.add_feature(cfeature.OCEAN, facecolor='#000000', zorder=1)
ax3.add_geometries([geom], crs=ccrs.PlateCarree(), facecolor='none', edgecolor='red', lw=3)
plt.text(-2320000, 1500000, 'D', weight='bold')

ticks = np.arange(0,                                                         \
                  species_richness.sel(organism='Reptiles').values.max()+10, \
                  10)

gl = ax3.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False)
gl.top_labels = False
gl.xlabel_style = {'rotation': 0}
gl.ylabel_style = {'rotation': 0}
gl.xlocator = mticker.FixedLocator([-120, -110, -100, -90, -80, -70])
gl.ylocator = mticker.FixedLocator([30, 40, 50])

cb_ax = fig.add_subplot(gs[7,0])
plt.colorbar(cb, orientation='horizontal', cax=cb_ax, ticks=ticks)

plt.draw()  #enable the use of ._lables()
for ea in gl._labels:
    if platform.platform() in ['macOS-10.16-x86_64-i386-64bit', 'Linux-3.10.0-1160.15.2.el7.x86_64-x86_64-with-glibc2.17']:
        if '120°W'==ea.artist.get_text():
            ea.artist.set_text("")
    else:
        if '120°W'==ea[2].get_text():
            ea[2].set_text("")
        if '70°W'==ea[2].get_text():
            ea[2].set_text("")
        
# AX4 - Mammals --------------------------------------------------------------

ax4 = fig.add_subplot(gs[8,0], projection=proj)
cb = ax4.pcolormesh(species_richness.x,
                    species_richness.y,
                    species_richness.sel(organism='Trees').values,
                    cmap=plt.cm.get_cmap('plasma'),
                    transform=proj,
                    rasterized=True)
ax4.set_title('Trees')
ax4.set_extent([-120, -70, 26, 50], crs=ccrs.PlateCarree())
ax4.add_feature(cfeature.COASTLINE, zorder=2, edgecolor='#000000')
ax4.add_feature(cfeature.LAKES, facecolor='#000000')
ax4.add_feature(cfeature.OCEAN, facecolor='#000000', zorder=1)
ax4.add_geometries([geom], crs=ccrs.PlateCarree(), facecolor='none', edgecolor='red', lw=3)
plt.text(-2320000, 1500000, 'E', weight='bold')

ticks = np.arange(0,                                                      \
                  species_richness.sel(organism='Trees').values.max()+20, \
                  20)

gl = ax4.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False)
gl.top_labels = False
gl.xlabel_style = {'rotation': 0}
gl.ylabel_style = {'rotation': 0}
gl.xlocator = mticker.FixedLocator([-120, -110, -100, -90, -80, -70])
gl.ylocator = mticker.FixedLocator([30, 40, 50])

cb_ax = fig.add_subplot(gs[9,0])
cbar = plt.colorbar(cb, orientation='horizontal', cax=cb_ax, ticks=ticks)
cbar.set_label('Number of Species')

plt.draw()  #enable the use of ._lables()
for ea in gl._labels:
    if platform.platform() in ['macOS-10.16-x86_64-i386-64bit', 'Linux-3.10.0-1160.15.2.el7.x86_64-x86_64-with-glibc2.17']:
        if '120°W'==ea.artist.get_text():
            ea.artist.set_text("")
    else:
        if '120°W'==ea[2].get_text():
            ea[2].set_text("")
        if '70°W'==ea[2].get_text():
            ea[2].set_text("")

plt.savefig('figures/richness_maps.png', bbox_inches='tight', dpi=300)
plt.savefig('figures/richness_maps.pdf', bbox_inches='tight', dpi=300)
