# from modules import nearest_neighbor
# from modules import model_skill
import pandas as pd
import numpy as np
import sys
sys.path.append('modules')
from nearest_neighbor import nearest_neighbor
from model_skill import model_skill

def parallel_skill(model, model_lat, model_lon, reps, ni, proxy, model_name):
    
    # Preparing the data for random sampling
    # Temporary list to hold grid_cell_ids
    id_list = list()
    
    # Loop through and find the grid cell where each proxy site belongs
    for ii, row in proxy.iterrows():
        bounding_grid_cell = nearest_neighbor(point_lat=row.lat_reproject, 
                                              point_lon=row.long_reproject,
                                              model_lat=model_lat,
                                              model_lon=model_lon,
                                              rectangle=False)
        id_list.append(np.array([bounding_grid_cell['index']['lat'], bounding_grid_cell['index']['lon']]))

    # Unique grid cells to be sampled from
    unique_grid_cells = np.unique(np.array(id_list), axis=0)
    
    # This long string returns the row index of unique_grid_cells that this proxy
    # site belongs to that will be used as a unique id that is a integer rathern than 
    # an array. Here's what the code does:
    #
    # 1. Compare the grid_cell_id of this row to all unqiue grid cells --> boolean array of unique_grid_cells shape
    # 2. Apply the all function along this new boolean array --> boolean array of unique_grid_cells nrow shape
    # 3. Find the index of the one True from (2) --> index as int
    # 4. Create lambda function from 1-3
    # 5. Apply lamdba function along rows of proxy data frame and assign results to proxy
    #
    # Forgive me future David, I didnt want to come up with 5 new variable names
    
    # Attach single number id for each grid cell to proxy data frame because
    # pandas does not like numpy arrays for its various functions
    proxy['grid_cell_id'] = id_list
    proxy['unique_grid_cells_row'] = proxy.apply(lambda x: np.where(np.apply_along_axis(np.all, 1, x.iloc[10] == unique_grid_cells) == True)[0].item(), axis=1)  
    
    # Now that we have a way to group sites based on the cells they belong to, 
    # the resampling procedure is looped within the while loop below.
    
    skill_ens = list()  # Ensemble list
    count = 0 # Resetting count with each loop
    while reps > count:
        count += 1
        # print(count)
        # Model skill is spatially biased, so given this we are taking a
        # resampling approach. Resample one site from each grid cell.
        
        random_sample = pd.DataFrame().reindex_like(proxy.iloc[0:unique_grid_cells.shape[0],0:12])
        for key in range(0, unique_grid_cells.shape[0]):
            random_sample.iloc[key] = proxy[proxy.unique_grid_cells_row == key].sample(n=1)[0:12]

        # Calculate skill for a specific season and model
        val = model_skill(mi=model,
                          mi_lat=model_lat,
                          mi_lon=model_lon,
                          oi=random_sample,
                          ni=ni)  # temporary value
        skill_ens.append(val)
    # Precipitation throws a lot of nans given the large error
    # this removes them.
    skill_ens = np.array(skill_ens)
    skill_ens = skill_ens[~np.isnan(skill_ens)]
    return({'skill_ens':skill_ens.tolist(), 'model_name':model_name})

if __name__ == '__main__':
    
    import xarray as xr
    import pandas as pd
    import numpy as np
    import sys
    sys.path.append('modules')
    from nearest_neighbor import nearest_neighbor
    from model_skill import model_skill
    import geopandas
    
    # Model
    model = xr.open_dataset('model_temperature_anomaly_50km_ena/TRACE_LORENZ_tas_anom_50km_ena.nc')
    
    proxy = pd.read_csv('proxy_data/tave_anomalies.csv')
    gdgt = pd.read_csv('proxy_data/gdgt_mat.csv')
    proxy = pd.concat([proxy, gdgt])
    
    # Drop NAs
    proxy = proxy.dropna(axis=0)
    
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
    
    # Drop geometry column
    proxy_gdf.drop('geometry', axis=1, inplace=True)
        
    # Temporary list to hold grid_cell_ids
    id_list = list()
    
    # Loop through and find the grid cell where each proxy site belongs
    for ii, row in proxy_gdf.iterrows():
        bounding_grid_cell = nearest_neighbor(point_lat=row.lat_reproject, 
                                              point_lon=row.long_reproject,
                                              model_lat=model.lat.values,
                                              model_lon=model.lon.values,
                                              rectangle=False)
        id_list.append(np.array([bounding_grid_cell['index']['lat'], bounding_grid_cell['index']['lon']]))

    # Unique grid cells to be sampled from
    unique_grid_cells = np.unique(np.array(id_list), axis=0)
    
    # This long string returns the row index of unique_grid_cells that this proxy
    # site belongs to that will be used as a unique id that is a integer rathern than 
    # an array. Here's what the code does:
    #
    # 1. Compare the grid_cell_id of this row to all unqiue grid cells --> boolean array of unique_grid_cells shape
    # 2. Apply the all function along this new boolean array --> boolean array of unique_grid_cells nrow shape
    # 3. Find the index of the one True from (2) --> index as int
    # 4. Create lambda function from 1-3
    # 5. Apply lamdba function along rows of proxy data frame and assign results to proxy
    #
    # Forgive me future David, I didnt want to come up with 5 new variable names
    
    # Attach single number id for each grid cell to proxy data frame because
    # pandas does not like numpy arrays for its various functions
    proxy_gdf['grid_cell_id'] = id_list
    proxy_gdf['unique_grid_cells_row'] = proxy_gdf.apply(lambda x: np.where(np.apply_along_axis(np.all, 1, x.iloc[10] == unique_grid_cells) == True)[0].item(), axis=1)  
    
    # Now that we have a way to group sites based on the cells they belong to, 
    # the resampling procedure is looped here:
    proxy_gdf[proxy_gdf.unique_grid_cells_row == 0].sample(n=1)

    random_sample = pd.DataFrame().reindex_like(proxy_gdf.iloc[0:unique_grid_cells.shape[0],0:12])
    for key in range(0, unique_grid_cells.shape[0]):
        random_sample.iloc[key] = proxy_gdf[proxy_gdf.unique_grid_cells_row == key].sample(n=1)[0:12]
    random_sample

    # Test function
    parallel_skill(model=model['tas'].values,
                   model_lat=model.lat.values,
                   model_lon=model.lon.values,
                   reps=10,
                   ni=0,
                   proxy=proxy_gdf,
                   model_name='x')
