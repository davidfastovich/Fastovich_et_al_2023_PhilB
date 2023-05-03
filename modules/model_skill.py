import numpy as np
# from modules import bilinear_interpolation

import sys
sys.path.append('modules')
from bilinear_interpolation import bilinear_interpolation

# def model_skill(mi, mi_lat, mi_lon, oi, ni=0):

#     # Formula from Hargreaves, J.C., Annan, J.D., Ohgaito, R., Paul,
#     # A., Abe-Ouchi, A., 2013. Skill and reliability of climate model
#     # ensembles at the Last Glacial Maximum and midHolocene.
#     # Clim. Past 9, 811e823.

#     # Model Skill is defined as:
#     # SS = 1 - \sqrt{\frac{\Sigma(m_{i} - o_{i})^{2} - \Sigma(e_{i}^{2})}{\Sigma(n_{i} - o_{i})^{2} - \Sigma(e_{i}^{2})}}
#     #
#     # where
#     #
#     # m_i = xarray of simulated results (model output).
#     # n_i = xarray of references results (model output). Set as 0 here assuming there is no change between climate states
#     # o_i = Pandas dataframe of observations (proxy record) containing
#     #       list of observation uncertainty (1 standard deviation as defined in Hargreaves et al., 2011)

#     oi = oi.dropna(axis=0, how='any') # Drop NaN values
#     proxy_sim_diff = [] # list to hold m_i - o_i
#     proxy_ref_diff = [] # list to hold n_i - o_i

#     # Loop through oi and take differences
#     for ii, row in oi.iterrows():

#         # Switching to numpy calculations for speed
#         lat_diff = (mi_lat - row['lat']) ** 2
#         lat_loc = np.where(lat_diff == np.min(lat_diff))

#         lon_diff = (mi_lon - row['long']) ** 2
#         lon_loc = np.where(lon_diff == np.min(lon_diff))
#         sim_value = mi[lat_loc, lon_loc]

#         # sim_value = mi.sel(lon=row['long'], method='nearest').sel(lat=row['lat'], method='nearest').values # Selecting nearest cell in the simulation
#         proxy_sim_diff.append((sim_value - row['mean'])**2)
#         proxy_ref_diff.append((ni - row['mean'])**2)

#     # Calculate skill
#     skill = 1 - ((sum(proxy_sim_diff) - sum(oi['sd']**2))/(sum(proxy_ref_diff) - sum(oi['sd']**2)))**(1/2)

#     # Skill = 1 -> perfect model performance
#     # Skill = 0 -> performance no better than reference state of no change
#     # Skill = -1 -> performance worse than reference state of no change

#     return(skill)

def model_skill(mi, mi_lat, mi_lon, oi, ni=0):

    # Formula from Hargreaves, J.C., Annan, J.D., Ohgaito, R., Paul,
    # A., Abe-Ouchi, A., 2013. Skill and reliability of climate model
    # ensembles at the Last Glacial Maximum and midHolocene.
    # Clim. Past 9, 811e823.

    # Model Skill is defined as:
    # SS = 1 - \sqrt{\frac{\Sigma(m_{i} - o_{i})^{2} - \Sigma(e_{i}^{2})}{\Sigma(n_{i} - o_{i})^{2} - \Sigma(e_{i}^{2})}}
    #
    # where
    #
    # m_i = xarray of simulated results (model output).
    # n_i = xarray of references results (model output). Set as 0 here assuming there is no change between climate states
    # o_i = Pandas dataframe of observations (proxy record) containing
    #       list of observation uncertainty (1 standard deviation as defined in Hargreaves et al., 2011)

    oi = oi.dropna(axis=0, how='any') # Drop NaN values
    proxy_sim_diff = [] # list to hold m_i - o_i
    proxy_ref_diff = [] # list to hold n_i - o_i

    # Loop through oi and take differences
    for ii, row in oi.iterrows():

        # Switching to numpy calculations for speed
        sim_value = bilinear_interpolation(new_lat=row['lat_reproject'], new_lon=row['long_reproject'], model_lat=mi_lat, model_lon=mi_lon, climate_var=mi.transpose())

        # sim_value = mi.sel(lon=row['long'], method='nearest').sel(lat=row['lat'], method='nearest').values # Selecting nearest cell in the simulation
        proxy_sim_diff.append((sim_value - row['mean'])**2)
        proxy_ref_diff.append((ni - row['mean'])**2)

    # Calculate skill
    skill = 1 - ((sum(proxy_sim_diff) - sum(oi['sd']**2))/(sum(proxy_ref_diff) - sum(oi['sd']**2)))**(1/2)

    # Skill = 1 -> perfect model performance
    # Skill = 0 -> performance no better than reference state of no change
    # Skill = -1 -> performance worse than reference state of no change

    return(skill)

if __name__ == '__main__':
    
    import xarray as xr
    import pandas as pd
    from modules import nearest_neighbor
    
    # Model
    model = xr.open_dataset('/media/david/DAV/models/trace_mwf_seasonal/b30.00_4kaDVTn.cam2.ncrcat.djf.nc')
    
    # Proxy
    proxy = pd.read_csv('pollen_reconstruction/tave/tave_anomalies.csv')
    gdgt = pd.read_csv('final_figures/data/gdgt_mat.csv')
    proxy = pd.concat([proxy, gdgt])
    
    # Changing longitude to [0, 360] since that is what the models use
    proxy['long'] = proxy['long'] + 360
    
    # Drop NAs
    proxy = proxy.dropna(axis=0)
    
    # Temporary list to hold grid_cell_ids
    id_list = list()
    
    # Loop through and find the grid cell where each proxy site belongs
    for ii, row in proxy.iterrows():
        bounding_grid_cell = nearest_neighbor(point_lat=row.lat, 
                                              point_lon=row.long,
                                              model_lat=model.lat.values,
                                              model_lon=model.lon.values)
        grid_cell_id = np.array([list(map(lambda x:x['index']['lat'], bounding_grid_cell)), 
                                list(map(lambda x:x['index']['lon'], bounding_grid_cell))]).ravel()
        grid_cell_id.sort()
        id_list.append(grid_cell_id.tolist())

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
    proxy['unique_grid_cells_row'] = proxy.apply(lambda x: np.where(np.apply_along_axis(np.all, 1, x.iloc[8] == unique_grid_cells) == True)[0].item(), axis=1)  
    
    # Now that we have a way to group sites based on the cells they belong to, 
    # the resampling procedure is looped here:
    proxy[proxy.unique_grid_cells_row == 0].sample(n=1)

    random_sample = pd.DataFrame().reindex_like(proxy.iloc[0:unique_grid_cells.shape[0],0:8])
    for key in range(0, unique_grid_cells.shape[0]):
        random_sample.iloc[key] = proxy[proxy.unique_grid_cells_row == key].sample(n=1).squeeze()[0:8]
    random_sample
    
    # Pass randome sample into model skill
    model_skill(mi=model['TS'].isel(time=0).values, mi_lat=model.lat.values, mi_lon=model.lon.values, oi=random_sample, ni=0)
