import numpy as np
from scipy import spatial

def nearest_neighbor(point_lat, point_lon, model_lat, model_lon, rectangle=True):
    '''
    This function algorithmically returns the grid cell points bounding a given
    spatial point. See source code for the conditional statement that makes this
    work. If rectangle is set to false the fuction only returns the nearest
    neighbor. Calculations are performed using scipy.spatial.KDTree.
    
    Nearest neighbor identification is modeled after:
        https://kbkb-wx-python.blogspot.com/2016/08/find-nearest-latitude-and-longitude.html
        https://github.com/blaylockbk/pyBKB_v3/blob/master/demo/KDTree_nearest_neighbor.ipynb
    
    Parameters
    ----------
    point_lat : int or float
        Latitude value of point data.
    point_lon : int or float
        Longitude value of point data.
    model_lat : array
        Array of latitude values the model.
    model_lon: array
        Array of latitude values the model.

    Returns
    -------
    An array of the four nearest neighbor points and their indices, or the
    nearest neighbor.
    '''
    
    # Create a grid corresponding to the model coordinates
    grid_lon, grid_lat = np.meshgrid(model_lon, model_lat)
    
    # KDTree requires a (nlon, nlat) array where each row represents all 
    # possible combinatinos of longitude to latitude
    tree = spatial.KDTree(np.column_stack([grid_lon.ravel(), grid_lat.ravel()]))
    
    # KDTree.query returns the distance and index of the nearest neighbor relative
    # to the input data
    dist_to_neighbor, idx_neighbor = tree.query([point_lon, point_lat])
    
    # Translate the KDTree.query index to the model lat/lon indices
    nearest_lon_index = np.where(tree.data[idx_neighbor][0] == model_lon)[0].item()
    nearest_lat_index = np.where(tree.data[idx_neighbor][1] == model_lat)[0].item()
    
    # Get nearest neighbor as the bases for the conditional statement
    nearest_lon = tree.data[idx_neighbor][0]
    nearest_lat = tree.data[idx_neighbor][1]
    
    
    # Condtional statement to return the four bounding grid points
    if nearest_lat >= point_lat and nearest_lon <= point_lon:
        
        # Neighbor 1 = same latitude one line of longitude to the east
        neighbor_1 = {'coords': {'lat': model_lat[nearest_lat_index], 'lon':model_lon[nearest_lon_index+1]},
                      'index': {'lat': nearest_lat_index, 'lon': nearest_lon_index+1}}
        
        # Neighbor 2 = one line of latitude of latitude south, same longitude
        neighbor_2 = {'coords': {'lat': model_lat[nearest_lat_index-1], 'lon': model_lon[nearest_lon_index]},
                      'index': {'lon': nearest_lon_index, 'lat': nearest_lat_index-1}}
        
        # Neighbor 3 = one line of latitude of latitude south, one line of longitude to the east
        neighbor_3 = {'coords': {'lat': model_lat[nearest_lat_index-1], 'lon': model_lon[nearest_lon_index+1]},
                      'index': {'lon': nearest_lon_index+1, 'lat': nearest_lat_index-1}}
        
        # Four nearest neighbors, including the point previously calculated
        nearest_neighbors = [{'coords':{'lat': nearest_lat, 'lon': nearest_lon}, 'index':{'lat': nearest_lat_index, 'lon': nearest_lon_index}}, 
                             neighbor_1,
                             neighbor_2, 
                             neighbor_3]
    
    elif nearest_lat >= point_lat and nearest_lon >= point_lon:
        
        # Neighbor 1 = same latitude one line of longitude to the west
        neighbor_1 = {'coords': {'lat': model_lat[nearest_lat_index], 'lon': model_lon[nearest_lon_index-1]},
                      'index': {'lon': nearest_lon_index-1, 'lat': nearest_lat_index}}
        
        # Neighbor 2 = one line of latitude of latitude south, same longitude
        neighbor_2 = {'coords': {'lat': model_lat[nearest_lat_index-1], 'lon': model_lon[nearest_lon_index]},
                      'index': {'lon': nearest_lon_index, 'lat': nearest_lat_index-1}}
        
        # Neighbor 3 = one line of latitude of latitude south, one line of longitude to the west
        neighbor_3 = {'coords': {'lat': model_lat[nearest_lat_index-1], 'lon': model_lon[nearest_lon_index-1]},
                      'index': {'lon': nearest_lon_index-1, 'lat': nearest_lat_index-1}}
        
        # Four nearest neighbors, including the point previously calculated
        nearest_neighbors = [{'coords':{'lat': nearest_lat, 'lon': nearest_lon}, 'index':{'lat': nearest_lat_index, 'lon': nearest_lon_index}}, 
                             neighbor_1,
                             neighbor_2, 
                             neighbor_3]
        
    elif nearest_lat <= point_lat and nearest_lon >= point_lon:
    
        # Neighbor 1 = same latitude one line of longitude to the west
        neighbor_1 = {'coords': {'lat': model_lat[nearest_lat_index], 'lon': model_lon[nearest_lon_index-1]},
                      'index': {'lon': nearest_lon_index-1, 'lat': nearest_lat_index}}
        
        # Neighbor 2 = one line of latitude of latitude north, same longitude
        neighbor_2 = {'coords': {'lat': model_lat[nearest_lat_index+1], 'lon': model_lon[nearest_lon_index]},
                      'index': {'lon': nearest_lon_index, 'lat': nearest_lat_index+1}}
        
        # Neighbor 3 = one line of latitude of latitude north, one line of longitude to the west
        neighbor_3 = {'coords': {'lat': model_lat[nearest_lat_index+1], 'lon': model_lon[nearest_lon_index-1]},
                      'index': {'lon': nearest_lon_index-1, 'lat': nearest_lat_index+1}}
        
        # Four nearest neighbors, including the point previously calculated
        nearest_neighbors = [{'coords':{'lat': nearest_lat, 'lon': nearest_lon}, 'index':{'lat': nearest_lat_index, 'lon': nearest_lon_index}}, 
                             neighbor_1,
                             neighbor_2, 
                             neighbor_3]
        
    elif nearest_lat <= point_lat and nearest_lon <= point_lon:
        
        # Neighbor 1 = same latitude one line of longitude to the east
        neighbor_1 = {'coords': {'lat': model_lat[nearest_lat_index], 'lon': model_lon[nearest_lon_index+1]},
                      'index': {'lon': nearest_lon_index+1, 'lat': nearest_lat_index}}
        
        # Neighbor 2 = one line of latitude of latitude north, same longitude
        neighbor_2 = {'coords': {'lat': model_lat[nearest_lat_index+1], 'lon': model_lon[nearest_lon_index]},
                      'index': {'lon': nearest_lon_index, 'lat': nearest_lat_index+1}}
        
        # Neighbor 3 = one line of latitude of latitude north, one line of longitude to the east
        neighbor_3 = {'coords': {'lat': model_lat[nearest_lat_index+1], 'lon': model_lon[nearest_lon_index+1]},
                      'index': {'lon': nearest_lon_index+1, 'lat': nearest_lat_index+1}}
        
        # Four nearest neighbors, including the point previously calculated
        nearest_neighbors = [{'coords':{'lat': nearest_lat, 'lon': nearest_lon}, 'index':{'lat': nearest_lat_index, 'lon': nearest_lon_index}}, 
                             neighbor_1,
                             neighbor_2, 
                             neighbor_3]
    
    if rectangle == True:
        return(nearest_neighbors)
    else:
        return({'coords':{'lat':nearest_lat, 'lon':nearest_lon}, 'index':{'lat':nearest_lat_index, 'lon':nearest_lon_index}})


if __name__ == '__main__':
    
    import numpy.random as random
    
    def rectangle_test(x, y, points):
        '''
        Tests if input points form a rectange and if test point lies within
        the rectangle.
    
        The four points are a list of four couplets:  (x, y).
        The four points can be in any order.  They should form a rectangle.
        
        Source for test is: https://stackoverflow.com/a/8662355
        
        Parameters
        ----------
        x : int or float
            Test point to check if its within the rectangle.
        y : int or float
            Test point to check if its within the rectangle.
        points : array
            Points defining the rectangle.

        Raises
        ------
        ValueError
            points do not form a rectangle.
            
        ValueError
            (x, y) not within the rectangle

        Returns
        -------
        None.
        
        Example
        -------
        >>> rectangle_test(12, 5.5,[(10, 4),
        ...                         (20, 4),
        ...                         (10, 6),
        ...                         (20, 6])
        '''
    
        points = sorted(points)               # order points by x, then by y
        (x1, y1), (_x1, y2), (x2, _y1), (_x2, _y2) = points
    
        if x1 != _x1 or x2 != _x2 or y1 != _y1 or y2 != _y2:
            raise ValueError('points do not form a rectangle')
        if not x1 <= x <= x2 or not y1 <= y <= y2:
            raise ValueError('(x, y) not within the rectangle')
    
        return(True)
    
    # Test nearest neighbor function with while loop
    count = 0
    results = []
    while count < 10000:
        # Simulate point data
        point_lat = random.uniform(low=-90, high=90, size=1)
        point_lon = random.uniform(low=0, high=360, size=1)
        
        # Simulate climate model coordinates
        model_lat = np.linspace(start=-90, stop=90, num=10)
        model_lon = np.linspace(start=0, stop=360, num=20)
        
        # Make grid for plotting purposes
        grid_lon, grid_lat = np.meshgrid(model_lon, model_lat)
        
        # Call function
        sim_neighbors = nearest_neighbor(point_lat=point_lat, 
                                         point_lon=point_lon, 
                                         model_lat=model_lat, 
                                         model_lon=model_lon)
        # Test function results
        results.append(rectangle_test(x=point_lon, y=point_lat, points=[(sim_neighbors[0]['coords']['lon'], sim_neighbors[0]['coords']['lat']),
                                                                        (sim_neighbors[1]['coords']['lon'], sim_neighbors[1]['coords']['lat']),
                                                                        (sim_neighbors[2]['coords']['lon'], sim_neighbors[2]['coords']['lat']),
                                                                        (sim_neighbors[3]['coords']['lon'], sim_neighbors[3]['coords']['lat'])]))
        count += 1
        
    np.all(results)
