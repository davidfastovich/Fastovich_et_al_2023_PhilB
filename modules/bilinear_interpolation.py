from scipy import interpolate
import xarray as xr
import numpy as np

def bilinear_interpolation(new_lat, new_lon, model_lat, model_lon, climate_var):
    '''
    Bilinearly interpolate gridded climate data to a point using 
    scipy.interpolate.RectBivariateSpline. Setting the bivariate spline 
    parameters to:
        
        kx = 1
        ky = 1
        s = 0
    
    reduces the spline to a bilinear interpolation.
    
    This only works with source data that is on a rectengular grid. The results 
    are confirmed to be accurate by comparing output with manually calculated
    bilinearly interpolated simulated data using a test function that solves
    the formula for bilinear interpolation from wikipedia.
    
    Parameters
    ----------
    new_lat : int, float, or array
        Latitude value where the climate data will be interpolated.
    new_lon : int, float, or array
        Longitude value where the climate data will be interpolated.
    model_lat : array
        Array of latitude values the model. Must be in ascending order.
    model_lon: array
        Array of latitude values the model. Must be in ascending order
    climate_var: array
        Climate variable of shape (model_lat, model_lon) that will be
        interpolated.

    Returns
    -------
    xarray of interpolated climate value(s).

    '''

    # Build model of the coordinate-climate data
    model = interpolate.RectBivariateSpline(x=model_lon, 
                                            y=model_lat, 
                                            z=climate_var,
                                            kx=1,
                                            ky=1,
                                            s=0)
    
    # Use model to predict climate data at the new coordinates
    return(model.ev(xi=new_lon, yi=new_lat))

def gridded_bilinear_interpolation(new_lat, new_lon, model_lat, model_lon, climate_var, varname):
    '''
    Bilinearly interpolate climate data to a point using 
    scipy.interpolate.RectBivariateSpline. Setting the bivariate spline 
    parameters to:
        
        kx = 1
        ky = 1
        s = 0
    
    reduces the spline to a bilinear interpolation.
    
    This only works with source data that is on a rectengular grid. The results 
    are confirmed to be accurate by comparing output with manually calculated
    bilinearly interpolated simulated data using a test function that solves
    the formula for bilinear interpolation from wikipedia.
    
    Parameters
    ----------
    new_lat : array
        Latitude value where the climate data will be interpolated.
    new_lon :array
        Longitude value where the climate data will be interpolated.
    model_lat : array
        Array of latitude values the model. Must be in ascending order.
    model_lon: array
        Array of latitude values the model. Must be in ascending order
    climate_var: array
        Climate variable of shape (model_lat, model_lon) that will be
        interpolated.

    Returns
    -------
    Array of interpolated climate value(s).

    '''

    # Build model of the coordinate-climate data
    model = interpolate.RectBivariateSpline(x=model_lon, 
                                            y=model_lat, 
                                            z=climate_var,
                                            kx=1,
                                            ky=1,
                                            s=0)
    
    # Before evaluating the model at new coordinates new_lat and new_lon
    # have to be reshaped into a flat, 2 column array
    grid_lat, grid_lon = np.meshgrid(new_lat, new_lon)
    flat_coords = np.column_stack([grid_lon.ravel(), grid_lat.ravel()])
    
    # Use model to predict climate data at the new coordinates
    interp_clim = model.ev(xi=flat_coords[:,0], yi=flat_coords[:,1])
    
    # Reshaping the interpolated data back to the correct lat-lon grid
    interp_clim = interp_clim.reshape(grid_lon.shape).transpose()
    
    # Prepare xarray Dataset as the return product
    # SAVE AS DATA ARRAY
    ds = xr.Dataset(
        data_vars=dict(
            var=(["lat", "lon"], interp_clim)
        ),
        coords=dict(
            lon=(["lon"], new_lon),
            lat=(["lat"], new_lat),
        ),
    )
    ds = ds.rename({'var':varname})

    return(ds)

def bilinear_interpolation_manual(x, y, points):
    '''Interpolate (x,y) from values associated with four points.

    The four points are a list of four triplets:  (x, y, value).
    The four points can be in any order.  They should form a rectangle.

        >>> bilinear_interpolation(12, 5.5,
        ...                        [(10, 4, 100),
        ...                         (20, 4, 200),
        ...                         (10, 6, 150),
        ...                         (20, 6, 300)])
        165.0

    '''
    # See formula at:  http://en.wikipedia.org/wiki/Bilinear_interpolation

    points = sorted(points)               # order points by x, then by y
    (x1, y1, q11), (_x1, y2, q12), (x2, _y1, q21), (_x2, _y2, q22) = points

    if x1 != _x1 or x2 != _x2 or y1 != _y1 or y2 != _y2:
        raise ValueError('points do not form a rectangle')
    if not x1 <= x <= x2 or not y1 <= y <= y2:
        raise ValueError('(x, y) not within the rectangle')

    return (q11 * (x2 - x) * (y2 - y) +
            q21 * (x - x1) * (y2 - y) +
            q12 * (x2 - x) * (y - y1) +
            q22 * (x - x1) * (y - y1)
           ) / ((x2 - x1) * (y2 - y1) + 0.0)


if __name__ == '__main__':
    
    import numpy.random as random
    from modules import nearest_neighbor
    
    ########################################################################
    # TEST IF RESULTS ARE CONSISTENT WITH FORMULA FOR BILINEAR INTERPOLATION
    ########################################################################
    
    comparison = []
    count = 0
    while count < 10000:
    
        # Simulate point data where temperature will be interpolated
        new_lat = random.uniform(low=-90, high=90, size=1)
        new_lon = random.uniform(low=0, high=360, size=1)
        
        # Simulate climate model coordinates
        model_lat = np.linspace(start=-90, stop=90, num=10)
        model_lon = np.linspace(start=0, stop=360, num=20)
        
        # Create synthetic temperature data
        climate_var = random.uniform(low=-30, high=30, size=(model_lon.shape[0], model_lat.shape[0]))
        
        # Interpolation results
        clim_interp = bilinear_interpolation(new_lat=new_lat, 
                                             new_lon=new_lon, 
                                             model_lat=model_lat, 
                                             model_lon=model_lon,
                                             climate_var=climate_var)
        
        # Test output against a function that solves the bilinear interpolation system of equations
        nearest_neighbors = nearest_neighbor(point_lat=new_lat, 
                                             point_lon=new_lon, 
                                             model_lat=model_lat, 
                                             model_lon=model_lon)
        
        # Coerce nearest neighbors into proper format for test function
        climate_var_t = climate_var.transpose()
        points = [(nearest_neighbors[0]['coords']['lon'], nearest_neighbors[0]['coords']['lat'], climate_var_t[nearest_neighbors[0]['index']['lat'], nearest_neighbors[0]['index']['lon']]),
                  (nearest_neighbors[1]['coords']['lon'], nearest_neighbors[1]['coords']['lat'], climate_var_t[nearest_neighbors[1]['index']['lat'], nearest_neighbors[1]['index']['lon']]),
                  (nearest_neighbors[2]['coords']['lon'], nearest_neighbors[2]['coords']['lat'], climate_var_t[nearest_neighbors[2]['index']['lat'], nearest_neighbors[2]['index']['lon']]),
                  (nearest_neighbors[3]['coords']['lon'], nearest_neighbors[3]['coords']['lat'], climate_var_t[nearest_neighbors[3]['index']['lat'], nearest_neighbors[3]['index']['lon']])]
        
        # Manual interpolation results
        clim_interp_manual = bilinear_interpolation_manual(x=new_lon, y=new_lat, points=points)
        
        # Compare calculations
        comparison.append(np.round(clim_interp, 6) == np.round(clim_interp_manual, 6))
        count += 1
    
    # Test if manual solving and scipy aggree - they do
    np.all(comparison)
    
    #######################################################
    # TEST IF MODEL INPUT PRODUCES SANE INTERPOLATED VALUES
    #######################################################
    
    # Model
    model = xr.open_dataset('/media/david/DAV/models/trace_mwf_seasonal/b30.00_4kaDVTn.cam2.ncrcat.djf.nc')
    model = xr.open_dataset('/media/david/DAV/models/hadcm3_ivanovic_2018_paleo/hs1/tecad.temp2m.monthly.nc')
    model = model.sortby('latitude')
    
    # Select new_lat and new_lon close to a grid cell
    new_lat = 55
    new_lon = 281.2
    bilinear_interpolation(new_lat=new_lat, 
                           new_lon=new_lon, 
                           model_lat=model.latitude.values, 
                           model_lon=model.longitude.values, 
                           climate_var=model.isel(ht=0).isel(t=0)['temp'].values.transpose())
    
    model.isel(ht=0).isel(t=0).sel(latitude=55).sel(longitude=281.25)['temp'].values
