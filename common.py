import netCDF4
import numpy as np

def get_lat_lon(netcdf_path):
    """read lat and lon from a model data netCDF file

    The lat and lon netCDF variables must be named "lat" and "lon".

    ARGS:
    netcdf_path (str): full path to the data file

    RETURNS:
    Two 2D numpy arrays of longitudes and latitudes for each model grid cell.
    """
    nc = netCDF4.Dataset(netcdf_path)
    lat = nc.variables['lat'][:]
    lon = nc.variables['lon'][:]
    lat_grid, lon_grid = np.meshgrid(lat, lon)
    return((lat_grid, lon_grid))
