import netCDF4
import numpy as np
import pandas as pd
from stem_pytools.domain import find_nearest_stem_xy

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

def assign_grid_cell_to_sites(sites, lon_grid, lat_grid):
    """find the nearest model grid cell to a pandas DataFrame of NOAA locations

    The DataFrame must contain columns named "Latitude" and
    "Longitude".  Typically the lon_grid and lat_grid arguments will
    come from get_lat_lon()

    ARGS:
    sites (pandas.DataFrame): DataFrame containing, at minimum, site
        latitude and longitude coordinates
    lon_grid (array-like): 2D numpy array of model cell longitudes
    lat_grid (array-like): 2D numpy array of model cell latitudes

    RETURNS:
    sites, with columns "X" and "Y" added.
    """
    sites['X'], sites['Y'] = find_nearest_stem_xy(sites.Longitude, sites.Latitude,
                                            lon_grid, lat_grid)
    return(sites)

def gather_sites_data():
    """main function: parse NOAA sites table and assign model cell X, Y coords

    wrapper function for assign_grid_cell_to_sites(), get_lat_lon()

    gets model cell latitude and longitude coordinates from
    /global/cscratch1/sd/twhilton/global_sandbox/anthro.nc and NOAA
    site data from ./noaa_flask_sites.dat

    """
    anthro_data = '/global/cscratch1/sd/twhilton/global_sandbox/anthro.nc'
    sites = pd.read_csv('./noaa_flask_sites.dat', sep='\t')
    lat_grid, lon_grid = get_lat_lon(anthro_data)
    sites = assign_grid_cell_to_sites(sites, lon_grid, lat_grid)
    return(sites)
