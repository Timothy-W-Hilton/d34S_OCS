"""calculate annual total OCS fluxes from Jim's model runs
"""

import netCDF4
import numpy as np
import os.path
import calendar
import pandas as pd

from stem_pytools.domain import calc_grid_area
from common import get_lat_lon

# define paths to data files as module globals
rootdir = os.path.join('/', 'global', 'cscratch1', 'sd', 'jstineci',
                       'geos_fluxes_for_visual')
ocean_dir = os.path.join(rootdir, 'OceanTotal', 'posterior-12month-kfl')
anthro_dir = os.path.join(rootdir, 'anthro_v3', 'v3_anthro', '2015')

def get_secs_per_month(year):
    """return 12-element array containing seconds in each month of given year
    """
    secs_per_day = 24 * 60 * 60
    return(np.array([calendar.monthrange(year, mon)[1]
                     for mon in range(1, 13)])
           * secs_per_day)

def get_area_all_gridcells(lon, lat, d_lon = 2.5, d_lat = 2.0):
    """calculate area in km^2 of every cell in two 2d arrays of lats, lons

    assumes lat, lon specify the center of the cell
    """
    area = np.zeros(lat.shape)
    for i in range(lat.shape[0]):
        for j in range(lat.shape[1]):
            area[i, j]  = calc_grid_area(lon[i, j] + (d_lon / 2.0),
                                         lon[i, j] - (d_lon / 2.0),
                                         lat[i, j] + (d_lat / 2.0),
                                         lat[i, j] - (d_lat / 2.0))
    area = area * 1e-6  # convert m^2 to km^2
    print("*LAUGH TEST* total Earth area (km^2): {:0.2e}".format(area.sum()))
    return(area)

def get_fluxes(dir, year, fluxes_dim=(12, 91, 144)):
    """read a year's worth of monthly OCS fluxes and calculate annual total

    RETURNS:
    pandas dataframe containing annual total ocean and anthropogenic
    OCS fluxes in gG S
    """
    fnames = [os.path.join(dir, "{:02d}.nc".format(mon))
                           for mon in range(1, 13)]
    fluxes = np.zeros(fluxes_dim)
    for this_mon in range(12):
        nc = netCDF4.Dataset(fnames[this_mon])
        fluxes[this_mon, ...] = nc.variables['COS_Flux'][...]
        nc.close()
    # lat, lon don't change so read once, outside the loop
    lat_grid, lon_grid = get_lat_lon(fnames[this_mon])
    areas = get_area_all_gridcells(lon_grid.transpose(), lat_grid.transpose())

    gG_per_kg = 1e-6
    secs_per_month = get_secs_per_month(year)[:, np.newaxis, np.newaxis]
    total = (secs_per_month * fluxes * areas * gG_per_kg).sum(axis=0)
    return(total)


if __name__ == "__main__":
    years = range(2004, 2008)
    df = pd.DataFrame({'anthro': None, 'ocean':None}, index=years)
    for this_year in years:
        df.anthro[this_year] = get_fluxes(anthro_dir, 2015).sum()
        this_ocean_dir = os.path.join(ocean_dir, "{:04d}".format(this_year))
        df.ocean[this_year] = get_fluxes(this_ocean_dir, this_year).sum()

    # fluxes are in kG S/Km^2/S -- convert kg S to gG S
