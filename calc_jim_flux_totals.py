"""calculate annual total OCS fluxes from Jim's model runs
"""

import netCDF4
import numpy as np
import os.path
import calendar
import pandas as pd

# define paths to data files as module globals
rootdir = os.path.join('/', 'global', 'cscratch1', 'sd', 'jstineci',
                       'geos_fluxes_for_visual')
ocean_dir = os.path.join(rootdir, 'OceanTotal', 'posterior')
anthro_dir = os.path.join(rootdir, 'anthro_v3', 'v3_anthro', '2015')

def get_secs_per_month(year):
    """return 12-element array containing seconds in each month of given year
    """
    secs_per_day = 24 * 60 * 60
    return(np.array([calendar.monthrange(year, mon)[1]
                     for mon in range(1, 13)])
           * secs_per_day)

def get_fluxes(dir, year, fluxes_dim=(12, 91, 144)):
    """read a year's worth of monthly OCS fluxes and calculate annual total
    """
    fnames = [os.path.join(dir, "{:02d}.nc".format(mon))
                           for mon in range(1, 13)]
    fluxes = np.zeros(fluxes_dim)
    for this_mon in range(12):
        nc = netCDF4.Dataset(fnames[this_mon])
        fluxes[this_mon, ...] = nc.variables['COS_Flux'][...]
        nc.close()
    secs_per_month = get_secs_per_month(year)
    total = (secs_per_month[:, np.newaxis, np.newaxis] * fluxes).sum(axis=0)
    return(total)


if __name__ == "__main__":
    years = range(2004, 2008)
    df = pd.DataFrame({'anthro': None, 'ocean':None}, index=years)
    for this_year in years:
        df.anthro[this_year] = get_fluxes(anthro_dir, this_year).sum()
        this_ocean_dir = os.path.join(ocean_dir, "{:04d}".format(this_year))
        df.ocean[this_year] = get_fluxes(this_ocean_dir, this_year).sum()

    # fluxes are in kG S/Km^2/S -- convert kg S to gG S
