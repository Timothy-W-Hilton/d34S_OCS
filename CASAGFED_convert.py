import numpy as np
import xarray as xr
import pandas as pd
import calendar
from stem_pytools.domain import calc_grid_area

def calc_cos_plant_uptake(GEE, LRU, CO2, COS):
    """
    calculate plant COS uptake from CO2 gross ecosystem
    productivity, leaf relative uptake, atmospheric CO2
    concentration, and atmospheric COS concentration.
    INPUT PARAMETERS:
    GEE: np.ndarray of gross primary production, kg C m-2 s-1
    LRU: leaf relative uptake
        (umol CO2 m-2 s-1 (ppm CO2)-1) /
        (pmol COS m-2 s-1 (ppt COS)-1)
    CO2: atmospheric CO2 concentration (ppm)
    COS: atmospheric COS concentration (ppt)

    RETURNS:
    plant COS flux (mol COS m-2 yr-1)

    NOTES:
    LRU, CO2, and COS may be numpy arrays of any dimensions that
    are broadcastable to self.GPP_124x124.shape.  This permits
    flexibility in specifying each as an invariant scalar, a
    spatially varying but temporally invariant 2 dimensional
    array, or a temporally and spatially varying 3 or 4
    dimensional array.
    """
    #define some constants for unit conversion
    g_per_kg = 1e3   #grams per kilogram
    molC_per_gC = 1.0 / 12.011   #moles carbon per gram carbon
    umol_per_mol = 1e6    #micromoles per mole
    mol_per_pmol = 1e-12
    #calculate the COS plant flux in pmol m-2 s-1
    f_COS_plant = (GEE * LRU * (COS/CO2) *
                   g_per_kg * molC_per_gC * umol_per_mol * mol_per_pmol)

    return(f_COS_plant)


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
    DEBUG_FLAG = True
    if DEBUG_FLAG:
        print("*LAUGH TEST* total Earth "
              "area (km^2): {:0.2e}".format(area.sum()))
    return(area)


def secs_per_month():
    days_per_month = np.array([calendar.monthrange(2012, m)[1]
                               for m in range(1, 13)])
    secs_per_day = 60 * 60 * 24
    return(days_per_month * secs_per_day)


def casa_aggregate(gee):
    dates = pd.date_range('2015-01-01',
                          periods=gee['t-3hr'].size,
                          freq='3H')
    # lons, lats = np.meshgrid(gee.lon, gee.lat)
    # gee = gee.assign(area=(('lat', 'lon'),
    #                        get_area_all_gridcells(lons, lats, 1.5, 1)))
    gee['GEE'] = gee['GEE'].assign_attrs(units="kgC m-2 s-1")
    gee = gee.assign(fOCS_s=xr.DataArray(
        data=calc_cos_plant_uptake(gee.GEE.values, 1.61, 1.1, 1.0),
        dims=('t-3hr', 'lat', 'lon'),
        name='fOCS_s'))
    gee['fOCS_s'] = gee['fOCS_s'].assign_attrs(units="pmol m-2 s-1")
    secs_per_3_hrs = 60 * 60 * 3
    gee['GEE'] = gee['GEE'] * secs_per_3_hrs
    gee = gee.assign(date=(('t-3hr'), dates),
                     month=(('t-3hr'), dates.month))
    gee = gee.groupby(gee['month']).sum()
    gee['fOCS'] = gee['fOCS_s'].assign_attrs(units="pmol m-2 mon-1")
    gee['GEE'] = gee['GEE'].assign_attrs(units="kgC m-2 mon-1")
    gee.assign_attrs(description=('OCS plant flux calculated from '
                                  'CASA-GFED CO2 GPP, LRU = 1.61, '
                                  'OCS/CO2 = 1.1'))
    return(gee)

# def CASA_GPP_2_monthly_fOCS():
#     """
#     """
if __name__ == '__main__':
    fname = '/project/projectdirs/m2319/Data/CASA/GEE.3hrly.1x1.25.2015.nc'
    gee = xr.open_dataset(fname)
    gee = casa_aggregate(gee)
    gee.to_netcdf('./CASAGFED_fOCS.nc')
