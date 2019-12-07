import xarray as xr
import pandas as pd
import numpy as np
import re
import os
import fortranformat
import idwnn_kdtree_interp as interp
import geoviews as gv
import xesmf as xe
from geoviews.operation.regrid import weighted_regrid
import calendar

from stem_pytools.domain import calc_grid_area
from holoviews.operation.datashader import regrid

DEBUG_FLAG = True
SECS_PER_MINUTE = 60
MINS_PER_HOUR = 60
HOURS_PER_DAY = 24
MONTHS_PER_YEAR = 12
SECS_PER_DAY = SECS_PER_MINUTE * MINS_PER_HOUR * HOURS_PER_DAY
DAYS_PER_MONTH = np.array([calendar.monthrange(2012, m)[1]
                           for m in range(1, 13)])
SECS_PER_MONTH = SECS_PER_DAY * DAYS_PER_MONTH

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


def parse_kettle_flux(fname):
    f = open(fname, 'r')
    found_format = False
    while not(found_format):
        this_line = f.readline()
        if "format" in this_line:
            found_format = True
            discard, fmt = this_line.split('=')
    reader = fortranformat.FortranRecordReader(fmt)
    f.readline()  # discard blank line
    remaining_content = f.readlines()
    data = np.array([reader.read(f) for f in remaining_content])
    lon = data[:, 0]
    lat = data[:, 1]
    mask = data[:, 2]
    data_reshaped = data[:, 3:].reshape(36, 72, 12)
    ds = xr.DataArray(data_reshaped,
                      coords={'month': range(12),
                              'lon': lon[:72],
                              'lat': lat[::72]},
                      dims=['lat', 'lon', 'month'],
                      name='ocean_flux')
    ds.name = 'ocean_flux'
    return(ds)


def regrid_fcos(fcos, new_lon, new_lat):
    """
    Regrid Kettle plant COS flux to new latitude, longiutde

    INPUT PARAMETERS
    fcos; pandas DataFrame: the Kettle plant fluxes.  This would
        generally be the output of parse_kettle_plant_flux
    fname_STEM_coords; string: full path to STEM topo file
    """
    old_lon, old_lat = np.meshgrid(fcos['lon'].values,
                                   fcos['lat'].values)
    fcos_regridded = interp.spherical_idw_kdtree_interp(
        old_lon, old_lat,
        new_lon, new_lat,
        data=np.moveaxis(fcos.data, 2, 0),
        n_nbr=1)
    return(fcos_regridded)


def get_kettle_data(data_array_name, fname, lon, lat):
    """parse & regrid Kettle flux data and place into xarray.DataArray
    """
    flux_kettle_grid = parse_kettle_flux(
        os.path.join('/', 'Users', 'tim', 'work', 'Data',
                     'kettle_fluxes', fname))
    grid_lon, grid_lat = np.meshgrid(lon, lat)
    flux_geoschem_grid = regrid_fcos(flux_kettle_grid, grid_lon, grid_lat)
    flux_da = xr.DataArray(flux_geoschem_grid,
                           coords={'month': range(12),
                                   'longitude': lon,
                                   'latitude': lat},
                           dims=['month', 'latitude', 'longitude'],
                           name=data_array_name)
    return(flux_da)


def format_GEOSChem_dataset(this_dataset):
    """promote longitude, latitude, lev from variable to coordinate
    variable; create a new coordinate variable tstep
    """
    this_dataset = this_dataset.assign_coords(
        longitude=this_dataset.longitude[0],
        latitude=this_dataset.latitude[0],
        lev=this_dataset.lev[0],
        tstep=xr.Variable(dims='tid',
                          data=range(this_dataset.dims['tid'])))
    return(this_dataset)


def get_CASAGFED_plant_ocs(new_lon1d, new_lat1d):
    """read CASA-GFED plant OCS flux
    """
    fOCS = xr.open_dataset('/Users/tim/work/Data/CASA_GFED/CASAGFED_fOCS.nc')
    grid = xe.util.grid_2d(new_lon1d.min() - 1.25,
                           new_lon1d.max() + 1.25,
                           2.5,
                           new_lat1d.min() - 1.0,
                           new_lat1d.max() + 1.0,
                           2.0)
    target = gv.Dataset(grid, kdims=['lon', 'lat'])

    # fOCS_new = xr.Dataset(
    #     data_vars= {'plant_flux': (('month', 'lat', 'lon'),
    #                                 fOCS['fOCS'])},
    #     coords= {'lon': (('lon'), fOCS['lon']),
    #              'lat': (('lat'), fOCS['lat']),
    #              'month': (('month'), range(1, 13))})
    #regrid adapted from
    images = gv.Dataset(fOCS).to(gv.Image, ['lon' ,'lat'])
    fOCS_lowres_quadmesh = weighted_regrid(images,
                                           target=target,
                                           streams=[],
                                           reuse_weights=True)
    # fOCS_lowres = [fOCS_lowres_quadmesh.data[(m, )].data for m in range(1, 13)]
    fOCS_lowres = xr.concat([fOCS_lowres_quadmesh.data[(m, )].data
                             for m in range(1, 13)],
                            dim='month')
    return(fOCS, fOCS_lowres)

def get_andrew_antho_cos(new_lon1d, new_lat1d):
    fname = '/Users/tim/work/Data/Anthro_COS/Total_Anth._COS1980-2012_v3.nc'
    anth_ocs_high_res = xr.open_dataset(fname)
    anth_ocs_high_res = anth_ocs_high_res.rename({'ul_longitude': 'lon',
                                                  'ul_latitude': 'lat'})
    anth_ocs_high_res = anth_ocs_high_res.rename(
        {'Total_Anth._COS': 'anthro_flux',
         'lon': 'longitude',
         'lat': 'latitude'})
    anth_ocs_raw = anth_ocs_high_res.assign_coords(
        longitude=anth_ocs_high_res.longitude,
        latitude=anth_ocs_high_res.latitude,
        years=anth_ocs_high_res.latitude)
    # This (below) dataset format works for plotting for reasons I do
    # not yet understand
    anth_highres = xr.Dataset(
        data_vars= {'anthro_flux': (('year', 'latitude', 'longitude'),
                                    anth_ocs_raw['anthro_flux'])},
        coords= {'longitude': (('longitude'),
                               anth_ocs_raw['longitude']),
                 'latitude': (('latitude'),
                              anth_ocs_raw['latitude']),
                 'year': (('year'),
                          range(anth_ocs_raw['t'].size)),
                 'month': (('month'),
                           range(1, 13))})
    #regrid adapted from
    #http://geoviews.org/user_guide/Resampling_Grids.html
    gvds = gv.Dataset(anth_highres.sel(year=32),
                      kdims=['longitude', 'latitude'])
    grid = xe.util.grid_2d(new_lon1d.min() - 1.25,
                           new_lon1d.max() + 1.25,
                           2.5,
                           new_lat1d.min() - 1.0,
                           new_lat1d.max() + 1.0,
                           2.0)
    target = gv.Dataset(grid, kdims=['lon', 'lat'])
    anth_lowres = weighted_regrid(gvds, target=target, streams=[]).data
    anth_lowres = anth_lowres.assign_coords(month=range(1, 13))
    anth_lowres['anthro_flux'] = anth_lowres['anthro_flux'].assign_attrs(
        units='pmol m-2 s-1')

    #flux (pmol / m2/ s)
    focs_permonth = xr.DataArray(
        data=np.tile(anth_lowres['anthro_flux'].values / MONTHS_PER_YEAR,
                     reps = (MONTHS_PER_YEAR, 1, 1)),
        dims=(('month', 'y', 'x')))
    #flux (pmol / m2/ month)
    focs_permonth = focs_permonth * SECS_PER_MONTH[:, np.newaxis, np.newaxis]
    #flux (pmol / gridcell / month)
    area_per_gridcell = get_area_all_gridcells(target.data['lon'].data,
                                               target.data['lat'].data)
    focs_permonth = focs_permonth * area_per_gridcell[np.newaxis, ...]
    focs_permonth = focs_permonth.assign_attrs(units='pmol / gridcell / month')
    anth_lowres = anth_lowres.assign(anthro_flux_permonth=focs_permonth)

    return(anth_highres, anth_lowres)

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
    if DEBUG_FLAG:
        print("*LAUGH TEST* total Earth "
              "area (km^2): {:0.2e}".format(area.sum()))
    return(area)


def main():
    global_conc_dir = os.path.join('/', 'Users', 'tim', 'work',
                                   'Data', 'Jim_global')
    ocean_conc, plant_conc = (xr.open_mfdataset(
        os.path.join(global_conc_dir, this_dir, 'out_monthly', '*.nc'),
        concat_dim='tid',
        combine='nested') for this_dir in ['ocean_only', 'plant_only'])
    ocean_conc = format_GEOSChem_dataset(ocean_conc)
    plant_conc = format_GEOSChem_dataset(plant_conc)

    lon_geoschem = ocean_conc['longitude'].values
    lat_geoschem = ocean_conc['latitude'].values

    anthro_conc = xr.open_dataset(os.path.join(global_conc_dir, 'anthro.nc'))
    # remove singleton dimemsion time.  the timestamp is still present
    # in tstep
    anthro_conc = anthro_conc.sel(time=0.0).drop('time')
    # make tstep a coordinate
    anthro_conc = anthro_conc.assign_coords(
        tstep=range(anthro_conc.dims['tstep']))
    anthro_conc.rename({'lat': 'latitude',
                        'lon': 'longitude'})

    da_of = get_kettle_data('ocean_flux',
                            'ocean_cos.dat',
                            lon_geoschem,
                            lat_geoschem)
    anth_highres, anth_lowres = get_andrew_antho_cos(
        lon_geoschem,
        lat_geoschem)
    plant_highres, plant_lowres = get_CASAGFED_plant_ocs(
        lon_geoschem,
        lat_geoschem)

    # combine all three flux and concentration datasets into a single
    # xarray dataset
    n_months = 48
    # repeat the 12-month flux time series 4 times to match the
    # 48-month concentration time series
    fanthro_48month = xr.concat([anth_lowres['anthro_flux_permonth']] * 4,
                                dim='month')
    focean_48month = xr.concat([da_of] * 4, dim='month')
    fplant_48month = xr.concat([plant_lowres['fOCS']] * 4, dim='month')
    OCS_all = xr.Dataset(
        data_vars= {'anthro_flux': (('month', 'lat', 'lon'), fanthro_48month),
                    'ocean_flux': (('month', 'lat', 'lon'), focean_48month),
                    'plant_flux': (('month', 'lat', 'lon'), fplant_48month),
                    'anthro_conc': (('month', 'lat', 'lon'),
                                    anthro_conc.sel(lev=1).drop('lev')['COS']),
                    'ocean_conc': (('month', 'lat', 'lon'),
                                   ocean_conc.sel(zid=1).drop('lev')['OCS']),
                    'plant_conc': (('month', 'lat', 'lon'),
                                   plant_conc.sel(zid=1).drop('lev')['OCS'])},
        coords = {'lon': (('lon'), lon_geoschem),
                  'lat': (('lat'), lat_geoschem),
                  'month': (('month'), range(n_months))})

    return(OCS_all)

if __name__ == "__main__":

    # todo calculate ocean flux per gridcell per month
    # todo make sure units are consistent
    OCS_all = main()
