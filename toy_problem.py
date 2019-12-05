import xarray as xr
import pandas as pd
import numpy as np
import re
import os
import fortranformat
import idwnn_kdtree_interp as interp

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


def format_dataset(this_dataset):
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


def main():
    global_conc_dir = os.path.join('/', 'Users', 'tim', 'work',
                                   'Data', 'Jim_global')
    ocean_conc, plant_conc = (xr.open_mfdataset(
        os.path.join(global_conc_dir, this_dir, 'out_monthly', '*.nc'),
        concat_dim='tid',
        combine='nested') for this_dir in ['ocean_only', 'plant_only'])
    ocean_conc = format_dataset(ocean_conc)
    plant_conc = format_dataset(plant_conc)

    anthro_conc = xr.open_dataset(os.path.join(global_conc_dir, 'anthro.nc'))
    # remove singleton dimemsion time.  the timestamp is still present
    # in tstep
    anthro_conc = anthro_conc.sel(time=0.0).drop('time')
    # make tstep a coordinate
    anthro_conc = anthro_conc.assign_coords(
        tstep=range(anthro_conc.dims['tstep']))
    anthro_conc.rename({'lat': 'latitude',
                        'lon': 'longitude'})

    ocean_flux_kettle_grid = parse_kettle_flux(
        os.path.join('/', 'Users', 'tim', 'work', 'Data',
                     'kettle_fluxes', 'ocean_cos.dat'))
    grid_lon, grid_lat = np.meshgrid(ocean_conc['longitude'].values,
                                     ocean_conc['latitude'].values)
    ocean_flux = regrid_fcos(ocean_flux_kettle_grid, grid_lon, grid_lat)

    da_of = xr.DataArray(ocean_flux,
                         coords={'month': range(12),
                                 'longitude': ocean_conc.longitude.values,
                                 'latitude': ocean_conc.latitude.values},
                         dims=['month', 'latitude', 'longitude'],
                         name='ocean_flux')

    return(da_of, ocean_flux, ocean_flux_kettle_grid,
           anthro_conc, ocean_conc, plant_conc)


if __name__ == "__main__":

    (da_of, ocean_flux, ocean_flux_kettle_grid,
     anthro_conc, ocean_conc, plant_conc) = main()
