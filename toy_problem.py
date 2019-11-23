import xarray as xr
import pandas as pd
import numpy as np
import re
import os
import fortranformat

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


def latlon_string_2_num(s):
    """
    parse Kettle longitude, latitude coordinate strings to numeric arrays

    Kettle provides lon and lat as strings in the format "-179.5-89.5"
    -- two decimal coordinates preceded by a minus sign (W longitude
    or S latitude) or space (E longitude or N latitude).
    """
    #regular expression to match the latitude and longitude strings
    latlon_re = r'([- ]*\d+\.\d+)'
    lon, lat = re.findall(latlon_re, s)
    return(float(lon), float(lat))


def parse_kettle_plant_flux(fname):
    fcos = pd.read_csv(fname, header=8)
    lonlat = fcos.MAX.apply(latlon_string_2_num)
    #lonlat is a pandas series of tuples.  convert to a data frame and
    #concatenate it to fcos
    lonlat = pd.DataFrame(lonlat.tolist(),
                          columns=['lon', 'lat'])
    fcos = pd.concat((fcos, lonlat), axis=1)

    return(fcos)

if __name__ == "__main__":

    global_conc_dir = os.path.join('/', 'Users', 'tim', 'work',
                                   'Data', 'Jim_global')
    ocean_conc, plant_conc = (xr.open_mfdataset(
        os.path.join(global_conc_dir, this_dir, 'out_monthly', '*.nc'),
        concat_dim='tid',
        combine='nested') for this_dir in ['ocean_only', 'plant_only'])
    anthro_conc = xr.open_dataset(os.path.join(global_conc_dir, 'anthro.nc'))

    foo = parse_kettle_flux('/Users/tim/work/Data/kettle_fluxes/ocean_cos.dat')
    #kettle_focean = parse_kettle_plant_flux('/Users/tim/work/Data/kettle_fluxes/ocean_cos.dat')
