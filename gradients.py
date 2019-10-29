import numpy as np
import netCDF4
import warnings
import xarray as xr


def get_ocean_production(infile='./ocean.nc'):
    """read ocean [OCS] into 48x191x90 array
    """
    ocean = xr.open_dataset(infile)
    ocean = ocean.assign_coords(time=range(ocean.year.size * ocean.month.size),
                                xid=ocean.xid,
                                yid=ocean.yid,
                                zid=ocean.zid,
                                year=ocean.year)
    ocean['OCS_y_m'] = ocean['OCS']
    ocean['OCS'] = xr.Variable(('time', 'zid', 'yid', 'xid'),
                               np.full((ocean.time.size,
                                        ocean.zid.size,
                                        ocean.yid.size,
                                        ocean.xid.size),
                                       fill_value=np.nan))
    # read OCS from file.  It has 4 years of data with years, months
    # as separate axes
    # copy the data into a new array with a continuous time axis
    # OCS_y_m = np.rollaxis(OCS_y_m, 3, 1)
    # OCS_y_m = np.rollaxis(OCS_y_m, 3, 2)
    idx = 0
    for y in range(ocean.OCS_y_m.year.size):
        for m in range(ocean.OCS_y_m.month.size):
            ocean.OCS[idx, ...] = ocean.OCS_y_m[y, m, ...]
            idx = idx + 1
    return(ocean)


def get_ocean_anomaly(infile='./ocean.nc', OCS=None):
    """read ocean [OCS] into 48x191x90 array
    """
    if infile is not None:
        OCS = get_ocean_production(infile)
    if OCS is None:
        raise(ValueError(("OCS was not provided and "
                          "could not be read from {}".format(infile))))
    # subtract time-varying global spatial mean OCS to get anomalies
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore",
                                category=RuntimeWarning,
                                message="Mean of empty slice")
        OCS_anom = OCS['OCS'].mean(skipna=True, dim=('xid', 'yid', 'zid'))
    return(OCS_anom.values)


def get_anthro_production(infile='./anthro.nc'):
    """read anthro [OCS] into 48x191x90 array
    """
    nc = netCDF4.Dataset(infile)
    OCS = nc.variables['COS'][...].squeeze()
    nc.close()
    return(OCS)


def get_anthro_anomaly(infile='./anthro.nc'):
    """read anthro [OCS] into 48x191x90 array
    """
    if infile is not None:
        OCS = get_anthro_production(infile)
    if OCS is None:
        raise(ValueError(("OCS was not provided and "
                          "could not be read from {}".format(infile))))
    # subtract time-varying global spatial mean OCS to get anomalies
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore",
                                category=RuntimeWarning,
                                message="Mean of empty slice")
        OCS_anom = OCS - np.apply_over_axes(np.nanmean, OCS, (1, 2, 3))
    return(OCS_anom)


if __name__ == "__main__":
    ocs_anthro = get_anthro_production()
    ocs_ocean = get_ocean_production()
    ocs_anthro_anomaly = get_anthro_anomaly()
    ocs_ocean_anomaly = get_ocean_anomaly()
