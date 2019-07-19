import numpy as np
import netCDF4
import warnings

def get_ocean_production(infile='./ocean.nc'):
    """read ocean [OCS] into 48x191x90 array
    """
    nc = netCDF4.Dataset(infile)
    # read OCS from file.  It has 4 years of data with years, months
    # as separate axes
    OCS_y_m = nc.variables['OCS'][...].squeeze()
    # copy the data into a new array with a continuous time axis
    # OCS_y_m = np.rollaxis(OCS_y_m, 3, 1)
    # OCS_y_m = np.rollaxis(OCS_y_m, 3, 2)
    OCS = np.full((48, 47, 91, 144), fill_value=np.nan)
    idx = 0
    for y in range(OCS_y_m.shape[0]):
        for m in range(OCS_y_m.shape[1]):
            OCS[idx, ...] = OCS_y_m[y, m, ...]
            idx = idx + 1
    nc.close()
    return(OCS)

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
        OCS_anom = OCS - np.apply_over_axes(np.nanmean, OCS, (1, 2, 3))
    import pdb; pdb.set_trace()
    return(OCS_anom)

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
