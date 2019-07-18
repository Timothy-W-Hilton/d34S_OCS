import numpy as np
import netCDF4
import warnings


def get_ocean_anomaly(infile='./ocean.nc'):
    """read anthro [OCS] into 48x191x90 array
    """
    nc = netCDF4.Dataset(infile)
    OCS = nc.variables['OCS'][...].squeeze()
    OCS = OCS.reshape( 48, 144, 91, 47)
    OCS = np.rollaxis(OCS, 3, 1)
    OCS = np.rollaxis(OCS, 3, 2)
    # subtract time-varying global spatial mean OCS to get anomalies
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore",
                                category=RuntimeWarning,
                                message="Mean of empty slice")
        OCS = OCS - np.apply_over_axes(np.nanmean, OCS, (1, 2, 3))
    nc.close()
    return(OCS)

def get_anthro_anomaly(infile='./anthro.nc'):
    """read anthro [OCS] into 48x191x90 array
    """
    nc = netCDF4.Dataset(infile)
    OCS = nc.variables['COS'][...].squeeze()
    # subtract time-varying global spatial mean OCS to get anomalies
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore",
                                category=RuntimeWarning,
                                message="Mean of empty slice")
        OCS = OCS - np.apply_over_axes(np.nanmean, OCS, (1, 2, 3))
    nc.close()
    return(OCS)

if __name__ == "__main__":

    ocs_anthro = get_anthro_anomaly()
    ocs_ocean = get_ocean_anomaly()
