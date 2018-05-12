import numpy as np
import netCDF4

def get_ocean_anomaly(infile='./ocean.nc'):
    """read anthro [OCS] into 48x191x90 array
    """
    nc = netCDF4.Dataset(infile)
    OCS = nc.variables['OCS'][...].squeeze()
    OCS = OCS.reshape( 48, 144, 91, 47)
    OCS = np.rollaxis(OCS, 3, 1)
    OCS = np.rollaxis(OCS, 3, 2)
    OCS = OCS - OCS.mean() # TODO: need apply_over_axis or apply_along_axis
    nc.close()
    return(OCS)

def get_anthro_anomaly(infile='./anthro.nc'):
    """read anthro [OCS] into 48x191x90 array
    """
    nc = netCDF4.Dataset(infile)
    OCS = nc.variables['COS'][...].squeeze()
    OCS = OCS - OCS.mean()
    nc.close()
    return(OCS)

if __name__ == "__main__":

    ocs_anthro = get_anthro_anomaly()
    ocs_ocean = get_ocean_anomaly()
