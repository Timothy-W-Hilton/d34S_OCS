import numpy as np
import netCDF4

def parse_ocean(infile='./ocean.nc'):
    """read anthro [OCS] into 48x191x90 array
    """
    nc = netCDF4.Dataset(infile)
    OCS = nc.variables['OCS'][...].squeeze()
    OCS = OCS.reshape( 48, 144, 91, 47)
    OCS = np.rollaxis(OCS, 3, 1)
    OCS = np.rollaxis(OCS, 3, 2)
    nc.close()
    return(OCS)

def parse_anthro(infile='./anthro.nc'):
    """read anthro [OCS] into 48x191x90 array
    """
    nc = netCDF4.Dataset(infile)
    OCS = nc.variables['COS'][...].squeeze()
    nc.close()
    return(OCS)

if __name__ == "__main__":

    ocs_anthro = parse_anthro()
    ocs_ocean = parse_ocean()
