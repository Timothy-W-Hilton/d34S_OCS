import numpy as np
import os.path
from glob import glob

def avg_anthro(p, workdir):
    """average daily files of hourly [OCS] into monthly means

    Anthro root directory should contain subdirectories labeled by
    year; each annual subdirectory should contain monthly
    subdirectories 01..12.

    ARGS:
    p (str): full path to anthro root directory
    """
    annual_dirs = sorted(glob(os.path.join(p, "20[01]*")))
    print('annual dirs: {}'.format(annual_dirs))
    for this_year in annual_dirs:
        for this_month in np.arange(1, 13):
            this_wildcard = os.path.join(this_year,
                                         "{:02d}".format(this_month),
                                         "*.nc")
            files = sorted(glob(this_wildcard))
            print(this_year, this_month, files)


if __name__ == "__main__":
    """
    """
    anthro_dir = os.path.join('/', 'global', 'cscratch1', 'sd',
                              'jstineci', 'seasonal_amplitude',
                              'anthro_flux_only')
    avg_anthro(anthro_dir, None)
