import numpy as np
import os
import os.path
from glob import glob
from nco import Nco
import calendar

def avg_anthro(p, workdir):
    """average daily files of hourly [OCS] into monthly means

    Anthro root directory should contain subdirectories labeled by
    year; each annual subdirectory should contain monthly
    subdirectories 01..12.

    ARGS:
    p (str): full path to anthro root directory
    """
    nco = Nco()
    annual_dirs = sorted(glob(os.path.join(p, "20[01]*")))
    all_monthly_files = []
    for this_dir in annual_dirs:
        this_year = os.path.basename(this_dir)
        for this_month in np.arange(1, 13):
            this_wildcard = os.path.join(this_dir,
                                         "{:02d}".format(this_month),
                                         "*.nc")
            ndays = calendar.monthrange(int(this_year), this_month)[1]
            files = sorted(glob(this_wildcard))
            # ignore incomplete months
            if len(files) == ndays:
                new_fname = os.path.join(workdir,
                                         "{}-{:02d}_anthro_avg.nc".format(
                                             this_year, this_month))
                print('creating {}... '.format(new_fname), end='')
                #average hourly data on time dimension
                nco.nces(input=files, output=new_fname, options=['-d time,0'])
                #add variables for year and month
                nco.ncap2(input=new_fname,
                          output=new_fname,
                          options=['-s year[time]={}'.format(this_year)])
                nco.ncap2(input=new_fname,
                          output=new_fname,
                          options=['-s month[time]={}'.format(this_month)])
                # 'time' variable, for time of day, now contains 0 at
                # all timestep and is meaningless, so remove it
                nco.ncks(input=new_fname,
                         output=new_fname,
                         options=['-x -v time'])
                print('done')
                all_monthly_files.append(new_fname)
    nco.ncecat(input=all_monthly_files,
               output=os.path.join(out_dir, "anthro.nc"),
               options=['-u tstep'])
    print('created anthro.nc')
    print('removing individual monthly files')
    for this_file in all_monthly_files:
        os.remove(this_file)

if __name__ == "__main__":
    """
    """
    anthro_dir = os.path.join('/', 'global', 'cscratch1', 'sd',
                              'jstineci', 'seasonal_amplitude',
                              'anthro_flux_only')
    out_dir = os.path.join('/', 'global', 'cscratch1', 'sd',
                          'twhilton', 'global_sandbox')
    avg_anthro(anthro_dir, out_dir)
