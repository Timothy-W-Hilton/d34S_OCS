import numpy as np
import os
import os.path
from glob import glob
from nco import Nco
import calendar
import tempfile

def concat_ocean_OCS(p, out_dir):
    """concatenate monthly ocean fluxes to a single file

    ARGS:
    p (str): full path to anthro root directory
    out_dir (str): full path to directory in which to place the result
    """
    nco = Nco()
    all_annual_files = sorted(glob(os.path.join(p, "COS*_monthly.nc")))
    all_tmp_files = []
    for this_annual_file in all_annual_files:
        fd, tmp = tempfile.mkstemp(suffix='.nc')
        this_year = os.path.basename(this_annual_file)[3:7]
        nco.ncap2(input=this_annual_file,
                  output=tmp,
                  options=['-s year[tid]={}'.format(this_year)])
        # rename variable 'tid' (timestep id) to 'month'
        nco.ncrename(input=tmp,
                     output=tmp,
                     options=['-d tid,month'])
        nco.ncrename(input=tmp,
                     output=tmp,
                     options=['-v mon,month'])

        all_tmp_files.append(tmp)
    outfile = os.path.join(out_dir, "ocean.nc")
    print('concatenating ', all_tmp_files)
    nco.ncecat(input=all_tmp_files,
               output=outfile,
               options=['-u year'])
    print("created {}".format(outfile))
    for this_tmp in all_tmp_files:
        os.remove(this_tmp)

def avg_anthro(p, out_dir):
    """average daily files of hourly [OCS] into monthly means

    Creates a single netCDF file of monthly mean global [OCS]

    Anthro root directory should contain subdirectories labeled by
    year; each annual subdirectory should contain monthly
    subdirectories 01..12.

    ARGS:
    p (str): full path to anthro root directory
    out_dir (str): full path to directory in which to place the result
    """
    nco = Nco()
    annual_dirs = sorted(glob(os.path.join(p, "20[01]*")))
    all_annual_files = []
    for this_dir in annual_dirs:
        all_monthly_files = []
        this_year = os.path.basename(this_dir)
        for this_month in np.arange(1, 13):
            this_wildcard = os.path.join(this_dir,
                                         "{:02d}".format(this_month),
                                         "*.nc")
            ndays = calendar.monthrange(int(this_year), this_month)[1]
            files = sorted(glob(this_wildcard))
            # ignore incomplete months
            if len(files) == ndays:
                new_fname = os.path.join(out_dir,
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
        this_annual_file = os.path.join(out_dir, "anthro_{}.nc".format(this_year))
        nco.ncecat(input=all_monthly_files,
                   output=this_annual_file,
                   options=['-u tstep'])
        all_annual_files.append(this_annual_file)
    # now concatenate annual files together
    nco.ncecat(input=all_annual_files,
               output=os.path.join(outdir, 'anthro.nc'),
               options=['-u year'])
    print('created anthro.nc')
    print('removing individual monthly files')
    for this_file in all_monthly_files:
        os.remove(this_file)
    for this_file in all_annual_files:
        os.remove(this_file)

if __name__ == "__main__":
    """
    """
    anthro_dir = os.path.join('/', 'global', 'cscratch1', 'sd',
                              'jstineci', 'seasonal_amplitude',
                              'anthro_flux_only')
    ocean_dir = os.path.join('/', 'global', 'cscratch1', 'sd',
                             'jstineci', 'seasonal_amplitude',
                             'ocean_flux_only', 'out_monthly')
    out_dir = os.path.join('/', 'global', 'cscratch1', 'sd',
                          'twhilton', 'global_sandbox')
    avg_anthro(anthro_dir, out_dir)
    #concat_ocean_OCS(ocean_dir, out_dir)
