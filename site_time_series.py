import matplotlib.pyplot as plt
import gradients
from noaa_sites import gather_sites_data

TO_PPT = 1e12

def get_site_xy(code):
    """return site info for one specified site
    """
    sites = gather_sites_data()
    # TODO: X and Y are switched?
    return(sites[sites.Code == code].Y.values,
           sites[sites.Code == code].X.values)



if __name__ == "__main__":
    ocs_anthro = gradients.get_anthro_anomaly()
    ocs_ocean = gradients.get_ocean_anomaly()

    all_sites = gather_sites_data()
    pacific_gradient = list(all_sites[(all_sites.Longitude > 150.0) |
                              (all_sites.Longitude < -150.0)].Code)
    atlantic_gradient = ['ZEP', 'ICE', 'MHD', 'AZR', 'IZO',
                    'ASC', 'HBA', 'NMB', 'CPT']
    indian_gradient = ['SEY', 'CRZ', 'SYO']

    for k, this_data in {'ocean': ocs_ocean, 'anthro':ocs_anthro}.items():
        plt.figure()
        for this_site in indian_gradient:
            this_x, this_y = get_site_xy(this_site)
            # show S hemisphere sites with dashed lines, N hemisphere with solid
            if (all_sites[all_sites.Code == this_site].Latitude < 0.0).all():
                linestyle='dashed'
            else:
                linestyle='solid'
            plt.plot(this_data.data[:, 0, this_x, this_y] * TO_PPT,
                     linestyle=linestyle,
                     label=this_site)
        plt.gca().set_title('{} OCS anomaly (from global mean)'.format(k))
        plt.gca().set_xlabel('month')
        plt.gca().set_ylabel('[OCS] (ppt)')
        plt.legend()
