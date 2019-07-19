import numpy as np
import matplotlib.pyplot as plt
import gradients
from noaa_sites import gather_sites_data

TO_PPT = 1e12
SFC_Z = 0  # Z level for the surface in [OCS] arrays

def delta_to_R(delta, R_ref=0.0422):
    """convert a isotope delta value to an abundance ratio (R)
    """
    R = R_ref * ((delta / 1000.0) + 1.0)
    return(R)

def epsilon_to_uptake_heavy(epsilon, R, U_light):
    """calculate uptake flux for heavy isotope

    calculate uptake flux for heavy isotope from isotopic
    fractionation epsilon, abundance ratio, uptake flux for light
    isotope
    """
    U_heavy = U_light * (1.0 + (epsilon_U / 1000.0)) * R


""" ocean adjustment factor

        %             anthro ocean_COS ocean_CS2 ocean_DMS ocean_missing ocean_post ocean_total
		2004  246.052   155.474   155.474   155.474       5.60375    2460.63     1085.17
		2005  246.052   154.926   154.926   154.926       5.58844    2454.14        1082

- as per 7 May 2018 email from Jim, these fluxes need to be scaled:

		Check the weighting for the components -- should be 0.65x for
        the Kettle OCS / CS2 / DMS, 75x for LUKAI, and 0.147x for
        posterior_12month_kfl. I come up with 967 after applying that
        scaling.

- applying the scaling, the total is

		0.65*(155.474 + 155.474 + 155.474) + (75 * 5.60375) + (0.147 * 2460.63)
		= 1085

- so, to estimate OCS concentrations for various ocean flux estimates use these scaling factors for the concentrations:
  - Lennartz: (345 / 1085) = 0.31797235023041476
  - Launois best guess: (813 / 1085) = 0.7493087557603687
  - Launois high estimate: (3997 / 1085) = 3.6838709677419357
"""
# [OCS] scaling factors for low, medium, high ocean OCS fluxes
ocean_flux = {'Lennartz': (345.0 / 1085.0),
              'Launois_best': (813.0 / 1085.0),
              'Launois_high': (3997.0 / 1085.0)}

# uncertainty in the detection equipment
instrument_uncert = 0.5

# d34S fractionation values assumed in text of BSF proposal (pp 7-8)
fractionation = {'ocean': 19,
                 'anthro': 3,
                 'plant': -5,
                 'soil': -3,
                 'oxidation': 8}
# d34S fractionation uncertainty ranges assumed in text of BSF
# proposal (pp 7-8)
fractionation_uncert  = {'ocean': 0,
                         'anthro': 2,
                         'plant': 0,
                         'soil': 1,
                         'oxidation': 0}

def get_site_xy(code):
    """return site info for one specified site
    """
    sites = gather_sites_data()
    # TODO: X and Y are switched?
    return(sites[sites.Code == code].Y.values,
           sites[sites.Code == code].X.values)


def calc_34S_concentration(ocs, permil_34S):
    """calculate concentration of 34S for a specified fractionation

    ocs (array-like): OCS concentration (ppt)
    permil_34S (float): 34S fractionation
    """
    return(ocs * (permil_34S / 1000.0))

if __name__ == "__main__":
    ocs_anthro = gradients.get_anthro_anomaly()
    ocs_ocean = gradients.get_ocean_anomaly()

    all_sites = gather_sites_data()
    pacific_gradient = list(all_sites[(all_sites.Longitude > 150.0) |
                              (all_sites.Longitude < -150.0)].Code)
    atlantic_gradient = ['ZEP', 'ICE', 'MHD', 'AZR', 'IZO',
                    'ASC', 'HBA', 'NMB', 'CPT']
    indian_gradient = ['SEY', 'CRZ', 'SYO']

    plt.figure()
    for k, this_ocean_flux in ocean_flux.items():
        for this_site in pacific_gradient:  # ('CRZ', ):
            this_x, this_y = get_site_xy(this_site)
            # show S hemisphere sites with dashed lines, N hemisphere with solid
            if (all_sites[all_sites.Code == this_site].Latitude < 0.0).all():
                linestyle='dashed'
            else:
                linestyle='solid'
            ocean_d34S = calc_34S_concentration(
                ocs=ocs_ocean.data[:, SFC_Z, this_x, this_y],
                permil_34S=fractionation['ocean']) * TO_PPT
            anthro_d34S = calc_34S_concentration(
                ocs=ocs_anthro.data[:, SFC_Z, this_x, this_y],
                permil_34S=fractionation['anthro']) * TO_PPT
            hi_ocean_flux_d34S = ocean_d34S * ocean_flux['Launois_high']
            med_ocean_flux_d34S = ocean_d34S * ocean_flux['Launois_best']
            low_ocean_flux_d34S = ocean_d34S * ocean_flux['Lennartz']
            plt.errorbar(np.arange(med_ocean_flux_d34S.size),
                         (med_ocean_flux_d34S / anthro_d34S).squeeze(),
                         yerr=np.array([low_ocean_flux_d34S / anthro_d34S,
                                        hi_ocean_flux_d34S / anthro_d34S]).squeeze(),
                         linestyle=linestyle,
                         label=this_site)
        plt.gca().set_title(r'ratio ocean:anthro $\delta$34S OCS anomaly (from global mean)')
        plt.gca().set_xlabel('month')
        plt.gca().set_ylabel(r'$\delta$34S [OCS] (ppt)')
        plt.legend()
