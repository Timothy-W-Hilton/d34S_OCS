"""TWH port of Alon's Matlab script demo to python"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

n_years = 20
n_months = 12

R_ref = 0.0422
R_atm = R_ref * (3 / 1000 + 1)  #lets start with atmospheric d34S of 3 permil, and later do a spin up.
COS32 = 500 # This is the amount in the atmosphere. for simplicity I
             # work here with concetrations in ppt instead of mass
COS34 = COS32*R_atm  # This is the amount of COS with 34S
epsilon_U = -5 #assuming fractionation of -5 permil in uptake by plants
#R_P = R_ref*(20/1000+1) # Assuming 20 permil for production from the ocean
#R_P = R_ref*(3/1000+1) # Assuming 3 permil for anthropogenic production

# Assuming 8 permil for total production from both sources.  TWH: This line
# rearranges the definition of isotope fraction (delta) to solve for
# R_P.
R_P=R_ref*(8/1000+1)

#assuming arbitrary production for each month (TWH: production is
#oceans + anthro, I think)
P_COS32 = np.array([50, 110, 130, 140, 150, 200, 250, 300, 250, 200, 150, 50])
#assuming arbitrary uptake for each month (TWH: uptake is from plants)
U_COS32 = np.array([25, 50, 60, 70, 100, 300, 500, 400, 200, 150, 100, 25])
Dt = 0.1  # time step

P_COS34 = P_COS32 * R_P
d34S = np.zeros((n_years, n_months))
COS = np.zeros((n_years, n_months))
for y in range(n_years): # years loop
    for m in range(n_months): # months loop
        R_atm = COS34 / COS32
        U_COS34 = U_COS32[m] * (1 + epsilon_U / 1000) * R_atm
        COS32 = COS32 + P_COS32[m] * Dt - U_COS32[m] * Dt  # COS 32 iteration
        COS34 = COS34 + P_COS34[m] * Dt - U_COS34 * Dt  # COS 34 iteration
        d34S[y, m] = ((COS34 / COS32) / R_ref - 1) * 1000  #calculate delta-34S
        COS[y, m] = COS32  # save the COS concetration for later


with PdfPages('multipage_pdf.pdf') as pdf:
    for y in range(n_years): # years loop
        fig, ax = plt.subplots(2, 1)
        ax[0].plot(range(1, n_months + 1), COS[y, :])
        ax[0].set_ylabel('[COS]')
        ax[0].set_ylim(COS.min(), COS.max())
        ax[1].plot(range(1, n_months + 1), d34S[y, :])
        ax[1].set_ylabel('d34S')
        ax[1].set_ylim(d34S.min(), d34S.max())
        fig.suptitle('year {:02d}'.format(y + 1))
        pdf.savefig()
        plt.close()
