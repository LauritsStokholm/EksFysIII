# # # # # # # # # # # # # # Loading Libraries # # # # # # # # # # # # # # # # 
# Used for basic commands
import os
import glob

# Data Analysis and basic mathematics
import pandas as pd
import numpy as np
import scipy as sp
from scipy.stats import norm
from scipy.optimize import curve_fit
from scipy import optimize as opt
from scipy import stats 
from scipy import signal
import matplotlib.pyplot as plt
plt.style.use('seaborn-deep')

# Matplotlib compatible TeX
params = {'legend.fontsize'     : '20',
          'axes.labelsize'      : '20',
          'axes.titlesize'      : '20',
          'xtick.labelsize'     : '20',
          'ytick.labelsize'     : '20',
          'legend.numpoints'    : 1,
          'text.latex.preamble' : [r'\usepackage{siunitx}',
                                   r'\usepackage{amsmath}'],
          #'axes.spines.right'   : False,
          #'axes.spines.top'     : False,
          'figure.figsize'      : [8.5, 6.375],
          'legend.frameon'      : False
          }

plt.rcParams.update(params)
plt.rc('text',usetex =True)
plt.rc('font', **{'family' : "sans-serif"})

# This is for awesome zoom
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset


# Preperations opgaverne
def prep():
    # Find the expectation relation between the scattering angle and the energy of
    # the scattered photon for Compton scattering of photons from 137Cs

    # Constants
    # # Electron mass
    m_e = 0.511 # MeV (mc2)
    # # Electron radius (classical)
    r0 = 2.818 # fm
    # # Fine structure
    alpha = 1/137
    # # Energy
    E_Cs137 = 662 * 10**(-3) # MeV


    # Scattering angle vs energy
    theta_rad = np.arange(0, np.pi, 0.01)
    theta_deg = theta_rad * (180/np.pi)

    # Compton scattering formula
    # # E1 is gamma energy before scattering
    # # m is the electron mass  (mc2)
    # # theta is scattering angle

    E_2 = lambda E_1, theta: E_1 / (1 + (E_1 / m_e) * (1-np.cos(theta)))

    plt.figure()
    plt.title("Preperation")
    plt.xlabel("Scattering angle")
    plt.ylabel(r"Scattering energy \ [$\si{\mega\electronvolt}$]")
    plt.grid()
    plt.plot(theta_deg, E_2(E_Cs137, theta_rad), label="Theoretical energy")
    plt.legend()


    # Calculate the variation of the cross section as a function of the
    # scattering angle for Compton scattering of photons from Cs137.

    dsdo = lambda theta: \
             (r0**2) * (1 / (1 + alpha*(1 - np.cos(theta) ) ) )**3 \
             * ( (1 + np.cos(theta) ) / 2) \
             * ( 1 + alpha**2 * ( 1 - np.cos(theta) )**2 / ( (1+
                 (np.cos(theta))**2) * (1 + alpha*(1 - np.cos(theta)))))

    plt.figure()
    plt.plot(theta_deg, dsdo(theta_rad), label="Theoretical compton cross section")
    plt.grid()
    plt.xlabel("Scattering angle")
    plt.ylabel("Cross section")
    plt.legend()

    return


# Experimental
# # # # # # # # # # # # # # Importing Data # # # # # # # # # # # # # # # # # #
""" Here we import the relevant data, and insert paths, directories in lists"""


