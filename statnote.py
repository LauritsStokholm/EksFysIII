# # # # # # # # # # # # # # # Loading Libraries # # # # # # # # # # # # # # # #
# Usual mathematics
import os
import pandas as pd
import numpy as np
import scipy as sp
from scipy.optimize import curve_fit
from scipy import optimize as opt
from scipy import stats
from scipy import signal
import matplotlib.pyplot as plt
plt.style.use('seaborn-deep')


# MatploLib running TeX
params = {'legend.fontsize'     : '20',
          'axes.labelsize'      : '20',
          'axes.titlesize'      : '20',
          'xtick.labelsize'     : '20',
          'ytick.labelsize'     : '20',
          'legend.numpoints'    : 1,
          'text.latex.preamble' : [r'\usepackage{siunitx}',
                                   r'\usepackage{amsmath}'],
          'axes.spines.right'   : False,
          'axes.spines.top'     : False,
          'figure.figsize'      : [8.5, 6.375],
          'legend.frameon'      : False
          }

plt.rcParams.update(params)
plt.rc('text',usetex =True)
#plt.rc('font', **{'family' : "sans-serif"})


# # # # # # # # # # # # # # # Starting code # # # # # # # # # # # # # # # # # #
# # Experiment 1: Gamma spectroscopy
def prep1():
    return
# # Experiment 2: X-rays
def prep2():
    return
# # Experiment 3: Compton Scattering
def prep3():
    return

# # Experiment 4: Rutherford Scattering
def prep4():
    # Conversion factors
    sr = (np.pi/180)**2 # Radian to Steradian

# # # Excercise A)
# # # #  Calculate the differential cross section for scattering 400 keV
# # # #  protons on gold and carbon

    # Defining variables 
    # Protonic number (#protons/#atoms)
    Z_Au = 79
    Z_C  = 6
    Z_p  = 1

    # Accelerated Energy (MeV)
    E = 400 * 10**(-3)

    # Conversion factor (radian to degrees)
    alpha_radToDeg = 360/(2*np.pi)

    # Domain of Angles (Theoretical plot for all angles)
    theta_rad = np.arange(0.01, np.pi, 0.01)
    theta_deg = theta_rad * alpha_radToDeg

    # Define cross section as y-values for plot (mb/sr)
    y_Au = RutherfordcrossSection(Z_Au, Z_p, E, theta_rad)
    y_C  = RutherfordcrossSection(Z_C, Z_p, E, theta_rad)

    # Plotting the distribution
    plt.figure()
    plt.semilogy(theta_deg, y_Au , label="Au")
    plt.semilogy(theta_deg, y_C, label="C")
    plt.title('Theoretical plot')
    plt.xlabel(r'Scattering angle $[\si{\degree}]$')
    plt.ylabel(r'Cross section $[\si{\milli\barn\per\steradian}]$')
    plt.grid()
    plt.xticks(np.arange(0, 190, 15))
    plt.legend()
    plt.savefig('Preperation_graph')

    plt.show()

# # # Excercise B)
# # # # Use this cross section to predict the count rate in a detector
# # # # of solid angle (domega) 1% of 4pi positioned at 160 degrees
# # # # when shooting a 1nA beam of protons on a target of 200(Angstrom) gold

    # Defining variables
    # Solid angle: 1percent of 4pi
    domega = (4*np.pi/100)  # steradians

    # Theta (positioned at 160 degrees)
    theta_0 = (8/9) * np.pi # radians

    # Conversion factor
    alpha_coulombToprotons = 6.24150913 * 10**(18) # elementary charges / C

    # Beam current (particles per sec)
    current_incomming = 1 * 10**(-9) # (C/s)
    N_incident = current_incomming * alpha_coulombToprotons # (#protons/S)

    # #  Thickness of target
    dx = 200 * 10**(-10) # (m)

    # Calculating target density (n)
    # n = (rho * N_a * A)/M

    # Density of target (g/cm3)
    rho_Au = 19.30
    rho_C  = 3.539

    # Atomic mass (g/mol)
    M_Au = 197.96655
    M_C  = 12.0107

    # Atomic number (#nucleons/#particles)


    # Target density (kg/m3)
    n_Au = atomdensity(rho_Au, M_Au, Z_Au)
    n_C  = atomdensity(rho_C,  M_C,  Z_C)

    # Cross Section (evaluated at 160 degrees) (mb/sr)
    # mb to m2 (barn = 10**(-28)m2) ==> mb = 10**(-3)b = 10**(-31)m2
    crossS_Au = RutherfordcrossSection(Z_p, Z_Au, E, theta_0)  * 10**(-31)
    crossS_C  = RutherfordcrossSection(Z_p, Z_C,  E, theta_0)  * 10**(-31)

    # Countrates
    N_scattering_Au = countrate(N_incident, n_Au, dx, crossS_Au, domega)
    N_scattering_C  = countrate(N_incident, n_C,  dx, crossS_C,  domega)

    print("The rate for Au, and C is respectively:")
    print(N_scattering_Au)
    print(N_scattering_C)
    print("That Au is bigger than C is likely, as Au is heavier")


    # Excercise C)
    # Calculate the energy of the scattered proton as a function of scattering
    # angle


    return()

# Numerically value of Cross Section (mb/sr)
def RutherfordcrossSection(Z1, Z2, E, theta):
    dsdo = (1.296*((Z1 * Z2)/(E))**(2)) *(1/(np.sin(theta/2))**4)
    return dsdo

# Returns the countrate by detector (domega) at angle (in crossS) after
# Rutherford scattering of N_incidents with n_target (density) of thickness dx
# (#protons #nucleon / s  ) 
def countrate(N_incident, n_target, dx,  crossS, domega):
    # N_incident = current of proton (#protons/s)
    # n_target = density of target   (#nucleons/m3)
    # dx = target thickness          (m)
    # CrossS = crossSection          (m2/sr)
    # domega = Solid angle.          (sr)
    dN = N_incident * n_target * dx * domega * crossS
    return dN


# Returns the number of Nucleons per meters cubed
def atomdensity(rho, M, A):
    # Defining variables:
    # # rho = density  (g/cm3)
    # # N_a = Avogadros number
    # # M   = molecular weight (g/mol)
    # # A   = Atomic number ( Z + N ) (#nucleons/#particles)
    N_a = 6.022 * 10**(23) # (#particles/mol)
    N = rho*N_a*A/M # (#nucleons/cm3)

    # Unit conversion (#nucleus / cm3) = (#nucleus / m3)
    # g/cm3 = 10**(-3) kg / (10**(-2))**3 m3 = 10**(3) kg/m3
    N = N * 10**(3)
    return N


# # # # # # # # # # # # # # # Output # # # # # # # # # # # # # # # # # # # # #

#prep1()
#prep2()
#prep3()
prep4()

# # # # # # # # # # # # # # # Skrald # # # # # # # # # # # # # # # # # # # # #

# Conversion factors
#st = (np.pi/180)**2 # Radian to Steradian
#br = 10**(28)       # m^2 to barn
#mbr = br * 10**(-3) # millibarn (rightfull unit)
