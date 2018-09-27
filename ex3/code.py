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
# This is for awesome zoom
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
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



# Preperations opgaverne
def prep():
    """Preperation tasks"""
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

# Current working directory
data_path     = os.path.join(os.getcwd(), 'data/')
data_path_cal = os.path.join(os.getcwd(), 'data/calibration/test/')

# List of files in ./Data
data_files     = os.listdir(data_path)
data_files_cal = os.listdir(data_path_cal)

# Filtering files
data     = [x for x in data_files if 'txt' in x]
data_cal = [x for x in data_files_cal if 'txt' in x]

# Paths for chosen datafiles
data_dir     = [os.path.join(data_path, x) for x in data]
data_dir_cal = [os.path.join(data_path_cal, x) for x in data_cal]
data_dir_cal.sort()
d1 = [os.path.join(data_path, 'calibration/Cs137.csv')]


# # # # # # # # # # # # # # Fitting Functions  # # # # # # # # # # # # # # # # 
""" Here we define the functions for which we will fit the parameters of"""
def gaussian(x, A, mu, std):
    """ Gaussian (3 parameters)"""
    return A * np.exp(-(x - mu)**2 / (2 * std**2))

def linear(x, a, b):
    """ Linear (2 parameters)"""
    return a*x + b

# # # # # # # # # # # # # # DataFrame Import # # # # # # # # # # # # # # # # # 
""" Here we import given list of directories into dataframes """
def Cs137(directory, threshold):
    lower = threshold[0]
    upper = threshold[1]

    df = pd.DataFrame()
    column_names = []
    for i in range(len(directory)):
        file_name = os.path.basename(directory[i])

        # Reading csv
        name1 = "Time"
        name2 = "Energy"

        column_names.append(name1)
        column_names.append(name2)
        df1 = pd.read_csv(directory[i], delimiter=',', names=[name1, name2])
        df1 = df1.loc[lower < df1[name2]]
        df1 = df1.loc[df1[name2] < upper]
        df = pd.concat([df, df1], axis=1)
        return df, column_names



def Preperation_Dataframe(directory, threshold):
    lower = threshold[0]
    upper = threshold[1]

    # Preparing dataframe structure
    df = pd.DataFrame()

    # Itterating over all files (saving each column name)
    column_names = []
    for i in range(len(directory)):
        # File name
        file_name = os.path.basename(directory[i])

        # Reading csv
        name1 = "Time {}".format(file_name)
        name2 = "Energy {}".format(file_name)

        column_names.append(name1)
        column_names.append(name2)

        df1 = pd.read_csv(directory[i], delimiter=' ', skiprows=[0,1,2,3,4], usecols=[0,1], names=[name1, name2])
        df1 = df1.loc[lower < df1[name2]]
        df1 = df1.loc[df1[name2] < upper]

        ## Concatenate
        df = pd.concat([df, df1], axis=1)
        df.dropna(axis=0, inplace=True)
    # We make Dataframes and return file name
    return df, column_names


# Output
#prep()

def cal1():
    """ Calibration work"""
    # Dataframe
    df , filename = Cs137(d1, [100, 7000])
    time     = filename[0]
    channels = filename[1]

    c_Cs = np.array(df[channels])
    t_Cs = np.array(df[time])

    x = np.arange(min(c_Cs), max(c_Cs), 0.1)

    n_bins = 100
    pdf = sp.stats.gaussian_kde(c_Cs)

    plt.figure()
    plt.plot(x, pdf(x) * np.size(x), label='PDF')
    plt.hist(c_Cs, bins=n_bins, normed=False)
    return()






    #df, file_names = Preperation_Dataframe(data_dir_cal, [100, 7000])
    #t_Co1 = df[file_names[0]]
    #c_Co1 = df[file_names[1]]
    #t_Co2 = df[file_names[6]]
    #c_Co2 = df[file_names[7]]

    #t_Cs1 = df[file_names[2]]
    #c_Cs1 = df[file_names[3]]
    #t_Cs2 = df[file_names[4]]
    #c_Cs2 = df[file_names[5]]

    ## Unpack values
    #pdf_Co1 = sp.stats.gaussian_kde(c_Co1)
    #pdf_Co2 = sp.stats.gaussian_kde(c_Co2)
    #pdf_Cs1 = sp.stats.gaussian_kde(c_Cs1)
    #pdf_Cs2 = sp.stats.gaussian_kde(c_Cs2)

    #print(pdf_Co1(t_Co1)[0])



    return()

def Calibration():
    df, column_names = Preperation_Dataframe(data_dir_cal, [200, 1000])

    # Assigning names
    t_Cs = column_names[0]
    c_Cs = column_names[1]
    t_BGO1 = column_names[2]
    c_BGO1 = column_names[3]

    t_Co = column_names[4]
    c_Co = column_names[5]
    t_BGO2 = column_names[6]
    c_BGO2 = column_names[7]


    Co = np.array([df[t_Co], df[c_Co]])
    Cs = np.array([df[t_Cs], df[c_Cs]])
    BGO = np.array([[df[t_BGO1], df[c_BGO1]], [df[t_BGO2], df[c_BGO2]]])

    # For calibration, we are not time interested
    x_Co = np.arange(min(Co[1]), max(Co[1]), 0.1)
    x_Cs = np.arange(min(Cs[1]), max(Cs[1]), 0.1)
    x_BGO1 = np.arange(min(BGO[0][1]), max(BGO[0][1]), 0.1)
    x_BGO2 = np.arange(min(BGO[1][1]), max(BGO[1][1]), 0.1)

    n_bins = 1000
    dx = np.ptp(x_Cs) / n_bins
    pdf_Co = sp.stats.gaussian_kde(Co[1])
    pdf_Cs = sp.stats.gaussian_kde(Cs[1])
    pdf_BGO1 = sp.stats.gaussian_kde(BGO[0][1])
    pdf_BGO2 = sp.stats.gaussian_kde(BGO[1][1])

    y_Co = pdf_Co(x_Co)
    y_Cs = pdf_Cs(x_Cs)
    y_BGO1 = pdf_BGO1(x_BGO1)
    y_BGO2 = pdf_BGO2(x_BGO2)

    plt.figure()
    plt.title("Cs137 Calibration")
    plt.grid()
    plt.xlabel("Channel Number")
    plt.ylabel("Counts")
    plt.plot(x_Cs, y_Cs * dx * np.size(x_Cs), label="PDF Cs")
    plt.hist(Cs[1], bins=n_bins, normed=False)#, log=True)
    plt.legend()
#    plt.savefig("test0.jpg")

    plt.figure()
    plt.title("Co Calibration")
    plt.grid()
    plt.xlabel("Channel Number")
    plt.ylabel("Counts")
    plt.plot(x_Co, y_Co * dx * np.size(x_Co), label="PDF Co")
    plt.hist(Co[1], bins=n_bins, normed=False)#, log=True)
    plt.legend()
#    plt.savefig("test1.jpg")

    fig, axs = plt.subplots(2, 1, sharex=True)
    fig.subplots_adjust(hspace=0)

    # Plot each graphc
    axs[0].plot(x_BGO1, y_BGO1 * dx, label="PDF BGO1")
    axs[0].hist(BGO[0][1], bins=n_bins, normed=False)

    axs[1].plot(x_BGO2, y_BGO2 * dx, label="PDF BGO2")
    axs[1].hist(BGO[1][1], bins=n_bins, normed=False)




    # Maximal values and corresponding bins numbers
    print(sp.signal.argrelmax(y_Co, axis=0, order=100))
    print(sp.signal.argrelmax(y_Cs, axis=0, order=100))





    return

#cal1()
#Calibration()
#x1 =np.array(df_cal[file_names_cal[2]])
#y1 =np.array(df_cal[file_names_cal[1]])
#x2 sd=np.array(df_cal[file_names_cal[1]])
#y2 =np.array(df_cal[file_names_cal[1]])




#df_cal[file_names_cal[1]].plot()
#plt.plot(df_cal[file_names_cal[0]], df_cal[file_names_cal[1]], '--')

#x = np.array(df_cal[file_names_cal[2]])
#y = np.array(df_cal[file_names_cal[3]])
#plt.plot(x, y)
#
#
Calibration()
plt.show()
