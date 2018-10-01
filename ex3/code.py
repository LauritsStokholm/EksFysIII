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
    plt.savefig("prep1")


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
    plt.savefig("prep2")

    return


# Experimental
# # # # # # # # # # # # # # Importing Data # # # # # # # # # # # # # # # # # #
""" Here we import the relevant data, and insert paths, directories in lists"""

# Current working directory
data_path     = os.path.join(os.getcwd(), 'data/')
data_path_cal = os.path.join(os.getcwd(), 'data/calibration/test/')
tree_path     = os.path.join(os.getcwd(), 'data/Treebuilder/')

# List of files in ./Data
data_files     = os.listdir(data_path)
data_files_cal = os.listdir(data_path_cal)
tree_files     = os.listdir(tree_path)

# Filtering files
data     = [x for x in data_files if 'txt' in x]
data_cal = [x for x in data_files_cal if 'txt' in x]
tree_data = [x for x in tree_files if 'txt' in x]

# Paths for chosen datafiles
data_dir     = [os.path.join(data_path, x) for x in data]
data_dir_cal = [os.path.join(data_path_cal, x) for x in data_cal]
data_dir.sort()
data_dir_cal.sort()
d1 = [os.path.join(data_path, 'calibration/Cs137.csv')]
tree_dir = [os.path.join(tree_path, x) for x in tree_data]
tree_dir.sort()


# # # # # # # # # # # # # # Fitting Functions  # # # # # # # # # # # # # # # # 
""" Here we define the functions for which we will fit the parameters of"""
def gaussian(x, A, mu, std):
    """ Gaussian (3 parameters)"""
    return A * np.exp(-(x - mu)**2 / (2 * std**2))

def linear(x, k):
    """ Linear (1 parameters)"""
    return k*x

def Cs137(directory):
    lower, upper = [300, 7000]

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
        df1 = df1.loc[(lower < df1[name2]) & (df1[name2] < upper)]
        df = pd.concat([df, df1], axis=1)
        return df, column_names

def Preperation_Dataframe(directory, TreeBuilder):
    lower_limit, upper_limit = [200, 1000]

    # Preparing dataframe structure
    df = pd.DataFrame()

    # Itterating over all files (saving each column name)
    column_names = []
    for i in range(len(directory)):
        # File name
        file_name = os.path.basename(directory[i])

        # Reading csv
        if TreeBuilder==False:
            name1 = "Time {}".format(file_name)
            name2 = "Energy {}".format(file_name)
            column_names.extend((name1, name2))
            df1 = pd.read_csv(directory[i], delimiter=' ', skiprows=[0,1,2,3,4], usecols=[0,1], names=[name1, name2])

        elif TreeBuilder==True:
            name1 = "t_BGO {}".format(file_name)
            name2 = "c_BGO {}".format(file_name)
            name3 = "t_NaI {}".format(file_name)
            name4 = "c_NaI {}".format(file_name)
            column_names.extend((name1, name2, name3, name4))
            df1 = pd.read_csv(directory[i], delimiter=' ', usecols=[0,1,2,3], names=[name1,name2,name3,name4])

        ## Concatenate
        df = pd.concat([df, df1], axis=1)
        df.dropna(axis=0, inplace=True)

    # We make Dataframes and return file name
    return df, column_names

def cal1():
    """ Calibration work"""
    # Dataframe
    df, filename   = Cs137(d1)
    time, channels = filename

    c_Cs = np.array(df[channels])
    t_Cs = np.array(df[time])

    x = np.arange(min(c_Cs), max(c_Cs), 0.1)

    # Bin width
    n_bins = 100
    dx = np.ptp(x) / n_bins

    pdf = sp.stats.gaussian_kde(c_Cs)

    y = pdf(x) * n_bins * dx

    plt.figure()
    plt.plot(x, y, label='PDF')
    plt.hist(c_Cs, bins=n_bins, normed=False)
    return()

def Calibration():
    # Calling DataFrame factory
    df, column_names = Preperation_Dataframe(data_dir_cal)

    # 000 is BGO and 001 is NaI
    # Assigning names
    t_BGO1, c_BGO1, t_Co, c_Co, t_BGO2, c_BGO2, t_Cs, c_Cs = column_names

    # Assigning values
    Co = np.array([df[t_Co], df[c_Co]])
    Cs = np.array([df[t_Cs], df[c_Cs]])
    BGO = np.array([[df[t_BGO1], df[c_BGO1]], [df[t_BGO2], df[c_BGO2]]])

    # For calibration, we are not time interested

    # Spline (HISTOGRAM)
    x_Co_hist   = np.arange(min(Co[1]),     max(Co[1]),     0.1)
    x_Cs_hist   = np.arange(min(Cs[1]),     max(Cs[1]),     0.1)
    x_BGO1_hist = np.arange(min(BGO[0][1]), max(BGO[0][1]), 0.1)
    x_BGO2_hist = np.arange(min(BGO[1][1]), max(BGO[1][1]), 0.1)

    # Number of bins
    n_bins_BGO1 = int(max(BGO[0][1]) - min(BGO[0][1]))
    n_bins_BGO2 = int(max(BGO[1][1]) - min(BGO[1][1]))
    n_bins_Co   = int(max(Co[1])     - min(Co[1]))
    n_bins_Cs   = int(max(Cs[1])     - min(Cs[1]))

    # Bin width
    dx_Cs   = np.ptp(x_Cs_hist)   / n_bins_Cs
    dx_Co   = np.ptp(x_Co_hist)   / n_bins_Co
    dx_BGO1 = np.ptp(x_BGO1_hist) / n_bins_BGO1
    dx_BGO2 = np.ptp(x_BGO2_hist) / n_bins_BGO2

    # Logics
    # Lower limits
    Co_filter   = Co[1][np.where(Co[1] > 850)]
    Cs_filter   = Cs[1][np.where(Cs[1] > 500)]
    BGO_filter1 = BGO[0][1][np.where(BGO[0][1] > 500)]
    BGO_filter2 = BGO[1][1][np.where(BGO[1][1] > 500)]

    # Upper limits
    Co_filter   = Co_filter[np.where(Co_filter < 1200)]
    Cs_filter   = Cs_filter[np.where(Cs_filter < 600)]
    BGO_filter1 = BGO_filter1[np.where(BGO_filter1 < 700)]
    BGO_filter2 = BGO_filter2[np.where(BGO_filter2 < 700)]

    # Spline
    x_Co   = np.arange(min(Co_filter),   max(Co_filter),   0.1)
    x_Cs   = np.arange(min(Cs_filter),   max(Cs_filter),   0.1)
    x_BGO1 = np.arange(min(BGO_filter1), max(BGO_filter1), 0.1)
    x_BGO2 = np.arange(min(BGO_filter2), max(BGO_filter2), 0.1)

    # Kernel Distribution Estimate (functions)
    pdf_Co = sp.stats.gaussian_kde(Co_filter)
    pdf_Cs = sp.stats.gaussian_kde(Cs_filter)
    pdf_BGO1 = sp.stats.gaussian_kde(BGO_filter1)
    pdf_BGO2 = sp.stats.gaussian_kde(BGO_filter2)

    # Evaluate functions on interval (Non normalized)
    y_Cs = pdf_Cs(x_Cs) * dx_Cs * np.size(Cs_filter)
    y_Co = pdf_Co(x_Co) * dx_Co * np.size(Co_filter)
    y_BGO1 = pdf_BGO1(x_BGO1) * dx_BGO1 * np.size(BGO_filter1)
    y_BGO2 = pdf_BGO2(x_BGO2) * dx_BGO2 * np.size(BGO_filter2)

    plt.figure()
    plt.title("Cs137 and Co Calibration")
    plt.grid()
    plt.xlabel("Channel Number")
    plt.ylabel("Counts")
    plt.xlim(200, 1300)
    plt.plot(x_Cs, y_Cs, label="PDF Cs")
    plt.hist(Cs[1], bins=n_bins_Cs, normed=False)#, log=True)
    plt.plot(x_Co, y_Co, label="PDF Co")
    plt.hist(Co[1], bins=n_bins_Co, normed=False)#, log=True)
    plt.legend()
    #plt.savefig("test0.jpg")

    fig, axs = plt.subplots(2, 1, sharex=True)
    fig.subplots_adjust(hspace=0)

    # Plot each graphc
    axs[0].plot(x_BGO1, y_BGO1, label="PDF BGO1")
    axs[0].hist(BGO[0][1], bins=n_bins_BGO1, normed=False)
    axs[0].set_xlim(300, 800)

    axs[1].plot(x_BGO2, y_BGO2, label="PDF BGO2")
    axs[1].hist(BGO[1][1], bins=n_bins_BGO2, normed=False)

    axs[1].set_xlim(300, 800)

#    plt.savefig("test1.jpg")

    # Maximal values and corresponding bins numbers
    Co_peak_index   = sp.signal.argrelmax(y_Co,   axis=0, order=250)
    Cs_peak_index   = sp.signal.argrelmax(y_Cs,   axis=0, order=700)
    BGO1_peak_index = sp.signal.argrelmax(y_BGO1, axis=0, order=400)
    BGO2_peak_index = sp.signal.argrelmax(y_BGO2, axis=0, order=400)

    channels_NaI = np.hstack([0, x_Cs[Cs_peak_index], x_Co[Co_peak_index]])
    channels_BGO = np.hstack([0, x_BGO1[BGO1_peak_index], x_BGO2[BGO2_peak_index]])
    energies_NaI = np.array([0, 661.661, 1183.238, 1332.513]) # keV
    energies_BGO = np.array([0, 661.661, 661.661])

    print("The peak channels are")
    print(channels_NaI)
    print(channels_BGO)

    popt_NaI, pcov_NaI = curve_fit(linear, channels_NaI, energies_NaI)
    perr_NaI = np.sqrt(np.diag(pcov_NaI))
    popt_BGO, pcov_BGO = curve_fit(linear, channels_BGO, energies_BGO)
    perr_BGO = np.sqrt(np.diag(pcov_BGO))

    k_NaI = popt_NaI[0]
    k_NaI_error = perr_NaI[0]
    k_BGO = popt_BGO[0]
    k_BGO_error = perr_BGO[0]

    plt.figure()
    plt.grid()
    plt.title("NaI Linear fit")
    plt.xlabel("Channel numbers")
    plt.ylabel("Energy")
    plt.scatter(channels_NaI, energies_NaI)
    plt.plot(channels_NaI, linear(channels_NaI, k_NaI), label="NaI Linear fit")
    plt.xlim(0, 1200)
    plt.savefig("calibration_fit_NaI")

    print("""Linear fit parameters are k={:.2f} with an error of {:.2f}""".format(k_NaI, k_NaI_error))

    plt.figure()
    plt.grid()
    plt.title("BGO Linear fit")
    plt.xlabel("Channel numbers")
    plt.ylabel("Energy")
    plt.scatter(channels_BGO, energies_BGO)
    plt.plot(channels_BGO, linear(channels_BGO, k_BGO), label="BGO Linear fit")
    plt.xlim(0, 1200)
    plt.savefig("calibration_fit_BGO")

    print("""Linear fit parameters are k={:.2f} with an error of {:.2f}""".format(k_BGO, k_BGO_error))


    # Values:
    # Bin numbers (Cs, Co1, Co2)
    # [552.2, 958.2, 1089.7]

    # Linear fit values
    # a = 1.26 +- 0.03
    # b = -30.09 +- 26.31

    return()


def data_analysis():
    Limit_bottomBGO = 100
    Limit_topBGO    = 200
    Limit_bottomNaI = 580
    Limit_topNaI    = 670

    # Calibration Values
    k_BGO=1.12
    k_NaI=1.22

    TreeBuilder = False
    df, filename = Preperation_Dataframe(data_dir, TreeBuilder)
    t_BGO_40, c_BGO_40, t_NaI_40, c_NaI_40 = filename

    np_BGO_40 = np.array([df[t_BGO_40], df[c_BGO_40]])
    np_NaI_40 = np.array([df[t_NaI_40], df[c_NaI_40]])

    # Logic
    filtered1 = np_BGO_40[1][np.where((np_BGO_40[1] <= 1200)  & (np_BGO_40[1] > 10))]
    filtered2 = np_NaI_40[1][np.where((np_NaI_40[1] <= 1200)  & (np_NaI_40[1] > 10))]

    e1 = k_BGO * filtered1
    e2 = k_NaI * filtered2

    filtered_cloud1  = (Limit_bottomBGO < e1) & (e1 < Limit_topBGO)
    filtered_cloud2  = (Limit_bottomNaI < e2) & (e2 < Limit_topNaI)

    rc_e1 = e1[np.where(filtered_cloud1)]
    rc_e2 = e2[np.where(filtered_cloud2)]


#    filtered_cloud = (Limit_bottomNaI < e2) & (Limit_topNaI > e2) & (Limit_bottomBGO < e1) & (Limit_topBGO > e1)
#    rc_e1 = e1(filtered_cloud)
#    rd_e2 = e2(filtered_cloud)

    filtered_time1 = np_BGO_40[0][np.where(filtered_cloud1)]
    filtered_time2 = np_NaI_40[0][np.where(filtered_cloud2)]

    Time_difference = (np_BGO_40[0] - np_NaI_40[0])
    #filtered_Time_difference = filtered_time1 - filtered_time2
    #filtered_time_difference = Time_difference[np.where[(filtered_cloud)]

    plt.figure()
    plt.hist(Time_difference,100)
    plt.xlabel("Time difference (ns)")
    plt.ylabel("Coun number")
    plt.title("Time difference BGO-NaI (unfiltered)")

def TreeBuilder():
    Limit_bottomBGO = 100
    Limit_topBGO    = 200
    Limit_bottomNaI = 580
    Limit_topNaI    = 670

    # Calibration Values
    k_BGO=1.12
    k_NaI=1.22

    TreeBuilder = True
    df2, column_names = Preperation_Dataframe(tree_dir, TreeBuilder)

    # Manual labour BGO
    t_BGO_100, c_BGO_100, t_NaI_100, c_NaI_100 = column_names[0:4]
    t_BGO_120, c_BGO_120, t_NaI_120, c_NaI_120 = column_names[4:8]
    t_BGO_20, c_BGO_20, t_NaI_20, c_NaI_20 = column_names[8:12]
    t_BGO_25, c_BGO_25, t_NaI_25, c_NaI_25 = column_names[12:16]
    t_BGO_40, c_BGO_40, t_NaI_40, c_NaI_40 = column_names[16:20]
    t_BGO_60, c_BGO_60, t_NaI_60, c_NaI_60 = column_names[20:24]
    t_BGO_80, c_BGO_80, t_NaI_80, c_NaI_80 = column_names[24:28]

    np_BGO_100 = np.array([df2[t_BGO_100], df2[c_BGO_100]])
    np_BGO_120 = np.array([df2[t_BGO_120], df2[c_BGO_120]])
    np_BGO_20 = np.array([df2[t_BGO_20], df2[c_BGO_20]])
    np_BGO_25 = np.array([df2[t_BGO_25], df2[c_BGO_25]])
    np_BGO_40 = np.array([df2[t_BGO_40], df2[c_BGO_40]])
    np_BGO_60 = np.array([df2[t_BGO_60], df2[c_BGO_60]])
    np_BGO_80 = np.array([df2[t_BGO_80], df2[c_BGO_80]])

    np_BGO = np.array([np_BGO_20, np_BGO_25, np_BGO_40, np_BGO_60, np_BGO_80, np_BGO_100, np_BGO_120])

    # Manual labour NaI
    t_NaI_100, c_NaI_100, t_NaI_100, c_NaI_100 = column_names[0:4]
    t_NaI_120, c_NaI_120, t_NaI_120, c_NaI_120 = column_names[4:8]
    t_NaI_20, c_NaI_20, t_NaI_20, c_NaI_20 = column_names[8:12]
    t_NaI_25, c_NaI_25, t_NaI_25, c_NaI_25 = column_names[12:16]
    t_NaI_40, c_NaI_40, t_NaI_40, c_NaI_40 = column_names[16:20]
    t_NaI_60, c_NaI_60, t_NaI_60, c_NaI_60 = column_names[20:24]
    t_NaI_80, c_NaI_80, t_NaI_80, c_NaI_80 = column_names[24:28]

    np_NaI_100 = np.array([df2[t_NaI_100], df2[c_NaI_100]])
    np_NaI_120 = np.array([df2[t_NaI_120], df2[c_NaI_120]])
    np_NaI_20 = np.array([df2[t_NaI_20], df2[c_NaI_20]])
    np_NaI_25 = np.array([df2[t_NaI_25], df2[c_NaI_25]])
    np_NaI_40 = np.array([df2[t_NaI_40], df2[c_NaI_40]])
    np_NaI_60 = np.array([df2[t_NaI_60], df2[c_NaI_60]])
    np_NaI_80 = np.array([df2[t_NaI_80], df2[c_NaI_80]])

    np_NaI = np.array([np_NaI_20, np_NaI_25, np_NaI_40, np_NaI_60, np_NaI_80, np_NaI_100, np_NaI_120])

    x = np.arange(7)
   # Logic
    E1 = np.zeros_like(x)
    E2 = np.zeros_like(x)
    print(E1)
    for i in range(len(np_NaI)):
        print(i)
        logic = (np_BGO[i][1] <= 1200) & (np_BGO[i][1] > 10) & (np_NaI[i][1] <= 1200) & (np_NaI[i][1] > 10)
        filtered1 = np_BGO[i][1][np.where(logic)]
        filtered2 = np_NaI[i][1][np.where(logic)]
    
        e1 = k_BGO * filtered1
        e2 = k_NaI * filtered2
        E1.itemset(i, np.average(e1))
        E2.itemset(i, np.average(e2))
    
        filtered_cloud  = (Limit_bottomBGO < e1) & (e1 < Limit_topBGO) & (Limit_bottomNaI < e2) & (e2 < Limit_topNaI)
    
        rc_e1 = e1[np.where(filtered_cloud)]
        rc_e2 = e2[np.where(filtered_cloud)]
        print(rc_e1)
        print(rc_e2)
    
        Time_difference = np_BGO[i][0] - np_NaI[i][0]
        filtered_time_difference = Time_difference[np.where(filtered_cloud)]
    
        plt.figure()
        plt.hist(Time_difference,100)
        plt.xlabel("Time difference (ns)")
        plt.ylabel("Coun number")
        plt.title("Time difference BGO-NaI (unfiltered)")
    
        plt.figure()
        plt.hist2d(e2, e1, bins=100)
        plt.xlabel("NaI Energy values")
        plt.ylabel("BGO Energy values")
        plt.title("Energy vs. Energy plot")
        plt.savefig("cloud_{}".format(i))
    
        plt.figure()
        plt.hist(filtered_time_difference, 100)
        plt.xlabel("Time difference (ns")
        plt.ylabel("Count number")
        plt.title("Time difference BGO-NaI, filtered")

    # Constants
    # # Electron mass
    m_e = 0.511 # MeV (mc2)
    # # Electron radius (classical)
    r0 = 2.818 # fm
    # # Fine structure
    alpha = 1/137
    # # Energy
    E_Cs137 = 662 * 10**(-3) # MeV

    print(E1)
    print(E2)

    # Scattering angle vs energy
    theta_rad = np.arange(0, np.pi, 0.01)
    theta_deg = theta_rad * (180/np.pi)

    angle = np.array([20, 25, 40, 60, 80, 100, 120])


    ## Compton scattering formula
    ## # E1 is gamma energy before scattering
    ## # m is the electron mass  (mc2)
    ## # theta is scattering angle

    E_2 = lambda E_1, theta: E_1 / (1 + (E_1 / m_e) * (1-np.cos(theta)))

    plt.figure()
    plt.title("Preperation")
    plt.xlabel("Scattering angle")
    plt.ylabel(r"Scattering energy \ [$\si{\mega\electronvolt}$]")
    plt.grid()
    plt.plot(theta_deg, E_2(E_Cs137, theta_rad), label="Theoretical energy")
    plt.plot(angle, E2*10**(-3), label="Data")
    plt.legend()
    plt.savefig("energy")

        
    return()
    


# Når filtreret, skal det valgte område vises i et histogram, så vi kan se, at data er normalfordelt.

#prep()
#cal1()
TreeBuilder()
#Calibration()
#data_analysis()
plt.show()
