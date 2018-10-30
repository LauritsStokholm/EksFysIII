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


# Input function
def fileIO(filetype, fileid):
    """ Treebuilding """
    lstdirs = [] 
    for (dirs, paths, files) in os.walk("."):
        for file in files:
            if (file.endswith(filetype)) & (all(id in file for id in fileid)):
                lstdirs.append(dirs+"/"+file)
    return lstdirs


# Relevant data directories listed
lstdirs_Al = fileIO("txt", "Al")
lstdirs_Br = fileIO("txt", "Brass")
lstdirs_Al.sort()
lstdirs_Br.sort()
print(lstdirs_Al)
print(lstdirs_Br)


# Preparing dataframe
def Preperation_Dataframe(directory):

    # Preparing dataframe structure
    df = pd.DataFrame()

    # Itterating over all files (saving each column name)
    column_names = []
    for i in range(len(directory)):
        # File name
        file_name = os.path.basename(directory[i])

        # Reading csv
        #name1 = "Time {}".format(file_name)
        name2 = "Channel {}".format(file_name)
        column_names.extend([name2])
        df1 = pd.read_csv(directory[i], delimiter=' ', skiprows=[0,1,2,3,4], usecols=[1], names=[name2])

        ## Concatenate
        df = pd.concat([df, df1], axis=1)
        df.dropna(axis=0, inplace=True)

    # We make Dataframes and return file name
    return df, column_names

df_Al, column_names_Al = Preperation_Dataframe(lstdirs_Al)

print(column_names_Al)

np_Al = np.array([np.array(df_Al.index), np.array(df_Al[column_names_Al[0]])])
plt.figure()
plt.plot(np_Al[0], np_Al[1])



#np_Al2 = np.array([np.array(df_Al.index), np.array(df_Al[column_names_Al[1]])])
#plt.figure()
#plt.plot(np_Al2[0], np_Al2[1])

plt.show()
#df_Br, column_names_Br = Preperation_Dataframe(lstdirs_Br)



#def Calibration():
##    # Calling DataFrame factory
##    df, column_names = Preperation_Dataframe(data_dir_cal)
##
##    # 000 is BGO and 001 is NaI
##    # Assigning names
##    t_BGO1, c_BGO1, t_Co, c_Co, t_BGO2, c_BGO2, t_Cs, c_Cs = column_names
##
##    # Assigning values
##    Co = np.array([df[t_Co], df[c_Co]])
##    Cs = np.array([df[t_Cs], df[c_Cs]])
##    BGO = np.array([[df[t_BGO1], df[c_BGO1]], [df[t_BGO2], df[c_BGO2]]])
##
##    # For calibration, we are not time interested
##
##    # Spline (HISTOGRAM)
##    x_Co_hist   = np.arange(min(Co[1]),     max(Co[1]),     0.1)
##    x_Cs_hist   = np.arange(min(Cs[1]),     max(Cs[1]),     0.1)
##    x_BGO1_hist = np.arange(min(BGO[0][1]), max(BGO[0][1]), 0.1)
##    x_BGO2_hist = np.arange(min(BGO[1][1]), max(BGO[1][1]), 0.1)
##
##    # Number of bins
##    n_bins_BGO1 = int(max(BGO[0][1]) - min(BGO[0][1]))
##    n_bins_BGO2 = int(max(BGO[1][1]) - min(BGO[1][1]))
##    n_bins_Co   = int(max(Co[1])     - min(Co[1]))
##    n_bins_Cs   = int(max(Cs[1])     - min(Cs[1]))
##
##    # Bin width
##    dx_Cs   = np.ptp(x_Cs_hist)   / n_bins_Cs
##    dx_Co   = np.ptp(x_Co_hist)   / n_bins_Co
##    dx_BGO1 = np.ptp(x_BGO1_hist) / n_bins_BGO1
##    dx_BGO2 = np.ptp(x_BGO2_hist) / n_bins_BGO2
##
##    # Logics
##    # Lower limits
##    Co_filter   = Co[1][np.where(Co[1] > 850)]
##    Cs_filter   = Cs[1][np.where(Cs[1] > 500)]
##    BGO_filter1 = BGO[0][1][np.where(BGO[0][1] > 500)]
##    BGO_filter2 = BGO[1][1][np.where(BGO[1][1] > 500)]
##
##    # Upper limits
##    Co_filter   = Co_filter[np.where(Co_filter < 1200)]
##    Cs_filter   = Cs_filter[np.where(Cs_filter < 600)]
##    BGO_filter1 = BGO_filter1[np.where(BGO_filter1 < 700)]
##    BGO_filter2 = BGO_filter2[np.where(BGO_filter2 < 700)]
##
##    # Spline
##    x_Co   = np.arange(min(Co_filter),   max(Co_filter),   0.1)
##    x_Cs   = np.arange(min(Cs_filter),   max(Cs_filter),   0.1)
##    x_BGO1 = np.arange(min(BGO_filter1), max(BGO_filter1), 0.1)
##    x_BGO2 = np.arange(min(BGO_filter2), max(BGO_filter2), 0.1)
##
##    # Kernel Distribution Estimate (functions)
##    pdf_Co = sp.stats.gaussian_kde(Co_filter)
##    pdf_Cs = sp.stats.gaussian_kde(Cs_filter)
##    pdf_BGO1 = sp.stats.gaussian_kde(BGO_filter1)
##    pdf_BGO2 = sp.stats.gaussian_kde(BGO_filter2)
##
##    # Evaluate functions on interval (Non normalized)
##    y_Cs = pdf_Cs(x_Cs) * dx_Cs * np.size(Cs_filter)
##    y_Co = pdf_Co(x_Co) * dx_Co * np.size(Co_filter)
##    y_BGO1 = pdf_BGO1(x_BGO1) * dx_BGO1 * np.size(BGO_filter1)
##    y_BGO2 = pdf_BGO2(x_BGO2) * dx_BGO2 * np.size(BGO_filter2)
##
##    plt.figure()
##    plt.title("Cs137 and Co Calibration")
##    plt.grid()
##    plt.xlabel("Channel Number")
##    plt.ylabel("Counts")
##    plt.xlim(200, 1300)
##    plt.plot(x_Cs, y_Cs, label="PDF Cs")
##    plt.hist(Cs[1], bins=n_bins_Cs, normed=False)#, log=True)
##    plt.plot(x_Co, y_Co, label="PDF Co")
##    plt.hist(Co[1], bins=n_bins_Co, normed=False)#, log=True)
##    plt.legend()
##    #plt.savefig("test0.jpg")
##
##    fig, axs = plt.subplots(2, 1, sharex=True)
##    fig.subplots_adjust(hspace=0)
##
##    # Plot each graphc
##    axs[0].plot(x_BGO1, y_BGO1, label="PDF BGO1")
##    axs[0].hist(BGO[0][1], bins=n_bins_BGO1, normed=False)
##    axs[0].set_xlim(300, 800)
##
##    axs[1].plot(x_BGO2, y_BGO2, label="PDF BGO2")
##    axs[1].hist(BGO[1][1], bins=n_bins_BGO2, normed=False)
##
##    axs[1].set_xlim(300, 800)
##
###    plt.savefig("test1.jpg")
##
##    # Maximal values and corresponding bins numbers
##    Co_peak_index   = sp.signal.argrelmax(y_Co,   axis=0, order=250)
##    Cs_peak_index   = sp.signal.argrelmax(y_Cs,   axis=0, order=700)
##    BGO1_peak_index = sp.signal.argrelmax(y_BGO1, axis=0, order=400)
##    BGO2_peak_index = sp.signal.argrelmax(y_BGO2, axis=0, order=400)
##
##    channels_NaI = np.hstack([0, x_Cs[Cs_peak_index], x_Co[Co_peak_index]])
##    channels_BGO = np.hstack([0, x_BGO1[BGO1_peak_index], x_BGO2[BGO2_peak_index]])
##    energies_NaI = np.array([0, 661.661, 1183.238, 1332.513]) # keV
##    energies_BGO = np.array([0, 661.661, 661.661])
##
##    print("The peak channels are")
##    print(channels_NaI)
##    print(channels_BGO)
##
##    popt_NaI, pcov_NaI = curve_fit(linear, channels_NaI, energies_NaI)
##    perr_NaI = np.sqrt(np.diag(pcov_NaI))
##    popt_BGO, pcov_BGO = curve_fit(linear, channels_BGO, energies_BGO)
##    perr_BGO = np.sqrt(np.diag(pcov_BGO))
##
##    k_NaI = popt_NaI[0]
##    k_NaI_error = perr_NaI[0]
##    k_BGO = popt_BGO[0]
##    k_BGO_error = perr_BGO[0]
##
##    plt.figure()
##    plt.grid()
##    plt.title("NaI Linear fit")
##    plt.xlabel("Channel numbers")
##    plt.ylabel("Energy")
##    plt.scatter(channels_NaI, energies_NaI)
##    plt.plot(channels_NaI, linear(channels_NaI, k_NaI), label="NaI Linear fit")
##    plt.xlim(0, 1200)
##    plt.savefig("calibration_fit_NaI")
##
##    print("""Linear fit parameters are k={:.2f} with an error of {:.2f}""".format(k_NaI, k_NaI_error))
##
##    plt.figure()
##    plt.grid()
##    plt.title("BGO Linear fit")
##    plt.xlabel("Channel numbers")
##    plt.ylabel("Energy")
##    plt.scatter(channels_BGO, energies_BGO)
##    plt.plot(channels_BGO, linear(channels_BGO, k_BGO), label="BGO Linear fit")
##    plt.xlim(0, 1200)
##    plt.savefig("calibration_fit_BGO")
##
##    print("""Linear fit parameters are k={:.2f} with an error of {:.2f}""".format(k_BGO, k_BGO_error))
##
##
##
##
##
##Aldirs = fileIO(".txt", ("Al", "ch001"))
##Brdirs = fileIO(".txt", ("Br", "ch001"))





