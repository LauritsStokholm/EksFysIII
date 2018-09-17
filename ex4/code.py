# # # # # # # # # # # # # # # Loading Libraries # # # # # # # # # # # # # # # # 
# Used for basic commands
import os

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
          'axes.spines.right'   : False,
          'axes.spines.top'     : False,
          'figure.figsize'      : [8.5, 6.375],
          'legend.frameon'      : False
          }

plt.rcParams.update(params)
plt.rc('text',usetex =True)
plt.rc('font', **{'family' : "sans-serif"})

# # # # # # # # # # # # # # # # Importing Data # # # # # # # # # # # # # # # # 

# Current working directory
data_path = os.path.join(os.getcwd(), 'data/')

# List of files in ./Data
data_files = os.listdir(data_path)

# Filtering files
data_cal = [x for x in data_files if 'cal' in x and 'csv' in x]
data     = [x for x in data_files if 'AuC' in x and 'csv' in x]

# Paths for chosen datafiles
data_cal_dir = [os.path.join(data_path, x) for x in data_cal]
data_dir     = [os.path.join(data_path, x) for x in data]


# # # # # # # # # # # # # # # # Calibration # # # # # # # # # # # # # # # # # # 
""" In this calibration we determine alpha and k0, such we can transform
channel numbers to energy by the E = alpha(k - k0), where k is the channel
number and k0 is the channelnumber for zero amplitude (extrapolation)"""

# Defining functions for fitting later
# Gaussian
def gaussian(x, amplitude, mu, std):
    return amplitude * np.exp(-(x - mu)**2 / (2 * std**2))

# Linear
def linear(a, b, x):
    return a*x + b

# Energy
def ENERGY():
    Ei = 400 # kev
    theta_deg = 160 # degrees
    theta = 160 * (180 / np.pi)
    mAu = 196.9665690
    mH = 1.008

    E = lambda mb, mt: (((mb*np.cos(theta) + np.sqrt(mt**2 - mb**2 *
            (np.sin(theta))**2)) / (mb+mt))**2) * Ei

    E1 = E(mH, mAu)
    E2 = E(2*mH, mAu)

    return E1, E2



# Data imported into dataframes
def Dataframe(directory):
    """ Welcome to the Datafram factory. """
    # Preparing dataframe structure
    df = pd.DataFrame()
    
    # Itterating over all files
    file_names = []
    for i in range(len(directory)):
        # file name
        file_names.append(os.path.basename(directory[i]))
    
        # Reading csv
        column_names = [file_names[i]]
        df1 = pd.read_csv(directory[i], delimiter=',', index_col=0, names=column_names)
    
        # Concatenating
        df = pd.concat([df, df1[column_names[0]]], axis=1)

    # We make Dataframes and return file name
    return df, file_names

def Gaussian_fit(df, file_names, k0_switch):
    """ This is a Gaussian fit factory"""
    # Gaussian mean and std saved for later
    gaussian_mean = []
    gaussian_std  = []
    
    if k0_switch == True:
        for name in file_names:
            # Choosing non zero values
            df_data = df[name].loc[df[name].nonzero()]
        
            # Transforming to array (Might be unnecesarry) [index, count]
            np_data = np.array([np.array(df_data.index), df_data.values])
            x = np_data[0] # Indices
            y = np_data[1] # Count number
            
            # Scipy Optimization after Gaussian (guessing start values)
            mean = sum(x * y) / sum(y)
            sigma = np.sqrt(sum(y * (x - mean)**2) / sum(y))
            popt,pcov = curve_fit(gaussian, x, y, p0=[max(y), mean, sigma])
            
            gaussian_mean.append(popt[1])
            gaussian_std.append(popt[2])
            
            
            # Plotting the data + fit
            plt.title("Calibration")
            plt.xlabel("Channel \#")
            plt.ylabel("Count \#")
            plt.grid()
            
                # Fit (smooth x and fitted parameters)
            gauss_x = np.arange(x[0], x[-1], 0.0001)
            gauss_y = gaussian(gauss_x, *popt)
   
            plt.plot(gauss_x, gauss_y, label='fit')
            plt.plot(x,y, '.', label='data')
            plt.legend(loc=1)

        return np.array(gaussian_mean), np.array(gaussian_std)
    else:
        df_data = df

        # Transforming to array (Might be unnecesarry) [index, count]
        np_data = np.array([np.array(df_data.index), df_data.values])
        x = np_data[0] # Indices
        y = np_data[1] # Count number
        
        # Scipy Optimization after Gaussian (guessing start values)
        mean = sum(x * y) / sum(y)
        sigma = np.sqrt(sum(y * (x - mean)**2) / sum(y))
        popt,pcov = curve_fit(gaussian, x, y, p0=[max(y), mean, sigma])
        
        gaussian_mean.append(popt[1])
        gaussian_std.append(popt[2])
        
        
        # Plotting the data + fit
        plt.title("Calibration")
        plt.xlabel("Channel \#")
        plt.ylabel("Count \#")
        plt.grid()
        
        # Fit (smooth x and fitted parameters)
        gauss_x = np.arange(x[0], x[-1], 0.0001)
        gauss_y = gaussian(gauss_x, *popt)
   
        plt.plot(gauss_x, gauss_y, label='fit')
        plt.plot(x,y, '.', label='data')
        plt.legend(loc=1)
        return np.array(gaussian_mean), np.array(gaussian_std)


def plotting(x, y, x_error, y_error, k0_switch):
    x_error = np.ones(len(x)) * x_error
    if k0_switch == True:
        # Weighing data
        weights = (1 / y_error)**2
        
        # Curve fitting data for linear fit, with weights
        popt, pcov = curve_fit(linear, x, y, sigma = weights)
        #popt, pcov = curve_fit(linear, x, y)
        
        plt.figure()
        plt.title("Calibration")
        plt.xlabel("Amplitude")
        plt.ylabel("Mean value of Gaussian Fit")
        plt.errorbar(x, y, xerr=x_error, yerr=y_error, fmt="o", label="Data")
        plt.plot(x, linear(x, *popt), label='Linear fit')
        plt.legend()
        plt.grid()

        # Determining parameters
        print("Linear fit parameters are a = {:.2f} and b = {:.2f}".format(popt[1],
            popt[0]))
        return popt
    else:
        popt, pcov = curve_fit(linear, x, y)

        plt.figure()
        plt.title("Calibration")
        plt.xlabel("Amplitude")
        plt.ylabel("Mean value of Gaussian fit")
        plt.plot(x, y, label="Data")
        plt.plot(x, linea(x, *popt), label="Linear fit")
        plt.legend()
        plt.grid()

        # Determine parameters
        print("Linear fit parameters are a = {:.2f} and b ={:.2f}".format(popt[1], popt[0]))
        return popt




# # # # # # # # # # # # # # # Output # # # # # # # # # # # # # # # # # # # # # 

# Determining k0 by changing amplitude (data_cal_dir)
def k0():
    # Determine k0 or alpha?
    k0_switch = True

    # Dataframe (data) + file_names (for itteration)
    df, file_names = Dataframe(data_cal_dir)
    
    # This is manual labour (extracting meassured amplitude from names)
    file_names.sort()
    Amps = np.array([1.62, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0])

    # Gaussian uncertainties and determining mean value of channel bin
    gaussian_mean, gaussian_std = Gaussian_fit(df, file_names, k0_switch)
    
    # Linear fit (parameters)
    # Meassured Amplitude is x and the mean channel number is y
    x = Amps
    y = gaussian_mean
    # Amplitude was meassured to lowest digit (0.1)
    x_error = np.array([0.1])
    # Mean channel number was meassured with gaussian uncertainty
    y_error = gaussian_std
    print(x_error)
    print(y_error)
    print(y_error**2)

#    b, a = plotting(Amps, gaussian_mean_k0, gaussian_std_k0, k0_switch)
    b, a = plotting(x, y, x_error, y_error, k0_switch)
    k0 = b
    print("k0 is {:.2f}".format(k0))

    return k0

# Determining alpha by changing matter; H to H2 (data_dir)
def alpha():
    # Determine k0 or alpha?
    k0_switch = False
    
    # Determine Energies
    EH, EH2 = ENERGY()

    # Dataframe (data) + file_names (for itteration)
    df_alpha, file_names_alpha = Dataframe(data_dir)

    name1 = file_names_alpha[0]
    name2 = file_names_alpha[1]

    # Manual labour part 2
    df_alpha2 = df_alpha.loc[100:500]

    # 3 Different peaks (for gaussians)
    df_alpha3 = df_alpha[name1].loc[100:260]
    df_alpha4 = df_alpha[name2].loc[260:367]
    df_alpha5 = df_alpha[name2].loc[367:500]

    plt.figure()
    df_alpha3.plot()
    df_alpha4.plot()
    df_alpha5.plot()

   # Gaussian uncertainties and determining mean value of channel bin
    gaussian_mean_alpha1, gaussian_std_alpha1 = Gaussian_fit(df_alpha3,
            name1, k0_switch)
    gaussian_mean_alpha2, gaussian_std_alpha2 = Gaussian_fit(df_alpha4,
            name2, k0_switch)
    gaussian_mean_alpha3, gaussian_std_alpha3 = Gaussian_fit(df_alpha5,
            [name1], k0_switch)
#
#
#    print(gaussian_mean_alpha1, gaussian_std_alpha1)
#    print(gaussian_mean_alpha2, gaussian_std_alpha2)
#    print(gaussian_mean_alpha3, gaussian_std_alpha3)
#
#
#    x = [gaussian_mean_alpha1, gaussian_mean_alpha3]
#    y = [EH, EH2]
#    x_error = [gaussian_std_alpha1, gaussian_std_alpha3]
#
#    # Linear fit (parameters) (chosen right values
#    a, b = plotting(x, y, x_error, k0_switch)
#
#    alpha = a
#    
#    print("alpha is {:.2f}".format(alpha))
    return


#k0()
alpha()

plt.show()



# # # # # # # # # # # # # # # Data Analysis # # # # # # # # # # # # # # # # # # 
def E():
    return alpha * (k - k0)




