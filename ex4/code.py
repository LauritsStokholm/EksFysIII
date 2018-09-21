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
    Ei = 350 # kev
    theta_deg = 160 # degrees
    theta = 160 * (np.pi / 180)
    mAu = 196.9665690
    mH = 1.008

    E = lambda mb, mt: (((mb*np.cos(theta) + np.sqrt(mt**2 - mb**2 *
            (np.sin(theta))**2)) / (mb+mt))**2) * Ei

    E1 = E(mH, mAu)
    E2 = E(2*mH, mAu)

    return np.array([E1, E2])



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

def Gaussian_fit(df, file_names, k0_switch, *alpha):
    """ This is a Gaussian fit factory"""
    if k0_switch == True:
        # Gaussian mean and std saved for later
        gaussian_mean = []
        gaussian_std  = []

         # Preparing figure
        fig, ax = plt.subplots()
        plt.title("Calibration")
        plt.xlabel("Channel \#")
        plt.ylabel("Number of counts")
        plt.grid()
        
        # Preparing zoom
        zoom_factor = 2
        
        # Make a double of figure untop of old
        axins = zoomed_inset_axes(ax, zoom_factor, loc=2) 
        
        # Setting limitation of double (minimal)
        x1, x2, y1, y2 = 220, 260, 0, 5000
        axins.set_xlim(x1, x2)
        axins.set_ylim(y1, y2)
        
        # Removing ticks on double
        plt.xticks(visible=False)
        plt.yticks(visible=False)
        
        # fc (fill colour) and ec (line colour)
        # loc1 and loc2 is the corners to connect
        mark_inset(ax, axins, loc1=3, loc2=4, fc="none", ec="0.3")

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
            print(sigma)
            print(mean)
            
            # Gaussian fit parameters
            popt,pcov = curve_fit(gaussian, x, y, p0=[max(y), mean, sigma])
            
            gaussian_mean.append(popt[1])
            gaussian_std.append(popt[2])
            
            
            # Plotting the data + fit
            # Fit (smooth x and fitted parameters)
            gauss_x = np.arange(x[0], x[-1], 0.0001)
            gauss_y = gaussian(gauss_x, *popt)
   
            ax.plot(gauss_x, gauss_y, label='Gaussian fit')
            ax.plot(x,y, 'o', label='data')
            if name == "cal2_0.csv":
                axins.plot(gauss_x, gauss_y)
                axins.plot(x, y, 'o')
        handles, labels = ax.get_legend_handles_labels()
        fig.legend(handles[:2], labels[:2], loc=1, ncol=2, borderaxespad=0, frameon=False)

        plt.savefig('gaussian_fit')

        return np.array([gaussian_mean, gaussian_std])
    else:
        fig, ax, axins = alpha
        # Gaussian mean and std saved for later
        gaussian_mean = []
        gaussian_std  = []

        df_data = df

        # Transforming to array (Might be unnecesarry) [index, count]
        np_data = np.array([np.array(df_data.index), df_data.values])
        x = np_data[0] # Indices
        y = np_data[1] # Count number
        
        # Scipy Optimization after Gaussian (guessing start values)
        mean = sum(x * y) / sum(y)
        sigma = np.sqrt(sum(y * (x - mean)**2) / sum(y))

        # Gaussian fit parameters
        popt,pcov = curve_fit(gaussian, x, y, p0=[max(y), mean, sigma])
        
        gaussian_mean.append(popt[1])
        gaussian_std.append(popt[2])
        
        # Plotting the data + fit
        
        # Fit (smooth x and fitted parameters)
        gauss_x = np.arange(x[0], x[-1], 0.0001)
        gauss_y = gaussian(gauss_x, *popt)

   
        ax.plot(gauss_x, gauss_y, label='Gaussian fit')
        ax.plot(x, y, 'o', label='data')

        if int(mean) == 206:
            axins.plot(gauss_x, gauss_y)
            axins.plot(x, y, 'o')

        return np.array([gaussian_mean, gaussian_std])


def plotting(x, y, x_error, y_error, k0_switch, fig, ax, axins):
    x_error = np.ones(len(x)) * x_error
    if k0_switch == True:
        # Weighing data
        weights = (1 / y_error)**2

        # Curve fitting data for linear fit, with weights
        popt, pcov = curve_fit(linear, x, y, sigma = weights)
    else:
        popt = curve_fit(linear, x, y)[0]
        
    # Plotting for zoom
    axins.errorbar(x, linear(x, *popt))
    axins.errorbar(x, y, xerr=x_error, yerr=y_error, fmt="o")
        
    # Plotting whole scale
    ax.errorbar(x, y, xerr=x_error, yerr=y_error, fmt="o", label="Data")
    ax.plot(x, linear(x, *popt), label='Linear fit')

    # Determining parameters
    print("Linear fit parameters are a = {:.2f} and b = {:.2f}".format(popt[1],
    popt[0]))
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
    gauss = Gaussian_fit(df, file_names, k0_switch)
    gaussian_mean = gauss[0]
    gaussian_std  = gauss[1]
    
    # Linear fit (parameters)
    # Meassured Amplitude is x and the mean channel number is y
    x = Amps
    y = gaussian_mean
    # Amplitude was meassured to lowest digit (0.1)
    x_error = np.array([0.1])
    # Mean channel number was meassured with gaussian uncertainty
    y_error = gaussian_std

    # Preparing figure
    fig, ax = plt.subplots()
    plt.title("Calibration")
    plt.xlabel("Amplitude")
    plt.ylabel("Mean value of Gaussian Fit")
    plt.grid()

    # Preparing zoom
    zoom_factor = 5

    # Make a double of figure untop of old
    axins = zoomed_inset_axes(ax, zoom_factor, loc=2) 

    # Setting limitation of double (minimal)
    x1, x2, y1, y2 = 2.8, 3.2, 338, 378
    axins.set_xlim(x1, x2)
    axins.set_ylim(y1, y2)

    # Removing ticks on double
    plt.xticks(visible=False)
    plt.yticks(visible=False)

    # fc (fill colour) and ec (line colour)
    # loc1 and loc2 is the corners to connect
    mark_inset(ax, axins, loc1=1, loc2=3, fc="none", ec="0.3")
    b, a = plotting(x, y, x_error, y_error, k0_switch, fig, ax, axins)

    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels, loc=1, ncol=2, borderaxespad=0, frameon=False)

    #plt.savefig("k0_plotting")
    k0 = b
    print("k0 is {:.2f}".format(k0))

    return k0

# Determining alpha by changing matter; H to H2 (data_dir)
def alpha():
    # Determine k0 or alpha?
    k0_switch = False
    
    # Determine Energies
    y = np.array([ENERGY()[1], ENERGY()[0]])
    print(y)

    # Dataframe (data) + file_names (for itteration)
    df, file_names = Dataframe(data_dir)

    name1 = file_names[0]
    name2 = file_names[1]

    # Manual labour part 2
    df1 = df.loc[100:500]

    # 3 Different peaks (for gaussians)
    df2 = df[name1].loc[100:260]
    df3 = df[name2].loc[260:367]
    df4 = df[name2].loc[367:500]


   # Gaussian uncertainties and determining mean value of channel bin


   # Figure of the gaussians
    # Preparing figure
    fig, ax = plt.subplots()
    plt.title("Calibration")
    plt.xlabel("Channel \#")
    plt.ylabel("Number of counts")
    plt.grid()
    
    # Preparing zoom
    zoom_factor = 2
    
    # Make a double of figure untop of old
    axins = zoomed_inset_axes(ax, zoom_factor, loc=2) 
    
    # Setting limitation of double (minimal)
    x1, x2, y1, y2 = 165, 255, 0, 42
    axins.set_xlim(x1, x2)
    axins.set_ylim(y1, y2)
    
    # Removing ticks on double
    plt.xticks(visible=False)
    plt.yticks(visible=False)
    
    # fc (fill colour) and ec (line colour)
    # loc1 and loc2 is the corners to connect
    mark_inset(ax, axins, loc1=1, loc2=3, fc="none", ec="0.3")

    gauss1 = Gaussian_fit(df2, name1, k0_switch, fig, ax, axins)
    gauss2 = Gaussian_fit(df3, name2, k0_switch, fig, ax, axins)
    gauss3 = Gaussian_fit(df4, name1, k0_switch, fig, ax, axins)

    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles[:2], labels[:2], loc=1, ncol=2, borderaxespad=0, frameon=False)

    plt.savefig('gaussian_fit2')

    x = np.concatenate((gauss1[0], gauss3[0]), axis=0)
    x_error = [gauss1[1], gauss3[1]]
    y_error = [0, 0] # Theoretical?


    # Preparing figure
    fig, ax = plt.subplots()
    plt.title("Calibration")
    plt.xlabel("Energy")
    plt.ylabel("Mean value of Gaussian Fit")
    plt.grid()
    
    # Preparing zoom
    zoom_factor = 5
    
    # Make a double of figure untop of old
    axins = zoomed_inset_axes(ax, zoom_factor, loc=2) 
    
    # Setting limitation of double (minimal)
    x1, x2, y1, y2 = 195, 220, 336.2, 336.5
    axins.set_xlim(x1, x2)
    axins.set_ylim(y1, y2)
    
    # Removing ticks on double
    plt.xticks(visible=False)
    plt.yticks(visible=False)
    
    # fc (fill colour) and ec (line colour)
    # loc1 and loc2 is the corners to connect
    mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.3")
    b, a = plotting(x, y, x_error, y_error, k0_switch, fig, ax, axins)
    
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels, loc=1, ncol=2, borderaxespad=0, frameon=False)
    
    #plt.savefig("k0_plotting")



    # Linear fit (parameters) (chosen right values
    a, b = plotting(x, y, x_error, y_error, k0_switch, fig, ax, axins)
    plt.savefig("alpha_plotting")

    alpha = a
    print("alpha is {:.2f}".format(alpha))
    return alpha


k0 = k0()
alpha = alpha()




# # # # # # # # # # # # # # # Data Analysis # # # # # # # # # # # # # # # # # # 
def E(alpha, k0, k):
    return alpha * (k - k0)

#k = np.arange(0, 1000)
#plt.figure()
#plt.plot(k, E(alpha, k0, k))

plt.show()
