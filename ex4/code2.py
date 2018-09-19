# # # # # # # # # # # # # # Loading Libraries # # # # # # # # # # # # # # # # 
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

# # # # # # # # # # # # # # Importing Data # # # # # # # # # # # # # # # # # #
""" Here we import the relevant data, and insert paths, directories in lists"""

# Current working directory
data_path = os.path.join(os.getcwd(), 'data/')
data_angles_path = os.path.join(os.getcwd(), 'data/Vinkler/')

# List of files in ./Data
data_files = os.listdir(data_path)
data_angles_files = os.listdir(data_angles_path)

# Filtering files
data_cal = [x for x in data_files if 'cal' in x and 'csv' in x]
data     = [x for x in data_files if 'AuC' in x and 'csv' in x]
data_angles = [x for x in data_angles_files  if 'csv' in x]

# Paths for chosen datafiles
data_cal_dir = [os.path.join(data_path, x) for x in data_cal]
data_dir     = [os.path.join(data_path, x) for x in data]
data_angles_dir = [os.path.join(data_angles_path, x) for x in data_angles]

# # # # # # # # # # # # # # Fitting Functions  # # # # # # # # # # # # # # # # 
""" Here we define the functions for which we will fit the parameters of"""

# Gaussians
""" A: Amplitude;    mu: Mean value;    std: Standard Deviation """

# # Single Gaussian
def gaussian(x, A, mu, std):
    return A * np.exp(-(x - mu)**2 / (2 * std**2))

# # Double Gaussian
def doublegaussian(x, A1, mu1, std1, A2, mu2, std2):
    return(
          A1 * np.exp(-(x - mu1)**2 / (2 * std1**2))
        + A2 * np.exp(-(x - mu2)**2 / (2 * std2**2)))

# Linear functions
# For 2 parameters (fitting k0)
def linear(x, a, b):
    return a*x + b

# For 1 parameter (having k0)
def linear2(x, a):
    global k0
    return a*x + k0

# # # # # # # # # # # # # # Energy for H+ and H2+ # # # # # # # # # # # # # # # 
def atomic_energy():
    # Incident energy (accelerator)
    Ei = 350 # keV
    # Meassured degree (theoretically, but identical to used experimental
    # setup)
    theta_deg = 160 # degrees
    theta = theta_deg * (np.pi/180) # radians
    # Target mass
    mAu = 196.9665690

    # Quanta of mass (assume m(H+)=m(H) and m(H2+)=2m(H)
    mH = 1.008

    # Formula for energy (derrived)
    E = lambda mb, mt: (((mb*np.cos(theta) + np.sqrt(mt**2 - mb**2 *
            (np.sin(theta))**2)) / (mb+mt))**2) * Ei

    # Energies
    E1 = E(mH, mAu)   #H+
    E2 = E(2*mH, mAu) #H2+

    return E1, E2

# # # # # # # # # # # # # # Energy-Calibration formula # # # # # # # # # # # # 
def energy_calibration(alpha, k0, k):
    return alpha * (k - k0)

# # # # # # # # # # # # # # DataFrame Import # # # # # # # # # # # # # # # # # 
""" Here we import given list of directories into dataframes """

def Dataframe(directory):

    # Preparing dataframe structure
    df = pd.DataFrame()

    # Itterating over all files
    file_names = []
    for i in range(len(directory)):
        # File name
        file_names.append(os.path.basename(directory[i]))

        # Reading csv
        column_names = [file_names[i]]
        df1 = pd.read_csv(directory[i], delimiter=',', index_col=0, names=column_names)

        # Concatenating
        df = pd.concat([df, df1[column_names[0]]], axis=1)

    # We make Dataframes and return file name
    return df, file_names

# # # # # # # # # # # # # # Gaussian fit # # # # # # # # # # # # # # # # # # #
def gaussian_fit(df, file_names, k0_switch, fig_ax_axins):
    fig, ax, axins = fig_ax_axins

    # Gaussian parameters with errors saved for later
    gaussian_amplitude, gaussian_amplitude_error = [], []
    gaussian_mean,      gaussian_mean_error      = [], []
    gaussian_std,       gaussian_std_error       = [], []


    # k0
    if k0_switch == True:
        for name in file_names:
            # Choosing non-zero values
            df_data = df[name].loc[df[name].nonzero()]

            # Transforming to np.array (Might be unnecesarry)
            np_data = np.array([np.array(df_data.index), df_data.values])
            bins = np_data[0]
            counts = np_data[1]

            # Scipy optimization (Gaussian)
            # # Guessing start values of parameters
            mean_bin = sum(bins * counts) / sum(counts)
            sigma = np.sqrt(sum(counts - (bins - mean_bin)**2) / sum(counts))
            amplitude = max(counts)


            # # Fit parameters
            popt, pcov = curve_fit(gaussian, bins, counts, p0=[amplitude,\
                mean_bin, sigma])
            # # Uncertainties of parameters
            perr = np.sqrt(np.diag(pcov))

            gaussian_amplitude.append(popt[0])
            gaussian_mean.append(popt[1])
            gaussian_std.append(popt[2])

            gaussian_amplitude_error.append(perr[0])
            gaussian_mean_error.append(perr[1])
            gaussian_std_error.append(perr[2])

            # Plotting data + fit
            gauss_x = np.arange(bins[0], bins[-1], 0.001)
            gauss_y = gaussian(gauss_x, *popt)

            ax.plot(gauss_x, gauss_y, label="Gaussian fit")
            ax.plot(bins, counts, 'o', label="Data")

            # Choosing 1 of the many to zoom in on
            if name == "cal2_0.csv":
                axins.plot(gauss_x, gauss_y, color="crimson")
                axins.plot(bins, counts, 'o', color="plum")


    #alpha
    else:
        # Defining dataframe
        df_data = df

        # Transforming to array (Might be unnecesarry)
        np_data = np.array([np.array(df_data.index), df_data.values])
        bins   = np_data[0]
        counts = np_data[1]

        # Scipy optimization (Gaussian)
        # # Guessing start values of parameters
        mean_bin = sum(bins * counts) / sum(counts)
        sigma = np.sqrt(sum(counts * (bins - mean_bin)**2) / sum(counts))
        amplitude = max(counts)

        # # Fit parameters
        popt, pcov = curve_fit(gaussian, bins, counts, p0=[amplitude,
            mean_bin, sigma])
        # # Uncertainties of parameters
        perr = np.sqrt(np.diag(pcov))

        gaussian_amplitude.append(popt[0])
        gaussian_mean.append(popt[1])
        gaussian_std.append(popt[2])

        gaussian_amplitude_error.append(perr[0])
        gaussian_mean_error.append(perr[1])
        gaussian_std_error.append(perr[2])

        # Plotting data + fit
        gauss_x = np.arange(bins[0], bins[-1], 0.001)
        gauss_y = gaussian(gauss_x, *popt)

        ax.plot(gauss_x, gauss_y, label="Gaussian fit")
        ax.plot(bins, counts, 'o', label="Data")

        if int(mean_bin) == 206:
            axins.plot(gauss_x, gauss_y)
            axins.plot(bins, counts, 'o')

    # The product of this function is all parameters with error
    np_amplitudes = np.array([gaussian_amplitude, gaussian_amplitude_error])
    np_mean = np.array([gaussian_mean, gaussian_mean_error])
    np_std  = np.array([gaussian_std, gaussian_std_error])

    return np_amplitudes, np_mean, np_std

# # # # # # # # # # # # # # Linear fit # # # # # # # # # # # # # # # # # # # #
def linear_fit(x, y, x_error, y_error, k0_switch, fig_ax_axins, *k0):
    fig, ax, axins = fig_ax_axins
    if k0_switch == True:
        # Weighing data
        weights = (1/y_error)**2

        # Fit parameter (with weights) (Linear)
        popt, pcov = curve_fit(linear, x, y, sigma = weights)
        perr = np.sqrt(np.diag(pcov))

        # Plotting for zoom
        axins.errorbar(x, linear(x, *popt))
        axins.errorbar(x, y, xerr=x_error, yerr=y_error, fmt="o")
    
        # Plotting 
        ax.errorbar(x, y, xerr=x_error, yerr=y_error, fmt="o", label="Data")
        ax.plot(x, linear(x, *popt), label="Linear fit")

        # The parameters are
        print("""Linear fit parameters are a={:.2f} with an error of {:.2f}
           and b={:.2f} with an error of {:.2f}""".format(popt[0], perr[0],
               popt[1], perr[1]))

    else:
        popt, pcov = curve_fit(linear2, x, y)
        print(popt)
        #popt, pcov = curve_fit(lambda x, a: linear(x, a, b), x, b)
        perr = np.sqrt(np.diag(pcov))

        # Plotting for zoom
        axins.errorbar(x, linear2(x, *popt))
        axins.errorbar(x, y, xerr=x_error, yerr=y_error, fmt="o")
    
        # Plotting 
        ax.errorbar(x, y, xerr=x_error, yerr=y_error, fmt="o", label="Data")
        ax.plot(x, linear2(x, *popt), label="Linear fit")

        # The parameters are
        print("""Linear fit parameters are a={:.2f} with an error of
                {:.2f}""".format(popt[0], perr[0]))

    return popt, perr

# # # # # # # # # # # # # # Zoom plot # # # # # # # # # # # # # # # # # # # # #
def zoomplot(figure_settings, zoom_settings, *colorz):
    """ Zooming in on a part of the graph """

    # Title, xlabel, ylabel are basic figure settings
    title, xlabel, ylabel = figure_settings

    # zoom_factor = strength
    # zoom_plot_loc = 0 (best),
    # 1(upper right), 2(upper left), 3(lower left) , 4(lower right)
    # x, y (min/max) are the intervals (zoomed square)

    # Extracting values
    # # zoom_factor: Strength;      zoom_plot_loc = location (0, 1, 2, 3, 4)
    # # x, y: square of zoom;       corner = connecting corners (1, 2, 3, 4)
    zoom_factor, zoom_plot_loc, \
    xmin, xmax, ymin, ymax, corner1, corner2 = zoom_settings

    # Preparing figure
    fig, ax = plt.subplots()
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.grid()

    # Zoom properties
    zoom_factor = zoom_factor

    # Make a double figure (on top of each other)
    axins = zoomed_inset_axes(ax, zoom_factor, loc=zoom_plot_loc)

    # Setting the limitations of double (area of zoom)
    axins.set_xlim(xmin, xmax)
    axins.set_ylim(ymin, ymax)

    # Removing ticks on double
    plt.xticks(visible = False)
    plt.yticks(visible = False)

    # fc = Fill coulour,    ec = line colour
    mark_inset(ax, axins, loc1=corner1, loc2=corner2, fc="none", ec="0.5")

    return fig, ax, axins

# # # # # # # # # # # # # # Calibration # # # # # # # # # # # # # # # # # # 
""" In this calibration we determine alpha and k0, such we can transform
channel numbers to energy by the E = alpha(k - k0), where k is the channel
number and k0 is the channelnumber for zero amplitude (extrapolation)"""

def k0():
    """ Here we fit a Gaussian for each plot, use the mean bin and plot as a
    function of the meassured energies (from a pulser). This is fitted
    linearly, to obtain the zero-amplitude bin (b-value). The inclince is not
    used, as it is determined from alpha()"""

    # Determine k0 or alpha?
    k0_switch = True

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Gaussian
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    # Dataframe (data) + file_names (for itteration)
    df, file_names = Dataframe(data_cal_dir)

    # This is manual labour (extracting meassured amplitude from names)
    file_names.sort() # Makine sure in the right order
    energy_input = np.array([1.62, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0])
    # Meassured to lowest digit (uncertainty)
    energy_input_error = np.ones(np.size(energy_input)) * 0.1

    # Plot settings
    title = "Calibration"
    xlabel = "Channel Numbers"
    ylabel = "Counts numbers"

    figure_settings = [title, xlabel, ylabel]

    zoom_factor = 2
    zoom_plot_loc = 2
    xmin, xmax, ymin, ymax = [220, 260, 0, 5000]
    corner1, corner2 = [3, 4]

    zoom_settings = [zoom_factor, zoom_plot_loc, xmin, xmax, ymin, ymax,
            corner1, corner2]

    # Calling zoom plot
    fig, ax, axins = zoomplot(figure_settings, zoom_settings)
    fig_ax_axins = [fig, ax, axins]

    # Gaussian uncertainties and determining mean value of channel bin
    # Each are numpy arrays with uncertainties on index [1]
    amplitude, mean, std = gaussian_fit(df, file_names, k0_switch,
            fig_ax_axins)

    # Last details of figure
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles[:2], labels[:2], loc=1, ncol=2, borderaxespad=0,
            frameon=False)
    #plt.savefig(Gaussian_fit)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Linear
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    # Linear fit (parameters)
    # Preparing figure
    # Figure settings
    title = "Calibration"
    xlabel = "Amplitude"
    ylabel = "Mean value of Gaussian Fit"
    figure_settings = [title, xlabel, ylabel]

    # Zoom settings
    zoom_factor = 5
    zoom_plot_loc = 2
    xmin, xmax, ymin, ymax = [2.8, 3.2, 338, 378]
    corner1, corner2 = [3, 4]
    zoom_settings = [zoom_factor, zoom_plot_loc, xmin, xmax, ymin, ymax,
            corner1, corner2]

    fig, ax, axins = zoomplot(figure_settings, zoom_settings)
    fig_ax_axins = [fig, ax, axins]

    popt, perr = linear_fit(energy_input, mean[0], energy_input_error, mean[1],
            k0_switch, fig_ax_axins)
    a = np.array([popt[0], perr[0]])
    b = np.array([popt[1], perr[1]])

    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels, loc=1, ncol=2, borderaxespad=0, frameon=False)

    #plt.savefig("k0_plotting")
    # We fit for (zero-amplitude/b) value, and obtain alpha in other function
    k0_val = b[0]
    k0_error = b[1]
    k0 = np.array([k0_val, k0_error])
    print("k0 is {:.2f} with an error of {:.2f}".format(k0_val, k0_error))

    return k0

def alpha(k0):
    """ Now we obtain the incline. Same procedure, but only 2 points. One from
    H+ another from H2+"""
    # Determine k0 or alpha?
    k0_switch = False

    # Determine Energies
    E = atomic_energy()

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Gaussian
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    df, file_names = Dataframe(data_dir)

    # We have 2 files (H+ and H2+)
    name1 = file_names[0] # H2+
    name2 = file_names[1] # H+

    # Manual labour part 2
    df1 = df.loc[100:500]

    # 3 Different peaks (for gaussians)
    df2 = df[name1].loc[100:260] # H2+ (GOLD)
    df3 = df[name2].loc[260:367] # H+  (CARBON)
    df4 = df[name2].loc[367:500] # H+  (GOLD)

   # Gaussian uncertainties and determining mean value of channel bin
   # Figure of the gaussians
    # Preparing figure
    title = "Calibration"
    xlabel = "Channel Number"
    ylabel = "Number of counts"
    figure_settings = [title, xlabel, ylabel]

    # Preparing zoom
    zoom_factor = 2
    zoom_plot_loc = 2
    xmin, xmax, ymin, ymax = [165, 255, 0, 42]
    corner1, corner2 = [1, 3]
    zoom_settings = [zoom_factor, zoom_plot_loc, xmin, xmax, ymin, ymax, corner1,
            corner2]

    # Calling function
    fig, ax, axins = zoomplot(figure_settings, zoom_settings)
    fig_ax_axins = [fig, ax, axins]

    # Those of interests are mean1 and mean3 
    amplitude1, mean1, std1 = gaussian_fit(df2, name1, k0_switch, fig_ax_axins)
    amplitude2, mean2, std2 = gaussian_fit(df3, name2, k0_switch, fig_ax_axins)
    amplitude3, mean3, std3 = gaussian_fit(df4, name1, k0_switch, fig_ax_axins)

    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles[:2], labels[:2], loc=1, ncol=2, borderaxespad=0, frameon=False)

    #plt.savefig('gaussian_fit2')

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Linear
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    # Fitting mean values for gold (1 and 3)
    x = np.concatenate((mean1[0], mean3[0]), axis=0)
    x_error = np.concatenate((mean1[1], mean3[1]), axis=0)

    # Energies for H2+ and H+ (Respectively)
    y = [E[1], E[0]]
    
    y_error = [0.2, 0.2] # Theoretical? No, we meassured energy


    # Preparing figure
    title = "Calibration"
    xlabel = "Energy"
    ylabel = "Mean value of Gaussian Fit"
    figure_settings = [title, xlabel, ylabel]

    ## Preparing zoom
    zoom_factor = 11
    zoom_plot_loc = 2
    xmin, xmax, ymin, ymax = [16, 18, 336.15, 336.6]
    corner1, corner2 = [3, 4]
    zoom_settings = [zoom_factor, zoom_plot_loc, xmin, xmax, ymin, ymax,
            corner1, corner2]

    # Calling zoomfunction
    fig, ax, axins = zoomplot(figure_settings, zoom_settings)
    fig_ax_axins = [fig, ax, axins]

    popt, perr = linear_fit(x, y, x_error, y_error, k0_switch, fig_ax_axins, k0)
    # Only fitting 1 parameter = the incline (alpha)
    a_val   = popt[0]
    a_error = perr[0]

    a = np.array([a_val, a_error])

    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels, loc=1, ncol=2, borderaxespad=0, frameon=False)

    #plt.savefig("alpha_plotting")

    alpha = a
    print("alpha is {:.2f} with an uncertainty of {:.2f}".format(alpha[0],
        alpha[1]))

    return alpha




# # # # # # # # # # # # # # Data Analysis # # # # # # # # # # # # # # # # # # #




# # # # # # # # # # # # # # Function calls # # # # # # # # # # # # # # # # # #

# Determine k0
k0 = k0()
k0_val, k0_error = [k0[0], k0[1]]

# Determine alpha
alpha = alpha(k0_val)
alpha_val, alpha_error = [alpha[0], alpha[1]]





plt.show()



