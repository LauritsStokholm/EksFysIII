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
data_thickness_path = os.path.join(os.getcwd(), 'data/thickness/')

# List of files in ./Data
data_files = os.listdir(data_path)
data_angles_files = os.listdir(data_angles_path)
data_thickness_files = os.listdir(data_thickness_path)

# Filtering files
data_cal = [x for x in data_files if 'cal' in x and 'csv' in x]
data     = [x for x in data_files if 'AuC' in x and 'csv' in x]
data_angles = [x for x in data_angles_files  if 'csv' in x]
data_thickness = [x for x in data_thickness_files  if 'csv' in x]

# Paths for chosen datafiles
data_cal_dir = [os.path.join(data_path, x) for x in data_cal]
data_dir     = [os.path.join(data_path, x) for x in data]
data_angles_dir = [os.path.join(data_angles_path, x) for x in data_angles]
data_thickness_dir = [os.path.join(data_thickness_path, x) for x in
        data_thickness]

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
    # E = alpha * (k-k0)
    global k0
    return a * (x - k0[0]) # k0 is a matrix with uncertainty, pick out value

# # # # # # # # # # # # # # Energy for H+ and H2+ # # # # # # # # # # # # # # # 
def atomic_energy(switch):
    # Incident energy (accelerator)
    Ei = 350  # keV

    # Target mass
    mAu = 196.9665690
    mC =  12.0107

    # Quanta of mass (assume m(H+)=m(H) and m(H2+)=2m(H)
    mH = 1.008

    # Meassured degree (theoretically, but identical to used experimental
    # setup)
    theta_deg = 160 # degrees
    theta = theta_deg * (np.pi/180) # radians

    # Formula for energy (derrived)
    # E2 = k**2 * E1
    E = lambda mb, mt: (((mb*np.cos(theta) + np.sqrt(mt**2 - (mb**2
        *(np.sin(theta))**2))) / (mb+mt))**2) * Ei
    # Meassuring on gold?
    if switch == True:
        E1 = E(mH, mAu)   #H+
        E2 = E(2*mH, mAu)/2 # H2+
        return E1, E2
    elif switch == False:
        E1 = E(mH, mC)     # H+
        E2 = E(2*mH, mC)/2 # H2+
        return E1 # Only this is seen (middle of the 3 Gaussians)
    else:
        print("atomic_energy has a false switch")

# # # # # # # # # # # # # # Energy-Calibration formula # # # # # # # # # # # # 
def energy_calibration(alpha, k0, k):
    return alpha*(k-k0)

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
                axins.plot(gauss_x, gauss_y, color="blue")
                axins.plot(bins, counts, 'o', color="red")


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
        perr = np.sqrt(np.diag(pcov))

        # Plotting for zoom
        axins.errorbar(x, linear2(x, *popt))
        axins.errorbar(x, y, xerr=x_error, yerr=y_error, fmt="o")
    
        # Plotting 
        ax.errorbar(x, y, xerr=x_error, yerr=y_error, fmt="o", label="Data")
        ax.plot(x, linear2(x, *popt), label="Linear fit")

        # The parameters are
        print("""Linear fit parameters are a={:.2f} with an error of
                {:.5f}""".format(popt[0], perr[0]))

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
    plt.savefig("gaussian_fit")

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

    plt.savefig("k0_plotting")
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
    E_Au = np.array(atomic_energy(True))
    E_C = np.array(atomic_energy(False))

    # Combining in right order
    E = np.array([E_Au[1], E_C, E_Au[0]])

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

    plt.savefig('gaussian_fit2')

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Linear
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    # Fitting mean values for gold (1 and 3)

    x = np.concatenate(([mean1[0], mean2[0], mean3[0]]), axis=0)
    x_error = np.concatenate(([mean1[1], mean2[1], mean3[1]]), axis=0)

    # Energies for H2+ and H+ (Respectively)
    y = E
    y_error = np.ones(np.size(y)) * 0.02 # Theoretical? No, we meassured energy

    # Preparing figure
    title = "Calibration"
    xlabel = "Mean value of Gaussian Fit"
    ylabel = r"Energy [\si{\kilo\electronvolt}]"
    figure_settings = [title, xlabel, ylabel]

    ## Preparing zoom
    zoom_factor = 5
    zoom_plot_loc = 2
    xmin, xmax, ymin, ymax = [x[0]-1, x[0]+1, y[0]-1, y[0]+1]
    corner1, corner2 = [3, 4]
    zoom_settings = [zoom_factor, zoom_plot_loc, xmin, xmax, ymin, ymax,
            corner1, corner2]

    # Calling zoomfunction
    fig, ax, axins = zoomplot(figure_settings, zoom_settings)
    fig_ax_axins = [fig, ax, axins]

    popt, perr = linear_fit(x, y, x_error, y_error, k0_switch, fig_ax_axins, k0)
    #ax.errorbar(mean2[0], E_C, xerr=mean2[1], yerr=0.2, fmt='o')

    # Only fitting 1 parameter = the incline (alpha)
    a_val   = popt[0]
    a_error = perr[0]

    a = np.array([a_val, a_error])

    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels, loc=1, ncol=2, borderaxespad=0, frameon=False)

    plt.savefig("alpha_plotting")

    alpha = a
    print("alpha is {:.2f} with an uncertainty of {:.5f}".format(alpha[0],
        alpha[1]))

    return alpha


def channel_to_energy(alpha, k0, channel):
    # Packing out values
    alpha_val = alpha[0]
    alpha_error = alpha[1]
    k0_val = k0[0]
    k0_error = k0[1]

    bin_val = channel[0]
    bin_error = channel[1]

    # Propagation of error
    # E = alpha * (k - k0)
    dE_dk0    = - alpha_val
    dE_dalpha = + alpha_val
    dE_dk     = + alpha_val
    

    # k0_error
    # alpha_error
    # Uncertainty of Energy (error of propagation)
    E_error = np.sqrt(((dE_dk0) * k0_error)**2 + ((dE_dalpha) * alpha_error)**2
            + ((dE_dk) * bin_error)**2)

    E_val = alpha_val * (bin_val - k0_val)

    E = np.array([E_val, E_error])

    print(E)
    return E

# # # # # # # # # # # # # # Data Analysis # # # # # # # # # # # # # # # # # # #
def doublegaussian_fit(plot_switch):
    """ Here we obtain all parameters of both Gold and Carbon gaussians"""
    # To avoid several function calls:
    k0 = 1.36
    k0_error = 0.15
    alpha = 0.80
    alpha_error = 0.00473

    # Dataframe (data) + file_names (for itteration)
    df, file_names = Dataframe(data_angles_dir)

    # Manual labour pt. 3
    file_names.sort()

    print(file_names)
    file_counts = np.array([6867, 9675, 9638, 47449, 48048, 45934,
        32847, 42593, 33412, 11769, 12881,16538, 16209])

    detector_angles = np.array([160, 150, 140, 40, 50, 130, 120, 60, 70,
        90, 100, 110, 80])

    print("Sizes:")
    print(np.size(file_counts))
    print(np.size(detector_angles))

    # Defining double gaussian parameters

    # Carbon
    gaussian_amplitude1       = []
    gaussian_mean1            = []
    gaussian_std1             = []
    gaussian_amplitude_error1 = []
    gaussian_mean_error1      = []
    gaussian_std_error1       = []
   
    # Gold
    gaussian_amplitude2       = []
    gaussian_mean2            = []
    gaussian_std2             = []
    gaussian_amplitude_error2 = []
    gaussian_mean_error2      = []
    gaussian_std_error2       = []
    
    # Targets
    dN_Au = []
    dN_C  = []

    for name in file_names:
        # Choosing non-zero values
        df_data = df[name].loc[df[name].nonzero()]
        df_data = df_data.loc[100:]

        # Transforming to np.array (Might be unnecesarry)
        np_data = np.array([np.array(df_data.index), df_data.values])
        bins = np_data[0]
        counts = np_data[1]

        # Scipy optimization (Double Gaussian)
        # # Guessing start values of parameters
        if "90" in name:
        # Carbon
            amplitude1 = 320
            mean_bin1 = 130
            sigma1 = 100
        # Gold
            amplitude2 = 284
            mean_bin2 = 430
            sigma = 25
        elif "100" in name:
        # Carbon
            amplitude1 = 75
            mean1_bin1 = 300
            sigma1 = 100

        # Gold
            amplitude2 = max(counts)
            mean_bin2 = sum(bins * counts) / sum(counts)
            sigma2 = np.sqrt(sum(counts + (bins - mean_bin2)**2) / sum(counts))

        else:
            amplitude1 = 50
            mean_bin1 = 200
            sigma1    = 10

            # Gold
            amplitude2 = max(counts)
            mean_bin2 = sum(bins * counts) / sum(counts)
            sigma2 = np.sqrt(sum(counts + (bins - mean_bin2)**2) / sum(counts))


        # # Fit parameters
        popt, pcov = curve_fit(doublegaussian, bins, counts, p0=[amplitude1,\
            mean_bin1, sigma1, amplitude2, mean_bin2, sigma2])
        # # Uncertainties of parameters
        perr = np.sqrt(np.diag(pcov))

        # Gold
        if "100det" in name:
            gaussian_amplitude2.append(popt[0])
            gaussian_mean2.append(popt[1])
            gaussian_std2.append(popt[2])

            gaussian_amplitude_error2.append(perr[0])
            gaussian_mean_error2.append(perr[1])
            gaussian_std_error2.append(perr[2])

        # Carbon
            gaussian_amplitude1.append(popt[3])
            gaussian_mean1.append(popt[4])
            gaussian_std1.append(popt[5])

            gaussian_amplitude_error1.append(perr[3])
            gaussian_mean_error1.append(perr[4])
            gaussian_std_error1.append(perr[5])
        else: 
            gaussian_amplitude1.append(popt[0])
            gaussian_mean1.append(popt[1])
            gaussian_std1.append(popt[2])

            gaussian_amplitude_error1.append(perr[0])
            gaussian_mean_error1.append(perr[1])
            gaussian_std_error1.append(perr[2])

        # Carbon
            gaussian_amplitude2.append(popt[3])
            gaussian_mean2.append(popt[4])
            gaussian_std2.append(popt[5])

            gaussian_amplitude_error2.append(perr[3])
            gaussian_mean_error2.append(perr[4])
            gaussian_std_error2.append(perr[5])



        if plot_switch == True:
            # Plotting data + fit
            plt.figure()
            plt.title("Double Gaussian fit")
            plt.xlabel("Channel number")
            plt.ylabel("Number of counts")
            plt.grid()

            gauss_x = np.arange(bins[0], bins[-1], 0.001)
            
            doublegauss_y = doublegaussian(gauss_x, *popt)
            gauss_C_y     = gaussian(gauss_x, *popt[:3])
            gauss_Au_y    = gaussian(gauss_x, *popt[3:])
    
            plt.plot(gauss_x, doublegauss_y, label="Double gaussian fit")
            plt.plot(gauss_x, gauss_C_y,     label="Single gaussian fit (C)")
            plt.plot(gauss_x, gauss_Au_y,    label="Single gaussian fit (Au)")
            plt.plot(bins, counts,           'o', label="Data")

            plt.legend()


        # Last thing is to calculate the dN for each target

        # Each value of gaussian evaluated in interval (histogram-like)
        x = np.arange(bins[0], bins[-1], 1)
        y_C  = gaussian(x, gaussian_amplitude1, gaussian_mean1, gaussian_std1)
        y_Au  = gaussian(x, gaussian_amplitude2, gaussian_mean2, gaussian_std2)

        plt.figure()
        plt.plot(x, y_C)
        plt.plot(x, y_Au)

        dN_C.append(np.sum(y_C))
        dN_Au.append(np.sum(y_Au))



    # Parameters and uncertainties
    # Carbon
    np_amplitudes1 = np.array([gaussian_amplitude1, gaussian_amplitude_error1])
    np_mean1       = np.array([gaussian_mean1, gaussian_mean_error1])
    np_std1        = np.array([gaussian_std1, gaussian_std_error1])

    # Gold
    np_amplitudes2 = np.array([gaussian_amplitude2, gaussian_amplitude_error2])
    np_mean2       = np.array([gaussian_mean2, gaussian_mean_error2])
    np_std2        = np.array([gaussian_std2, gaussian_std_error2])

    dN_C = np.array(dN_C)
    dN_Au = np.array(dN_Au)

    print(np.size(dN_C))
    print(np.size(dN_Au))

    return(detector_angles, np_amplitudes1, np_mean1, np_std1, np_amplitudes2,
np_mean2, np_std2, dN_C, dN_Au)

def energy_angle_plot(theta, E_Au, E_C):
    """ Plotting energy as a function of scattering angle (detector)"""
    detector_angles = theta

    # Meassured error
    theta_error = np.ones(np.size(detector_angles)) * 0.1

    E_Au_val = E_Au[0]
    E_Au_error = E_Au[1]

    E_C_val = E_C[0]
    E_C_error = E_C[1]

    plt.figure()
    plt.grid()
    plt.xlabel("Detector Angles")
    plt.ylabel(r"Energy $[\si{\kilo\electronvolt}]$")
    plt.errorbar(theta, E_Au_val, xerr=theta_error, yerr=E_Au_error, fmt="o", label="Au")
    plt.errorbar(theta, E_C_val,  xerr=theta_error, yerr=E_C_error,  fmt="o", label="C")
    plt.legend()


    return

def RutherfordCrossSection(dN, target_switch):

# # # # # # # # # # # # # # Define constants  # # # # # # # # # # # # # # # # #
    # Atomic number
    A_Au = 197
    A_C  = 12
    A_p  = 1

    # Density of elements
    rho_Au = 19.30 * 10**(3) # kg/m3
    rho_C  = 3.50 * 10**(3) # kg/m3
    
    # Molar masses
    u = 1.660539 * 10**(-27) # atomic mass units
    M_Au = 197.96655 * u # kg
    M_C  = 12.0107  * u # kg
    sigma_M_Au   = 0.000004 # Google
    sigma_M_C    = 0.0008 # Google
# # # # # # # # # # # # # # Meassured values  # # # # # # # # # # # # # # # # #

    file_counts = np.array([6867, 9675, 9638, 47449, 48048, 45934,
        32847, 42593, 33412, 11769, 12881, 16538, 16209])

    # Counts (10**11 counts per coulomb) (READ OFF MACHINE)

    # Conversion factor
    beta = 6.24150913 * 10**(18) # e/C
    N = file_counts * 10**(-11) * beta # [#protoner]
    # [counts] * [C/counts] * [e/C]

    if target_switch == True:
        print("Gold")
        # Target thickness (read off whiteboard)
        # 25 Angstrom
        dx = 25 * 10**(-10) # [m]
        M = M_Au
        sigma_M = sigma_M_Au
        n = rho_Au / (M_Au)
        #n = 5.7 * 10**(28)
        print(n)
    elif target_switch == False: 
        print("Carbon")
        # Target thickness (read off whiteboard)
        # 200 Angstrom
        dx = 200 * 10**(-10) # [m]
        M = M_C
        sigma_M = sigma_M_C
        n = rho_C / (M_C)
        # [kg/m3] * [mol/kg] * [#atomer/mol] * [#protomer/#atomer]
    else: print("Wrong switch")
    
    # Solid angle
    # Distance from detector to target
    r = 49 * 10**(-3) # m
    # Diameter of detector
    D = 1.88 * 10**(-3) # m
    dA = np.pi * (D/2)**2 # [m2]

    domega = dA / r**2 # per definition

    # Rutherford Cross section
    dsdo = (1 / (N*n*dx)) * (dN/domega) # [#protoner] [m2/sr]
    # (([#protoner]/[sr])  / ([#protoner][#protoner/m3][m])))
    dsdo = dsdo * (10**(31)) # [#protoner] [mb/sr]
    print(dsdo)


# # # # # # # # # # # # # # Error of Propagation # # # # # # # # # # # # # # # #

    # Poisson (detectors)
    # sigma_N
    sigma_N = np.sqrt(N) 
    sigma_dN = np.sqrt(dN)

    print("Poisssssssssson")
    print(sigma_N)
    print(sigma_dN)

    # Target thickness (Assumed???)
    sigma_dx = 10**(-10) # 10 Angstrom ?

    # Densities (as found on google at lowest digit)
    sigma_rho = 0.001 
    # sigma_m is defined in if/else

    # # Error of propagation (For other calculated values)
    # # # sigma_n # # # 
    # # # n = rho/M

    dn_drho = 1 / M
    dn_dM  = -rho_Au / M**2

    sigma_n = np.sqrt((dn_drho * sigma_rho)**2 + (dn_dM * sigma_M)**2)

    # # # sigma_do
    # do = (1/r**2) * (pi * (D/2)**2)

    do_dr = -2* (1/r**3) * (dA)
    do_dD = (np.pi * D) / (r**2)

    # Read off whiteboard in lab
    sigma_r =  1 * 10**(-3)
    sigma_D =  0.001 * 10**(-3)

    sigma_domega = np.sqrt( (do_dr * sigma_r)**2 + (do_dD * sigma_D)**2)


    # Now for the real parameter error of propagation
    # choose dsdo = y
    dy_dN = - (1/N**2) * (dN/(n*dx*domega))
    dy_dn = - (1/n**2) * (dN/(N*dx*domega))
    dy_ddx = - (1/dx**2) * (dN/(n*N*domega))
    dy_ddomega = - (1/domega**2) * (dN/(n*N*dx))
    dy_ddN = 1 / (N*n*dx*domega)

    sigma_y = np.sqrt( (dy_dN * sigma_N)**2 + (dy_dn * sigma_n)**2 + (dy_ddx *
                sigma_dx)**2 + (dy_ddomega * sigma_domega)**2 + (dy_ddN *
                    sigma_dN)**2)

    dsdo = np.array([dsdo, sigma_y])

    return dsdo

def RutherfordcrossSection_theoretical(Z1, Z2, theta):
    E = 350 * 10**(-3) # MeV
    return (1.296*((Z1 * Z2)/(E))**(2)) * (1/(np.sin(theta/2))**4)

def harryplotter(theta, dsdo_C, dsdo_Au):

    # Meassured error
    theta_error = np.ones(np.size(theta)) * 0.1

    # Unpacking values
    dsdo_C_val = dsdo_C[0]
    dsdo_C_error = dsdo_C[1]
    dsdo_Au_val = dsdo_Au[0]
    dsdo_Au_error =dsdo_Au[1]

    print("ERRORSSSW")
    print(dsdo_Au_error)
    print(dsdo_C_error)

    z1 = 1
    zAu = 79
    zC  = 6

    #### plotting
    # Theoretical
    theta_rad =  np.arange(0.3491, np.pi, 0.01)
    theta_deg = theta_rad * (180/np.pi)

    dsdo_Au_the = RutherfordcrossSection_theoretical(z1, zAu, theta_rad)
    dsdo_C_the  = RutherfordcrossSection_theoretical(z1, zC, theta_rad)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_yscale("log")
    # Experimental data
    ax.errorbar(theta, dsdo_Au_val, xerr=theta_error, yerr=dsdo_Au_error,
    fmt="o", label="Data Au")
    ax.errorbar(theta, dsdo_C_val , xerr=theta_error, yerr=dsdo_C_error, fmt='o', label="Data C")
    # Theoretical plots

    ax.plot(theta_deg, dsdo_Au_the, label="Theoretical Au")
    ax.plot(theta_deg, dsdo_C_the,  label="Theoretical C")

    plt.title("Rutherford Crosssection")
    plt.xlabel(r'Scattering Angle $[\si{\degree}]$')
    plt.ylabel(r'Cross Section $[\si{\milli\barn\per\steradian}]$')
    plt.grid()
    plt.xticks(np.arange(0, 190, 15))
    plt.legend()
    #plt.savefig('rutherford_graph')

    return

def target_thickness():
    # Dataframe (data) + file_names (for itteration)
    df, file_names = Dataframe(data_thickness_dir)
    for name in file_names:
        print(name)
        # Choosing non-zero values
        df_data = df[name].loc[df[name].nonzero()]
        df_data = df_data.loc[100:]

        # Transforming to np.array (Might be unnecesarry)
        np_data = np.array([np.array(df_data.index), df_data.values])
        bins = np_data[0]
        counts = np_data[1]

    return


# # # # # # # # # # # # # # Function calls # # # # # # # # # # # # # # # # # #
# To avoid several function calls:
k0_val = 1.36
k0_error = 0.15
alpha_val = 0.80
alpha_error = 0.00473

alpha = np.array([alpha_val, alpha_error])
k0    = np.array([k0_val, k0_error])

# Determine k0
#k0 = k0()
#k0_val, k0_error = [k0[0], k0[1]]

# Determine alpha
#alpha = alpha(k0_val)
#alpha_val, alpha_error = [alpha[0], alpha[1]]

# Determing gaussian parameters
plot_switch = False
detector_angles, A1, mean1, std1, A2, mean2, std2, dN_C, dN_Au =\
doublegaussian_fit(plot_switch)

# Convertion of Channel numbers to Energies (calibration)
E_C  = channel_to_energy(alpha, k0, mean1) # [Val, Error]
E_Au = channel_to_energy(alpha, k0, mean2) # [Val, Error]

E_C_val    = E_C[0]
E_C_error  = E_C[1]
E_Au_val   = E_Au[0]
E_Au_error = E_Au[1]

print("The energies of carbon is:")
print(E_C[0])
print("The energies of Gold is:")
print(E_Au[0])

# Plotting energy angle dependency

energy_angle_plot(detector_angles, E_Au, E_C)

# Rutherford cross sections
dsdo_C  = RutherfordCrossSection(dN_C , False)
dsdo_Au = RutherfordCrossSection(dN_Au, True)

harryplotter(detector_angles, dsdo_C, dsdo_Au)










plt.show()



