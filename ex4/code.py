# # # # # # # # # # # # # # # Loading Libraries # # # # # # # # # # # # # # # # 
# Used for basic commands
import os

# Data Analysis and basic mathematics
import pandas as pd
import numpy as np
import scipy as sp
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
data     = [x for x in data_files if 'data' in x and 'csv' in x]

# Paths for chosen datafiles
data_cal_dir = [os.path.join(data_path, x) for x in data_cal]
data_dir     = [os.path.join(data_path, x) for x in data]

# # # # # # # # # # # # # # # # Calibration # # # # # # # # # # # # # # # # # # 
# Preparing dataframe structure
df = pd.DataFrame()

file_names = []
# Itterating over all files
for i in range(len(data_cal_dir)):
    # file name
    file_names.append(os.path.basename(data_cal_dir[i]))

    # Reading csv
    column_names = [file_names[i]]
    df1 = pd.read_csv(data_cal_dir[i], delimiter=',', index_col=0, names=column_names)

    # Plot Calibration
    df = pd.concat([df, df1[column_names[0]]], axis=1)

# Plotting non zero values
for name in file_names:
    df[name].loc[df[name].nonzero()].plot()
    plt.grid()
    plt.title("Calibration")
    plt.xlabel("Channel \#")
    #plt.legend()



plt.show()



## Making two dataframes, one with integer indices and the other with lambda
## integers. (We need integer index for scipy, but lambda is better for plot)
#df_intidx = pd.concat([df1[column_names[0]], df], axis=1)
#df_lamidx = df_intidx.copy()
#df_lamidx.set_index(df_lamidx[column_names[0]], inplace=True)
#df_lamidx.drop([column_names[0]], axis=1, inplace=True)
#
#
#
## # # # # # # # Part A ~ Light bulbs and the solar spectrum # # # # # # # # # #
#
## Determining the temperature
#def Temperature(df):
#    ''' This function will determine the temperature of the given spectra '''
#
#    # Determine peak of wavelength
#    lambda_peak = df.idxmax()
#    print(lambda_peak)
#
#    # Determine the temperature
#    constant = 2.898 * 10**-3 # m*K
#    temperature = constant / lambda_peak
#    return temperature
#
## Sorting relevant files and calculating temperature
#T = []
#for column in df_lamidx:
#    if any(s for s in df_lamidx if column in ('solar1.csv', 'solar2.csv',
#        'solar3.csv', 'notseethrough.csv', 'seethrough.csv', 'curly.csv')):
#        T.append(column)
#        T.append(Temperature(df_lamidx[column]))
#
#print('The temperature is given by:')
#print(T)
#print(np.average(np.array([T[1], T[-1], T[-3]])))
#
#
## Setting better indices for plot
#plot_index = pd.Series(df_lamidx.index)
#df_lamidx.set_index(plot_index.apply(lambda x: x*10**9), inplace=True)
#
## Compare results with solar spectra
#i = 0
#for item in data_lamps:
#    # the file name
#    lamp = os.path.basename(item)
#
#    fig = plt.figure()
#    ax = fig.add_subplot(111)
#    ax.set_title('Solar spectrum')
#    df_lamidx['solar1.csv'].plot()
#    df_lamidx['solar2.csv'].plot()
#    df_lamidx['solar3.csv'].plot()
#    df_lamidx[lamp].plot()
#    ax.grid()
#    ax.legend()
#    ax.set_xlabel(r'Wavelength [\si{\nano\meter}]')
#    ax.set_ylabel(r'Intensity')
#    plt.savefig("SolarComparison"+str(i))
#
#    # Do not run out of memory
#    #plt.close(fig)
#    i += 1
#
#
## Determine elements in solar spectrum
#
## Finding extrema
#x = np.array(df_intidx['solar2.csv'], dtype=np.float)
#order_val = 40
#list_of_peaks = sp.signal.argrelextrema(x, np.less,
#        order=order_val)[0].tolist()
#list_of_peaks = [x for x in list_of_peaks if 295<x and x<1500]
#
#
#x_val = list_of_peaks
#y_val = df_intidx['solar2.csv'][list_of_peaks]
#labels = ['Fe', 'Fe/Ca', 'Mg/Fe', 'Na', 'O_2', 'O_2', 'O_2', 'O_2']
#
## Making plot with labels at scatterplot
#plt.figure()
#plt.title('Fraunhofer lines')
#plt.scatter(x_val, y_val, marker='o', color='red')
#
## For name in labels
#for label, x, y in zip(labels, x_val, y_val):
#    plt.annotate(label,
#		xy=(x, y), xytext=(-20, 20),
#        textcoords='offset points', ha='right', va='bottom',
#        bbox=dict(boxstyle='round,pad=0.5', fc='red', alpha=0.5),
#        arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'))
#df_intidx['solar2.csv'].plot()
#plt.grid()
#plt.legend()
#plt.xlabel(r"Wavelength [\si{\nano\meter}]")
#plt.ylabel(r"Intensity")
#
#print('Interesting wavelengths are')
#print(df_intidx.wavelength[list_of_peaks] * 10**9)
#plt.savefig('Fraunhofer')
#
#
## # # # # # # # # # Part B ~ Spectral Lamps and Lasers # # # # # # # # # # # #
#
#
#
## # # # # # # # # # Part C ~ Absorbance of three cuvette# # # # # # # # # # # #
##
#for item in data_absorbance:
#    cuvette = os.path.basename(item)
#    print(item)
#    print(cuvette)
#    plt.figure()
#    df_lamidx[cuvette].plot()
#    plt.grid()
#    plt.title('Spectrum of cuvette ' + str(cuvette[0]))
#    plt.legend()
#    plt.xlabel(r'Wavelength [\si{\nano\meter}]')
#    plt.ylabel(r'Intensity')
#    plt.xlim(180, 650)
#    plt.savefig('abs'+cuvette[0])
#
##'plt.show()
##data_C = []
##for item in data_absorbance:
##    data_C.append(os.path.basename(item))
##print(data_C)
#
##for item in data_C:
##    plt.figure()
##    df_lamidx[item].plot()
##    plt.grid()
##    plt.legend()
##    plt.xlabel(r'Wavelength [\si{\nano\meter}]')
##    plt.ylabel(r'Intensity')
##    plt.xlim(180, 650)
#
#
