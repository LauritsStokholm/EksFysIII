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
#plt.rc('text',usetex =True)
plt.rc('font', **{'family' : "sans-serif"})



""" Experimental setup 2: Michelson Interferometer"""


# # # # # # # # # # # # # # # # Importing Data # # # # # # # # # # # # # # # # 

# Current working directory
data_path = os.path.join(os.getcwd(), 'Data')

# List of files in ./Data
data_files = os.listdir(data_path)

# Picking out .csv files (datafiles)
data = [x for x in data_files if '.csv' in x]

# Paths for chosen datafiles
data_dir = [os.path.join(data_path, x) for x in data]


# # # # # # # # # # # # # # # # Define Analysis # # # # # # # # # # # # # # # # 

def myFunc(data_dir):
    ''' This is the data analysis per datafile, defined for determining the
    Piezo Electric constant, C '''

    # Hello data:
    print('The current data is ' +
            str(os.path.basename(data_dir)))

    # Reading csv
    df = pd.read_csv(data_dir, skiprows=[0, 1, 2], 
                     names=['Time', 'ChannelA', 'ChannelB'])

    # First plot: Raw data
    plt.figure()
    plt.title('Raw Data')
    df.ChannelA.plot()
    df.ChannelB.plot()
    plt.xlabel('Index')
    plt.ylabel('Voltage [V]')
    plt.legend()
    plt.grid()
    plt.savefig('raw')

    # Unit Conversion
    df.ChannelA = df.ChannelA.apply(lambda x: x*10**-3)
    #df.CHannelB = df.ChannelB.apply(lambda x: x*10**-3)
    df.Time = df.Time.apply(lambda x: x*10**-3)

    # Determining index of extremas
    min_index = df.ChannelB.idxmin()
    max_index = df.ChannelB.idxmax()

    # Relations between indices (for appropriate plotting)
    low_index = min(min_index, max_index)
    high_index = max(min_index, max_index)
    
    #print('The extremas are')
    #print('Minima ' + str(min_index))
    #print('Maxima ' + str(max_index))
    
    # Filtering data to first interval
    v_piezo   = df.ChannelB[low_index : high_index+1]
    intensity = df.ChannelA[low_index : high_index+1]
    
    # Smoothening intensity data by method
    df['intensity_hat'] = pd.Series()
    df['intensity_hat'][low_index : high_index+1]=\
    pd.Series(sp.signal.savgol_filter(intensity, 51, 3))
    df.fillna(0, inplace=True)
    intensity_hat = df.intensity_hat[low_index : high_index+1]

    
    # Determining the wavelenght of V piezo

    # Defining orders for each data file
    # (how many points to consider on each side of extrema)

    if '9_3' in data_dir:
        order_vals = 1000
    elif '1_0' in data_dir:
        order_vals = 1200
    elif '3_8' in data_dir:
        order_vals = 1000
    elif '2_4' in data_dir:
        order_vals = 500
    elif '1_7' in data_dir:
        order_vals = 1000
    else:
        order_vals = 500
    
    # Determining relative extrema
    maximas = sp.signal.argrelextrema(np.array(intensity_hat),\
                                      np.greater, order=order_vals)
    # Returns tuple, but we have only 1d array, so picking relevant data
    maximas = maximas[0]


    print('this is maximas ' + str(maximas))

    if '9_3' in data_dir:
        print('9_3 way')
        max_index1_hat = maximas[1]
        max_index2_hat = maximas[2]

    elif '2_4' in data_dir:
        print('2_4 way')
        max_index1_hat = maximas[0]
        max_index2_hat = maximas[2]
    elif '3_8' in data_dir:
        max_index1_hat = maximas[2]
        max_index2_hat = maximas[3]
    elif '7_6' in data_dir:
        max_index1_hat = maximas[1]
        max_index2_hat = maximas[2]
    elif '4_5' in data_dir:
        max_index1_hat = maximas[1]
        max_index2_hat = maximas[2]
    elif '5_2' in data_dir:
        max_index1_hat = maximas[1]
        max_index2_hat = maximas[2]
    elif '6_8' in data_dir:
        max_index1_hat = maximas[1]
        max_index2_hat = maximas[2]
    else:
        max_index1_hat = maximas[0]
        max_index2_hat = maximas[1]
    
    #print('this is the relevant indices')
    #print(max_index1_hat + low_index)
    #print(max_index2_hat + low_index)
    
    # Visualisation of chosen indices
    idx1 = max_index1_hat + low_index
    idx2 = max_index2_hat + low_index
    
    plt.figure()
    plt.title('Data for ' + r"%s" %(os.path.basename(data_dir)))
    plt.xlabel('Voltage [V]')
    plt.ylabel('Channel A [V]')
    intensity.plot()
    intensity_hat[max_index1_hat:max_index2_hat].plot()
    plt.plot([idx1, idx2],
             [intensity_hat[idx1], intensity_hat[idx2]],
             marker='o', markersize=10, linestyle=' ')
    plt.grid()
    
    # Corresponding piezo voltages
    # dataframe corresponding index
    V1 = v_piezo[idx1]
    V2 = v_piezo[idx2]

    
    # Difference of piezo voltages (wavelength so to say)
    # and laser wavelength (We used a red HeNe-Laser)
    lambda_signal = abs(V2 - V1)*10
    lambda_laser = 633 * 10**-9
    
    C = lambda_laser / (2 * lambda_signal)
    plt.savefig('plot' + str(os.path.basename(data_dir) + '.png'))
    plt.show()
    return(C)
    

# Piezo constants
C = [x for x in range(len(data_dir))]
C = []

for i in range(len(data_dir)):
    if any(s in data_dir[i] for s in ('1_0', '1_7', '2_4')):
        pass
    else:
        C.append(myFunc(data_dir[i]))

print(C)

C_mean = np.average(np.array(C))
C_std = np.std(C)
print('The average of all determined Piezo constants are ' + str(C_mean))
print('Uncertainty is ' + str(C_std))
