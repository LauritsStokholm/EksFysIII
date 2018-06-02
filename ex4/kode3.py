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
plt.rc('font', **{'family' : "sans-serif"})



""" Experimental setup 3: Mach Zehnder Interferometer"""


# Importing Data

# Current working directory
data_path = os.path.join(os.getcwd(), 'DataII')

# List of files in ./Data
data_files = os.listdir(data_path)

# Picking out .csv files (datafiles)
data = [x for x in data_files if '.csv' in x]

# Paths for chosen datafiles
data_dir = [os.path.join(data_path, x) for x in data]

# Picking out a file to work on
print(len(data_dir))
data_file = data_dir[0]


# Data Analysis
''' This is the data analysis per datafile, defined for determining the
material constant'''

# Reading csv
df = pd.read_csv(data_file, skiprows=[0, 1, 2], names=['Time', 'ChannelA'])

## Slicing data to seperate air flow in/out
df['Data1'] = df.ChannelA.where(df.ChannelA.index <
        df.ChannelA.size/2)
df['Data2'] = df.ChannelA.where(df.ChannelA.index >
        df.ChannelA.size/2)

plt.figure()
plt.title('Raw data')
df.Data1.plot()
df.Data2.plot()
plt.xlabel('Index')
plt.ylabel('Voltage [V]')
plt.grid()
plt.savefig('MZ.png')

m_in  = [22, 23, 23, 23, 24, 23, 22, 22, 21, 21, 21, 24, 24, 23, 25]
m_out = [24, 25, 25, 25, 25, 25, 25, 25, 23, 24, 25, 26, 26, 24, 26]

m_in_avg = np.mean(m_in)
m_out_avg = np.mean(m_out)

m = (m_in_avg + m_out_avg) / 2
print(m)

lambda_laser = 633 * 10**-9
length_glass = 6.00 * 10**-2

del_n = m*lambda_laser/length_glass
print(del_n)
plt.show()
