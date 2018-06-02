# Preamble
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['legend.numpoints']=1
plt.rcParams.update({'font.size': 18})
plt.rc('text',usetex=True)
plt.rc('font',family='serif')



# Funktion
def f(x):
    return x**2
x = np.arange(0, 10, 0.01)


plt.figure()
plt.plot(x, f(x))
plt.xlabel('x--values')
plt.ylabel('y--values')
plt.title('Vores figur')


plt.show()

