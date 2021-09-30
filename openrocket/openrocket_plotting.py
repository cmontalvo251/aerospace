import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('In_Class_OpenRocket_CLEAN.csv',delimiter=',')

time = data[:,0]
altitude = data[:,1]
aoa = data[:,4]

plt.plot(time,altitude)

plt.show()