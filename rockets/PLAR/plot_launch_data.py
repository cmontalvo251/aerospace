import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


data = np.loadtxt('Launch_t_axyz_p_rH_Te_T.txt')

time = data[:,0]
l1 = np.where(time > 1280)[0][0]
l2 = np.where(time > 1500)[0][0]
time = time[l1:l2]

pressure = data[:,4]
pressure = pressure[l1:l2]

plt.figure()
plt.plot(time,pressure)

##CONVERT TO ALTITUDE
plt.figure()
pressure_pascals = pressure * 100.0
print(pressure_pascals[0])
altitude = (1.0-pow((pressure_pascals/101325.0),1.0/5.25588))/(2.2557*pow(10,-5.0)) - 71.0
plt.plot(time,altitude)

###Velocity
velocity = (altitude[1:] - altitude[0:-1])/(time[1:] - time[0:-1])
plt.figure()
plt.plot(time[0:-1],velocity)

##Acceleration
acceleration = (velocity[1:] - velocity[0:-1])/(time[1:-1] - time[0:-2])
plt.figure()
plt.plot(time[0:-2],acceleration)


plt.figure()
plt.plot(time[1:]-time[0:-1])

plt.show()