import numpy as np
import matplotlib.pyplot as plt

h = 1000.0 #altitude at apogee in feet
theta = 10.0 #launch angle from vertical in degrees
g = 32.2 #gravity on Earth in Eng units

##Compute initial velocity
v = np.sqrt(2*h*g)/np.cos(theta*np.pi/180.0)
print("Delta V Required (ft/s) = ",v)
##Compute time to apogee
ta = v*np.cos(theta*np.pi/180.0)/g
print('Time to Apogee (sec) = ',ta)
##Compute Trajectory
t = np.linspace(0,2*ta,1000)
x = v*np.sin(theta*np.pi/180.0)*t
y = v*np.cos(theta*np.pi/180.0)*t - 0.5*g*t**2
plt.plot(x,y)
plt.xlabel('X (ft)')
plt.ylabel('Y (ft)')
plt.grid()
plt.show()