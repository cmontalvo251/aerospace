#!/usr/bin/python

import numpy as np
import os
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


#Change directory to Satellite Files
os.chdir(os.path.expanduser("~/Dropbox/Code/FORTRAN/Satellites"))

#Run make rebuild
os.system('make rebuild')

#Run Satellite Program
os.system('./Satellites')

#1st run data
test = np.loadtxt('SatelliteResults.txt')

#Store run data
time = test[:,0]
x1 = test[:,1]
y1 = test[:,2]
z1 = test[:,3]
phi1 = test[:,4]
theta1 = test[:,5]
psi1 = test[:,6]
u1 = test[:,7]
v1 = test[:,8]
w1 = test[:,9]
x2 = test[:,13]
y2 = test[:,14]
z2 = test[:,15]
phi2 = test[:,16]
theta2 = test[:,17]
psi2 = test[:,18]
u2 = test[:,19]
v2 = test[:,20]
w2 = test[:,21]

#TURN ON INTERACTIVE MODE
plt.ion() 

# 3d Plot
fig = plt.figure('3-D')
ax = fig.add_subplot(111,projection='3d')
ax.plot_wireframe(x1,y1,z1, color = 'black', linestyle = 'solid', label = "CubeSat 1")
ax.plot_wireframe(x2,y2,z2, color = 'black', linestyle = 'dashed', label = 'CubeSat 2')
ax.invert_yaxis()
ax.invert_zaxis()
plt.title('3-d')
plt.legend()

# X Plot
plt.figure('X')
plt.plot(time,x1, 'k', label = "CubeSat 1")
plt.plot(time,x2, 'k--', label = "CubeSat 2")
plt.title('X (m)')
plt.xlabel('Time (sec)')
plt.ylabel('X (m)') 
plt.legend(loc = 2)

# Y Plot
plt.figure('Y')
plt.plot(time,y1, 'k', label = 'CubeSat 1')
plt.plot(time,y2, 'k--', label = 'CubeSat 2')
plt.title('Y (m)')
plt.xlabel('Time (sec)')
plt.ylabel('Y (m)')
plt.legend(loc = 2)

# Z Plot
plt.figure('Z')
plt.plot(time,z1, 'k', label = 'CubeSat 1')
plt.plot(time,z2, 'k--', label = 'CubeSat 2')
plt.title('Z (m)')
plt.xlabel('Time (sec)')
plt.ylabel('Z (m)')
plt.legend(loc = 2)

# U Plot
plt.figure('U')
plt.plot(time,u1, 'k', label = 'CubeSat 1')
plt.plot(time,u2, 'k--', label = 'CubeSat 2')
plt.title('U (m/s)')
plt.xlabel('Time (sec)')
plt.ylabel('U (m/s)')
plt.legend(loc = 2)


# V Plot
plt.figure('V')
plt.plot(time,v1, 'k', label = 'CubeSat 1')
plt.plot(time,v2, 'k--', label = 'CubeSat 2')
plt.title('V (m/s)')
plt.xlabel('Time (sec)')
plt.ylabel('V (m/s)')
plt.legend(loc = 2)


# W Plot
plt.figure('W')
plt.plot(time,w1, 'k', label = 'CubeSat 1')
plt.plot(time,w2, 'k--', label = 'CubeSat 2')
plt.title('W (m/s)')
plt.xlabel('Time (sec)')
plt.ylabel('W (m/s)')
plt.legend(loc = 2)

#Phi Plot
plt.figure('Phi')
plt.plot(time,phi1, 'k', label = 'CubeSat 1')
plt.plot(time,phi2, 'k--', label = 'CubeSat 2')
plt.title('Phi')
plt.xlabel('Time (sec)')
plt.ylabel('Phi')
plt.legend()


#Theta Plot
plt.figure('Theta')
#plt.subplot(322)
plt.plot(time,theta1, 'k', label = 'CubeSat 1')
plt.plot(time,theta2, 'k--', label = 'CubeSat 2')
plt.title('Theta')
plt.xlabel('Time (sec)')
plt.ylabel('Theta')
plt.legend()


#Psi Plot
plt.figure('Psi')
plt.plot(time,psi1, 'k', label = 'CubeSat 1')
plt.plot(time,psi2, 'k--', label = 'CubeSat 2')
plt.title('Psi')
plt.xlabel('Time (sec)')
plt.ylabel('Psi')
plt.legend(loc = 5)

# Print all plots to screen
plt.show()

#Hang Code until user hits ENTER
raw_input('Hit ENTER to close all windows')
plt.close('all')
