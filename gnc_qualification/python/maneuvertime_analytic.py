#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 12 15:51:27 2021

@author: carolinefranklin
"""

import numpy as np
import matplotlib.pyplot as plt


#variables,values from RWP500 RWS
power = 6.0 #(W)
f = np.linspace(0,0.5,100) #factor is arbitrary and did not come from GNC spreadsheet
max_torque = 0.025 #(Nm)
angleinput = input('Enter angle:')
angle = float(angleinput)
final_angle = angle*np.pi/180.0

#x axis plot
inertiax = 0.385 #inertia from ABEX data file
angular_accelerationx = max_torque/inertiax #angular acceleration
t1x = np.sqrt(2*final_angle*f/angular_accelerationx) #time of first maneuver
energyx = 2*power*t1x #energy (Wh)
timex = 2*t1x + final_angle*(1.0-2.0*f)/(angular_accelerationx*t1x) #total time
plt.plot(timex,energyx) 


#y axis plot
inertiay = 0.431 #inertia from ABEX data file
angular_accelerationy = max_torque/inertiay 
t1y = np.sqrt(2*final_angle*f/angular_accelerationy)
energyy = 2*power*t1y
timey = 2*t1y + final_angle*(1.0-2.0*f)/(angular_accelerationy*t1y)
plt.plot(timey,energyy) 


#z axis plot
inertiaz = 0.296 #inertia from ABEX data file
angular_accelerationz = max_torque/inertiaz
t1z = np.sqrt(2*final_angle*f/angular_accelerationz)
energyz = 2*power*t1z
timez = 2*t1z + final_angle*(1.0-2.0*f)/(angular_accelerationz*t1z)
plt.plot(timez,energyz)

plt.xlabel('Time (s)')
plt.ylabel('Energy (W-h)')
plt.grid('both')
plt.title('Power Required to Move {}°'.format(angle))
plt.legend(['X axis','Y axis','Z axis'])
plt.show()