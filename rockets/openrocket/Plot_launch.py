# -*- coding: utf-8 -*-
"""
Created on Sat Jun  4 12:22:51 2022

@author: Justin

LAUNCH DATA
CPB + MS8607 
"""

import numpy as np
import matplotlib.pyplot as plt


#Load Data
data2 = np.loadtxt('Launch_t_axyz_p_rH_Te_T.txt')


#Assign variables to data
time2 = data2[:,0]
ax = data2[:,1]
ay = data2[:,2]
az = data2[:,3]
p = data2[:,4]
rH = data2[:,5]
Te = data2[:,6]
T2 = data2[:,7]


#Set launch time
tlaunch2 = 1294


#Adjust data to launch time
ax = ax[time2>tlaunch2]
ay = ay[time2>tlaunch2]
az = az[time2>tlaunch2]
p = p[time2>tlaunch2]
rH = rH[time2>tlaunch2]
Te = Te[time2>tlaunch2]
T2 = T2[time2>tlaunch2]
time2 = time2[time2>tlaunch2]


## Set launch land time
tland = 1550
tland2 = 1494


#Adjust data to land time
ax = ax[time2<tland2]
ay = ay[time2<tland2]
az = az[time2<tland2]
p = p[time2<tland2]
rH = rH[time2<tland2]
Te = Te[time2<tland2]
T2 = T2[time2<tland2]
time2 = time2[time2<tland2]
time2 -= time2[0]


####PLOT RAW DATA


#Constants
R = 8.31432 #Universal gas constant (Nm)/(mol*K)
g_0 = 9.80665 #gravity accel constant m/s**2
p = p #pressure data hPa
p_0 = p[0] #initial pressure hPa
T_0 = Te[0]+273.15 #initial temp K
z_0 = 0 #initial alt m
L_b = -0.0065 #standard temp lapse rate K/m
M = 0.0289644 #molar mass of air kg/mol


## Apply international barometric formula to get altitude
## SEE THEORY
alt = z_0+(T_0/L_b)*((p/p_0)**((-R*L_b)/(g_0*M))-1)


#Close all figs
plt.close('all')


###PLOT RAW Temperature2
plt.figure()
plt.plot(time2,T2, 'b-')
plt.xlabel('Time (seconds)')
plt.ylabel('Temperature (Celcius)')
plt.title('Raw Temperature 2')
plt.grid()


###PLOT RAW External Temperature
plt.figure()
plt.plot(time2,Te, 'b-')
plt.xlabel('Time (seconds)')
plt.ylabel('Temperature (Celcius)')
plt.title('Raw External Temperature')
plt.grid()


###PLOT RAW Internal ACCEL2
plt.figure()
plt.plot(time2,ax, 'b-', label = 'X Axis')
plt.plot(time2,ay, 'r-', label = 'Y Axis')
plt.plot(time2,az, 'g-', label = 'Z Axis')
plt.xlabel('Time (seconds)')
plt.ylabel('Acceleration (m/s^2)')
plt.title('Raw Internal Acceleration2')
plt.grid()
plt.legend()


###PLOT RAW Relative Humidity
plt.figure()
plt.plot(time2,rH, 'b-')
plt.xlabel('Time (seconds)')
plt.ylabel('Relative Humidity (%)')
plt.title('Raw Relative Humidity')
plt.grid()

###PLOT RAW Pressure
plt.figure()
plt.plot(time2,p, 'b-')
plt.xlabel('Time (seconds)')
plt.ylabel('Pressure (hPa)')
plt.title('Raw Pressure')
plt.grid()



###PLOT altitude from pressure
plt.figure()
plt.plot(time2,alt,'b-')
plt.xlabel('Time (seconds)')
plt.ylabel('Altitude (m)')
plt.title('Altitude from Pressure Data (SI)')
plt.grid()


plt.figure()
plt.plot(time2,alt*3.28084,'b-')
plt.xlabel('Time (seconds)')
plt.ylabel('Altitude (ft)')
plt.title('Altitude from Pressure Data (English)')
plt.grid()


### Take the derivative of Altitude to get vertiacal velocity
velo = 0*alt
for i in range(0,len(alt)-1):
    velo[i+1] = ((alt[i+1]-alt[i])/(time2[i+1]-time2[i]))


### Plot vertical velocity
plt.figure()
plt.plot(time2,velo,'b-')
plt.xlabel('Time (seconds)')
plt.ylabel('Velocity (m/s)')
plt.title('Vertical Velocity from Pressure Data')
plt.grid()


### Take the derivative of vertiacal velocity to get vertiacal 
##acceleration
accel = 0*velo
for i in range(0,len(velo)-1):
    accel[i+1] = ((velo[i+1]-velo[i])/(time2[i+1]-time2[i]))


##Plot vertical acceleration
plt.figure()
plt.plot(time2,accel,'b-')
plt.xlabel('Time (seconds)')
plt.ylabel('Acceleration (m/s^2)')
plt.title('Vertical Acceleration from Pressure Data')
plt.grid()

