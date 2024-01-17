# -*- coding: utf-8 -*-
"""
Created on Fri Oct  6 10:06:13 2023

@author: colem
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as sci 

plt.close("all")

##Define Custom Parameters

##Planet CONSTANTS
Rplanet = 6357000 #meters
mplanet = 5.972e24 #kg
G = 6.6742*10**-11; #Gravitational Constant (SI Unit)

#Gravitational Acceleration Model
def gravity(z):
    global Rplanet,mplanet
    r = np.sqrt(z**2)
    if r< Rplanet:
        accel = 0.0
    else: 
        accel = G*mplanet/(r**3)*r 
        #print(accel) 
    return accel

def control_system(zdot,z):
    tp = zdot / 9.81 ##First order approximation
    zp = (z-Rplanet) + zdot * tp - 0.5*9.81*tp**2  ##First order approximation
    zc = 300
    error_altitude = zp - zc
    kp = 10.0/1000.0   ##THIS NEEDS TO BE ADJUSTED
    df = kp * error_altitude #Fin distance - variable from 0 - 1 0 = closed and 1 = fully deployed
    if df > 1:
        df = 1
    if df < 0:
        df = 0
    if zdot < 0:
        df = 0
    return df

#Second Order Diff Equations
def Derivatives(state,t):
    ##Globals
    z = state[0]
    velz = state[1]
    #compute zdot
    zdot = velz
    
    #mass
    mass = 1.0 ##LOOK UP MASS AFTER SECOND STAGE BURNOUT
    
    ###CONTROL SYSTEM######
    df = control_system(zdot,z)
    
    ###HIL####
    #Send the df command via serial to the fins
    
    ########FORCES############
    #Gravity
    gravityF = -gravity(z)*mass
    #Aerodynamics
    h = (z-Rplanet)/1000.0
    beta = 0.1354 #scale height
    rho = 1.225*np.exp(-beta*h)
    D = 0.102
    kappa = 0.5  ## NEED AERODYNAMIC MODELING FOR THIS 
    CD = 0.5 + df*kappa ##AND THIS
    aeroF = -np.pi/8.0*rho*zdot**2*D**2*CD
    #print(aeroF,CD,rho,z)
    if z < Rplanet:
        aeroF = 0.0 
    thrustF = 0
    Forces = gravityF + aeroF + thrustF
    
    #Compute zddot
    zddot = Forces/mass
    
    statedot = np.asarray([zdot,zddot])
    
    return statedot

####EVERYTHING BELOW HERE IS MAIN SCRIPT###

#Test Surface Gravity
print('Surface Gravity (m/s^2) = ', gravity(Rplanet))

#Initial Conditions
altitude_SSECO = 100.0 ##ALTITUDE AT SECOND STAGE CUTOFF
z0 = Rplanet+altitude_SSECO #m
##VELOCITY BELOW IS AT SECOND STAGE CUTOFF
velz0 = 100.0 #m/s)

##Time Window
dt = 0.1
tout = np.arange(0,15,dt)

##Numerical Integration Call
zout = 0*tout
velzout = 0*tout
dfout = 0*tout
zout[0] = z0
velzout[0] = velz0
for i in range(0,len(tout)-1):
    print('T = ',tout[i])
    state = np.array([zout[i],velzout[i]])
    dfout[i] = control_system(velzout[i],zout[i])
    t = tout[i]
    k1 = Derivatives(state,t)
    k2 = Derivatives(state + k1*dt/2.0,t+dt/2)
    k3 = Derivatives(state + k2*dt/2.0,t+dt/2)
    k4 = Derivatives(state + k3*dt,t+dt)
    phi = (1.0/6.0)*(k1 + 2*k2 + 2*k3 + k4)
    statenext = state + phi*dt
    zout[i+1] = statenext[0]
    velzout[i+1] = statenext[1]
    ###ADD A WAIT FUNCTION FOR HIL

##Rename Variables
altitude = zout - Rplanet
##Plot

#ALTITUDE
plt.plot(tout,altitude)
plt.xlabel('Time (sec)')
plt.ylabel('Altitude (m)')
plt.grid()

#VELOCITY
plt.figure()
plt.plot(tout,velzout)
plt.xlabel('Time (sec)')
plt.ylabel('Normal Speed (m/s)')
plt.grid()

##PLOT DF
plt.figure()
plt.plot(tout,dfout)
plt.xlabel('Time (sec)')
plt.ylabel('df (nd)')
plt.grid()

plt.show()

 
    
