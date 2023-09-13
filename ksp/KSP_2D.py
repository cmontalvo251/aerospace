#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May 28 07:22:57 2020

@author: carlos
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as I

###Closes all Figures
plt.close("all")

def atmosphere_model(altitude):
    global altx,deny
    return np.interp(altitude,altx,deny)    

def planet_parameters():
    ##https://wiki.kerbalspaceprogram.com/wiki/Kerbin
    ##A few things you can compute to compare to the wiki above
    #mu 
    #surface gravity
    #rotational velocity
    
    G = 6.6742*10**-11; #%%Gravitational constant
    Mkerbin = 5.2915158*10**22 #
    muKerbin = G*Mkerbin
    Rkerbin = 600000 #meters
    sidereal_period = 21549.425
    sidereal_angular_velocity = 2*np.pi/sidereal_period
    sidereal_rotational_velocity = sidereal_angular_velocity*Rkerbin
    surface_gravity = muKerbin*Rkerbin/Rkerbin**3
    return sidereal_rotational_velocity,muKerbin,surface_gravity,Rkerbin
   

#######THIS IS A FUNCTION CALLED A DEFINITION IN PYTHON
#z and t are input. so x = [z,V], x is a 2x1 array with position
# and velocity in it and t is just time.
# matrices are NxM arrays and vectors are Nx1 arrays. Arrays are
# just a way for the computer to handle multiple numbers in one variable
def Derivatives(state,t):

    x = state[0]
    z = state[1]
    velx = state[2]
    velz = state[3]
    mass = state[4]
    
    ###Kinematics
    xdot = velx
    zdot = velz
    
    ##Dynamics
    ## F = m*a
    
    #now let's compute Forces on rocket
    
    ##Gravitational Acceleration
    rSat = np.sqrt((x**2) +(z**2))
    #MNeed to get parameters of the planet
    sidereal_rotational_velocity,mu,surface_gravity,R = planet_parameters()
    if rSat < 1.0:
        gravx = 0.0
        gravz = 0.0
    else:
        gravx = -mu*x/(rSat**3)
        gravz = -mu*z/(rSat**3)

    ##Now let's do Aerodynamics 
    #https://wiki.kerbalspaceprogram.com/wiki/Atmosphere
    #https://wiki.kerbalspaceprogram.com/wiki/Kerbin#Atmosphere
    altitude = rSat-R
    if altitude < 0:
        rho = 1.0
    else:
        rho = atmosphere_model(altitude)
    #Need AeroGUI Mod to get these parameters
    #https://forum.kerbalspaceprogram.com/index.php?/topic/105524-105-aerogui-v30-14-nov/
    #https://www.youtube.com/watch?v=ASHRPo4sw80
    #A few problems here. First Cd is a function of Mach Number and Reynolds number
    #so......I think I'll just leave this off
    #Make sure that S is D^2
    qinf = -np.pi/8.0*rho*S*Cd/mass #This aero model was checked with ARL work in 2012
    aerox = qinf*abs(velx)*velx
    aeroz = qinf*abs(velz)*velz

    #And of course thrust
    if t > stage_1_time:
        thrustx = 0.0
        thrustz = 0.0
    else:
        if GNC:
            thrust = T1*1000.0
            theta = 90.0*altitude/apogee
            if theta > 90.0:
                theta = 90.0
            thrustx = thrust*np.cos(theta*np.pi/180.0)
            thrustz = thrust*np.sin(theta*np.pi/180.0)
        else:
            thrustx = T1*1000.0
            thrustz = 0.0 
                
    #Fire stage 2
    if t > stage_2_start and t < stage_2_end:
        thrust = T2*1000.0
        ##What angle?
        ##Tangent to the current x,z coordinate
        #First get the angle to the current coordinate
        thetaSat = np.arctan2(z,x)            
        ##Then compute thrust based on that angle
        thrustx = thrust*np.sin(thetaSat)
        thrustz = thrust*np.cos(thetaSat)
            
    if mass < 1.0:
        thrustx = 0.0
        thrustz = 0.0
        
    thrust = np.sqrt(thrustx**2 + thrustz**2)
    if abs(thrust) > 0:
        #print('thrust',thrust)
        #But when thrust is fired we lose mass
        ve = Isp*abs(surface_gravity)
        mdot = -thrust/ve
    else:
        mdot = 0.0
    
    ##Now we can put Newton's EOMs together
    xdbldot = thrustx/mass + gravx + aerox
    zdbldot = thrustz/mass + gravz + aeroz
    
    statedot = np.asarray([xdot,zdot,xdbldot,zdbldot,mdot])
    if altitude < 0:
        statedot = np.asarray([0,0,0,0,0])
    
    #make [xdot,xdbldot] an array
    return statedot

##############END OF FUNCTION SEPARATED BY TABS#######
    
#Read in Atmospher models
atm_model = np.loadtxt('kerbin_atmosphere.txt')
altx = atm_model[:,0]
deny = atm_model[:,3]
#And planet parameters
sidereal_rotational_velocity,mu,surface_gravity,R = planet_parameters()

"""
##Initial Conditions for hop Orbit
x0 = R
z0 = 0.
velx0 = 0.0
velz0 = 0.0
masstons = 2.4
T1 = 38.28
Isp = 140.
Cd = 0.4
S = 0.18 #This needs to be D^2
stage_1_time = 38.0
stage_2_start = -99
stage_2_end = -99
period = 100.0
GNC = 0
apogee = 10000.
"""


###Conditions for Sub orbital flight

x0 = R
z0 = 0.
velx0 = 0.0
velz0 = 0.0
masstons = 8.0
mass_end = 3.0
stage_1_time = 63.0
stage_2_start = -99
stage_2_end = -99
T1 = 205.16
Isp = 265.
Cd = 0.23
D = 1.25
S = np.pi*D**2/4.0
period = 1000.0
GNC = 0
apogee = 80000.

"""
###Let's make an orbit
desired_orbit_altitude_km = 70. #kilometers
r = R + desired_orbit_altitude_km*1000.
vorbit = np.sqrt(mu/r)
##Vorbit for 70 km is 2,295.9 m/s <- this is the delta v we 
#need to get to orbit
##Initial conditions in orbit
x0 = r
z0 = 0.
velx0 = 0.0
velz0 = vorbit
T1 = 0.0
Isp = 0.0
Cd = 0.0
S = 0.0
masstons = 5.3
stage_1_time = -99
stage_2_start = -99
stage_2_end = -99
#Orbit Time
semi_major = r
period = 2*np.pi/np.sqrt(mu)*semi_major**(3.0/2.0)
GNC = 0
"""

"""
###Two Stage Rocket
x0 = R
z0 = 0.0
velx0 = 0.0
velz0 = 0.0
T1 = 180.0
T2 = 160.0
Isp = 250.
Cd = 0.1
S = 0.01
masstons = 5.3+2.0
stage_1_time = 38.0
stage_2_start = 150.0
stage_2_end = stage_2_start + 35.0
period = 3000.0
GNC = 1
apogee = 70000.
"""

###################################################################
mass0 = masstons*2000/2.2
tout = np.linspace(0,period,100000)  #linspace(start,end,number of data points)
stateinitial = np.asarray([x0,z0,velx0,velz0,mass0])
stateout = I.odeint(Derivatives,stateinitial,tout) ##This is the ode toolbox from scipy (Scientific Python)

xout = stateout[:,0]
zout = stateout[:,1]
vxout = stateout[:,2]
vzout = stateout[:,3]
mout = stateout[:,4]

vout = np.sqrt(vxout**2+vzout**2)

plt.figure()
plt.plot(tout,vout)
plt.grid()
plt.xlabel('Time (sec)')
plt.ylabel('Velocity (m/s)')

plt.figure()
plt.plot(xout,zout)
thetakerbal = np.linspace(0,2*np.pi,1000)
xkerbal = R*np.cos(thetakerbal)
zkerbal = R*np.sin(thetakerbal)
plt.plot(xkerbal,zkerbal)
plt.grid()
plt.xlabel('X (m)')
plt.ylabel('Z (m)')

plt.figure()
plt.plot(tout,zout)
plt.grid()
plt.xlabel('Time (sec)')
plt.ylabel('Z (m)')

plt.figure()
plt.plot(tout,np.sqrt(xout**2+zout**2)-R)
plt.grid()
plt.xlabel('Time (sec)')
plt.ylabel('AGL (m)')

plt.figure()
plt.plot(tout,mout)
plt.grid()
plt.xlabel('Time (sec)')
plt.ylabel('Mass (kg)')

plt.show()

