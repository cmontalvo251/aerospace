#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May 28 07:22:57 2020

@author: carlos
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as I
import pandas as pd

#data = pd.read_csv('First_Open_Rocket_Design.csv')

###Closes all Figures
plt.close("all")
#
def atmosphere_model(altitude):
#    #global altx,deny
    return 1.225*np.exp(-altitude/1000.0*0.1354)

def planet_parameters():
    ##https://wiki.kerbalspaceprogram.com/wiki/Kerbin
    ##A few things you can compute to compare to the wiki above
    #mu 
    #surface gravity
    #rotational velocity
    
    G = 6.6742*10**-11; #%%Gravitational constant
    MEarth = 5.972*10**24
    muEarth = G*MEarth
    REarth = 6371000 #meters
    sidereal_period = 21549.425
    sidereal_angular_velocity = 2*np.pi/sidereal_period
    sidereal_rotational_velocity = sidereal_angular_velocity*REarth
    surface_gravity = muEarth*REarth/REarth**3
    return sidereal_rotational_velocity,muEarth,surface_gravity,REarth
   

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
    gravx = -mu*x/(rSat**3)
    gravz = -mu*z/(rSat**3)
    
    ##Now let's do Aerodynamics 
    #https://wiki.kerbalspaceprogram.com/wiki/Atmosphere
    #https://wiki.kerbalspaceprogram.com/wiki/Kerbin#Atmosphere
    altitude = rSat-R
    if altitude < 0:
        rho = atmosphere_model(0.0)
        gravx = -surface_gravity
        gravz = 0.0
    else:
        rho = atmosphere_model(altitude)
    #Need AeroGUI Mod to get these parameters
    #https://forum.kerbalspaceprogram.com/index.php?/topic/105524-105-aerogui-v30-14-nov/
    #https://www.youtube.com/watch?v=ASHRPo4sw80
    #A few problems here. First Cd is a function of Mach Number and Reynolds number
    #so......I think I'll just leave this off
    if t > t_parachute:
        Cdi = Cd_parachute
        Si = S_parachute
    else:
        Cdi = Cd
        Si = S
    qinf = -np.pi/8.0*rho*Si*Cdi/mass
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
            
    if mass < 200./1000.:
        print('Mass < 0')
        thrustx = 0.0
        thrustz = 0.0
        
    thrust = np.sqrt(thrustx**2 + thrustz**2)
    if abs(thrust) > 0:
        #print('thrust',thrust)
        #But when thrust is fired we lose mass
        ve = Isp*abs(surface_gravity)
        mdot = -thrust/ve
        #print('mdot = ',mdot)
    else:
        #print('Mdot = 0')
        mdot = 0.0
    
    ##Now we can put Newton's EOMs together
    xdbldot = thrustx/mass + gravx + aerox
    zdbldot = thrustz/mass + gravz + aeroz
    
    statedot = np.asarray([xdot,zdot,xdbldot,zdbldot,mdot])
    #if altitude < 0:
    #    statedot = np.asarray([0.,0.,0.,0.,0.])
    
    #make [xdot,xdbldot] an array
    #print('mass = ',mass,'time = ',t,'atitude = ',altitude)
    #print('gx,gz = ',gravx,gravz,' t = ',t)
    return statedot

##############END OF FUNCTION SEPARATED BY TABS#######
    
#Read in Atmospher models
#atm_model = np.loadtxt('kerbin_atmosphere.txt')
#altx = atm_model[:,0]
#deny = atm_model[:,3]
#And planet parameters
sidereal_rotational_velocity,mu,surface_gravity,R = planet_parameters()


##Initial Conditions for open rocket
x0 = R
z0 = 0.
velx0 = 0.0
velz0 = 0.0
mass0 = 1.6 #kg
T1 = 175/1000.0  
massflowrate = (168./1000.)/2.4
exit_velocity = T1*1000 / massflowrate
Isp = exit_velocity / 9.81
print('Isp = ',Isp)
Cd = 0.75
Cd_parachute = 13.0
t_parachute = 9.0
D = 7.62/100.
S = D**2
S_parachute = (45./100.0)**2
stage_1_time = 2.4
stage_2_start = -99
stage_2_end = -99
period = 75.0
GNC = 0
apogee = 10000.

###################################################################
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
plt.plot(tout,np.sqrt(xout**2+zout**2)-R,label='Python Simulation')
#plt.plot(data['time'],data['altitude'],label='OpenRocket Simulation')
plt.grid()
plt.legend()
plt.xlabel('Time (sec)')
plt.ylabel('AGL (m)')

plt.figure()
plt.plot(tout,mout,label='Python Simulation')
#plt.plot(data['time'],data['mass']/1000.0,label='OpenRocket Simulation')
plt.grid()
plt.xlabel('Time (sec)')
plt.ylabel('Mass (kg)')

plt.show()

