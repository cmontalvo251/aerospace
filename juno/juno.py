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

#Gemini Pro 3.1 wrote this function to make my atmosphere model more accurate
def get_atmosphere(altitude):
    """
    Returns atmospheric density (rho), pressure ratio, and speed of sound (a)
    based on Droo's standard atmosphere.
    """
    R_air = 287.058 # Specific gas constant for air J/(kg*K)
    gamma = 1.4 # Heat capacity ratio
    rho_sl = 1.225 # Sea level density kg/m^3
    T_sl = 288.15 # Sea level temperature K
    
    if altitude <= 0:
        return rho_sl, 1.0, np.sqrt(gamma * R_air * T_sl)
    if altitude >= 60000:
        return 0.0, 0.0, 0.0 # Vacuum
        
    if altitude < 11000:
        # Troposphere: Temperature drops linearly
        T = T_sl - 0.0065 * altitude
        p_ratio = (1 - (0.0065 * altitude) / T_sl)**5.2558
        # Ideal Gas Law density calculation relative to sea level
        rho = rho_sl * p_ratio * (T_sl / T)
    else:
        # Stratosphere: Temperature is constant in this simplified model
        T = 216.65 
        p_11k = 0.2233 
        p_ratio = p_11k * np.exp(-(altitude - 11000) / 6340)
        rho = rho_sl * p_ratio * (T_sl / T)
        
    speed_of_sound = np.sqrt(gamma * R_air * T)
    return rho, p_ratio, speed_of_sound


#Gemini Pro 3.1 also wrote this function
def get_drag_coefficient(velocity, speed_of_sound):
    """
    Interpolates the drag coefficient based on Mach number 
    to simulate the transonic wave drag wall.
    """
    if speed_of_sound == 0:
        return 0.0 # Vacuum, no drag
        
    mach_number = velocity / speed_of_sound
    
    # Define the Mach vs Cd curve (The "Transonic Drag Wall")
    mach_points = [0.0,  0.8,  1.05, 1.2,  2.0,  4.0,  10.0]
    cd_points   = [0.35, 0.36, 0.85, 0.80, 0.45, 0.35, 0.30]
    
    # np.interp mathematically connects the dots between the arrays above
    cd = np.interp(mach_number, mach_points, cd_points)
    return cd

def planet_parameters():
    ##https://junoneworigins.fandom.com/wiki/Juno_System
    ##A few things you can compute to compare to the wiki above
    #mu 
    #surface gravity
    #rotational velocity
    
    G = 6.6742*10**-11; #%%Gravitational constant
    Mdroo = 5.286*10**22 #
    mudroo = G*Mdroo
    Rdroo = 600000 #meters
    sidereal_period = 21549.425
    sidereal_angular_velocity = 2*np.pi/sidereal_period
    sidereal_rotational_velocity = sidereal_angular_velocity*Rdroo
    surface_gravity = mudroo*Rdroo/Rdroo**3
    return sidereal_rotational_velocity,mudroo,surface_gravity,Rdroo

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
    altitude = rSat-R
    rho, p_ratio, speed_of_sound = get_atmosphere(altitude)

    ##Thrust ratio based on altitude
    T_ALT = T1 + (TVAC-T1)*(1-p_ratio)
    ISP_ALT = Isp + (IspVAC-Isp)*(1-p_ratio)
    
    #A few problems here. First Cd is a function of Mach Number and Reynolds number
    #As Such I had Gemini Pro 3.1 write a nice function to interpolate drag
    velocity = np.sqrt(velx**2 + velz**2)
    Cd = get_drag_coefficient(velocity,speed_of_sound)
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
            thrust = T_ALT*1000.0
            theta = 90.0*altitude/apogee
            if theta > 90.0:
                theta = 90.0
            thrustx = thrust*np.cos(theta*np.pi/180.0)
            thrustz = thrust*np.sin(theta*np.pi/180.0)
        else:
            thrustx = T_ALT*1000.0
            thrustz = 0.0 
                
    #Fire stage 2
    if t > stage_2_start and t < stage_2_end:
        thrust = T2*1000.0
        ##What angle?
        ##Tangent to the current x,z coordinate
        #First get the angle to the current coordinate
        if GNC:
            thetaSat = np.arctan2(z,x)
        else:
            thetaSat = np.pi/2.0
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
        ve = ISP_ALT*abs(surface_gravity)
        mdot = -thrust/ve
    else:
        mdot = 0.0
        
    ##If we are in between
    if t > stage_1_time and t < stage_2_start:
        mdot = mass_stage_1*2000/2.2 / (stage_1_time - stage_2_start)
    
    ##Now we can put Newton's EOMs together
    xdbldot = thrustx/mass + gravx + aerox
    zdbldot = thrustz/mass + gravz + aeroz
    
    statedot = np.asarray([xdot,zdot,xdbldot,zdbldot,mdot])
    if altitude < 0:
        statedot = np.asarray([0,0,0,0,0])
    
    #make [xdot,xdbldot] an array
    return statedot

##############END OF FUNCTION SEPARATED BY TABS#######
    
#Set up planet parameters
sidereal_rotational_velocity,mu,surface_gravity,R = planet_parameters()

##Initial Conditions for suborbital Orbit
x0 = R
z0 = 0.
velx0 = 0.0
velz0 = 0.0
masstons = 11.8
T1 = 475.0 ##Thrust ASL - At sea level
TVAC = 539. ##Thrust in Vacuum
Isp = 226. #This is ISP at sea level dv = Isp*ln(MR)
IspVAC = 256. #This is ISP in vacuum
D = 2.01 #meters using width and depth
S = D**2 #This needs to be D^2
stage_1_time = 41.
stage_2_start = -99
stage_2_end = -99
period = 2000.0
GNC = 0
apogee = 10000.


###Conditions for Sub orbital flight (Two Stage)
"""
x0 = R
z0 = 0.
velx0 = 0.0
velz0 = 0.0
masstons = 13.5
mass_stage_1 = 2.2 #this is the mass of the first stage
stage_1_time = 58.0
stage_2_start = stage_1_time + 1.0
stage_2_end = stage_2_start + 58.0
T1 = 167.97
T2 = 167.97
Isp = 250.0
Cd = 0.5
D = 1.25
S = D**2
period = 1000.0
GNC = 0
apogee = 80000.
"""
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
mass0 = masstons*1000 #CONVERT to kg via metric tons
tout = np.linspace(0,period,1000000)  #linspace(start,end,number of data points)
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

altitude = (np.sqrt(xout**2+zout**2)-R)

plt.figure()
plt.plot(tout,altitude/1000.)
plt.grid()
plt.xlabel('Time (sec)')
plt.ylabel('AGL (km)')

rho = []
rho0,_,_ = get_atmosphere(0.0) ##This will always be 1.225
for a in altitude:
    rhoi,_,_ = get_atmosphere(a)
    rho.append(rhoi)
rho = np.asarray(rho)
##Thrust ratio based on altitude
ratio = rho/rho0 ##This number will be between 0 and 1, start at 1 and go to 0
inverse_ratio = 1.0 - ratio ##This number will start at 0 and go to 1
T_ALT = T1*ratio + TVAC*inverse_ratio
ISP_ALT = Isp*ratio + IspVAC*inverse_ratio
plt.figure()
plt.plot(tout,T_ALT)
plt.grid()
plt.xlabel('Time (sec)')
plt.ylabel('Thrust (kN)')

plt.figure()
plt.plot(tout,ISP_ALT)
plt.grid()
plt.xlabel('Time (sec)')
plt.ylabel('Isp (sec)')

plt.figure()
plt.plot(tout,mout/1000.)
plt.grid()
plt.xlabel('Time (sec)')
plt.ylabel('Mass (metric tons)')

###Compute the delta V
dv = Isp*surface_gravity*np.log(mass0/mout[-1]) #in m/s
print('Delta V = ',dv,'m/s')

plt.show()

