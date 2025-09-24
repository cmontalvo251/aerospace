#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May 28 07:22:57 2020

@author: carlos
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as I
import pandas as pd #Try tools manage packages and search for pandas if that doesn't work use the shell and type 'pip3 install pandas'

data = pd.read_csv('fall_2024.csv') ##Use the plot/export tool and export your openrocket data to a csv. Then make sure to change the first line to match the data[] keys in this script

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
mass0 = 2.065 #kg - you get this from the rocket design tab in openrocket
T1 = 99.6/1000.0 #Newtons converted to kN - double click your motor in Motor and Configurations and grab the average thrust
massflowrate = (121./1000.)/2.3 #Take the difference betweeen the initial motor mass in g and final motor mass in grams, divide by 1000 to get to kg and then divide by the burn time
exit_velocity = T1*1000 / massflowrate
Isp = exit_velocity / 9.81
print('Isp = ',Isp)
Cd = 0.61 #open your spreadsheet and grab the average of your drag coefficient over the length of the burn time
Cd_parachute = 0.8 #double click the parachute and grab the Cd of the parachute from the properties
t_parachute = 9.0+2.3 #GO back to your motor properties and check your ejection delay. This t_parachute is your ejection delay + your burn time
D = 8.89/100. #Double click your body tube and then check for your outer diameter in cm and divide by 100 to get to meters
S = D**2
S_parachute = ((91.4)/100.0)**2 ##This is the area of the parachute so I double click the parachute and look at the diamter
stage_1_time = 2.3 ##This is your single stage burn time which is also in your motor properties
stage_2_start = -99 ##YOu have a single stage rocket so leave the next
stage_2_end = -99 #two lines here at -99
period = 58.0 ##This is the time it takes for my rocket to ascend and descend.
GNC = 0 ##leave these alone
apogee = 10000. #leave these alone

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
plt.plot(data['time'],data['altitude'],label='OpenRocket Simulation')
plt.grid()
plt.legend()
plt.xlabel('Time (sec)')
plt.ylabel('AGL (m)')


plt.figure()
plt.plot(data['time'],data['vertical acceleration'],label='OpenRocket Simulation')
plt.grid()
plt.legend()
plt.xlabel('Time (sec)')
plt.ylabel('Vertical Acceleration (m/s^2)')

plt.figure()
plt.plot(tout,vxout,label='Python Simulation')
plt.plot(data['time'],data['vertical velocity'],label='OpenRocket Simulation')
plt.grid()
plt.legend()
plt.xlabel('Time (sec)')
plt.ylabel('Velocity (m/s)')

plt.figure()
plt.plot(tout,mout,label='Python Simulation')
plt.plot(data['time'],data['mass'],label='OpenRocket Simulation')
plt.grid()
plt.legend()
plt.xlabel('Time (sec)')
plt.ylabel('Mass (kg)')

plt.show()

