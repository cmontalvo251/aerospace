#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep  4 21:56:24 2021

@author: carlos
"""

####Import all the modules we need
import numpy as np ###numeric python
import matplotlib.pyplot as plt ###matlab style plotting
import scipy.integrate as sci ##integration toolbox
import time

plt.close("all")

##DEFINE SOME CONSTANT PARAMETERS
G = 6.6742*10**-11; #%%Gravitational constant (SI Unit)

###PLANET
###EARTH
Rplanet = 6357000.0 #meters
mplanet = 5.972e24 #kg
name = 'Earth'
###KERBIN
#Rplanet = 600000 #meters
#mplanet = 5.2915158*10**22
#name = 'Kerbin'

##PARAMETERS OF ROCKET
###Initial Conditions for single stage rocket
x0 = Rplanet 
z0 = 0.0
velz0 = 0.0
velx0 = 0.0
period = 300 #2*np.pi/np.sqrt(G*mplanet)*r0**(3.0/2.0)*1.5
mass0 = 13.47 #kg
average_thrust = 1348.0 #N
Isp1 = 265.61 #seconds
Isp2 = 400.0
tMECO = 3.84 #main engine cutoff time
tsep1 = 2.0 #length of time to remove 1st stage
mass1tons = 0.0
mass1 = mass1tons*1000
t2start = 100000.0
t2end = t2start + 17.5
D = (6.0/12.0)/3.28 #meters
CD = 0.4839 #nd

###Create a Aerodynamics Class
class Aerodynamics():
    def __init__(self,name):
        self.name = name
        if name == 'Kerbin':
            ###import the aero model for interpolation
            data = np.loadtxt('kerbin_aerodynamics.txt')
            #print(data)
            self.altitude = data[:,0]
            #print(self.altitude)
            self.density = data[:,3]
            self.rhos = self.density[0]
            self.beta = 0.0
        elif name == 'Earth':
            ##going to use the Earth aero model
            self.beta = 0.1354/1000.0 ##density constant
            self.rhos = 1.225 #kg/m^3
    
    def getDensity(self,altitude):
        if self.name == 'Kerbin':
            ###interpolate
            rho = np.interp(altitude,self.altitude,self.density)
        elif self.name == 'Earth':
            ###Use special equation
            rho = self.rhos*np.exp(-self.beta*altitude)
        return rho
            
            
##Create the Aeromodel variable which is an instance
##of the class Aerodynamics and I'm putting it up here
##so that it is a global
aeroModel = Aerodynamics(name)

##Gravitational Acceleration Model
def gravity(x,z):
    global Rplanet,mplanet
    
    r = np.sqrt(x**2 + z**2)
    
    if r < Rplanet:
        accelx = 0.0
        accelz = 0.0
    else:
        accelx = G*mplanet/(r**3)*x
        accelz = G*mplanet/(r**3)*z
        
    return np.asarray([accelx,accelz]),r

def propulsion(t):
    global average_thrust,Isp,tMECO
    ##Timing for thrusters
    if t < tMECO:
        #We are firing the main thruster
        theta = 0*np.pi/180.0
        thrustF = average_thrust
        ve = Isp1*9.81 #m/s
        mdot = -thrustF/ve
    if t > tMECO and t < (tMECO + tsep1):
        theta = 0.0
        thrustF = 0.0
        ## masslost = mass1
        mdot = -mass1/tsep1
    if t > (tMECO + tsep1):
        theta = 0.0
        thrustF = 0.0
        mdot = 0.0
    if t > (t2start) and t < (t2end):
        theta = 90.0*np.pi/180.0
        thrustF = max_thrust
        ve = Isp2*9.81 #m/s
        mdot = -thrustF/ve
    if t > t2end:
        theta = 0.0
        thrustF = 0.0
        mdot = 0.0
    
    thrustx = thrustF*np.cos(theta)
    thrustz = thrustF*np.sin(theta)
    
    return np.asarray([thrustx,thrustz]),mdot
    
###Equations of Motion
###F = m*a = m*zddot
## z is the altitude from the center of the planet along the north pole
### x  is the altitude from center along equator through Africa 
## this is in meter
## zdot is the velocity along z
## zddot is the acceleration along z
###Second Order Differential Equation
def Derivatives(state,t):
    global aeroModel,Rplanet
    #state vector
    x = state[0]
    z = state[1]
    velx = state[2]
    velz = state[3]
    mass = state[4]
    
    #Compute zdot - Kinematic Relationship
    zdot = velz
    xdot = velx
    
    ###Compute the Total Forces
    ###GRavity
    accel,r = gravity(x,z)
    gravityF = -accel*mass
    #print(gravityF)
    
    ###Aerodynamics
    altitude = r - Rplanet ##altitude above the surface
    rho = aeroModel.getDensity(altitude) ##air density
    V = np.sqrt(velz**2+velx**2)
    qinf = (np.pi/8.0)*rho*(D**2)*abs(V)
    aeroF = -qinf*CD*np.asarray([velx,velz])
    #print(aeroF)
    
    ###Thrust
    thrustF,mdot = propulsion(t)
    #print(thrustF)
    
    Forces = gravityF + aeroF + thrustF
    #print(Forces)
    #print('------')
    
    #time.sleep(5.0)
    
    #Compute Acceleration
    if mass > 0:
        ddot = Forces/mass
    else:
        ddot = 0.0
        mdot = 0.0
    
    #Compute the statedot
    statedot = np.asarray([xdot,zdot,ddot[0],ddot[1],mdot])
    
    return statedot

###########EVERYTHING BELOW HERE IS THE MAIN SCRIPT###

###Test Surface Gravity
print('Surface Gravity (m/s^2) = ',gravity(0,Rplanet))

###Plot the Air Density as a function of AGL
test_altitude = np.linspace(0,100000,100)
test_rho = aeroModel.getDensity(test_altitude)
plt.figure()
plt.plot(test_altitude,test_rho,'b-')
plt.xlabel('altitude (m)')
plt.ylabel('Air Density (kg/m^3')
plt.grid()

##Compute Exit Velocity
##Populate Initial Condition Vector
stateinitial = np.asarray([x0,z0,velx0,velz0,mass0])

##Time window
tout = np.linspace(0,period,1000)

###Numerical Integration Call
stateout = sci.odeint(Derivatives,stateinitial,tout)

###REname variables
xout = stateout[:,0]
zout = stateout[:,1]
altitude = np.sqrt(xout**2+zout**2) - Rplanet
velxout = stateout[:,2]
velzout = stateout[:,3]
velout = np.sqrt(velxout**2 + velzout**2)
massout = stateout[:,4]

###Plot

###ALTITUDE
plt.figure()
plt.plot(tout,altitude)
plt.xlabel('Time (sec)')
plt.ylabel('Altitude (m)')
plt.grid()

###VELOCITY
plt.figure()
plt.plot(tout,velout)
plt.xlabel('Time (sec)')
plt.ylabel('Total Speed (m/s)')
plt.grid()

###Mass
plt.figure()
plt.plot(tout,massout)
plt.xlabel('Time (sec)')
plt.ylabel('Mass (kg)')
plt.grid()

##2D Orbit
plt.figure()
plt.plot(xout,zout,'r-',label='Orbit')
plt.plot(xout[0],zout[0],'g*')
theta = np.linspace(0,2*np.pi,1000)
xplanet = Rplanet*np.sin(theta)
yplanet = Rplanet*np.cos(theta)
plt.plot(xplanet,yplanet,'b-',label='Planet')
plt.grid()
plt.legend()

plt.show()
