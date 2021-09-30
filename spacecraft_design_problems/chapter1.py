#!/usr/bin/python

#Import all the Modules
import numpy as np
import matplotlib.pyplot as plt
import plotting as P
from pdf import *
import copy as C
from Universe import *

##Close all Windows
plt.close("all")         
        
##Simulate Earth Sun and Moon
CB = UniverseParameters() ##class of standard universe parameters

#####################EARTH/MOON No Movement#########################################################
Earth = C.copy(CB.Earth)

#Then Let's create a Satellite
#Position of Satellite = 
altitude = 600000.
pos = np.asarray([Earth.r+altitude,0,0])
#Velocity of Satellite
vel = 1.2*np.asarray([0,Earth.CircularVelocity(altitude),0])
JagSat = Satellite(1.3,10./100.,pos,vel,'JagSat','red',altitude+Earth.r)

#Then Let's Kick off the Simulation by creating a SolarSystem class
EarthSatellite = SolarSystem([Earth,JagSat],'EarthSystem')

#And Then Simulate the System
#tfinal = 18000.
tfinal = 36000.
timestep = 1.
tnext = 100.
EarthSatellite.Simulate(tfinal,timestep,tnext,False)

##Finally Plot the Output of the Systems
print('Creating Plots')
pp = PDF(0,plt)

###This is equation 1.1-2
Earth.plotgravity(pp)
#What is distance from Moon to Earth in AU?
print('Moon to Earth (Earth Radii) = ',CB.Moon.toCenter/CB.Earth.r)
print('Moon to Earth (AU) = ',CB.Moon.toCenter/CB.Earth.toCenter)

###This is Essentially Plotting Equation 1.3-4
EarthSatellite.PlotSystem(pp,0)

#Plotting position and velocity to show effect of Angular Momentum
EarthSatellite.plotPositionVelocity(pp)

##Plotting Eqaution 1.4-2    
EarthSatellite.PlotMechanicalEnergy(pp)

##plotting 1.4-3
#and Equation 1.4-4 -- Remember that FPA is 0 for circular orbits
EarthSatellite.plotAngularMomentum(pp)

pp.close()