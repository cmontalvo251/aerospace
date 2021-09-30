#!/usr/bin/python

from Universe import *
from pdf import * 


##In this routine we will simulate a planet or satellite of our choosing.
CB = UniverseParameters() ##class of standard universe parameters

#####################EARTH/MOON No Movement#########################################################
Earth = C.copy(CB.Earth)

#Then Let's create a Satellite
#Position of Satellite = 
altitude = 600000.
pos = np.asarray([Earth.r+altitude,0,0])
#Velocity of Satellite
vel_norm = 1.2*Earth.CircularVelocity(altitude)
inclination = 0.2
vel = np.asarray([0.,np.cos(inclination)*vel_norm,np.sin(inclination)*vel_norm])
JagSat = Satellite(1.3,10./100.,pos,vel,'JagSat','red',altitude+Earth.r)

#Then Let's Kick off the Simulation by creating a SolarSystem class
EarthSatellite = SolarSystem([Earth,JagSat],'EarthSystem')

#And Then Simulate the System
#tfinal = 18000.
tfinal = 36000.
timestep = 1.
tnext = 100.
#This uses an RK4 Routine to simulate the satellites
#EarthSatellite.Simulate(tfinal,timestep,tnext,False)
#This uses the orbital elements to compute the orbit of the satellites
EarthSatellite.Orbit()

##Finally Plot the Output of the Systems
print('Creating Plots')
pp = PDF(0,plt)

###This is Essentially Plotting Equation 1.3-4
#This plots the solution to the RK4 routine
#EarthSatellite.PlotSystem(pp,0)
#This plots the solution to the orbital elements
EarthSatellite.PlotOrbit(pp,0)

pp.close()

