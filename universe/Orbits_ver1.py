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
        
##Let's Simulate the Two Body Problem from Chapter 1

CB = UniverseParameters() ##class of standard universe parameters

#######################EARTH + SATELLITE################################################3

#Let's the create the Earth Planet
Earth = Satellite(5.972e24,6357000,np.asarray([0.,0.,0.]),np.asarray([0.,0.,0.]),'Earth','blue',CB.Earth.toCenter)

#Then Let's create a Satellite
#Position of Satellite = 
altitude = 600000.
pos = np.asarray([Earth.r+altitude,0,0])
#Velocity of Satellite
vel = np.asarray([0,Earth.CircularVelocity(altitude),0])
JagSat = Satellite(1.3,10./100.,pos,vel,'JagSat','red',Earth.r+altitude)

#Then Let's Kick off the Simulation by creating a SolarSystem class
EarthSystem = SolarSystem([Earth,JagSat],'EarthSystem')

#And Then Simulate the System
tfinal = 5800.
timestep = 10.
tnext = 100.
EarthSystem.Simulate(tfinal,timestep,tnext,False)

#####################EARTH + MOON #########################################################

###Let's Create the Earth and Moon
#Earth will stay the same but we need to create a new one 
#Moon will need to be added
moon_distance = 238900*5280/3.28
moon_altitude = moon_distance - Earth.r
moon_vel = np.asarray([0,Earth.CircularVelocity(moon_altitude),0])
Moon = Satellite(7.34767309e22,1079*5280./3.28,np.asarray([moon_distance,0,0]),moon_vel,'Moon','gray',CB.Moon.toCenter)

#Create Earth Moon System
EarthMoon = SolarSystem([Earth,Moon],'EarthMoon')

#Simulate this system
tfinal = 2373496.5656129662
timestep = 100.
tnext = 1000.
EarthMoon.Simulate(tfinal,timestep,tnext,False)

###################EARTH + MOON + SATELLITE########################################

EarthMoonSat = SolarSystem([Earth,JagSat,Moon],'EarthMoonSat')
#Simulate this system
tfinal = 2373496.5656129662
timestep = 100.
tnext = 1000.
EarthMoonSat.Simulate(tfinal,timestep,tnext,True)

###################EARTH + SUN#############################
#Sun = Satellite(1.989e30,432169*5280./3.28,np.asarray([0,0,0]),np.asarray([0,0,0]),'Sun','yellow')
##Need to make another Earth that is moving
#Earth_Distance = (92.96e6)*5280/3.28
#pos = np.asarray([Earth_Distance,0,0])
#vel = np.asarray([0,Sun.CircularVelocity(Earth_Distance),0])
#Earth3 = Satellite(Earth.M,Earth.r,pos,vel,'Earth3','blue')

#Helios = SolarSystem([Sun,Earth3],'Helios')

##Simulate Helios
#tfinal = 31788686.529976141
#timestep = 100.
#tnext = 1000.
#Helios.Simulate(tfinal,timestep,tnext)

################EARTH + SUN + MOON########################
#In order to put the Moon in the proper orbit we need to compute Earth3
#print(Earth3.initial_pos)
#Moon_pos3 = Earth3.initial_pos + Moon.initial_pos
#Moon_vel3 = Earth3.initial_vel + Moon.initial_vel
#Moon3 = Satellite(Moon.M,Moon.r,Moon_pos3,Moon_vel3,'Moon3','gray')
#HeliosLuna = SolarSystem([Sun,Earth3,Moon3],'HeliosLuna')
##Simulate Helios Plus Luna
#tfinal = 31788686.529976141
#tfinal = 40000.
#timestep = 100.
#tnext = 1000.
#HeliosLuna.Simulate(tfinal,timestep,tnext) #This doesn't work but more than likely it's because of my initial conditions?

##Finally Plot the Output of the Systems
print('Creating Plots')
pp = PDF(0,plt)
EarthSystem.PlotSystem(pp,-1)
EarthMoon.PlotSystem(pp,-1)
EarthMoonSat.PlotSystem(pp,-1)
EarthMoonSat.PlotSystem(pp,0)
#Helios.PlotSystem(pp)
#HeliosLuna.PlotSystem(pp)
pp.close()