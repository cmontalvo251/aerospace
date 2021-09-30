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
Sun = Satellite(CB.Sun.M,CB.Sun.r,np.asarray([0,0,0]),np.asarray([0,0,0]),'Sun','yellow',0)
vnorm = Sun.CircularVelocity(CB.Earth.toCenter)

###Let's Create the Earth and Moon
Earth_distance = np.asarray([CB.Earth.toCenter,0.,0.])
Earth_vel = np.asarray([0,0.,0.])
Earth = Satellite(CB.Earth.M,CB.Earth.r,Earth_distance,Earth_vel,'Earth','blue',CB.Earth.toCenter)

moon_distance = np.asarray([CB.Earth.toCenter+CB.Moon.toCenter,0,0])
moon_vel = np.asarray([0,Earth.CircularVelocity(CB.Moon.toCenter),0])

Moon = Satellite(CB.Moon.M,CB.Moon.r,moon_distance,moon_vel,'Moon','gray',CB.Moon.toCenter)

#Create Earth Moon System
EarthMoon = SolarSystem([Earth,Moon],'EarthMoonInertial')

#Simulate this system
#tfinal = 31788686.529976141
#tfinal = 317886.86529976141
tfinal = 2432590.3757229168
timestep = 100.
tnext = 1000.
#EarthMoon.Simulate(tfinal,timestep,tnext,True)

###################EARTH/MOON With Orbital velocity of Earth#############################

###Let's Create the Earth and Moon
Earth_distance = np.asarray([CB.Earth.toCenter,0.,0.])
Earth_vel = np.asarray([0,vnorm,0.])
Earth = Satellite(CB.Earth.M,CB.Earth.r,Earth_distance,Earth_vel,'Earth','blue',CB.Earth.toCenter)

moon_distance = np.asarray([CB.Earth.toCenter+CB.Moon.toCenter,0,0])
moon_vel = np.asarray([0,Earth.CircularVelocity(CB.Moon.toCenter)+vnorm,0])

Moon = Satellite(CB.Moon.M,CB.Moon.r,moon_distance,moon_vel,'Moon','gray',CB.Moon.toCenter)

#Create Earth Moon System
EarthMoonOrbital = SolarSystem([Earth,Moon],'EarthMoonSunCentered')
#EarthMoonOrbital.Simulate(tfinal,timestep,tnext,False)

###EARTH MOON SUN
EarthMoonSun = SolarSystem([Sun,Earth,Moon],'EarthMoonSun')
EarthMoonSun.Simulate(tfinal,timestep,tnext,False)

##Finally Plot the Output of the Systems
print('Creating Plots')
pp = PDF(0,plt)
#EarthMoon.PlotSystem(pp,False)
#EarthMoonOrbital.PlotSystem(pp,False)
EarthMoonSun.PlotSystem(pp,0)
pp.close()