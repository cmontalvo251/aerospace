#!/usr/bin/python

import numpy as np

##constants. don't change
g = 32.2
pi2 = 2*np.pi

##This code will compute the inertia based on simple pendulum dynamics
#To do this experiment you need to hang your vehicle from a hinge point and
#compute the distance of the rope and the distance from the connection point
#on the rigid body to the center of mass of the vehicle. You need to make
#sure to connect the rope to a point on the vehicle that is not the
#center of mass otherwise the physics won't work right

#Used for Ixx and Iyy
lrope = 145./12. #feet
lconnection2cg = 13.5/12. #feet

#Used for Izz
#lrope = 124./12. #feet
#lconnection2cg = 24./12. #feet

#You also need to weigh the vehicle. Feel free to weigh in pounds and
#convert or in grams and convert

#weight_grams = 735.
#Conversion from grams
#mass = (weight_grams/1000.0)*2.2/g

#Conversion from lbf
weight_lbf = 50.3
mass = weight_lbf/g

##Need weight of the rope
weight_rope = 1.0 #lbf
mass_rope = weight_rope/g

##Then you need to swing this bitch and time one full period. That is,
##let go and time how long it takes to come back

#Period for Ixx
period = 3.61 ##seconds
#Period for Iyy
#period = 3.34 ##seconds
#Period for Izz
#period = 3.60 #seconds

##The inertia of the vehicle about the axis of rotation about the
##center of mass is equal to this equation. The derivation is done
##using sum(Moments) = dH/dt the pdf of the derivation is in the
##derivations folder assuming I pushed it otherwise Dr. C has it on
##his computer

##the length l used for these calculations is actually a bit more complex than I though.
##First let's compute the length from the hinge joint to the center of mass of the vehicle
##you are swinging
lvehicle_cg = (lrope + lconnection2cg) ##Combine the length of the rope and the length of the connection point to the cg of the vehicle
##The center of mass of the rope is just the length of the rope over 2. 
lrope_cg = lrope/2.

##the length of the entire system is just the center of mass formula
l = (mass_rope*lrope_cg + mass*lvehicle_cg)/(mass_rope+mass)

print('Length of swing (ft) = ',l)

Icg = ((period/pi2)**2)*mass*g*l - mass*l**2
print('Inertia (slugs-ft^2) = ',Icg)

###Let's compute a couple checks. 
##If you had an Icg = 0 the system would swing with a natural frequency of sqrt(g/L)
ideal_natural_freq = np.sqrt(g/l)
## wn = 2*pi*f
## f = wn/(2*pi)
ideal_freq = ideal_natural_freq/pi2
## T = 1/f
ideal_period = 1/ideal_freq

print('Experimental Period (sec) = ',period)
print('Ideal Period (sec) = ',ideal_period)

if (period < ideal_period):
	print('Your experimental period is less than the period with inertia = 0')
	print('This means that you measured the rope incorrectly or you measured')
	print('The period incorrectly')



