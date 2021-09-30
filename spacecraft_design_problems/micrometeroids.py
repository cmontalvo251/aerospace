import numpy as np

###Mass of meteroid
m = 0.1 #grams

###Equation for mass above 10e-6
log10Nt = -14.37 - 1.213*np.log10(m)
print(log10Nt)
###Compute Nt 
Nt = 10**log10Nt ##This is in units of particles per square meter per second
print('Nt = ',Nt)

##Ok so the space shuttle is how big?
S = 367.0 #square meter

###So Nt is really our flux so let's rename the variable
flux = Nt

##In order to figure out how many meteroids hit the space shuttle we need to figure
##out how many particles hit the space shuttle every second.
particles_per_second = flux * S ##( particles / m^2 / second * m^2 = particles / sec)

print('Particles Per Second = ',particles_per_second)

##So that's a really small number but that's fine. So how many seconds
##need to go by before the space shuttle gets hit by 1 meteroid of 0.1 g or greater?
time = 1.0/particles_per_second ##this is in untis of seconds

###That's also a big number so let's convert to years
years = (time / 84600.) / 365. ##84,600 seconds in a day and 365 days in a year ish

print('Years = ',years) ##So 1,266 years. Assuming our math is right.