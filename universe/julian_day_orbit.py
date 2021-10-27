#!/usr/bin/python

from pdf import *
from Universe import *
import matplotlib.pyplot as plt

#This will plot the position of the planets on a given Julian Day

#Here is the Julian Day
#julian_day = 2444240 ##this is jan 1 1980 #this was a leap year
#julian_day = 2445701 #this is jan 1 1984 #1984 was a leap year
julian_day = 2446796 ##this is jan 1 1987 -- really? yes really. You can
#double check on this website
#https://www.heavens-above.com/planets.aspx
#julian_day = 2446797
#julian_day = 2447162 ##this is jan 1 1988 #1988 was a leap year
#julian_day = 2451545. ##this is jan 1 2000
#julian_day = 2458120.-365 ##this is jan 1 2017
#julian_day = 2458120. ##this is jan 1 2018
#julian_day = 2458485 ##this is jan 1 2019
#julian_day = 2458850 ##this is jan 1 2020
#julian_day = 2459215 #this is jan 1 2021
#julian_day = 2459580  #this is jan 1 2022
#julian_day -= 15
#julian_day += 275 #October 2
# - 10 to get to the winter solstice from the next year
# + 79 to get to the spring equinox
# + 172 to get to the summer solstice
# + 245 to get to Sep 2nd
# + 265 to get to the fall equinox
# + 355 to get to the winter solstice
julian_day += 245 #Sep 2nd
#julian_day += 284
#julian_day += 49
#julian_day += 211
#julian_day += 182
#julian_day += 300
#julian_day += 365 - 30
#julian_day -= 10
#julian_day += 245
#julian_day += 355

#julian_day += 492*365

##Then compute all the planets using the JPL class
planets = JPL(julian_day)

#Finally compute all the orbits
planets.MilkyWay.Orbit()

##Finally Plot the Output of the Systems
print('Creating Plots')
pp = PDF(0,plt)

##Plot All planets
planets.MilkyWay.PlotOrbit(pp,-1)

##Only plot inner planets
planets.MilkyWay.numsatellites = 5
planets.MilkyWay.PlotOrbit(pp,-1)
pp.close()
sys.exit()

##Animation??
#planets.MilkyWay.numsatellites = 6 + 1 #The plus one is because of loops in Python. I know I probably need to fix the loop but whatever
#pa = PDF(1,plt)
###
#Inputs are the 
#julian_day = defined above
#day_skip = how many days to skip
#day_skip = 1.0
#num_skips = the number of skips (must be integer!!)
#num_skips = 365
#pause_time = how long to pause
#pause_time = 0.1
#planets.AnimateOrbits(pa,julian_day,day_skip,num_skips,pause_time)

##Use Mayavi if you're using Python3
#planets.MilkyWay.numsatellites = 10
#planets.MilkyWay.PlotMayavi()
