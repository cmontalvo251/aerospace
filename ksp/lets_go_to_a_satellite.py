import numpy as np
import matplotlib.pyplot as plt

####CONSTANTS

###MILY WAY GALAXY
'''
Radius_of_Planet = 6379000 #meters of Planet
Altitude_of_Parking_Orbit_Planet = 200000 #meters
G = 6.67e-11 #SI
mPlanet = 5.972e24 #kg
Radius_of_Satellite_from_Center_of_Satellite = 1738100 #meters of Satellite we intend to travel to (Moon, Mun, Minmus...)
Distance_bt_Planet_Center_2_Satellite_Center = 384467000 #meters (bt = between)
mSatellite = 0.07346e24
Altitude_of_Parking_Orbit_Satellite = 112000 #meter - This is what Apollo
'''

###KERMS SOLAR SYSTEM - Kerbin to Mun
Radius_of_Planet = 600000 #meters of Planet
Altitude_of_Parking_Orbit_Planet = 80000 #meters
G = 6.67e-11 #SI
mPlanet = 5.291e22 #kg
Radius_of_Satellite_from_Center_of_Satellite = 20000 #meters of Satellite we intend to travel to (Moon, Mun, Minmus...)
Distance_bt_Planet_Center_2_Satellite_Center = 12000000 #meters (bt = between)
mSatellite = 9.759e20
Altitude_of_Parking_Orbit_Satellite = 5000 #meter - This is what Apollo


###COMPUTATION

#Compute Mu
muPlanet = G*mPlanet
muSatellite = G*mSatellite

##Creating a circle for the planet
angle = np.linspace(0,2*np.pi,10000)
xPlanet = Radius_of_Planet*np.cos(angle)
yPlanet = Radius_of_Planet*np.sin(angle)
#Create a circle for satellite
xSatellite = Radius_of_Satellite_from_Center_of_Satellite*np.cos(angle) - Distance_bt_Planet_Center_2_Satellite_Center
ySatellite = Radius_of_Satellite_from_Center_of_Satellite*np.sin(angle)
xSatellite_Orbit = Distance_bt_Planet_Center_2_Satellite_Center*np.cos(angle)
ySatellite_Orbit = Distance_bt_Planet_Center_2_Satellite_Center*np.sin(angle)
##Compute the Orbital Elements of our Parking Orbit
a = Radius_of_Planet + Altitude_of_Parking_Orbit_Planet
e = 0
p = a
h = np.sqrt(p*muPlanet)
r = a
v_parking_orbit_planet = h/r
print('Velocity of Parking Orbit around Planet = ',v_parking_orbit_planet,' m/s')
xParkingPlanet = r*np.cos(angle)
yParkingPlanet = r*np.sin(angle)
###DELTAV TO GET TO PARKING ORBIT AROUND PLANET
DV1 = 2*v_parking_orbit_planet
print('DV of First Manuever (Surface to Parking Orbit) = ',DV1,' m/s')
###COMPUTE OUR TRANSFER ORBIT
ra = Distance_bt_Planet_Center_2_Satellite_Center + Radius_of_Satellite_from_Center_of_Satellite + Altitude_of_Parking_Orbit_Satellite
rp = a
aTIO = (ra+rp)/2
eTIO = (ra-rp)/(ra+rp)
print("Eccentricity of TIO Orbit = ",eTIO)
pTIO = aTIO*(1-eTIO**2)
hTIO = np.sqrt(pTIO*muPlanet)
va = hTIO/ra
vp = hTIO/rp
print("Velocity at Perigee = ",vp," m/s  Velocity at Apogee = ",va,' m/s')
DV2 = vp - v_parking_orbit_planet
print("DV of Trans Injection Orbit = ",DV2," m/s")
##Generate TIO Orbit
rTIO = pTIO/(1+eTIO*np.cos(angle))
xTIO = rTIO*np.cos(angle)
yTIO = rTIO*np.sin(angle)

###PLOTS

###OVERALL PLOT
plt.figure()
plt.plot(xPlanet,yPlanet,'b-')
plt.plot(xParkingPlanet,yParkingPlanet,'r-')
plt.plot(xSatellite,ySatellite,'k-')
plt.plot(xSatellite_Orbit,ySatellite_Orbit,'k-')
plt.plot(xTIO,yTIO,'g-')
plt.axis('equal')
plt.grid()

###ZOOMED IN AROUND PLANET PLOT
plt.figure()
plt.plot(xPlanet,yPlanet,'b-')
plt.plot(xParkingPlanet,yParkingPlanet,'r-')
plt.plot(xSatellite,ySatellite,'k-')
plt.plot(xSatellite_Orbit,ySatellite_Orbit,'k-')
plt.plot(xTIO,yTIO,'g-')
plt.xlim([-1.5*Radius_of_Planet,1.5*Radius_of_Planet])
plt.ylim([-1.5*Radius_of_Planet,1.5*Radius_of_Planet])
plt.grid()

###ZOOMED IN AROUND SATELLITE PLOT
plt.figure()
plt.plot(xPlanet,yPlanet,'b-')
plt.plot(xParkingPlanet,yParkingPlanet,'r-')
plt.plot(xSatellite,ySatellite,'k-')
plt.plot(xSatellite_Orbit,ySatellite_Orbit,'k-')
plt.plot(xTIO,yTIO,'g-')
var = Distance_bt_Planet_Center_2_Satellite_Center
rvar = Radius_of_Satellite_from_Center_of_Satellite
plt.xlim([-var-1.5*rvar,-var+1.5*rvar])
plt.ylim([-1.5*rvar,1.5*rvar])
plt.grid()

plt.show()