#!/usr/bin/python3
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import csv
import os
import sys
from gpscomms import GPSCOMMS

def convertTLE(LEADING_DIGITS,TRAILING_DIGITS,INSTR):
        DEC = INSTR.find('.')
        INSTRPRE = INSTR[0:DEC]
        INSTRPOST = INSTR[DEC+1:]
        if len(INSTRPOST) > TRAILING_DIGITS:
                INSTRPOST = INSTRPOST[0:TRAILING_DIGITS]
        else:
                while len(INSTRPOST) < TRAILING_DIGITS:
                        INSTRPOST+='0'
        while len(INSTRPRE) < LEADING_DIGITS:
                INSTRPRE = ' ' + INSTRPRE
        return INSTRPRE + '.' + INSTRPOST

##Make the GPS class
GPS = GPSCOMMS()

###TLE FILE FORMAT
#Line 0 - Satellite name
#Line 1 - Line Number, Catalog Number, Classification, Last two digits of launch year, launch number of the year, piece of the launch, last two digits of current year, current julian day, first derivative of mean motion, second derivative of mean motion, drag term, ephemeris type, element set number, checksum
#Line 2 - Line Number, Catalog Number, INC, RAAN, ECC, ARG, V0, Mean Motion, Revolution number at current TLE, checksum 
#The checksums for each line are calculated by adding all numerical digits on that line, including the line number. One is added to the checksum for each negative sign (-) on that line. All other non-digit characters are ignored.

###CREATE AN ORBIT AROUND THE EARTH
KEPLER = False
NAME = 'ISS' 
CATALOG = 'XXXXX'
CLASS = 'U'
LAUNCHYR = '22'
LAUNCHNO = 'XXX'
LAUNCHPIECE = 'A  '
YR = '22'
JULIANDAY = 165
JULIANHR = 9
JULIANMIN = 5
JULIANSEC = 0
BALLISTIC = -0.081629
MEANMOTION2 = ' 00000-0 ' 
RADIATIONPRESSURE = '-00000-0'
TLENO = 292
CHECKSUM = 7
REVNUMBER = '00000'
if KEPLER:
        print('Orbital Elements Given')
        #height_at_perigee_km = 600.0 #in kilometers
        #ECC = 0.5 #0 is circular orbit
        #INC = 56. #degrees
        #LAN = 30.0 #Longitude of the ascending node in degrees
        #ARG = 30. #argument of the periaps in degrees
        height_at_perigee_km = 412.4219041130347
        ECC =  0.0016186250257476963
        INC =  51.65411785007067
        LAN =  356.4565
        ARG =  75.9857377967462
        V0 = 298.2527
        ##The problem with Keplerian orbits is we don't have the time component. In order to do that we
        ##Need to get an initial state vector and then integrate the equations of motion
        #The trailing 0 just means give the coordinate at a true anomaly of 0
        print('ECC = ',ECC)
        print('INC = ',INC)
        print('LAN = ',LAN)
        print('ARG = ',ARG)
        print('Height At Perigee = ',height_at_perigee_km)
        x0,y0,z0,u0,v0,w0 = GPS.getStateVector(height_at_perigee_km,ECC,INC,LAN,ARG,V0)
else:
        print('State Vector Given')
        ##Sometimes though we are given a state vector
        #CIRCULAR 56 INC
        #x0 = 6.97814000e6
        #y0 = 0.0
        #z0 = 0.0
        #u0 = 0.0
        #v0 = 4.22627963e3
        #w0 = 6.26571722e3
        #CIRCULAR 56 INC ECC 0.3
        #x0 = 6.97814000e+06
        #y0 = 0.00000000e+00
        #z0 = 0.00000000e+00
        #u0 = 0.00000000e+00
        #v0 = 4.81870017e+03
        #w0 = 7.14401679e+03
        #CIRCULAR 56 INC ECC 0.5 LAN 30.0 ARG 30.0
        #x0 = 4258073.408143662
        #y0 = 4711293.537065895
        #z0 = 2892570.123274619
        #u0 = -6249.462870766547
        #v0 = 1567.9852296033864
        #w0 = 6645.79670240959
        #SAME ORBIT BUT A DIFFERNT POINT IN THE ORBIT
        #x0 = -11373129.869698685
        #y0 = -1582590.9250512433
        #z0 = 6398734.607771511
        #u0 = -3824.3727885792405
        #v0 = -3708.0441704344503
        #w0 = -1925.9564901311678
        #CUBESAT F1 FILE CONVERTED TO ORBITAL ELEMENTS AND THEN BACK TO STATE
        #x0 = 1388650.3862920234
        #y0 = 4181338.8491511224
        #z0 = 5167183.614583657
        #u0 = -7496.5015435881705
        #v0 = 689.9959793034844
        #w0 = 1456.2890215638377
        ##SAME ORBIT BUT RANDOM SPOT IN ORBIT
        #x0 = -3450261.6064182753
        #y0 = -3778046.447778934
        #z0 = -4497266.240788944
        #u0 = 6578.214958392202
        #v0 = -2159.4810810827257
        #w0 = -3238.621725744133
        #EXAMPLE CUBESAT F1 FILE from 2020
        x0 = 6639069.208831
        y0 = 627950.722629
        z0 = 1311013.452319
        u0 = -1606.655729
        v0 = 4713.714334
        w0 = 5821.859626
        #ANOTHER EXAMPLE SAT FROM 2004
        #x0 = 9686275.84
        #y0 = 41043317.21
        #z0 = 370926.32
        #u0 = -2976.641143
        #v0 = 699.873647
        #w0 = 315.983197
        print(x0,y0,z0,u0,v0,w0)
        ##And then we can get orbital elements from the state vector
        height_at_perigee_km,ECC,INC,LAN,ARG,V0 = GPS.getOrbitalElements(x0,y0,z0,u0,v0,w0)

##Now we integrate the EOMs with out state vecto
xsat_n,ysat_n,zsat_n,xdot_n,ydot_n,zdot_n,tsat = GPS.sixdof_orbit(x0,y0,z0,u0,v0,w0)

##USING THE ORBITAL ELEMENTS YOU CAN GET ANALYTIC ORBITS JUST USING KEPLERS EQUATIONS
xsat_a,ysat_a,zsat_a,xdot_a,ydot_a,zdot_a,nu_deg = GPS.kepler_orbit(height_at_perigee_km,ECC,INC,LAN,ARG)

##COMPUTE FRACTIONAL PORTION OF THE DAY
JULIANFRACTION = JULIANDAY + (JULIANHR + JULIANMIN/60. + JULIANSEC/3600.)/24.
JULIANSTR = str(JULIANFRACTION)
JULIANSTR = JULIANSTR[0:12]
##CONVERT BALLISTIC COEFFICIENT
BALLISTICSTR = str(BALLISTIC)
BALLISTICSTR = BALLISTICSTR[:1] + BALLISTICSTR[2:]
while len(BALLISTICSTR) < 10:
        BALLISTICSTR+='0'
#CONVERT RADIATION PRESSURE
RADSTR = str(RADIATIONPRESSURE)
RADSTR = RADSTR.replace('.','-')
##CONVERT INCLINATION
INCSTR = str(INC)
INCSTR = convertTLE(3,4,INCSTR)
#CONVERT THE LONGITUDE OF THE ASCENDING NODE
if LAN < 0:
        LAN+=360.
LANSTR = str(LAN)
LANSTR = convertTLE(3,4,LANSTR)
#CONVERT ECCENTRICITY
ECCSTR = str(ECC)
ECCSTR = convertTLE(0,7,ECCSTR)
ECCSTR = ECCSTR[2:]
#ARGUMENT OF PERIGEE
ARGSTR = str(ARG)
ARGSTR = convertTLE(3,4,ARGSTR)
#MEAN ANOMALY
V0STR = str(V0)
V0STR = convertTLE(3,4,V0STR)
##MEAN MOTION
MMSTR = str(GPS.MM)
MMSTR = convertTLE(2,8,MMSTR)

##PRINT THE TLE FORMAT FILE
#Check out this page - https://en.wikipedia.org/wiki/Two-line_element_set
print('============================')
print('TLE:')
print(NAME)
print('1 '+ CATALOG + CLASS + ' ' + LAUNCHYR + LAUNCHNO + LAUNCHPIECE + ' ' + YR + JULIANSTR + ' ' + BALLISTICSTR + MEANMOTION2 + RADSTR + ' 0 ' + ' ' + str(TLENO) + str(CHECKSUM))
print('2 '+ CATALOG + ' ' + INCSTR + ' ' + LANSTR + ' ' + ECCSTR + ' ' + ARGSTR + ' ' + V0STR + ' ' + MMSTR + REVNUMBER + str(CHECKSUM))
print('============================')

###PLOT 3D ORBIT
fig = plt.figure('3-D')
ax = fig.add_subplot(111,projection='3d')
ax.plot(xsat_a,ysat_a,zsat_a,color = 'red', linestyle = 'solid',label='Analytic Orbit')
ax.plot(xsat_n,ysat_n,zsat_n,color = 'blue', linestyle = 'solid',label='Numerical Orbit')
ax.plot(xsat_a[0:1],ysat_a[0:1],zsat_a[0:1],'m*',markersize = 20,label='Perigee')
ax.plot([x0,x0],[y0,y0],[z0,z0],'y*',markersize = 20,label='Initial Position')
plt.title('Satellite Orbit')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
GPS.plotEARTH(ax)
ax.legend()

##PLOT VELOCITY
fig = plt.figure()
plt.plot(nu_deg,xdot_a,label='X')
plt.plot(nu_deg,ydot_a,label='Y')
plt.plot(nu_deg,zdot_a,label='Z')
plt.grid()
plt.legend()
plt.xlabel('True Anamoly (deg)')
plt.ylabel('Velocity (m/s)')

fig = plt.figure()
plt.plot(tsat,xdot_n,label='X')
plt.plot(tsat,ydot_n,label='Y')
plt.plot(tsat,zdot_n,label='Z')
plt.grid()
plt.legend()
plt.xlabel('Time (sec)')
plt.ylabel('Velocity (m/s)')

#Plot delta
#plt.figure()
#plt.plot(time,np.sqrt((x0-xsat_n)**2 + (y0-ysat_n)**2 + (z0 - zsat_n)**2),'b*')
#plt.grid()

##Compute Distance from centroid of Earth
xsat_n = np.array(xsat_n)
ysat_n = np.array(ysat_n)
zsat_n = np.array(zsat_n)
tsat = np.array(tsat)
norm = np.sqrt(xsat_n**2 + ysat_n**2 + zsat_n**2) 

##Compute Latitude Longitude and (Geocentric)
print('Computing Geocentric Lat/Lon')
latitude,longitude = GPS.computeGeocentricLATLON(xsat_n,ysat_n,zsat_n,tsat)
#Convert to Radians because degrees suck
longitude_rad = longitude*np.pi/180.0
latitude_rad = latitude*np.pi/180.0

###Lat/Long (Geodetic)


###Swath Angle
print('Computing Swath Angle')
observer_angle_deg = 0.0 ##this will be higher if you can't see the horizon (this is in degrees)
sinrho = GPS.EarthRadius_ae/norm
sineta = np.cos(observer_angle_deg*np.pi/180.0)*sinrho
eta_rad = np.arcsin(sineta)
swath_angle_deg = 90.0 - eta_rad*180.0/np.pi - observer_angle_deg
swath_angle_rad = swath_angle_deg*np.pi/180.0

###Compute the Swath Vectors
yscalar = norm*np.tan(swath_angle_deg*np.pi/180.0)
latitudesw1 = 0*latitude
latitudesw2 = 0*latitude
longitudesw1 = 0*latitude
longitudesw2 = 0*latitude
print('Computing Swath Vectors and Swath Lat/Lons')
for n in range(0,len(xsat_n)):
    xyz = np.asarray([xsat_n[n],ysat_n[n],zsat_n[n]])
    vxyz = np.asarray([xdot_n[n],ydot_n[n],zdot_n[n]])
    hn = np.cross(xyz,vxyz)
    #h[:,n] = hn
    yn = yscalar[n]*hn/np.linalg.norm(hn)
    sw1 = xyz + yn
    sw2 = xyz - yn
    lat1,long1 = GPS.computeGeocentricLATLON(sw1[0],sw1[1],sw1[2],tsat[n]) #0 is here so that gps.py does not turn the Earth
    lat2,long2 = GPS.computeGeocentricLATLON(sw2[0],sw2[1],sw2[2],tsat[n]) #Eventually we need to turn the Earth
    latitudesw1[n] = lat1
    latitudesw2[n] = lat2
    longitudesw1[n] = long1
    longitudesw2[n] = long2
print('Done')

##Plot longitude vs time
#plt.figure()
#plt.plot(tsat,longitude)
#plt.show()
#sys.exit()

###Plot Ground Track
plt.figure()
n = 0
az_deg = np.linspace(0,180,100)
az_rad = az_deg * np.pi/180.0

nframes = 1.0
skip = int(len(xsat_n)/nframes)

##Location of Ground Station
#KMOB
lonG = -88.175228
latG = 30.691615
#MOB Shifted
#latG = 30.691615
#lonG = 25.0
#Origin
#lonG = 0.0
#latG = 0.0
#Random
#lonG = 100.0
#latG = 40.0

INSIDE = 0
if skip < 0:
	skip = 1
if skip > len(xsat_n):
	skip = len(xsat_n)

startTime = []
endTime = []
deltaTime = []

for ifov in range(0,len(xsat_n),skip):
	#Clear the axis
	plt.cla()

	##Plot Ground Track and Swatch Track
	plt.plot(longitude,latitude)
	plt.plot(longitudesw1,latitudesw1,'r-')
	plt.plot(longitudesw2,latitudesw2,'r-')

	##Plot Current Position
	plt.plot(longitude[ifov],latitude[ifov],'b*',markerSize=20)

	##Plot Location of Ground Station
	plt.plot(lonG,latG,'rs',markerSize=20)

	##Plot line from ground station to satellite
	plt.plot([longitude[ifov],lonG],[latitude[ifov],latG],'y-')

	##Compute the IFOV
	deltas = latitude_rad[ifov]
	deltas_deg = latitude[ifov]
	Ls = longitude_rad[ifov]
	Ls_deg = longitude[ifov]
	swath = swath_angle_rad[ifov]
	dt,dL = GPS.IFOV(swath,deltas,az_rad)
	dt_all = np.hstack((dt,dt[-1::-1]))
	Lt_all = np.hstack((Ls_deg-dL,Ls_deg+dL[-1::-1]))
	
	##Compute Azimuth to GND
	dlong_deg = lonG - Ls_deg
	dlat_deg  = latG - deltas_deg
	dlong_rad = dlong_deg*np.pi/180.0
	dlat_rad = dlat_deg*np.pi/180.0
	azG_rad = np.arctan2(dlong_rad,dlat_rad)
	azG_deg = azG_rad*180.0/np.pi

	##Compute the limit of lat/lon in that direction
	latFOV_G,dLG = GPS.IFOV(swath,deltas,azG_rad)
	longFOV_G1 = longitude[ifov] + dLG
	longFOV_G2 = longitude[ifov] - dLG

	###Compute the norm from sat to ground station
	r = GPS.EarthRadius_ae

	##Compute norm assuming flat earth with satellite as origin
	xG,yG,zG = GPS.LATLON2Cartesian(latG,lonG,r)
	xs,ys,zs = GPS.LATLON2Cartesian(latitude[ifov],longitude[ifov],r)
	sat2G = np.sqrt((xG-xs)**2 + (yG-ys)**2+(zG-zs)**2)
	###Compute norm from sat to IFOV
	xFOV1,yFOV1,zFOV1 = GPS.LATLON2Cartesian(latFOV_G,longFOV_G1,r)
	sat2FOV1 = np.sqrt((xs-xFOV1)**2 + (ys-yFOV1)**2+(zs-zFOV1)**2)
	xFOV2,yFOV2,zFOV2 = GPS.LATLON2Cartesian(latFOV_G,longFOV_G2,r)
	sat2FOV2 = np.sqrt((xs-xFOV2)**2 + (ys-yFOV2)**2+(zs-zFOV2)**2)

	##Plot the IFOV in black or green
	if sat2G < sat2FOV1 or sat2G < sat2FOV2:
		plt.plot(Lt_all,dt_all,'g-')
		plt.plot(longFOV_G1,latFOV_G,'gs',markerSize=20)
		plt.plot(longFOV_G2,latFOV_G,'gs',markerSize=20)
		"""
		print('g',xG,yG,zG)
		print('s',xs,ys,zs)
		print('fov1',xFOV1,yFOV1,zFOV1)
		print('lat/lon',latitude[ifov],longitude[ifov])
		print('lat/lonG',latG,longG)
		print('lat/lonFOV',latFOV_G,longFOV_G1)
		print('sat2G',sat2G/1000.0)
		print('sat2FOV1',sat2FOV1/1000.0)
		"""
		if INSIDE == 0:
			tstart = tsat[ifov]
			startTime.append(tstart)
			print('IN VIEW = ',tstart)
		INSIDE = 1
	else:
		plt.plot(Lt_all,dt_all,'k-')
		plt.plot(longFOV_G1,latFOV_G,'ks',markerSize=20)
		plt.plot(longFOV_G2,latFOV_G,'ks',markerSize=20)
		if INSIDE == 1:
			tend = tsat[ifov]
			endTime.append(tend)
			deltaTime.append((tend - tstart)/60.0)
			print('OUT OF VIEW = ',tend)
			print('Delta Time (min) = ',(tend - tstart)/60.0)
		INSIDE = 0
	plt.xlabel('Longitude (Deg)')
	plt.ylabel('Latitude (Deg)')
	plt.title(str(tsat[ifov]))
	plt.grid()
	plt.pause(0.1)



plt.show()
