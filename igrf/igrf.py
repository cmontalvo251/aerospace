import sys
import matplotlib.pyplot as plt
import numpy as np
import Earth_Orbit as orb
try:
	import pyIGRF as IGRF
except:
	print('pyIGRF not found')
	print('pip3 install pyIGRF to install')
	sys.exit()

def bfield_arr(lat,lon,alt,date):
	#This function takes latitude and longitude arrays in degrees
	#Altitude in units of km above ground level (AGL) and a date in a
	#int format like 2020 or 2021
	if len(lat) != len(lon) or len(lon) != len(alt):
		print('Array length not correct')
		return 0,0,0,0
	bx_arr = 0*lat
	by_arr = 0*lat
	bz_arr = 0*lat
	btotal_arr = 0*lat
	for i in range(0,len(lat)):
		##Make the call to the IGRF module (assuming date is a constant)
		result = IGRF.igrf_value(lat[i],lon[i],alt[i],date)
		##Parse out the result (mag field strength is in nT)
		bx_arr[i] = result[3]
		by_arr[i] = result[4]
		bz_arr[i] = result[5]
		btotal_arr[i] = result[6]
	return bx_arr,by_arr,bz_arr,btotal_arr

if __name__ == '__main__':
	print('Running Example Script')

	lat = 0.0 ##right at the equator
	lon = 0.0 ##right at the prime meridian
	alt = 100.0 #I'm assuming this is in kilometers above the surface of the Earth
	date = 2021 ##I'm assuming this is the date in a int format and just the year

	##Make the call to the IGRF module
	result = IGRF.igrf_value(lat,lon,alt,date)

	##Parse out the result (mag field strength is in nT)
	bx = result[3]
	by = result[4]
	bz = result[5]
	btotal = result[6]

	print('B-field = ',bx,by,bz,btotal)

	##Plot out B-field as a function of a specific orbit
	#Using custom Earth_Orbit module
	rp = 400.0 #km
	ra = 160000.0 #km
	orbit = orb.Earth_Orbit(ra,rp)

	#Make the plots of the orbit
	orbit.make_plots()

	##First Plot Altitude as a function of Nu
	plt.figure()
	plt.plot(orbit.nu,orbit.alt/1000.0)
	plt.xlabel('True Anamoly (rad)')
	plt.ylabel('Altitude (km)')
	plt.grid()

	##Now feed altitude vector plus latitude and longitude vectors to 
	##Custom function. Assume that orbit is planar (i = 0, around the equator)
	#multipling the altitude by 0 just makes an array of zeros that's the same length
	#as altitude
	lat_arr = 0*orbit.alt #latitude is constant
	lon_arr = orbit.nu*180.0/np.pi ##longitude is the same as True Anamoly for planar orbits
	date = 2021
	bx_arr,by_arr,bz_arr,btotal_arr = bfield_arr(lat_arr,lon_arr,orbit.alt/1000.0,date)

	##Now plot the result
	plt.figure()
	plt.plot(orbit.nu,bx_arr,'b-',label='X Component')
	plt.plot(orbit.nu,by_arr,'r-',label='Y Component')
	plt.plot(orbit.nu,bz_arr,'g-',label='Z Component')
	plt.plot(orbit.nu,btotal_arr,'m-',label='Total Component')
	plt.xlabel('True Anamoly (rad)')
	plt.ylabel('B-Field (nT)')
	plt.grid()
	plt.legend()

	plt.show()
