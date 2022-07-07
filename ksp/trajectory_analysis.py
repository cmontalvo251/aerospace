import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import sys
import os

################BEGIN SUBROUTINES############

def atmosphere_model(altitude):
    global altx,deny
    return np.interp(altitude,altx,deny)

def muplanet(mass):
    G = 6.6742*10**-11; #%%Gravitational constant
    return G*mass

def mun_parameters():
    a = 12000000 #meters
    e = 0
    i = 0
    v = 543.0
    m = 9.7599066*10**20 #kg
    mu = muplanet(m)
    return a,e,v

def minmus_parameters():
    a = 47000000
    e = 0
    i = 6 
    v = 274.0
    m = 2.6457580*10**19
    R = 60000
    mu = muplanet(m)
    return a,e,v,R

def kerbin_parameters():
    ##https://wiki.kerbalspaceprogram.com/wiki/Kerbin
    ##A few things you can compute to compare to the wiki above
    #mu 
    #surface gravity
    #rotational velocity
    m = 5.2915158*10**22 #
    mu = muplanet(m)
    R = 600000 #meters
    sidereal_period = 21549.425 ##seconds
    sidereal_angular_velocity = 2*np.pi/sidereal_period
    sidereal_rotational_velocity = sidereal_angular_velocity*R
    surface_gravity = mu*R/R**3
    return sidereal_rotational_velocity,mu,surface_gravity,R

def orbit_parameters(rah,rph,mu,R):
    ra = rah*1000 + R
    rp = rph*1000 + R
    a = (rp+ra)/2
    e = 1-(rp/a)
    p = a*(1-e**2)
    h = np.sqrt(mu*p)
    vp = h/rp
    va = h/ra
#    print(a,e,p,va,vp)
    return va,vp,ra,rp

def get_orbit(ra,rp):
    ##True Anamoly
    nu = np.linspace(0,2*np.pi,100)
    e = (ra - rp)/(ra+rp)
    a = (ra+rp)/2.0
    p = a*(1-e**2)
    r = p/(1+e*np.cos(nu))
    xp = r*np.cos(nu)
    yq = r*np.sin(nu)
    return xp,yq

############START PROGRAM HERE##########

###SETUP PLOTTING
pdfhandle = PdfPages('plots.pdf')

###First our parking orbit
apoaps_parking_altitude = 90.0 #km
periaps_parking_altitude = 90.0 #km 
vsideKerbin,muKerbin,gKerbin,RKerbin = kerbin_parameters()
v_aps,v_peri,a_aps,a_peri = orbit_parameters(apoaps_parking_altitude,periaps_parking_altitude,muKerbin,RKerbin)
print('Parking Orbit = ',a_aps,a_peri,v_aps,v_peri)

###PLOT KERBIN IN BLUE
plt.figure()
x_parking,y_parking = get_orbit(a_aps,a_peri)
plt.plot(x_parking,y_parking,'r-',label='Parking Orbit')
thetakerbal = np.linspace(0,2*np.pi,1000)
xkerbal = RKerbin*np.cos(thetakerbal)
zkerbal = RKerbin*np.sin(thetakerbal)
plt.plot(xkerbal,zkerbal,label='Kerbin')
plt.legend()
plt.grid()
plt.xlabel('X (m)')
plt.ylabel('Y (m)')
print('Creating Kerbin Parking Orbit Plot')
pdfhandle.savefig()

####END PROGRAM
pdfhandle.close()
if sys.platform == 'linux2' or sys.platform == 'linux':
    print('Opening PDF')
    os.system('evince plots.pdf &')
sys.exit()

##Transfer orbit
amun,emun,vmun = mun_parameters()
apoaps_transfer = amun
periaps_transfer = periaps_parking
v_aps_transfer,v_peri_transfer = orbit_parameters(apoaps_transfer,periaps_transfer,Mkerbin)
print('Transfer Orbit = ',apoaps_transfer,v_aps_transfer,periaps_transfer,v_peri_transfer)
dv_transfer2 = v_peri_transfer - v_peri
print('Transfer DV = ',dv_transfer2)

###Mun's Orbit
v_am,v_pm = orbit_parameters(amun,amun,Mkerbin)
print('Mun Orbit = ',amun,v_am,amun,v_pm)

##Velocity at Perigee
v_mun_perigee_relative = v_am - v_aps_transfer
print('Velocity Relative to Mun = ',v_mun_perigee_relative)

##Parking Orbit around the Mun
Rmun = 200000
mun_aps = 800000 + Rmun
mun_peri = 143000 + Rmun
v_aps_mun,v_peri_mun = orbit_parameters(mun_aps,mun_peri,mmun)
print('Mun Parking = ',mun_aps,v_aps_mun,mun_peri,v_peri_mun)

##Third DV Burn
dv3_mun = v_peri_mun - v_mun_perigee_relative
print('DV3 = ',dv3_mun)

##Fourth DV Burn - To escape from Mun
dv4_mun = dv3_mun

##Last burn to bring you back to a parking orbit
dv5_parking = dv_transfer2

###Overall DV
print('Overall DV = ',dv_transfer2+dV_parking+dv3_mun+dv4_mun+dv5_parking)