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

def circular_velocity(mu,a):
    return np.sqrt(mu/a)

def mun_parameters():
    a = 12000000 #meters
    e = 0
    i = 0
    R = 200000
    vside,muKerbin,gKerbin,RKerbin = kerbin_parameters()
    v = circular_velocity(muKerbin,a)
    m = 9.7599066*10**20 #kg
    mu = muplanet(m)
    return a,e,v,R,mu

def minmus_parameters():
    a = 47000000
    e = 0
    i = 6 
    vside,muKerbin,gKerbin,RKerbin = kerbin_parameters()
    v = circular_velocity(muKerbin,a)
    m = 2.6457580*10**19
    R = 60000
    mu = muplanet(m)
    return a,e,v,R,mu

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
    print('Kerbin Sidereal Rotational Velocity = ',sidereal_rotational_velocity)
    surface_gravity = mu*R/R**3
    return sidereal_rotational_velocity,mu,surface_gravity,R

def orbit_parameters(rah,rph,mu,R):
    print('Computing Orbital Elements for Orbit...')
    ra = rah*1000 + R
    rp = rph*1000 + R
    print('Apoaps = ',ra,' Periaps = ',rp)
    a = (rp+ra)/2
    print('Semi major axis = ',a)
    e = 1-(rp/a)
    print('Eccentricity = ',e)
    p = a*(1-e**2)
    h = np.sqrt(mu*p)
    vp = h/rp
    print('Velocity at Periaps = ',vp)
    va = h/ra
    print('Velocity at Apoaps = ',va)
#    print(a,e,p,va,vp)
    return va,vp,ra,rp

def get_orbit(ra,rp,mu):
    print('Comptuting Orbit with Apoaps = ',ra,' and Periaps = ',rp)
    ##True Anamoly
    nu = np.linspace(0,2*np.pi,100)
    ##Eccentricity
    e = (ra - rp)/(ra+rp)
    print('Eccentricity = ',e)
    ##Semi Major Axis
    a = (ra+rp)/2.0
    print('Semi Major Axis = ',a)
    ##Parameter
    p = a*(1-e**2)
    ##Angular Momentum
    h = np.sqrt(mu*p)
    ##Velocity at Peri
    vp = h/rp
    print('Velocity at Periaps = ',vp)
    ##Velocity at Apo 
    va = h/ra
    print('Velocity at Apoasp = ',va)
    r = p/(1+e*np.cos(nu))
    xp = r*np.cos(nu)
    yq = r*np.sin(nu)
    return xp,yq,vp,va

############START PROGRAM HERE##########

###SETUP PLOTTING
pdfhandle = PdfPages('plots.pdf')

###First our parking orbit
apoaps_parking_altitude = 90.0 #km
periaps_parking_altitude = 90.0 #km 
vsideKerbin,muKerbin,gKerbin,RKerbin = kerbin_parameters()
v_aps_parking,v_peri_parking,a_aps,a_peri = orbit_parameters(apoaps_parking_altitude,periaps_parking_altitude,muKerbin,RKerbin)
print('Parking Orbit = ',a_aps,a_peri,v_aps_parking,v_peri_parking)
dvKerbin = (v_aps_parking-vsideKerbin)/(2./3.)
print('DV = ',dvKerbin) ##The 2/3rd is because of the atmosphere

###PLOT KERBIN IN BLUE
plt.figure()
x_parking,y_parking,va,vp = get_orbit(a_aps,a_peri,muKerbin)
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

###TRANSFER ORBIT
aMun,eMun,vMun,RMun,muMun = mun_parameters()
x_transfer,y_transfer,v_peri_transfer,v_apo_transfer = get_orbit(aMun,a_peri,muKerbin)
dvtransfer = v_peri_transfer - v_aps_parking
print('DV for Transfer Orbit = ',dvtransfer)

###Mun's Orbit
x_Mun,y_Mun,v_peri_Mun,v_apo_Mun = get_orbit(aMun,aMun,muKerbin)

###Plot Tranfer orbit
plt.figure()
plt.plot(x_transfer/1000,y_transfer/1000,'m-',label='Transfer Orbit')
plt.plot(xkerbal/1000,zkerbal/1000,'b-',label='Kerbin')
plt.plot(x_parking/1000,y_parking/1000,'r-',label='Parking Orbit')
plt.plot(x_Mun/1000,y_Mun/1000,'k-',label='Mun Orbit')
plt.legend()
plt.grid()
plt.xlabel('X (km)')
plt.ylabel('Y (km)')
print('Creating Transfer Orbit Plot')
pdfhandle.savefig()

###Parking Orbit around Mun
apoaps_parking_altitude_Mun = 10.0 #km
periaps_parking_altitude_Mun = 10.0 #km
v_aps_parking_Mun,v_peri_parking_Mun,a_aps_Mun,a_peri_Mun = orbit_parameters(apoaps_parking_altitude_Mun,periaps_parking_altitude_Mun,muMun,RMun)
x_Mun_parking,y_Mun_parking,v_peri_parking_Mun,v_apo_parking_Mun = get_orbit(a_aps_Mun,a_peri_Mun,muMun)

plt.figure()
plt.plot(x_Mun_parking,y_Mun_parking,'k-',label='Parking Orbit')
xMun = RMun*np.cos(thetakerbal)
zMun = RMun*np.sin(thetakerbal)
plt.plot(xMun,zMun,color='grey',label='Mun')
plt.legend()
plt.grid()
print('Creating Mun Parking Orbit Plot')
pdfhandle.savefig()

##DV required for parking around Mun
vretro = v_peri_Mun - v_apo_transfer
print('Speed of Spacecraft in Retrograde = ',vretro)
dVMun = v_aps_parking_Mun - vretro
print('dV for Parking around Mun = ',dVMun)

###TOTAL DV
dV = dvKerbin + 2*dvtransfer + dVMun + 2*v_aps_parking_Mun
print('Total dV Required for Mission = ',dV)

####END PROGRAM
pdfhandle.close()
if sys.platform == 'linux2' or sys.platform == 'linux':
    print('Opening PDF')
    os.system('evince plots.pdf &')
sys.exit()