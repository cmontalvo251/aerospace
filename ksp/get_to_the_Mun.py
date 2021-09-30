import numpy as np
import matplotlib.pyplot as plt

def atmosphere_model(altitude):
    global altx,deny
    return np.interp(altitude,altx,deny)

def muplanet(mass):
    G = 6.6742*10**-11; #%%Gravitational constant
    return G*mass

def mun_parameters():
    a = 12000000 #meters
    e = 0
    vmun = 543.0
    mmun = 9.7599066*10**20 #kg
    G = 6.6742*10**-11; #%%Gravitational constant
    muMun = G*mmun
    return a,e,vmun

def kerbin_parameters():
    ##https://wiki.kerbalspaceprogram.com/wiki/Kerbin
    ##A few things you can compute to compare to the wiki above
    #mu 
    #surface gravity
    #rotational velocity
    
    G = 6.6742*10**-11; #%%Gravitational constant
    Mkerbin = 5.2915158*10**22 #
    muKerbin = G*Mkerbin
    Rkerbin = 600000 #meters
    sidereal_period = 21549.425
    sidereal_angular_velocity = 2*np.pi/sidereal_period
    sidereal_rotational_velocity = sidereal_angular_velocity*Rkerbin
    surface_gravity = muKerbin*Rkerbin/Rkerbin**3
    return sidereal_rotational_velocity,muKerbin,surface_gravity,Rkerbin

def orbit_parameters(ra,rp,mass):
    a = (rp+ra)/2
    e = 1-(rp/a)
    p = a*(1-e**2)
    mu = muplanet(mass)
    h = np.sqrt(mu*p)
    vp = h/rp
    va = h/ra
#    print(a,e,p,va,vp)
    return va,vp

##Mass of Kerbin
Mkerbin = 5.2915158*10**22 #

##Mass of Mun
mmun = 9.7599066*10**20 #kg

###First our parking orbit
apoaps_parking = 730000 #m
periaps_parking = 670000 #m 
dV_parking = 3200 #m/s
v_aps,v_peri = orbit_parameters(apoaps_parking,periaps_parking,Mkerbin)
print('Parking Orbit = ',apoaps_parking,v_aps,periaps_parking,v_peri)

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