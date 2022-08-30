#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import sys
import os

##In order to import this toolbox into a python script you need to 
##do the following. Copy the following lines of code below
# import sys
# sys.path.append('/home/carlos/Dropbox/BlackBox/plotting')
# from gps import *

# or 
# In order to get python to search for all of your lovely blackbox 
# python routines. Add this to your .bashrc file

# for d in /home/carlos/Dropbox/BlackBox/*/; do
# 	PYTHONPATH+=:$d
# done
# export PYTHONPATH

# In order to get this to work in Thonny you need to navigate to
# /home/user/.thonny/BundledPython36/lib/python3.6/site-packages and place
# a symbolic link here

# In Enthough you need to make symbolic links here
# /home/carlos/.local/share/canopy/edm/envs/User/lib/python2.7/site-packages

class GPSCOMMS():
    def __init__(self):
        self.NM2FT=6076.115485560000
        self.FT2M=0.3048
        self.G = 6.6742*10**-11; #%%Gravitational constant
        self.GPSVAL = 60.0*self.NM2FT*self.FT2M
        self.EarthRadius_ae = 6378140. ##meters mean eq radius
        self.MEarth = 5.97219e24 #kg
        self.muEarth = self.G*self.MEarth
        #23 hours 56 minutes 4 seconds
        self.solar_day_sec = 23.*3600.0 + 56.*60.0 + 4.0
        self.wEarth = 2*np.pi/self.solar_day_sec
        self.elements = False

    def NMEA_LAT_LON(self,IN):
        l = np.float(IN)/100.0;
        l_deg = np.floor(l);
        l_min = (l-l_deg)*100.0;
        out = l_deg + l_min/60.0;
        return out

    def NMEA_TIME(self,time_raw,units):
        #try splitting by :
        times = time_raw.split(':')
        if len(times) == 1:
            #split didnt work which means there are no colons
            #print time_raw
            hour = np.float(time_raw[0:2])
            minute = np.float(time_raw[2:4])
            sec = np.float(time_raw[4:6])
            msec = 0
            #print hour,minute,sec
        else:
            hour = np.float(times[0])
            minute = np.float(times[1])
            sec = np.float(times[2])
            try:
                msec = np.float(times[3])
            except:
                msec = 0;
        time = 0
        if units == 'sec':
            time = msec/1000 + sec + minute*60 + hour*3600
        elif units == 'hrs':
            time = hour + minute/60 + sec/3600 + (msec/1000)/3600
        else:
            print('Invalid Unit in NMEA_TIME. Returning 0')
        #print time
        return time

    def computeGeocentricLATLON(self,x,y,z,t):
        xprime = np.sqrt(x**2 + y**2)
        zprime = z
        latitude = np.arctan2(zprime,xprime)*180/np.pi
        longitude = np.asarray(np.arctan2(y,x)*180.0/np.pi)-180.0/np.pi*self.wEarth*t
        #longitude[longitude<0]+=360.0
        #norm = np.sqrt(x**2 + y**2 + z**2)
        #phi = np.arccos(self.zsat / rho)
        #the = np.arctan2(self.ysat,self.xsat);
        #self.latitude = 90 - phi*(180 /np.pi);
        #self.longitude = the*(180/np.pi);
        return latitude,longitude

    def IFOV(self,swath,deltas,az_rad):
        cosdtprime = np.cos(swath)*np.sin(deltas)+np.sin(swath)*np.cos(deltas)*np.cos(az_rad)
        dtprime = np.arccos(cosdtprime)
        dt = 90.0 - dtprime * 180.0/np.pi
        dt_rad = dt*np.pi/180.0
        cosdL = np.asarray((np.cos(swath)-np.sin(deltas)*np.sin(dt_rad))/(np.cos(deltas)*np.cos(dt_rad)))
        cosdL[abs(cosdL)>1.0] = 1.0
        dL = np.arccos(cosdL)*180/np.pi
        return dt,dL

    def LATLON2Cartesian(self,lat,lon,alt):
        z = alt*np.sin(lat*np.pi/180.0)
        xyprime = alt*np.cos(lat*np.pi/180.0)
        x = xyprime*np.cos(lon*np.pi/180.0)
        y = xyprime*np.sin(lon*np.pi/180.0)
        return x,y,z

    def convertXY2LATLON(self,xy,origin):
        x = xy[0]
        y = xy[1]
        Ox = origin[0]
        Oy = origin[1]
        lat = x/self.GPSVAL + Ox
        lon = y/(self.GPSVAL*np.cos(Ox*np.pi/180)) + Oy
        return np.asarray([lat,lon])

    def convertLATLON(self,lat_lon,origin):
        lat = lat_lon[0]
        lon = lat_lon[1]
        Ox = origin[0]
        Oy = origin[1]
        x = (lat-Ox)*self.GPSVAL
        y = (lon-Oy)*self.GPSVAL*np.cos(Ox*np.pi/180)
        return np.asarray([x,y])

    def plotEARTH(self,ax):
        u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
        xsph = np.cos(u)*np.sin(v)*self.EarthRadius_ae
        ysph = np.sin(u)*np.sin(v)*self.EarthRadius_ae
        zsph = np.cos(v)*self.EarthRadius_ae
        ax.plot_wireframe(xsph,ysph,zsph,color='green')
        return

    def Derivatives(self,statei):
        xi = statei[0]
        yi = statei[1]
        zi = statei[2]
        ui = statei[3]
        vi = statei[4]
        wi = statei[5]
        xdoti = ui
        ydoti = vi
        zdoti = wi
        r = np.sqrt(xi**2 + yi**2 + zi**2)
        xddoti = -self.muEarth*xi/r**3
        yddoti = -self.muEarth*yi/r**3
        zddoti = -self.muEarth*zi/r**3
        return np.asarray([xdoti,ydoti,zdoti,xddoti,yddoti,zddoti])

    def sixdof_orbit(self,x0,y0,z0,u0,v0,w0):
        #We want to simulate until we complete 1 orbit
        orbit = True

        #Setup vectors
        x = [x0]
        y = [y0]
        z = [z0]
        u = [u0]
        v = [v0]
        w = [w0]
        t = [0]
        xi = x0
        yi = y0
        zi = z0
        ui = u0
        vi = v0
        wi = w0
        ti = 0.
        statei = np.asarray([xi,yi,zi,ui,vi,wi])
        print('Initial State = ',statei)
        timestep = 1.0
        
        while orbit:
            ###RK4 Function Calls
            k1 = self.Derivatives(statei)
            k2 = self.Derivatives(statei+k1*timestep/2)
            k3 = self.Derivatives(statei+k2*timestep/2)
            k4 = self.Derivatives(statei+k3*timestep)
            phi = (1./6.)*(k1 + 2*k2 + 2*k3 + k4)
            statei = statei + phi*timestep
            ti = ti + timestep
            
            ##Need some way to determine if we've completed 1 orbit
            xi = statei[0]
            yi = statei[1]
            zi = statei[2]
            ui = statei[3]
            vi = statei[4]
            wi = statei[5]
            r = np.sqrt(xi**2 + yi**2 + zi**2)
            if np.sqrt((x0-xi)**2 + (y0-yi)**2 + (z0 - zi)**2) < 10000.0 and ti > 100:
                print('1 Orbit complete')
                orbit = False
            ##Need a check for inside the Earth
            if r < self.EarthRadius_ae:
                print('Inside Earth')
                orbit = False
            ##Finally we need a fail safe
            if ti > 20000:
                print('Maximum time of 20000 seconds has been reached')
                orbit = False
            #Append Vectors
            x.append(xi)
            y.append(yi)
            z.append(zi)
            u.append(ui)
            v.append(vi)
            w.append(wi)
            t.append(ti)
                
        return x,y,z,u,v,w,t

    def getOrbitalElements(self,x0,y0,z0,u0,v0,w0):
        height_at_perigee_km = 0.0
        ECC = 0.0
        INC = 0.0
        LAN = 0.0
        ARG = 0.0
        return height_at_perigee_km,ECC,INC,LAN,ARG

    def computeOrbitalElements(self,height_at_perogee_km,ECC,INC,LAN,ARG):
        if self.elements == False:
            ##PEROGEE IS GIVEN
            self.rp = height_at_perogee_km*1000 + self.EarthRadius_ae
            #SOLVE FOR APOGEE
            ##ECC = (ra - rp)/(ra + rp)
            self.ra = self.rp*(1+ECC)/(1-ECC)
            ##SEMI MAJOR AXIS
            self.a = (self.ra + self.rp)/2.0
            #inclination in radians
            self.i = INC*np.pi/180.0
            ###Longitude of the Ascending Node in radians
            self.W = LAN*np.pi/180.0
            #Argument of the periaps in radians
            self.w = ARG*np.pi/180.0
            ###The parameter
            self.p = self.a*(1-ECC**2)
            ##Coefficient of velocity
            self.coeff = np.sqrt(self.muEarth/self.p)
            self.elements = True
            ###Rotation matrix to Earth Centered Inertial Frame (ECI)
            self.TPI = np.asarray([[np.cos(self.W)*np.cos(self.w)-np.sin(self.W)*np.sin(self.w)*np.cos(self.i),-np.cos(self.W)*np.sin(self.w)-np.sin(self.W)*np.cos(self.w)*np.cos(self.i),np.sin(self.W)*np.sin(self.i)],
                 [np.sin(self.W)*np.cos(self.w)+np.cos(self.W)*np.sin(self.w)*np.cos(self.i),-np.sin(self.W)*np.sin(self.w)+np.cos(self.W)*np.cos(self.w)*np.cos(self.i),-np.cos(self.W)*np.sin(self.i)],
                 [np.sin(self.w)*np.sin(self.i),np.cos(self.w)*np.sin(self.i),np.cos(self.i)]])
        return
        
    def getStateVector(self,height_at_perogee_km,ECC,INC,LAN,ARG,nu):
        ##First compute the orbital elements if needed
        self.computeOrbitalElements(height_at_perogee_km,ECC,INC,LAN,ARG)

        #Then we compute the following:
        ##Radius
        r = self.p/(1+ECC*np.cos(nu))
        ##Orbital plane x and y
        xp = r*np.cos(nu)
        yq = r*np.sin(nu)
        zr = 0*xp
        ##Velocity
        up = self.coeff*(-np.sin(nu))
        vq = self.coeff*(ECC+np.cos(nu))
        wr = 0*up
        #Rotate to ECI frame
        xyzO = np.asarray([xp,yq,zr]) ##3x1 vector
        uvwO = np.asarray([up,vq,wr])
        xyzi = np.matmul((self.TPI),xyzO)
        uvwi = np.matmul((self.TPI),uvwO)
        xi = xyzi[0]
        yj = xyzi[1]
        zk = xyzi[2]
        ui = uvwi[0]
        vj = uvwi[1]
        wk = uvwi[2]
        
        return xi,yj,zk,ui,vj,wk

    def kepler_orbit(self,height_at_perogee_km,ECC,INC,LAN,ARG):
        ##First compute the orbital elements if needed
        self.computeOrbitalElements(height_at_perogee_km,ECC,INC,LAN,ARG)

        ##True Anamoly
        nu = np.linspace(0,2*np.pi,100)

        #Then run through the true anomaly vector and get the state vector
        xi = 0*nu
        yj = 0*nu
        zk = 0*nu
        ui = 0*nu
        vj = 0*nu
        wk = 0*nu
        for x in range(0,len(xi)):
            xii,yjj,zkk,uii,vjj,wkk = self.getStateVector(height_at_perogee_km,ECC,INC,LAN,ARG,nu[x])
            xi[x] = xii
            yj[x] = yjj
            zk[x] = zkk
            ui[x] = uii
            vj[x] = vjj
            wk[x] = wkk
        return xi,yj,zk,ui,vj,wk,nu*180.0/np.pi


    def HHMM_Format(self,time_vec,del_min):
        ##Assume this time_vec is in format HH.(HH/60)
        time_label = []
        last_time = time_vec[0]-del_min
        ctr = -1
        for time in time_vec:
            hour = int(np.floor(time))
            minute = int(np.round((time-hour)*60))
            #print time,hour,minute
            str_minute = str(minute)
            if (minute < 10):
                str_minute = '0' + str_minute
            ##Check and make sure enough time has passed
            str_time = str(hour)+':'+str_minute
            if (time-last_time)*60.0 >= del_min-1e-2:
                last_time = time
                #print str_time,ctr
                ctr+=1
                time_label.append(str_time)

        #print time_label
        xticks = np.linspace(time_vec[0],time_vec[-1],ctr+1)

        return time_label,xticks

# Copyright - Carlos Montalvo 2016
# You may freely distribute this file but please keep my name in here
# as the original owner
