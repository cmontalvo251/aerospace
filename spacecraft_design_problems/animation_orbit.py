import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

##Kerbin Parameters
G = 6.6742*10**-11; #%%Gravitational constant
Mkerbin = 5.2915158*10**22 #
#MMun = 9.7599066*10**20 #
muKerbin = G*Mkerbin
#muMun = G*MMun
Rkerbin = 600000. #meters
#RMun = 200000 #meters

##True Anamoly
nu = np.linspace(0,2*np.pi,100)
##Semi Major Axis of an 80 km parking orbit
alt_AGL = 80000.
a = Rkerbin + alt_AGL
#alt_AGL = 6000
#a = RMun + alt_AGL
##Eccentricity
e = 0.0
##inclination
i = 98.0*np.pi/180.0 ##Drew's random satellite he wants a just slightly over polar retrograde orbit
###Longitude of the Ascending Node
W = 0*45*np.pi/180.0
#Argument of the periaps
w = 0.

##plot in the orbital plane
###Phat and Qhat
p = a*(1-e**2)
r = p/(1+e*np.cos(nu))
xp = r*np.cos(nu)
yq = r*np.sin(nu)

theta = np.linspace(0,2*np.pi,100)
xkerbin = Rkerbin*np.cos(theta)
ykerbin = Rkerbin*np.sin(theta)

###Rotate to Kerbin Centered Inertial Frame (KCI)
zr = 0*xp

TPI = np.asarray([[np.cos(W)*np.cos(w)-np.sin(W)*np.sin(w)*np.cos(i),-np.cos(W)*np.sin(w)-np.sin(W)*np.cos(w)*np.cos(i),np.sin(W)*np.sin(i)],
                 [np.sin(W)*np.cos(w)+np.cos(W)*np.sin(w)*np.cos(i),-np.sin(W)*np.sin(w)+np.cos(W)*np.cos(w)*np.cos(i),-np.cos(W)*np.sin(i)],
                 [np.sin(w)*np.sin(i),np.cos(w)*np.sin(i),np.cos(i)]])

xi = 0*xp
yj = 0*yq
zk = 0*zr
for x in range(0,len(xp)):
  xyzO = np.asarray([xp[x],yq[x],zr[x]]) ##3x1 vector
  xyzi = np.matmul((TPI),xyzO)
  xi[x] = xyzi[0]
  yj[x] = xyzi[1]
  zk[x] = xyzi[2]
  
##Convert xyz to lat/lon/alt
rho = np.sqrt((pow(xi, 2) + pow(yj, 2) + pow(zk, 2)));
phi = np.arccos(zk/rho)
the = np.arctan2(yj , xi)
lat = 90 - phi*(180.0/ np.pi);
lon = the*(180.0/np.pi);
h = rho-Rkerbin;

#fig = plt.figure()
#plt.plot(nu,lat)
#plt.plot(nu,lon)
#plt.show()

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
xsph = np.cos(u)*np.sin(v)
ysph = np.sin(u)*np.sin(v)
zsph = np.cos(v)
for i in range(0,len(xi)):
    ax.plot_wireframe(Rkerbin*xsph,Rkerbin*ysph,Rkerbin*zsph,color='blue')
    ax.text(9, 0, 0, str(lat[i]), color='red')
    ax.text(100,0,0,str(lon[i]),color='green')
    ax.plot(xi,yj,zk,'r-')
    ax.scatter(xi[i],yj[i],zk[i],'r',marker='s')
    plt.pause(0.01)
    plt.cla()
    
