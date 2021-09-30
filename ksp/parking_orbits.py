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
i = 0*98.0*np.pi/180.0 ##Drew's random satellite he wants a just slightly over polar retrograde orbit
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

plt.plot(xp,yq,'r-')
plt.plot(xp[0],yq[0],'r*',markerSize=20)
theta = np.linspace(0,2*np.pi,100)
xkerbin = Rkerbin*np.cos(theta)
ykerbin = Rkerbin*np.sin(theta)
plt.plot(xkerbin,ykerbin,'b-')
#xMun = RMun*np.cos(theta)
#yMun = RMun*np.sin(theta)
#plt.plot(xMun,yMun,'b-')
plt.axis('equal')
plt.title('Orbital Plane')

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
  
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(xi,yj,zk,'r-')
ax.scatter(xi[0],yj[0],zk[0],'r*',s=20)
u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
xsph = np.cos(u)*np.sin(v)
ysph = np.sin(u)*np.sin(v)
zsph = np.cos(v)
ax.plot_wireframe(Rkerbin*xsph,Rkerbin*ysph,Rkerbin*zsph,color='blue')
#ax.plot_wireframe(RMun*xsph,RMun*ysph,RMun*zsph,color='grey')
#ax.axis('square')

###Now let's plot velocity
vx = np.sqrt(muKerbin/p)*(-np.sin(nu))
vy = np.sqrt(muKerbin/p)*(e+np.cos(nu))
#vx = np.sqrt(muMun/p)*(-np.sin(nu))
#vy = np.sqrt(muMun/p)*(e+np.cos(nu))
v = np.sqrt(vx**2 + vy**2)

plt.figure()
plt.plot(nu,vx,label='Vx')
plt.plot(nu,vy,label='Vy')
plt.plot(nu,v,label='V')
plt.grid()
plt.xlabel('True Anomaly (rad)')
plt.ylabel('Velocity (m/s)')
plt.legend()
    
plt.show()
