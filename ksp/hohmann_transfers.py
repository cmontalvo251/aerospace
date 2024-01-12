import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

##Kerbin Parameters
G = 6.6742*10**-11; #%%Gravitational constant
M = 5.972*10**24  #Earth
#Mkerbin = 5.2915158*10**22 #
#MMun = 9.7599066*10**20 #
#muKerbin = G*Mkerbin
#muMun = G*MMun
mu = G*M
#Rkerbin = 600000. #meters Kerbin
R = 6371000 #meters Earth
#RMun = 200000 #meters Mun

##True Anamoly
nu = np.linspace(0,2*np.pi,1000)
##Type in your apogee and perigee
peri_AGL = 500000.
rp = R + peri_AGL
apo_AGL = 5000000.
ra = R + apo_AGL
#ra = rp
a = (ra+rp)/2.
#alt_AGL = 6000
#a = RMun + alt_AGL
##Eccentricity
e = (ra - rp)/(ra+rp) 
print('Eccentricity = ',e)
##inclination
i = 45.0*np.pi/180.0 ##Drew's random satellite he wants a just slightly over polar retrograde orbit
###Longitude of the Ascending Node
W = 45*np.pi/180.0
#Argument of the periaps
w = 45.*np.pi/180.0

##plot in the orbital plane
###Phat and Qhat
p = a*(1-e**2)
r = p/(1+e*np.cos(nu))
xp = r*np.cos(nu)
yq = r*np.sin(nu)

plt.plot(xp,yq,'r-')
plt.plot(xp[0],yq[0],'r*',markersize=20)
theta = np.linspace(0,2*np.pi,100)
xplanet = R*np.cos(theta)
yplanet = R*np.sin(theta)
plt.plot(xplanet,yplanet,'b-')
#xMun = RMun*np.cos(theta)
#yMun = RMun*np.sin(theta)
#plt.plot(xMun,yMun,'b-')
#plt.axis('equal')
plt.title('Orbital Plane')
plt.xlabel('X')
plt.ylabel('Y')

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

u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
xsph = np.cos(u)*np.sin(v)
ysph = np.sin(u)*np.sin(v)
zsph = np.cos(v)
ax.plot_surface(R*xsph,R*ysph,R*zsph,color='blue',zorder=1)
ax.plot(xi,yj,zk,'r-',zorder=0)
ax.scatter(xi[0],yj[0],zk[0],'r*',s=20)
#ax.plot_wireframe(RMun*xsph,RMun*ysph,RMun*zsph,color='grey')
#ax.axis('square')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

###Now let's plot velocity
vx = np.sqrt(mu/p)*(-np.sin(nu))
vy = np.sqrt(mu/p)*(e+np.cos(nu))
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

