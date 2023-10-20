import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

##Kerbin Parameters
##Radius of Kerbin
radius_of_Kerbin = 600000 #meters
#Mu
mass_of_Kerbin = 5.2915158e22
G = 6.67e-11
mu = G*mass_of_Kerbin

##True Anamoly
nu = np.linspace(0,2*np.pi,100)
##Compute my Apoaps and Periaps
altitude_apoaps = 82.109*1000 #meters
altitude_periaps = 79.9871*1000
ra = radius_of_Kerbin + altitude_apoaps
rp = radius_of_Kerbin + altitude_periaps
##Semi Major Axis
a = (ra + rp)/2.0
##Eccentricity
e = (ra - rp)/(ra+rp)
print('e = ',e)
##inclination
i = 4.3*np.pi/180.0
###Longitude of the Ascending Node
W = 242.2*np.pi/180.0
#Argument of the periaps
w = 262.8*np.pi/180.0

###Phat and Qhat
p = a*(1-e**2)
r = p/(1+e*np.cos(nu))

##Angular Momentum
h = np.sqrt(p*mu)

##Velocity at Peri and Apo
vapo = h/ra
vperi = h/rp
print('Velocity = ',vapo,vperi)

##plot in the orbital plane
xp = r*np.cos(nu)
yq = r*np.sin(nu)

plt.plot(xp,yq)
plt.plot(xp[0],yq[0],'b*',markersize=20)
plt.plot(0,0,'ys',markersize=20)
plt.axis('equal')
plt.title('Orbital Plane')

###Rotate to Ecliptic Plane
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
ax.plot(xi,yj,zk)
ax.scatter(xi[0],yj[0],zk[0],'b*',s=20)
ax.scatter(0,0,0,'y*',s=100)
ax.set_title('Look at this sweet plot!!!!')
ax.set_xlim([-ra,ra])
ax.set_ylim([-ra,ra])
ax.set_zlim([-ra,ra])
    
plt.show()
