import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

##Kerbin Parameters

##True Anamoly
nu = np.linspace(0,2*np.pi,100)
##Semi Major Axis
a = 13599840256. #meters
##Eccentricity
e = 0.8
##inclination
i = 45*np.pi/180.0
###Longitude of the Ascending Node
W = 45*np.pi/180.0
#Argument of the periaps
w = 0.

##plot in the orbital plane
###Phat and Qhat
p = a*(1-e**2)
r = p/(1+e*np.cos(nu))
xp = r*np.cos(nu)
yq = r*np.sin(nu)

plt.plot(xp,yq)
plt.plot(xp[0],yq[0],'b*',markerSize=20)
plt.plot(0,0,'ys',markerSize=20)
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
ax.set_title('Look at this f#@$&*@#@ing sweet plot!!!!')
#ax.axis('equal')
    
plt.show()
