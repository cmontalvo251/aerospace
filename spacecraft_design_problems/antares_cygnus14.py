import numpy as np
from mayavi import mlab
import matplotlib.pyplot as plt

###Alright Let's plot an orbit around Kerbin using this new Mayavi Library
##Kerbin Parameters
G = 6.6742*10**-11; #%%Gravitational constant
Mkerbin = 5.2915158*10**22 #
MEarth = 5.972e24 #kg
#MMun = 9.7599066*10**20 #
muKerbin = G*Mkerbin
muEarth = G*MEarth
#muMun = G*MMun
Rkerbin = 600000. #meters
REarth = 6.371e6 #meters
#RMun = 200000 #meters

##True Anamoly
nu = np.linspace(0,2*np.pi,100)
##Semi Major Axis of an 80 km parking orbit
ra = REarth + 254.0*1000.
rp = REarth + 182.0*1000.
#ra = 12000000
a = (ra+rp)/2.
#alt_AGL = 6000
#a = RMun + alt_AGL
##Eccentricity
e = (ra - rp)/(ra+rp)
print('Eccentricity = ',e)
##inclination
i = 51.6*np.pi/180.0 ##Drew's random satellite he wants a just slightly over polar retrograde orbit
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
plt.plot(xp[0],yq[0],'r*',markersize=20)
theta = np.linspace(0,2*np.pi,100)
#xkerbin = Rkerbin*np.cos(theta)
#ykerbin = Rkerbin*np.sin(theta)
xearth = REarth*np.cos(theta)
yearth = REarth*np.sin(theta)
#plt.plot(xkerbin,ykerbin,'b-')
plt.plot(xearth,yearth,'b-')
#xMun = RMun*np.cos(theta)
#yMun = RMun*np.sin(theta)
#plt.plot(xMun,yMun,'b-')
plt.axis('equal')
plt.title('Orbital Plane')

###Rotate to Earth Centered Inertial Frame (ECI)
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


###Let's plot velocity
#vx = np.sqrt(muKerbin/p)*(-np.sin(nu))
#vy = np.sqrt(muKerbin/p)*(e+np.cos(nu))
#vx = np.sqrt(muMun/p)*(-np.sin(nu))
#vy = np.sqrt(muMun/p)*(e+np.cos(nu))
vx = np.sqrt(muEarth/p)*(-np.sin(nu))
vy = np.sqrt(muEarth/p)*(e+np.cos(nu))
v = np.sqrt(vx**2 + vy**2)

plt.figure()
plt.plot(nu,vx,label='Vx')
plt.plot(nu,vy,label='Vy')
plt.plot(nu,v,label='V')
plt.grid()
plt.xlabel('True Anomaly (rad)')
plt.ylabel('Velocity (m/s)')
plt.legend()

###From here we can plot the orbit in 3D using matplotlib.
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
xsph = np.cos(u)*np.sin(v)
ysph = np.sin(u)*np.sin(v)
zsph = np.cos(v)
ax.plot_surface(REarth*xsph,REarth*ysph,REarth*zsph,color='blue',zorder=1)
ax.plot(xi,yj,zk,'r-',zorder=0)
ax.scatter(xi[0],yj[0],zk[0],'r*',s=20)
#ax.plot_wireframe(RMun*xsph,RMun*ysph,RMun*zsph,color='grey')
#ax.axis('square')

plt.show()

##Unfortunately When you use matplotlib the orbit does not plot properly.
##So instead we will use this mayavi module
##First let's plot the orbit of the satellite
#dphi, dtheta = np.pi/250.0, np.pi/250.0
#theta0 = np.arange(0, 2*np.pi+dtheta*1.5, dtheta) #
#orbit_radius = 5
#xorb = orbit_radius * np.cos(theta0)
#yorb = orbit_radius * np.sin(theta0)
#zorb = 0 * np.cos(theta0)
#s = mlab.plot3d(xorb,yorb,zorb)
f = 100000
s=mlab.plot3d(xi/f,yj/f,zk/f,color=(1.0,0.0,0.0))

##Then let's plot a sphere the size of the satellite
rsat = 0.05*REarth # km?
s=mlab.mesh((rsat*xsph+xi[0])/f,(rsat*ysph+yj[0])/f,(rsat*zsph+zk[0])/f,color=(0.7,0.7,0.7))

##Finally let's plot kerbin in blue at the center of the orbit
s=mlab.mesh(REarth*xsph/f,REarth*ysph/f,REarth*zsph/f,color=(0.0,0.9,0.0))

##Then we set the proper view port
#mlab.view(azimuth=-0, elevation=45, distance=Rkerbin*10)
mlab.view(azimuth=-0, elevation=45, distance=5*ra/f)
mlab.orientation_axes()
#mlab.savefig(f'moonearth/Above{int(time*10):04d}.png')
mlab.show()
