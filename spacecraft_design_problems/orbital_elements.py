#!/usr/bin/python

import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
try:
    from pdf import *
    usepdf = 1
except:
    print('No module pdf')
    usepdf = 0

if usepdf:
    pp = PDF(1,plt)

##Position and velocity
#Problem from Book
r = np.asarray([-.707,0.707,0])
v = np.asarray([0,0.5,0.0])
    
#Perfect circular orbit
#r = np.asarray([1,0,0])
#v = np.asarray([0,1,0])
    
#45 degree inclined circular orbit
#r = np.asarray([1,0,0])
#v = np.asarray([0,0.707,0.707])

#Eccentric Planar Orbit
#r = np.asarray([1,0,0])
#v = np.asarray([0,1.3,0])

#Ecentric Orbit rotated in plane 45 degrees
#rp = np.sqrt(2)/2.0
#r = np.asarray([rp,rp,0])
#v = np.asarray([-1.3*rp,1.3*rp,0])

#Angular Momentum
h = np.cross(r,v)
hnorm = np.linalg.norm(h)
print('h = ',h)

#Line of Nodes
K = np.asarray([0,0,1])
n = np.cross(K,h)
print('n=',n)

##Eccentricity
mu = 1.0
vnorm = np.linalg.norm(v)
rnorm = np.linalg.norm(r)
e = 1.0/mu*((vnorm**2-mu/rnorm)*r - np.dot(r,v)*v)
print('e = ',e) #This points to the apoapsis point
enorm = np.linalg.norm(e)

##Parameter
p = hnorm**2/mu
print('p=',p)

#Inclination
inc = np.arccos(h[2]/hnorm)
print('inc = ',inc)

##Longitude of Ascending Node
nnorm = np.linalg.norm(n)
if abs(nnorm) < 1e-10:
    print('2-D Orbit, Longitude of Ascending Node is undefined')
    I = np.asarray([1,0,0])
    #Longitude of Periapsis - eccentricicty vector points towards periapsis point
    bigPI = np.arccos(np.dot(e,I)/enorm)
    #Since Longitude of ascending node and argument of the periapsis are undefined we simply let
    OHM = 0.0
    omega = bigPI #arument of periapsis and longitude of periapsis are the same in this case
else:
    OHM = np.arccos(n[0]/nnorm)
    omega = np.dot(n,e)/(nnorm*enorm)
    bigPI = OHM + omega
    
print('Longitude of Periapsis = ',bigPI)
    
##Plot orbital plane
nu = np.linspace(0,2*np.pi,1000)
rorbit = p/(1+enorm*np.cos(nu))
rx = rorbit*np.cos(nu)
ry = rorbit*np.sin(nu)
plt.plot(rx,ry)
plt.plot(0,0,'y*')
plt.title('Orbital Plane')
plt.grid()
plt.axis('equal')
if usepdf:
    pp.savefig()

##Now let's rotate the entire orbit to the ecliptic plane
xecl = (np.cos(omega)*np.cos(OHM) - np.sin(omega)*np.sin(OHM)*np.cos(inc)) * rx + (-np.sin(omega)*np.cos(OHM)-np.cos(omega)*np.sin(OHM)*np.cos(inc))*ry
yecl = (np.cos(omega)*np.sin(OHM) + np.sin(omega)*np.cos(OHM)*np.cos(inc)) * rx + (-np.sin(omega)*np.sin(OHM)+np.cos(omega)*np.cos(OHM)*np.cos(inc))*ry
zecl = (np.sin(omega)*np.sin(inc))*rx + (np.cos(omega)*np.sin(inc))*ry

fig = plt.figure('Ecliptic')
#ax = fig.add_subplot(111,projection='3d')
ax = Axes3D(fig)
ax.plot(xecl,yecl,zecl, color = 'blue', linestyle = 'solid')
plt.title('Ecliptic')
plt.grid()
#plt.axis('equal')
ax.plot([0.],[0.],[0.],'y*',markersize=5)
ax.plot([r[0]],[r[1]],[r[2]],'r*',markersize=5)
if usepdf:
    pp.savefig()

if usepdf:
    pp.close()
else:
    plt.show()
