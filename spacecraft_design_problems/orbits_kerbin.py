import numpy as np
import matplotlib.pyplot as plt

Rkerbin = 600000  #RADIUS of KERBIN
ro = Rkerbin + 80*1000.0 ##The radius of the orbit is the radius of kerhin + the altitude above the surface
v = np.linspace(0,2*np.pi,1000) ###The true anomaly of an orbit goes from 0 to 2pi
xp = np.cos(v)*Rkerbin ##Parametric coordinates.
yp = np.sin(v)*Rkerbin ##of a circle
xo = np.cos(v)*ro ##Go email your math teacher and said you used parametric eqns
yo = np.sin(v)*ro ##in college to do rocket science
plt.plot(xp,yp,'b-') ###Plot a circle that's the same size as the surface of kerbin
plt.plot(xo,yo,'r-') ###Plot a circle that's the same size as the orbit around kerbin
rcs = ro ###Compute the orbital velocity of the orbit by getting the circular radisu of the orbit
G = 6.67408e-11 ###Gravitational constant (Google)
M = 5.29e22 ##Mass of Kerbin (KSP wiki)
mu = G*M ###Mu is a variable in the Fundamentals of Astrodynamics Textbook (FAT)
VCS = np.sqrt(mu/rcs) ##This is equation 1.8-2 in the FAT textbook
print('VS = ',VCS) ##Print to stdout (standard out) the console
plt.show() ##Show the plot from above