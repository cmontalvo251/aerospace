import numpy as np

MMars = 0.6417e24 #kg
G = 6.67e-11 #units?

muMars = MMars*G

T = 3.802*3600.

a = (T*np.sqrt(muMars)/(2*np.pi))**(2./3.)

RMars = 3389000

rp = 500000 + RMars

e = 1-rp/a

p = rp*(1+e)

nu = np.pi/2.0
vx = np.sqrt(muMars/p)*(-np.sin(nu))
vy = np.sqrt(muMars/p)*(np.cos(nu))
vtheta = np.sqrt(vx**2 + vy**2)

di = 10*np.pi/180.0

delV = 2*vtheta*np.sin(di/2.0)

print(delV)