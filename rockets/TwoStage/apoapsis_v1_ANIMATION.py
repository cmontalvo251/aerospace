#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt


v1 = np.linspace(1,4000,1000)

sig = 0.0001
rp  = 3396200.0
mp = 6.4171e23
G = 6.67408e-11
mu  = G*mp

a = (rp/mu - rp*(np.cos(sig)**2)/mu)
b = (-rp*np.cos(sig)*np.sin(sig)/mu)

e = np.sqrt((v1**4)*a**2 - 2.0*(v1**2)*a + 1.0 + (v1**4)*b**2)

# print("ex = ",(v1**2)*a-1)
# print("ey = ",(v1**2)*b)
# print("e = ",e)

ra = (rp**2)*(v1**2)*(np.sin(sig)**2)/(mu*(1-e))

plt.plot(v1,(ra-rp)/1000.0)

plt.xlabel('v1 (m/s)')
plt.ylabel('ra-rp (km)')

plt.show()

#r = 
#v = 
