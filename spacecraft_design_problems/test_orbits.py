import matplotlib.pyplot as plt
import numpy as np

nu = np.linspace(0,2*np.pi,1000)
rk = 1.0
xk = rk*np.cos(nu)
yk = rk*np.sin(nu)
dm = 5
rm = 0.1
xm = rm*np.cos(nu)+dm
ym = rm*np.sin(nu)

rp = rk + 0.3
ra = dm + rm + 0.1
a = (ra+rp)/2.0
e = (ra-rp)/(ra+rp)
p = a*(1-e**2)
r = p/(1-e*np.cos(nu))
xe = r*np.cos(nu)
ye = r*np.sin(nu)

plt.plot(xk,yk)
plt.plot(xm,ym)
plt.plot(xe,ye)
plt.axis('equal')
plt.show()