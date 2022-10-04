import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as I

###Closes all Figures
plt.close("all")

#######THIS IS A FUNCTION CALLED A DEFINITION IN PYTHON
#z and t are input. so x = [z,V], x is a 2x1 array with position
# and velocity in it and t is just time.
# matrices are NxM arrays and vectors are Nx1 arrays. Arrays are
# just a way for the computer to handle multiple numbers in one variable
m = 4
def Derivatives(x,t):
    global m
    z = x[0]
    V = x[1]
    zdot = -V
    D = 0.5
    rho = 1.225
    Cx = 0.47
    g = 9.81
    qinf = -np.pi*rho*D**3*Cx/(8*m)
    Vdot = qinf/D*V**2 - g
    xdot = np.asarray([zdot,Vdot]) #np is numpy (Numeric Python) and asarray says
    #make [xdot,xdbldot] an array
    return xdot
##############END OF FUNCTION SEPARATED BY TABS#######
    
tout = np.linspace(0,5,100)  #linspace(start,end,number of data points)
t = 3.
Impulse = 100.
T = Impulse/t
a = T/m
z0 = -0.5*a*t**2
V0 = a*t
xinitial = np.asarray([z0,V0])
xout = I.odeint(Derivatives,xinitial,tout) ##This is the ode toolbox from scipy (Scientific Python)

zout = xout[:,0]
Vout = xout[:,1]
plt.figure()
plt.plot(tout,-zout)
plt.grid()
plt.xlabel('Time (sec)')
plt.ylabel('Altitude (m)')

plt.figure()
plt.plot(tout,Vout)
plt.grid()
plt.xlabel('Time (sec)')
plt.ylabel('Velocity (m/s)')

plt.show()
