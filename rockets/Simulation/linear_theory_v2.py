import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as I
import sys

###Closes all Figures
plt.close("all")

#######THIS IS A FUNCTION CALLED A DEFINITION IN PYTHON
#z and t are input. so x = [z,V], x is a 2x1 array with position
# and velocity in it and t is just time.
# matrices are NxM arrays and vectors are Nx1 arrays. Arrays are
# just a way for the computer to handle multiple numbers in one variable

#Initial Mass
T = 0.

def Derivatives(x,t):
    global T
    z = x[0]
    V = x[1]
    m = x[2]

    zdot = -V
    D = 0.0508
    rho = 1.225
    Cx = 0.47
    g = 9.81
    qinf = -np.pi*rho*D**3*Cx/(8*m)
    ##Create our thrust profile
    tburn0 = 3.45
    tburn1_start = 8.6
    mprop = 60./1000.0
    Impulse = 49.61
    if t < tburn0:
        T = Impulse/tburn0
        mdot = -(60./1000.)/tburn0
    elif t > tburn1_start and t < (tburn1_start + tburn0):
        T = Impulse/tburn0
        mdot = -(60./1000.)/tburn0
    else:
        T = 0
        mdot = 0.0
    Vdot = qinf/D*V**2 - g + T/m
    tdrop = 0.1
    if t > tburn0 and t < tburn0 + tdrop:
        mdot = -(113./1000.)/tdrop 
    xdot = np.asarray([zdot,Vdot,mdot]) #np is numpy (Numeric Python) and asarray says
    #make [xdot,xdbldot] an array
    return xdot
##############END OF FUNCTION SEPARATED BY TABS#######

z0 = 0.
V0 = 0.
m0 = ((517.798)/1000.)
tout = np.linspace(0,8.6+3.45+10,100000)
xinitial = np.asarray([z0,V0,m0])
xout = I.odeint(Derivatives,xinitial,tout) 

zout = xout[:,0]
Vout = xout[:,1]
mout = xout[:,2]
Tout = 0*zout
for i in range(0,len(tout)):
    Derivatives([zout[i],Vout[i],mout[i]],tout[i])
    Tout[i] = T

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

plt.figure()
plt.plot(tout,mout)
plt.grid()
plt.xlabel('Time (sec)')
plt.ylabel('Mass (kg)')

plt.figure()
plt.plot(tout,Tout)
plt.xlabel('Time (sec)')
plt.ylabel('Thrust (N)')
plt.grid()

plt.show()