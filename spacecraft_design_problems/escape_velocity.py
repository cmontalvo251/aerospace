#!/usr/bin/python

###Integrate an ordinary differential equation
#in MATLAB that's using the function ode45.
#in Python we're going to use the Scipy toolbox and odeint
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as I
import control as ctl
import scipy.signal as S

##Create a function

def Derivatives(state,t):
    r = state[0]
    v = state[1]
    rdot = v
    G = 6.67e-11 #Gravity Constant
    m = 5.9736e24 #earth mass
    vdot = -G*m/r**2

    return np.asarray([rdot,vdot])

##Main script
plt.close("all")

#integrate for 10 seconds
tout = np.linspace(0,100000,10000000)
r0 = 100000 ##Initial height
G = 6.67e-11 #Gravity Constant
m = 5.9736e24 #earth mass
v0 = np.sqrt(2*G*m/r0)  ##from wikipedia
xinitial = np.asarray([r0,v0])
stateout = I.odeint(Derivatives,xinitial,tout)

plt.plot(tout,stateout[:,0])
plt.grid()

plt.figure()
plt.plot(tout,stateout[:,1])

a = G*m/(stateout[:,0]**2)
v = stateout[:,1]
r = stateout[:,0]
odd = 0.5*np.sqrt(2*G*m)*v*pow(r,-1.5)

plt.figure()
plt.plot(tout,a-odd)

plt.show()

###Integrate 0 to inifinity Gm/r^2
#dt = tout[1]-tout[0]
#vsum = np.sum(a)*dt
# vint = 0
# for val in a:
#     vint = vint + val*dt

# vsumodd = np.sum(odd)*dt

# print r0,v0,vsum,vsumodd
