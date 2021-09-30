import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as I

###Closes all Figures
plt.close("all")

#######THIS IS A FUNCTION CALLED A DEFINITION IN PYTHON
#z and t are input. so z = [x,xdot], z is a 2x1 array with position
# and velocity in it and t is just time.
# matrices are NxM arrays and vectors are Nx1 arrays. Arrays are
# just a way for the computer to handle multiple numbers in one variable
def Derivatives(z,t):
    global mu
    x = z[0]
    y = z[1]
    xdot = z[2]
    ydot = z[3]
    ##This model is correct
    r = np.sqrt(x**2 + y**2)
    ax = -x*mu/r**3
    ay = -y*mu/r**3
    ##Below this is right
    xdbldot = ax
    ydbldot = ay
    zdot = np.asarray([xdot,ydot,xdbldot,ydbldot]) #np is numpy (Numeric Python) and asarray says
    #make [xdot,xdbldot] an array
    return zdot
##############END OF FUNCTION SEPARATED BY TABS#######

##Earth
numplanets = 3
for ctr in np.arange(0,numplanets):
    mu = 1.0
    av = [1.0,9.53667594,5.202887]
    ev = [0.0167,0.05386179,0.04838624]
    a = av[ctr]
    e = ev[ctr]
    p = a*(1-e**2)
    r = p/(1+e)
    rx = r
    ry = 0.0
    vx = 0.0
    vy = np.sqrt(mu/p)*(1+e)
    T = 2*np.pi/np.sqrt(mu)*a**(3./2.)
    tout = np.linspace(0,T,100)  #linspace(start,end,number of data points)
    zinitial = np.asarray([rx,ry,vx,vy])
    zout = I.odeint(Derivatives,zinitial,tout) ##This is the ode toolbox from scipy (Scientific Python)

    xout = zout[:,0]
    yout = zout[:,1]

    plt.plot(xout,yout)
plt.grid()
plt.xlabel('X')
plt.ylabel('Y')
plt.show()
