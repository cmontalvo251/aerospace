import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as I


def Derivatives(state,t):
    ### x - position along x axis - state[0]
    ### y - position along y axis - state[1]
    x = state[0]
    y = state[1]
    ### xdot - velocity along x axis - state[2]
    ### ydot - veloicty along y axis - state[3]
    mass = state[4]
    
    ###Kinematics
    xdot = state[2]
    ydot = state[3]
    
    ###Force Model
    ## Gravity Model
    G = 6.6742*10**-11; #%%Gravitational constant
    Mkerbin = 5.2915158*10**22 #
    muKerbin = G*Mkerbin
    r = np.sqrt(x**2 + y**2)
    gravx = -x*muKerbin*mass/r**3
    gravy = -y*muKerbin*mass/r**3
    
    ### Thrust Model
    if t < 51:    
        thrustx = 0*168410 ##Newtons
    else:
        thrustx = 0.0        
    thrusty = 0.0
    
    thrust = np.sqrt(thrustx**2 + thrusty**2)
    
    ## Aerodynamic Model
    dragx = 0.0
    dragy = 0.0
    
    ##Total Forces
    forcex = gravx + thrustx + dragx
    forcey = gravy + thrusty + dragy
    ##Dynamics
    ## Force = mass * accel -> accel = Force /mass
    xdbldot = forcex/mass
    ydbldot = forcey/mass
    
    ##Model Mass Flow Rate
    Isp = 251.0 #seconds
    surface_gravity = 9.81 #m/s^2
    exit_velocity = Isp*surface_gravity
    mdot = -thrust/exit_velocity
    
    ### statedot = [xdot,ydot,xdbldot,ydbldot,mdot]
    
    return np.asarray([xdot,ydot,xdbldot,ydbldot,mdot])

###PROPERTIES OF KERBIN
Rkerbin = 600000 #meters
sidereal_period = 21549.425
sidereal_angular_velocity = 2*np.pi/sidereal_period
sidereal_rotational_velocity = sidereal_angular_velocity*Rkerbin

tout = np.linspace(0,100*60.0,1000)
x0 = Rkerbin + 80000
y0 = 0.0
xdot0 = 0.0
ydot0 = sidereal_rotational_velocity + 2300
mass0 = (12.3*1000) #kilograms (assume a metric ton)
state_initial = np.asarray([x0,y0,xdot0,ydot0,mass0])
stateout = I.odeint(Derivatives,state_initial,tout) ##This is the ode toolbox from scipy (Scientific Python)

xout = stateout[:,0]
yout = stateout[:,1]
vxout = stateout[:,2]
vyout = stateout[:,3]
mout = stateout[:,4]

vout = np.sqrt(vxout**2+vyout**2)

plt.figure()
plt.plot(tout,vout)
plt.grid()
plt.xlabel('Time (sec)')
plt.ylabel('Velocity (m/s)')

plt.figure()
plt.plot(xout,yout)
thetakerbal = np.linspace(0,2*np.pi,1000)
xkerbal = Rkerbin*np.cos(thetakerbal)
ykerbal = Rkerbin*np.sin(thetakerbal)
plt.plot(xkerbal,ykerbal)
plt.grid()
plt.xlabel('X (m)')
plt.ylabel('Y (m)')

plt.figure()
plt.plot(tout,np.sqrt(xout**2+yout**2)-Rkerbin)
plt.grid()
plt.xlabel('Time (sec)')
plt.ylabel('AGL (m)')

plt.figure()
plt.plot(tout,mout)
plt.grid()
plt.xlabel('Time (sec)')
plt.ylabel('Mass (kg)')

plt.show()



