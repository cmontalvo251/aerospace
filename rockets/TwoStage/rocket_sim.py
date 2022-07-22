from control import *
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as I

plt.close("all")

############################################################################
#FUNCTIONS SECTION

def derivatives(state,t):

    x = state[0]
    y = state[1]
    xdot = state[2]
    ydot = state[3]

    r = np.sqrt(x**2+y**2)

    #1ST BURN
    if t <= 1.0:
        yaw = launch_ang*np.pi/180 #launch angle, rad
        xddot = V1*np.cos(yaw)
        yddot = V1*np.sin(yaw)
    #2ND BURN
    elif t >= t_a and t <= t_a+1.0:
        yaw = np.pi/2 + np.arccos(x/r) #launch angle, rad
        xddot = V2*np.cos(yaw)
        yddot = V2*np.sin(yaw)
    #COAST
    else:
        yaw = 0.0 #launch angle, rad
        T = 0.0 #thrust, N 
        xddot = T/m_r*np.cos(yaw) - (G*m_p/r**3)*x
        yddot = T/m_r*np.sin(yaw) - (G*m_p/r**3)*y 

    statedot = np.asarray([xdot,ydot,xddot,yddot])
    
    return statedot

#END FUNCTIONS SECTION
############################################################################



############################################################################
#VARIABLES SECTION

#PARAMETERS
m_r = 200.0 #mass of rocket, kg
r_p = 3396200 #radius of planet, m
m_p = 0.64171e24 #mass of planet, kg
G = 6.67408e-11 #gravitational constant, m^3/(kg*s^2)
mu = G*m_p #gravitational parameter

#EXAMPLE 1
V1 = 3625.5 #delta v1, m/s
V2 = 72.85 #delta v2, m/s
launch_ang = 90.0 #launch angle, deg
t_a = 3194.52 #apoapsis time, sec

#EXAMPLE 2
#V1 = 2000.0 #delta v1, m/s
#V2 = 2000.0 #delta v2, m/s
#launch_ang = 49.56 #launch angle, deg
#t_a = 474.74 #apoapsis time, sec

#EXAMPLE 3
#V1 = 1440.0 #delta v1, m/s
#V2 = 3400.0 #delta v2, m/s
#launch_ang = 0.001 #launch angle, deg
#t_a = 436.15 #apoapsis time, sec

#INITIAL CONDITIONS
x = r_p #x position, m
y = 0.0 #y position, m
xdot = 0.0 #x velocity, m/s
ydot = 0.0 #y velocity, m/s
state_init = np.asarray([x,y,xdot,ydot])

#TIME SPAN TO INTEGRATE OVER
t_init = 0.0 #start time
t_final = 10000.0 #end time
dt = 0.01 #time step
tspan = np.linspace(t_init,t_final,int(t_final/dt+1))

#END VARIABLES SECTION
############################################################################



############################################################################
#INTEGRATION SECTION

state_out = I.odeint(derivatives,state_init,tspan,hmax = dt)
x_out = state_out[:,0]
y_out = state_out[:,1]
xdot_out = state_out[:,2]
ydot_out = state_out[:,3]

#END INTEGRATION SECTION
############################################################################



############################################################################
#PLOTTING SECTION

print((np.max(np.sqrt(x_out**2+y_out**2))-r_p)/1000)
print(tspan[np.argmax(np.sqrt(x_out**2+y_out**2))])

#XY
plt.figure()
planet = plt.Circle((0,0), radius = r_p, fc='r')
plt.gca().add_patch(planet)
plt.plot(y_out,x_out)
plt.xlabel('Y (km)')
plt.ylabel('X (km)')
plt.gca().set_aspect("equal")
plt.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
plt.ticklabel_format(axis='x', style='sci', scilimits=(-2,2))
plt.grid()
"""
#X
plt.figure()
plt.plot(tspan,x_out)
plt.xlabel('Time (sec)')
plt.ylabel('X (m)')
plt.grid()

#Y
plt.figure()
plt.plot(tspan,y_out)
plt.xlabel('Time (sec)')
plt.ylabel('Y (m)')
plt.grid()

#XDOT
plt.figure()
plt.plot(tspan,xdot_out)
plt.xlabel('Time (sec)')
plt.ylabel('X Velocity (m/s)')
plt.grid()

#YDOT
plt.figure()
plt.plot(tspan,ydot_out)
plt.xlabel('Time (sec)')
plt.ylabel('Y Velocity (m/s)')
plt.grid()
"""
plt.show()

#END PLOTTING SECTION
############################################################################