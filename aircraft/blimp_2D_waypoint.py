import numpy as np
import matplotlib.pyplot as plt

##DERIVATIVES ROUTINE
def Derivatives(t,state,delta_steer,delta_throttle):
    ##EXTRACT OUR STATE VECTOR
    x = state[0]
    y = state[1]
    psi = state[2]
    xdot = state[3]
    ydot = state[4]
    psidot = state[5]
    
    ##Mass Properties
    mass = 1.0 #kg
    Inertia = 1.0 #kg-m^2
    
    ##Friction or Drag terms
    CR = 1.0
    CD = 1.0
    CV = 25.0
    
    ##CONTROL TERMS
    DS = 1.0/500.0
    DT = 0.025
    
    ##KINEMATICS
    uwinds = 0.0
    vwinds = 2.0
    u = xdot * np.cos(psi) + ydot * np.sin(psi) + uwinds
    v = -xdot * np.sin(psi) + ydot * np.cos(psi) + vwinds
    
    ##Control Forces
    Steer_Moment = DS * (delta_steer - 1500)
    Throttle_Force = DT * (delta_throttle - 1000)
    
    ##Total Forces and Moments
    Moments = Steer_Moment - CR*psidot
    ForceXbody = Throttle_Force - CD*u
    ForceYbody = -CV*v
    ForceX = ForceXbody * np.cos(psi) - ForceYbody * np.sin(psi)
    ForceY = ForceXbody * np.sin(psi) + ForceYbody * np.cos(psi)
    
    ##Dynamics
    ##Translation Dynamics
    xddot = ForceX / mass
    yddot = ForceY / mass
    
    ##Rotational Dynamics
    psiddot = Moments / Inertia
    
    ##RETURN THE DERIVATIVE OF OUR STATE VECTOR
    dstatedt = np.array([xdot,ydot,psidot,xddot,yddot,psiddot])
    return dstatedt

##Estimation routine
def Estimate(state,stateE):
    xest = stateE[0]
    yest = stateE[1]
    psi = state[2]
    xdot = state[3]
    ydot = state[4]
    u = xdot * np.cos(psi) + ydot * np.sin(psi)
    v = -xdot * np.sin(psi) + ydot * np.cos(psi)
    psiest = psi
    uest = 10
    vest = 0
    xdotest = uest * np.cos(psiest) - vest * np.sin(psiest)
    ydotest = uest * np.sin(psiest) + vest * np.cos(psiest)
    xest = xest + uest*np.cos(psiest)*dt
    yest = yest + uest*np.sin(psiest)*dt
    stateE = np.array([xest,yest,psiest,xdotest,ydotest])
    return stateE
    
###CONTROL ROUTINE
def Control(t,state,stateE):
    global error_int,XWPs,YWPs,WPctr, heading_bias
    #Extract State Vector
    x = state[0]
    y = state[1]
    heading_noise = 0.0*np.random.normal()
    psi = state[2] + 0*heading_bias + 0*heading_noise
    xdot = state[3]
    ydot = state[4]
    psidot = state[5]
    
    #Extract estimated states
    xest = stateE[0]
    yest = stateE[1]
    psiest = stateE[2]
    xdotest = stateE[3]
    ydotest = stateE[4]
    
    ##Kinematics
    u = xdot * np.cos(psi) + ydot * np.sin(psi) ##Assuming that you are getting GPS data to obtain velocity.
    uest = xdotest*np.cos(psiest) + ydotest*np.sin(psiest)
    
    ##Create a waypoint
    if WPctr >= len(XWPs):
        return 1500,1000
        #WPctr = 0
    XWP = XWPs[WPctr]
    YWP = YWPs[WPctr]
    dx = XWP - xest
    dy = YWP - yest
    norm = np.sqrt(dx**2 + dy**2)
    if norm < 1:
        WPctr +=1
        print('Arrive at WP: ',XWP,YWP,' T = ',t, 'WPctr = ',WPctr)

    psic = np.arctan2(dy,dx)
    #print(psic)
    
    ##Compute Heading Error using Non Wrapping Function
    #//%%%returns delta psi from a heading and a heading command without
    #//%worrying about wrapping issues
    #//%%%This computes delpsi = psi-psic
    #psic = 45*np.pi/180.
    spsi = np.sin(psiest);
    cpsi = np.cos(psiest);
    spsic = np.sin(psic);
    cpsic = np.cos(psic);
    delpsi = np.arctan2(spsi*cpsic-cpsi*spsic,cpsi*cpsic+spsi*spsic)
    
    ##Heading Angle Controller
    kpp = -1000/np.pi
    delta_steer = 1500.0 + kpp*delpsi
    
    #Throttle Controller
    ucommand = 10.0
    error_signal = uest - ucommand
    error_int += error_signal
    #print(error_int)
    kp = -120.0
    ki = -2.0
    delta_throttle = 1000 + kp*error_signal + ki*error_int
    
    ##Saturation Controller
    if delta_throttle < 1000:
        delta_throttle = 1000
    if delta_throttle > 2000:
        delta_throttle = 2000
        
    if delta_steer > 2000:
        delta_steer = 2000
    if delta_steer < 1000:
        delta_steer = 1000
    
    #print(delta_throttle)
    #sprint(delta_steer)
    
    return delta_steer,delta_throttle
    
##RUN SIMULATION
dt = 0.1
tout = np.arange(0,2000,dt)
xout = 0*tout
yout = 0*tout
xdotout = 0*tout
ydotout = 0*tout
psiout = 0*tout
psidotout = 0*tout
delta_throttle_out = 0*tout
delta_steer_out = 0*tout
ctr = 0
#Initial Conditions
x = 0.
xest = x
y = 0.
yest = y
xdot = 10.
xdotest = xdot
ydot = 0.
ydotest = ydot
psi = 0*45*np.pi/180.
psiest = psi
psidot = 0.
error_int = -200.0
state0 = np.array([x,y,psi,xdot,ydot,psidot])
stateE = np.array([xest,yest,psiest,xdotest,ydotest])
XWPs = np.array([0,10000])
YWPs = np.array([0,0])
WPctr = 0

##Initialize Bias
heading_bias = 0.0*np.random.normal()
print('Heading Bias = ',heading_bias)

for t in tout:
    ##Save States
    xout[ctr] = state0[0]
    yout[ctr] = state0[1]
    psiout[ctr] = state0[2]
    xdotout[ctr] = state0[3]
    ydotout[ctr] = state0[4]
    psidotout[ctr] = state0[5]
    #Run Control loop once per timestep
    delta_steer,delta_throttle = Control(t,state0,stateE)
    delta_steer_out[ctr] = delta_steer
    delta_throttle_out[ctr] = delta_throttle
    #Compute Estimated position (once per timestep)
    stateE = Estimate(state0,stateE)
    #RK4
    k1 = Derivatives(t,state0,delta_steer,delta_throttle)
    k2 = Derivatives(t+dt/2,state0+k1*dt/2,delta_steer,delta_throttle)
    k3 = Derivatives(t+dt/2,state0+k2*dt/2,delta_steer,delta_throttle)
    k4 = Derivatives(t+dt,state0+k3*dt,delta_steer,delta_throttle)
    phi = (1./6.)*(k1 + 2*k2 + 2*k3 + k4)
    state0 += phi*dt
    ctr+=1
    #print('T = ',t)

##PLOTS
plt.figure()
plt.plot(xout,yout)
plt.grid()
plt.xlabel('x')
plt.ylabel('y')

##Body Frame Velocity
uout = xdotout * np.cos(psiout) + ydotout * np.sin(psiout)
vout = -xdotout * np.sin(psiout) + ydotout * np.cos(psiout)
plt.figure()
plt.plot(tout,uout,label='U')
plt.plot(tout,vout,label='V')
plt.grid()
plt.legend()
plt.xlabel('t')
plt.ylabel('velocity')

plt.figure()
plt.plot(tout,psidotout*180/np.pi)
plt.grid()
plt.xlabel('t')
plt.ylabel('Yaw Rate (deg/s)')

plt.figure()
plt.plot(tout,psiout*180/np.pi)
plt.grid()
plt.xlabel('t')
plt.ylabel('Euler Yaw Angle (deg)')

plt.figure()
plt.plot(tout,delta_steer_out)
plt.grid()
plt.xlabel('t')
plt.ylabel('Delta Steer (us)')

plt.show()