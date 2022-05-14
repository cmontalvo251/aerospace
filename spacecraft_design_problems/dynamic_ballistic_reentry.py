import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as sci

###################ANALYTIC##############################
##Altitude
he = 122000.0 #entry altitude
ha = np.linspace(0,he,100000) #100 km
R = 6378000 #Radius of Earth in meters
Ve = 7500.0 #m/s
beta = 0.1354/1000.0 ##density constant  (1/km) = (1/m)*(1/1000)
rhos = 1.225 #kg/m^3
m = 1350.0 #kg
CD = 1.5 #unitless - non dimensional coefficient
S = 2.8 #m^2 - planform area of entry vehicle
gammae = -2.0*np.pi/180.0
Va = Ve*np.exp((1.0/(2*beta))*(rhos/(np.sin(gammae)))*(S*CD/m)*np.exp(-beta*ha))
#Compute Acceleration
drdt = Va*np.sin(gammae)
dhdt = drdt
dVdh = (Va[0:-2]-Va[1:-1])/(ha[0:-2]-ha[1:-1])
accela = dVdh*dhdt[0:-2]
amax = -(beta*Ve**2)/(2*np.exp(1))*np.sin(gammae)
print('Max Accel = ',amax,' Max G = ',amax/9.81)
##COmpute Altitude
h_analytic = Va*np.sin(gammae)
##############################################################

##########################NUMERICAL##########################

def Derivatives(state,t):
    ###State vector
    V = state[0]
    h = state[1]
    gamma = state[2]
    
    ###Density model
    rho = rhos*np.exp(-beta*h)
    
    ###Aerodynamic Model
    L = 0.
    D = 0.5*rho*V**2*S*CD
    
    ###Gravity Model
    gs = 9.81
    g = gs*(R/(R+h))**2
    r = R + h
    
    #Dynamics    
    dVdt = -D/m - g*np.sin(gamma)
    dgammadt = (L/m - (g-V**2/r)*np.cos(gamma))/V
    dhdt = V*np.sin(gamma)
    
    ##Return statedot
    statedot = np.asarray([dVdt,dhdt,dgammadt])
    return statedot

####So now we integrate in odeint
stateinitial = np.asarray([Ve,he,gammae])
tout = np.linspace(0,450,10000)
stateout = sci.odeint(Derivatives,stateinitial,tout)

Vnum = stateout[:,0]
hnum = stateout[:,1]
gammanum = stateout[:,2]

accelnum = (Vnum[0:-2]-Vnum[1:-1])/(tout[0:-2]-tout[1:-1])
#############################################################

plt.figure()
plt.plot(ha,Va,label='Analytic')
plt.plot(hnum,Vnum,label='Numerical')
plt.xlabel('Altitude (m)')
plt.ylabel('Velocity (m/s)')
plt.grid()
plt.legend()

plt.figure()
plt.plot(ha[0:-2],accela/9.81,label='Analytic')
plt.plot(hnum[0:-2],accelnum/9.81,label='Numerical')
plt.xlabel('Altitude (m)')
plt.ylabel('Gs')
plt.grid()
plt.legend()

plt.figure()
plt.plot(tout,gammanum,'m-')
plt.ylabel('Flight Path Angle (rad)')
plt.xlabel('Time (sec)')
plt.grid()

plt.figure()
plt.plot(tout,hnum)
plt.xlabel('Time (sec)')
plt.ylabel('Altitude (m)')
plt.grid()

plt.show()
