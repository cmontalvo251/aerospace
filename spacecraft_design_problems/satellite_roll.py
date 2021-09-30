from numpy import *
from matplotlib.pyplot import *
from scipy.integrate import *

close("all")

## I * phidbldot = M

#z is the array of states
#z = [phi,phidot]
#zdot = [phidot,phidbldot]

def Derivatives(z,t):
    global M
    #Inertia satellite
    I = 1.0
    I_RW = 0.3
    #States
    phi = z[0]
    phidot = z[1]
    rw_phi = z[2]
    rw_phidot = z[3]
    #Moment on the RW
    #if phi < 20*pi/180.:
    #   M = -1.0
    #if phi > 20*pi/180.0:
    #    M = 1.0
    phic = 0*pi/180.0
    phidotc = 0.0
    err = phic-phi
    errdot = phidotc - phidot
    if rw_phidot > 20:
        M = 0.0
    else:
        M = -1.0*err - 0.1*errdot
    ##Angular acceleration of Satellite
    phidbldot = -M/I
    ##Angular Acceleration of RW
    rw_phidbldot = M/I_RW
    zdot = np.asarray([phidot,phidbldot,rw_phidot,rw_phidbldot])
    return zdot

###MAIN ROUTINE
tout = linspace(0,20,10000) ##Simulate for 10 seconds
##Initial conditions 
zinitial = np.asarray([0,3,0,0])
zout = odeint(Derivatives,zinitial,tout)

#Extract control effort
Mout = 0*tout
for idx in range(0,len(tout)):
    Derivatives(zout[idx,:],tout[idx])
    Mout[idx] = M

plot(tout,zout[:,0],label='Angle Satellite')
plot(tout,zout[:,1],label='Rate Satellite')
grid()
xlabel('Time (sec)')
ylabel('State')
legend()

figure()
plot(tout,zout[:,2],label='Angle RW')
plot(tout,zout[:,3],label='Rate RW')
grid()
xlabel('Time (sec)')
ylabel('RW')
legend()

figure()
plot(tout,Mout)
grid()
xlabel('Time (sec)')
ylabel('Moment Applied to RW (N-m)')

show()