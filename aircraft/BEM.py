import numpy as np
import matplotlib.pyplot as plt

plt.close("all")

def BEM(R,c,w,angle_of_pitch,v_inf,rho):
    T = 0
    S = R*c
    AR = R**2/S
    ClaIdeal = 2.*np.pi
    Cla = ClaIdeal/(1+ClaIdeal/(np.pi*AR))
    ##Need to integrate from 0 to R*B
    ##Due to tip losses you have B = 0.95
    B = 0.95
    r = np.linspace(0.001,B*R,100)
    dr = r[1]-r[0]
    for ri in r:
        vi = ri*w
        V = np.sqrt((ri*w)**2 + v_inf**2)
        alfa_inf = np.arctan(v_inf/vi)
        alfai = angle_of_pitch - alfa_inf
        Cl = Cla*alfai
        Ti = 0.5*rho*V**2*c*dr*Cl
        T += Ti
    if T < 0:
        T = 0
    return T
    
def BEMQ(R,c,w,angle_of_pitch,v_inf,rho):
    Q = 0
    #S = R*c
    #AR = R**2/S
    #ClaIdeal = 2.*np.pi
    #Cla = ClaIdeal/(1+ClaIdeal/(np.pi*AR))
    #Cda = Cla**2/(np.pi*AR)
    Cda = 1.3
    Cd0 = 0.008
    ##Need to integrate from 0 to R*B
    ##Due to tip losses you have B = 0.95
    B = 0.95
    r = np.linspace(0.001,B*R,100)
    dr = r[1]-r[0]
    for ri in r:
        vi = ri*w
        V = np.sqrt((ri*w)**2 + v_inf**2)
        alfa_inf = np.arctan(v_inf/vi)
        alfai = angle_of_pitch - alfa_inf
        Cd = Cd0 + Cda*alfai
        Qi = 0.5*rho*V**2*c*ri*dr*Cd
        Q += Qi
    if Q < 0:
        Q = 0
    return Q

###In order to compute performance of a prop you need a few things
##Radius of Prop
#R = (29./2.)/12. ##This is a 6 inch prop or really 12 in diameter prop
#R = (14./2.)/12.
R_vec = np.linspace(1,10*12,100)/24.

T_vec = []

for R in R_vec:
    ###Revolutions per Second
    ##To get revolutions per second a reasonable approximation is to take the
    ##voltage of the battery times the KV of the motor
    #KV = 100. #These are the motors I have on my quad.
    #Volt = 22.2 #Again this is the voltage of my battery on my quad
    #N = (Volt*KV)/60.  ###This is revolutions per minute converted to revs per second
    #wIdeal = N*2*np.pi

    KV = 50.
    Volt = 100.
    N = (Volt*KV)/60.
    wIdeal = N*2*np.pi

    ##Need to get w actual
    #Ts = 1.0 #seconds - settling time
    #rmotor = (86.8/1000.)/3.28 ##feet
    #mass_of_motor = 0.267*2.2/32.2 #kg to slugs
    #mstator = 0.3*mass_of_motor ##30% of motor
    #Jmotor = mstator*rmotor**2/2.
    #cdrag = 4*Jmotor/Ts
    #Torque = wIdeal*cdrag
    
    ##Chord of Prop
    #c = 2./12. #inches to feet
    AR = 9.0
    c = R/AR
    #c = 1.0/12.
    
    ##Pitch of props is 4.5 to 5
    #pitch = 9.5/12. ###moves 4 inches for every revolution
    pitch = 9./12.
    angle_of_pitch = np.arctan(pitch/(2*np.pi*R))
    
    ##Density at sea-level
    rho = 0.00238 ##slugs/ft^3
    
    ###Alright so now we need DR for wIdeal
    #wdot = 100
    w = wIdeal*1.0
    #step = 10.
    #while abs(wdot) > 0.01:
    #    Q = BEMQ(R,c,w,angle_of_pitch,0,rho)
    #    wdot_new = Torque-c*w-Q
    #    if np.sign(wdot_new) != np.sign(wdot):
    #        step/=2.
    #    wdot = wdot_new
    #    if wdot > 0:
    #        w+=step
    #    if wdot < 0:
    #        w-=step
    #    print w,wdot
    
    #Inflow Velocity
    #v_inf_vec = np.linspace(0,100,100)
    v_inf = 0
    
    Number_of_Blades = 2.0

    #for v_inf in v_inf_vec:
    T = Number_of_Blades*BEM(R,c,w,angle_of_pitch,v_inf,rho)
    #    T_vec.append(T)
    
    T_vec.append(T)
    
    #Thrust
    #T0 = 1.46*2.2 ##This is from Tiger Motors

    #print(T0,T_vec[0])

#Plot Lift as a function of Vinf
#plt.plot(v_inf_vec,T_vec)
#plt.xlabel('Inflow Velocity (ft/s)')
plt.plot(R_vec,T_vec)
plt.xlabel('Radius (feet)')
plt.ylabel('Thrust (lbf)')
plt.grid()
plt.show()