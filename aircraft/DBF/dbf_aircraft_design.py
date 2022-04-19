# -*- coding: utf-8 -*-
###DBF Aircraft Design 2018-2019
import numpy as np
import matplotlib.pyplot as plt

plt.close("all")

def BEM(R,c,w,pitch,v_inf,rho):
    T = 0
    Cla = 2.*np.pi
    ##Need to integrate from 0 to R
    r = np.linspace(0.001,R,100)
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

###Wing Cube Loading (WCL= W/S^(3/2))
#Slow flyers and thermal gliders – under 4
#Trainers, park flyers, 3D – 5 to 7
#General sport and scale aerobatics – 7 to 10
#Sport and scale models – 10 to 13
#Warbirds and racers – 13 and over

#The chord length 
c_vec = np.linspace(0.2,1.5,10) #feet
#c = 0.7 #based on iteration above

#Wingspan is minimum 4 feet but we can change this
#b_vec = np.linspace(4,20,10)
#b = 4.0 #feet
b = 7.0 #Based on iterations above

##Variation on attack stores
#Na_vec = np.asarray([0,1,2,3,4])
Na = 0.

M_vec = []
ctr = 0
for c in c_vec:
#for b in b_vec:
    ctr+=1
    print(ctr)
    ##We also need to get airfoil parameters
    ###Assume main wing and tail airfoils are the same for 
    ###Simplicity
    Cl0 = 0.0 
    Cd0 = 0.008
    Cda = 1.35
    Cla = 5.8
    
    ###Need fuselage drag coefficients too
    ##If it looks good it flys good. So let's assume the diamter 
    #of the fuselage is 25 percent of the wingspan
    d = 0.10*b
    Cdf0 = 0.47 ##ANF
    Cdf2 = 3.32 ##These are taken from an Army Navy Finner (ANF)
    
    ##Area (ft^2)
    S = b*c
    
    #Need to know tail area as well
    st = 0.33 #I'm pretty sure most aircraft tails are a third the weight of main wing
    St = st*S
    
    ##Aspect Ra    b = 4.0 #feettio
    AR = b**2/S
    ARt = AR #Assume aspect ratio of tail is the same as main wing.
    
    #Compute wingspan and chord of tail just for kicks
    bt = np.sqrt(St*ARt)
    ct = St/bt
    
    ##Compute 3D lift
    CL0w = Cl0 / (1 + Cl0/(np.pi*AR))
    CD0w = Cd0 / (1 + Cd0/(np.pi*AR))
    CDAw = Cda / (1 + Cda/(np.pi*AR))
    CLAw = Cla / (1 + Cla/(np.pi*AR))
    CL0t = Cl0 / (1 + Cl0/(np.pi*ARt))
    CD0t = Cd0 / (1 + Cd0/(np.pi*ARt))
    CDAt = Cda / (1 + Cda/(np.pi*ARt))
    CLAt = Cla / (1 + Cla/(np.pi*ARt))
    
    ###Combine Coefficients to get total coeffs
    CL0 = CL0w + CL0t*St/S
    CD0 = CD0w + CD0t*St/S
    CDA = CDAw + CDAt*St/S
    CLA = CLAw + CDAt*St/S
    
    ##Compute Weight which is based on wing area and attack stores
    rhob = 0.244 ###this is density of balsa - need lbf/ft^2
    Ww = rhob*S ##weight of wing
    ss = 1.0 ##not sure on this one either. This is a weight factor due to wing
    Wb = 1.1*2.2 ##weight of 1 battery pack
    Wm = 0.2*2.2 ##one motor
    Ws = ss*Ww + Wb + Wm  ##weight of structure is a percentage of wing weight
    Wt = st*Ww ##weight of tail is percentage of wing weight
    Wa = (0.7*2.2) ##weight of attack store
    #Na = 0.0 #Number of attack stores
    Wp = Na*Wa
    W = Ws + Ww + Wt + Wp
    
    #Compute alfa as a function of flight speed
    #Assume Lift = Weight to get angle of attack
    V = np.linspace(10,100,1000)
    rho = 0.00238 ##slugs/ft^3
    alfa = (2*W/(rho*S*V**2)-CL0)/CLA
    #plt.figure()
    #plt.plot(V,alfa*180./np.pi)
    #plt.xlabel('V (ft/s)')
    #plt.ylabel('AoA (deg)')
    #plt.grid()
    
    ##Make sure Lift equals weight
    #L = 0.5*rho*V**2*S*(CL0 + CLA*alfa)
    #plt.figure()
    #plt.plot(V,L)
    
    ##Compute Thrust Required - Assume Thrust = Drag and Lift = Weight
    ##This includes drag from the fuselage as well
    #The easiest way to approximate fuselage drag is just to assume it 
    #is a rocket is diameter d then Fuselage Drag is
    Df = (np.pi*rho*(V**2)*(d**2)/8.0)*(Cdf0 + Cdf2*alfa**2)
    D = 0.5*rho*V**2*S*(CD0 + CDA*alfa**2)
    #plt.figure()
    #plt.plot(V,D,label='Thrust Required')
    #plt.xlabel('V (ft/s)')
    #plt.ylabel('T (lbf)')
    #plt.grid()
    
    ###In order to compute performance of a prop you need a few things
    R = (13./2.)/12. ##This is a 6 inch prop or really 12 in diameter prop
    cprop = 1.5/12.
    ##To get revolutions per second a reasonable approximation is to take the
    ##voltage of the battery times the KV of the motor
    KV = 600 #These are the motors I have on my quad.
    Volt = 14.8 #Again this is the voltage of my battery on my quad
    N = (Volt*KV)/60.  ###This is revolutions per minute converted to revs per second
    J = V/(N*2.*R) ###Advance ratio is just a linear thing with V since N and R are constants
    w = N*2*np.pi
    
    ##Pitch of props is 4.5 to 5
    pitch = 4./12. ##4 inches per revolution converted to feet
    angle_of_pitch = np.arctan(pitch/(2*np.pi*R))
    
    ##In order to get efficiency (npr) we need a few more things
    ##npr is a function of advance ratio
    #Need some blade element momentum theory here to get this
    ##P*npr = T*V where P is power available
    ##When V = 0, T is not zero (static thrust) but npr = 0 even though P is not zero
    #So 0 = 0
    #What I want to know though is what is Thrust as a function of V?
    #Alright according to Anderson - Power is constant with speed which makes sense
    #That means that P*npr/V = T
    #P is constant, V changes linearly and npr is a function of efficiency
    #Problem of course is that npr = f(V) so that's where we use BEM to get it
    #This also highlights another issue. If V=0 then P/0 = inf but that's obviously not true
    #What happens is that npr is also zero. Therefore we actually have
    # P * 0 / 0 = T and since 0/0 is undefined T can be whatever it wants
    # So we either need T = f(V) or npr = f(V). Apparently npr is obtained
    # Experimentally but like I said BEM should work.
    
    #How much thrust do you have?
    Tmax = []
    for v_inf in V:
        T = BEM(R,cprop,w,pitch,v_inf,rho)
        Tmax.append(T)
    #plt.plot(V,Tmax,label='Thrust Available')
    #plt.legend()
    
    ##So another idea is to take T = D and multiply both sides by V
    # then we have T*V = D*V = P
    # so if we multiply D*V we get power required rather than thrust required
    #P = D*V
    #plt.figure()
    #plt.plot(V,P)
    
    ##Solve for Max Speed
    #plt.figure()
    #plt.plot(V,abs(Tmax-T))
    Vs = V[np.where(abs(Tmax-D)<0.01)]
    if len(Vs) == 0:
        Vmax = 0
        print('Not enough thrust')
    else:
        Vmax = Vs[-1]
    print('Vmax=',Vmax)
    
    ###Now Compute Mission Score
    distance_per_lap = 2500. #feet
    if Vmax != 0:
        Our_Time = 3*distance_per_lap/Vmax
        number_of_laps = Na #Number of attack stores
    else:
        Our_Time = 600. #Max is 10 minutes
        number_of_laps = -2.
        
    M1 = 1.0
    Min_Time = 120. #2 minutes seems unreasonabl for min time but whatever
    M2 = 1.0 + Min_Time/Our_Time
    M3 = 2.0 + number_of_laps
    GND_M = 1.0
    
    ## Final Score
    M = M1 + M2 + M3 + GND_M
    
    print('Final Mission Score =',M)
    
    M_vec.append(M)


plt.figure()
plt.plot(c_vec,M_vec-np.min(M_vec))
#plt.plot(main_vec,M_vec-np.min(M_vec),label=variable[var_ctr])
plt.xlabel('Chord Length (ft)')
#plt.xlabel('Wingspan (ft)')
#plt.xlabel('% Change in Parameter')
#plt.xlabel('Number of Attack Stores')
plt.ylabel('Change in Mission Score')
plt.grid()
plt.show()
