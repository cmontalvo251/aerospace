import numpy as np
import sys
import matplotlib.pyplot as plt

plt.close("all")

#%%%Weight Estimate
weight_wing = 0.378;
weight_structure = 0.740;
weight_battery = 0.259;

weight = (weight_wing + weight_structure + weight_battery)*9.81

#%%%I get to pick the Thrust to Weight Ratio
TW = 1.2;

#%%%Calculate my Thrust
T = TW*weight

T_lbf = T/4.44

#%%%From here you can go online and find a motor 
#%%%That provides enough thrust. then you can compute the
#%%%weight more accurately given the motor/battery weight.

#%%%%Wing Loading - W/S
WS_ozsqft = 15;
WS = WS_ozsqft*4.44*(3.28**2)/16

#%%%Now we can compute what? Wing Area
S = weight/WS

#%%%Choose Aspect Ratio
AR = (58**2)/(5.08*9.0)

#%%%Can Compute Wingspan now
b = np.sqrt(AR*S)

#%%%Can compute chord now
c = S/b

#%%%But first. Need a Flight Speed Regime
V = np.arange(5,20,1);

#%%%Also need altitude
h = 0;
rho = 1.225;

#%%%Compute Re
mu = 1.81e-5;
#%format long g
Re = rho * V * c / mu;

#%%%Mach Number
a_inf = 330;
M = V/a_inf;
#%%%Pick 0.05

#%%%Matrix for XFLR5
[Re,M];

#%%%Pick an Airfoil - NACA 2414
#%%%Assume Re is 50,000 to 300,000
alpha0 = -2.1*np.pi/180.0;
alpha0_deg = alpha0*180.0/np.pi;
a0 = 6.21; #%%All from XFLR5

#%%%Find the Lift Curve Slope of the wing
e = 0.9;
a = a0/(1+a0/(np.pi*e*AR));
a_deg = a*np.pi/180;

#%%%Can find CL vs Alpha
alpha_deg = np.linspace(-15,15,100);
CL = a_deg*(alpha_deg-alpha0_deg);

plt.figure()
plt.plot(alpha_deg,CL)
plt.xlabel('Angle of Attack (deg)')
plt.ylabel('CL')
plt.grid()

#%%Now we can compute Angle of Attack as function of
#%%velocity
AoA = 2*weight/(rho*(V**2)*S*a) + alpha0;
AoA_deg = AoA*180/np.pi;

plt.figure()
plt.plot(V,AoA_deg)
plt.ylabel('Angle of Attack (deg)')
plt.xlabel('Velocity (m/s)')

#%%%So for a given flight speed
#%%%How much thrust do we need?
#%%%Well Let's plot drag
plt.figure()
Cd0 = 0.01; #%%%This is the drag coefficient at zero lift
Cdfit = 0.037;
Clfit = 1.05;
plt.plot(0,Cd0,'bx')
plt.plot(Clfit,Cdfit,'bx')
plt.plot(-Clfit,Cdfit,'bx')

#%%%Solve for k
k = (Cdfit - Cd0)/(Clfit**2)
Cl = np.linspace(-Clfit,Clfit,100);
Cd = Cd0 + k*Cl**2;
plt.plot(Cl,Cd,'r--')
plt.xlabel('Cl')
plt.ylabel('Cd')

#%%%Now let's get the Wing Drag Coefficients
CD0 = Cd0/(1+Cd0/(np.pi*e*AR))

#%%%Now we can compute CD
CD = CD0 + k*CL**2;

plt.figure()
plt.plot(CL,CD)
plt.xlabel('CL')
plt.ylabel('CD')

#%%%%D = 0.5*rho*V^2*S*CD
#%%Vary flight speed from 5:20 and we've already done that
#%%%Compute AoA which we've done on line 82
#%%Using that AoA compute Cl
CLflight = a*(AoA-alpha0);
plt.figure()
plt.plot(V,CLflight)
plt.xlabel('Velocity (m/s)')
plt.ylabel('CL')
#%%%Now let's compute Drag Coefficient
CDflight = CD0 + k*CLflight**2;
plt.figure()
plt.plot(V,CDflight)
plt.xlabel('Velocity (m/s)')
plt.ylabel('CD')
#%%%Now Let's compute total Drag
D = 0.5*rho*V**2*S*CDflight;
plt.figure()
plt.plot(V,D)
plt.xlabel('Velocity (m/s)')
plt.ylabel('Drag (N)')

#%%%Lift to Drag Ratio
plt.figure()
plt.plot(V,CLflight/CDflight)
plt.xlabel('Velocity (m/s)')
plt.ylabel('L/D')

plt.figure()
plt.plot(AoA_deg,CLflight/CDflight)
plt.xlabel('Angle of Attack (deg)')
plt.ylabel('L/D')

#%%%Empennage design
#%%%Assume percentage of main wing
bt = 0.5*b
ARt = AR
St = bt**2/ARt
ct = St/bt

#%%%%Compute Center of Mass based on weight 
#%%%of components
#%%%For this example we built the aircraft
#%%%and measured it in class
xcg = (4.0/12.0)/3.28

#%%%%Determine location of Empennage
lt = (18.0/12.0)/3.28
#%%%%Then compute the aerodynamic center 
#%%%of the entire aircraft. Why?
#%%%Because I want xacwb to be bigger
#%%%Than xcg
xacwb = (St*(c/4.0+lt)+S*c/4.0)/(St+S)

#%%%%Stability Margin
Sm_inches = (xacwb-xcg)*3.28*12.0


#%%%Let's also compute stall speed
Clmax = 1.5; #%%%This would come from XFLR5
CLmax = Clmax/(1+Clmax/(np.pi*e*AR))

Vstall = np.sqrt(2*weight/(rho*S*CLmax))
#%%%If this is too fast
#%%%Add flaps or increase camber
#%%%Increase the size of wing

plt.show()