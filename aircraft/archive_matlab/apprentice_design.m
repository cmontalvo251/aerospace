clear
clc
close all

format bank

%%%Weight Estimate
weight_wing = 0.378;
weight_structure = 0.740;
weight_battery = 0.259;

weight = (weight_wing + weight_structure + weight_battery)*9.81

%%%I get to pick the Thrust to Weight Ratio
TW = 1.2;

%%%Calculate my Thrust
T = TW*weight

T_lbf = T/4.44

%%%From here you can go online and find a motor 
%%%That provides enough thrust. then you can compute the
%%%weight more accurately given the motor/battery weight.

%%%%Wing Loading - W/S
WS_ozsqft = 15;
WS = WS_ozsqft*4.44*(3.28^2)/16

%%%Now we can compute what? Wing Area
S = weight/WS

%%%Choose Aspect Ratio
AR = (58^2)/(58*9)

%%%Can Compute Wingspan now
b = sqrt(AR*S)

%%%Can compute chord now
c = S/b

%%%But first. Need a Flight Speed Regime
V = 5:20;

%%%Also need altitude
h = 0;
rho = 1.225;

%%%Compute Re
mu = 1.81e-5;
%format long g
Re = rho * V * c / mu;

%%%Mach Number
a_inf = 330;
M = V/a_inf;
%%%Pick 0.05

%%%Matrix for XFLR5
[Re',M'];

%%%Pick an Airfoil - NACA 2414
%%%Assume Re is 50,000 to 300,000
alpha0 = -2.1*pi/180;
alpha0_deg = alpha0*180/pi;
a0 = 6.21; %%All from XFLR5

%%%Find the Lift Curve Slope of the wing
e = 0.9;
a = a0/(1+a0/(pi*e*AR));
a_deg = a*pi/180;

%%%Can find CL vs Alpha
alpha_deg = linspace(-15,15,100);
CL = a_deg*(alpha_deg-alpha0_deg);

figure()
plot(alpha_deg,CL)
xlabel('Angle of Attack (deg)')
ylabel('CL')

%%Now we can compute Angle of Attack as function of
%%velocity
AoA = 2*weight./(rho*V.^2*S*a) + alpha0;
AoA_deg = AoA*180/pi;

figure()
plot(V,AoA_deg)
ylabel('Angle of Attack (deg)')
xlabel('Velocity (m/s)')

%%%So for a given flight speed
%%%How much thrust do we need?
%%%Well Let's plot drag
figure()
Cd0 = 0.01; %%%This is the drag coefficient at zero lift
Cdfit = 0.037;
Clfit = 1.05;
plot(0,Cd0,'bx')
hold on
plot(Clfit,Cdfit,'bx')
plot(-Clfit,Cdfit,'bx')

%%%Solve for k
k = (Cdfit - Cd0)/(Clfit^2)
Cl = linspace(-Clfit,Clfit,100);
Cd = Cd0 + k*Cl.^2;
plot(Cl,Cd,'r--')
xlabel('Cl')
ylabel('Cd')

%%%Now let's get the Wing Drag Coefficients
CD0 = Cd0/(1+Cd0/(pi*e*AR))

%%%Now we can compute CD
CD = CD0 + k*CL.^2;

figure()
plot(CL,CD)
xlabel('CL')
ylabel('CD')

%%%%D = 0.5*rho*V^2*S*CD
%%Vary flight speed from 5:20 and we've already done that
%%%Compute AoA which we've done on line 82
%%Using that AoA compute Cl
CLflight = a*(AoA-alpha0);
figure()
plot(V,CLflight)
xlabel('Velocity (m/s)')
ylabel('CL')
%%%Now let's compute Drag Coefficient
CDflight = CD0 + k*CLflight.^2;
figure()
plot(V,CDflight)
xlabel('Velocity (m/s)')
ylabel('CD')
%%%Now Let's compute total Drag
D = 0.5*rho*V.^2.*S.*CDflight;
figure()
plot(V,D)
xlabel('Velocity (m/s)')
ylabel('Drag (N)')

%%%Lift to Drag Ratio
figure()
plot(V,CLflight./CDflight)
xlabel('Velocity (m/s)')
ylabel('L/D')

figure()
plot(AoA_deg,CLflight./CDflight)
xlabel('Angle of Attack (deg)')
ylabel('L/D')

%%%Empennage design
%%%Assume percentage of main wing
bt = 0.5*b
ARt = AR
St = bt^2/ARt
ct = St/bt

%%%%Compute Center of Mass based on weight 
%%%of components
%%%For this example we built the aircraft
%%%and measured it in class
xcg = (4/12)/3.28

%%%%Determine location of Empennage
lt = (18/12)/3.28
%%%%Then compute the aerodynamic center 
%%%of the entire aircraft. Why?
%%%Because I want xacwb to be bigger
%%%Than xcg
xacwb = (St*(c/4+lt)+S*c/4)/(St+S)

%%%%Stability Margin
Sm_inches = (xacwb-xcg)*3.28*12


%%%Let's also compute stall speed
Clmax = 1.5; %%%This would come from XFLR5
CLmax = Clmax/(1+Clmax/(pi*e*AR))

Vstall = sqrt(2*weight/(rho*S*CLmax))
%%%If this is too fast
%%%Add flaps or increase camber
%%%Increase the size of wing