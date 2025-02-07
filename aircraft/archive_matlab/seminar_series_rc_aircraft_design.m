clear
clc
close all

format bank

%%%%Wing Cube Loading - I get to decide this
WCL = 9;

%%%Weight Estimate (oz)
weight_servo = 3;
x_servo = 2; %%%in reference to the LE of the main wing
weight_battery = 5.0; %%
x_battery = -3;
weight_esc = 4; %%%  Power supply 
x_esc = -2;
weight_motor = 4; %% and propulsion
x_motor = -4;
weight_prop = 3; %% Thrust
x_prop = -4;
weight_receiver = 3;
x_receiver = 0;
weight_control_horns_rods = 2;
x_control_horns_rods = 0;
weight_fuselage = 16;
x_fuselage = 8;
weight = (weight_servo+weight_battery+weight_esc+weight_motor+weight_prop+weight_receiver+weight_control_horns_rods+weight_fuselage);
weight_lbf = weight/16;


disp(['Weight = ',num2str(weight),' oz'])

%%%One Design point is to pick a motor and figure out it's thrust
Thrust = 33; %%%Thrust tables online
TW = Thrust/weight

%%%This is another design point where you select the TW
%TW = 0.8;
%Thrust = TW*weight

%%%From here you can go online and find a motor
%%%That provides enough thrust. then you can compute the
%%%weight more accurately given the motor/battery weight.

%%%Now we can compute
%%% WCL = W/S^(3/2)
S = (weight/WCL)^(2/3)

%%%Choose Aspect Ratio
AR = 4.5

%%%Can Compute Wingspan now
b = sqrt(AR*S)

%%%Can compute chord now
c = S/b
c_SI = c/3.28;

%%%But first. Need a Flight Speed Regime
V_SI = linspace(5,20,100); %%%m/s
V = V_SI*3.28; %%%ft/s

%%%Density at altitude
rho_SI = 1.225; %%%kg/m^3
rho = 0.002378; %%%slugs/ft^3

%%%Compute Re
mu = 1.81e-5; %%%SI units
%format long g
Re = rho_SI * V_SI * c_SI / mu;

%%%Mach Number
a_inf = 330;
M = V_SI/a_inf;

figure()
plot(V_SI,M)
figure()
plot(V_SI,Re)

%%%MACH and RE NUMBER FOR XFLR5
%%% M = 0.05
%%% Re = 200000

%%%Pick an Airfoil
%%%Assume Re is 50,000 to 300,000

%%%PICK OUT ALL OF THESE VALUES FROM XFLR5
alpha0 = -5.0*pi/180; %%%ZERO LIFT AOA
alpha0_deg = alpha0*180/pi
Cla = 0.6/(-alpha0) %%%Lift curve slope

%%%Find the Lift Curve Slope of the wing
e = 0.9; %%%Efficiency
CLA = Cla/(1+Cla/(pi*e*AR))
CLA_deg = CLA*pi/180;

%%%Can find CL vs Alpha
alpha_deg = linspace(-15,15,1000);
alpha = alpha_deg*pi/180;
CL = CLA_deg*(alpha_deg-alpha0_deg);

figure()
plot(alpha_deg,CL)
xlabel('Angle of Attack (deg)')
ylabel('CL')

%%Now we can compute Angle of Attack as function of
%%velocity
%%% L = W = 0.5*rho*V^2*S*CL
%%% CL = 2*W/(rho*V^2*S)
CL_req = 2*weight_lbf./(rho*(V.^2)*S);
figure()
plot(V,CL_req)
xlabel('Velocity (ft/s)')
ylabel('CL (L=W)')

AoA_req_deg = CL_req/CLA_deg + alpha0_deg;
AoA_req = AoA_req_deg*pi/180;
figure()
plot(V,AoA_req_deg)
xlabel('Velocity (ft/s)')
ylabel('AoA (deg)')

%%%So for a given flight speed
%%%How much thrust do we need?
%%%Well Let's plot drag
Cd0 = 0.01; %%%This is the drag coefficient at zero lift
Cd2 = 0.037;

%%%Now let's get the Wing Drag Coefficients
CD0 = Cd0/(1+Cd0/(pi*e*AR))
CD2 = Cd2/(1+Cd2/(pi*e*AR))
CDf = 0.2; %%%Fuselage drag (get this from another program. SolidWorks, DATCOMM)

%%%Now we can compute CD
CD = CD0 + CD2*alpha.^2 + CDf;

figure()
plot(alpha_deg,CD)
xlabel('AoA (deg)')
ylabel('CD')

%%%%D = 0.5*rho*V^2*S*CD
%%Vary flight speed from 5:20 and we've already done that
%%%Compute AoA which we've done on line 82
%%Using that AoA compute Cl
Thrust_req = 0.5*rho*V.^2*S.*(CD0 + CD2*AoA_req.^2 + CDf);
figure()
plot(V,Thrust_req*16)
xlabel('Velocity (ft/s)')
ylabel('Thrust (oz)')

%%%Empennage design
%%%Assume percentage of main wing
bt = 0.5*b
ARt = AR
St = bt^2/ARt
ct = St/bt

xcg = ((weight_servo*x_servo+weight_battery*x_battery+weight_esc*x_esc+weight_motor*x_motor+weight_prop*x_prop+weight_receiver*x_receiver+weight_control_horns_rods*x_control_horns_rods+weight_fuselage*x_fuselage)/weight)/12

%%%%Then compute the aerodynamic center
xac = c/4

Sm = (xcg - xac)/c


