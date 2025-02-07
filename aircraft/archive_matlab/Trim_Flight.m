%%%%Lets Trim an Airplane, Yea? Yea.
clear
clc
close all

c = 0.3215;
cbar = c; %%%Rectangular wing
b = 2.04;
S = b*c;
AR = b^2/S;
e = 0.8;

ct = 0.12;
bt = 0.3;
St = ct*bt;
ARt = bt^2/St;
et = 0.4; %%%this is really bad. Let's see what happens
lt = 0.7;
VH = St*lt/(S*cbar);

AoAdeg = (-30:0.1:30);
AoA = AoAdeg*pi/180;

%%%Airfoil for main wing? 0012? Change to 2412
Cmac = 0.1; %%%Zero for symmetric and nozero for cambered
Clalpha = 6.0; %%for a flat plate but symmetric is an ok approximation %%%Made this up for 2412
alpha0 = -2*pi/180; %%%zero lift angle of attack for a cambered airfoil
Clalphat = 2*pi; %%%for a flat plate but symmetric is an ok approximation
CLalpha = Clalpha/(1+Clalpha/(pi*e*AR));
CLalphat = Clalphat/(1+Clalphat/(pi*et*ARt));

CL = CLalpha*(AoA-alpha0);
CLt = CLalphat*AoA;

xacW = c/4; %%%for a flat plate but is a good approximation for symmetric airfoils
xact = c/4 + lt;
xac = (xacW*S + xact*St)/(St + S)

CLwb = CL + St/S*CLt;

figure()
plot(AoAdeg,CL,'b-')
hold on
plot(AoAdeg,CLt,'r-')
plot(AoAdeg,CLwb,'g-')
legend('Main Wing','Tail','Wing-Body')
ylabel('Lift Coefficient')
xlabel('Angle of Attack (deg)')

%%%What's First? Pitch Moment Coefficient?
xcg = -0.12;
%%%Static Margin
xsmbar = xcg/cbar - xac/cbar;
%%%break into 3 components
%%%This is from the airfoil
Cm1 = Cmac*ones(length(AoA),1);
%%%This is from the wingbody which includes the tail
Cm2 = CLwb*xsmbar;
%%%Combining like terms created this which is pitch moment of the tail
Cm3 = -VH*CLt;
%%%Go back and make sure you understand where the derivation for this
%%%equation comes from
Cmde = -0.8; %%%Dr. C pulled this out of thin air
de = -5*pi/180; %%%Jesse just felt like it
Cm = Cmac + CLwb*xsmbar - VH*CLt + Cmde*de;

figure()
plot(AoAdeg,Cm1,'b-')
hold on
plot(AoAdeg,Cm2,'r-')
plot(AoAdeg,Cm3,'g-')
plot(AoAdeg,Cm,'k-')
xlabel('Angle of Attack (deg)')
ylabel('Pitch Moment Coefficient')
legend('Mac','LiftWB','Tail')

%%%%Weight and Balance
W = 5.6*9.81;

%%%%This is similar to equations 3.26 and 3.27 in your notes
%%%This is solving for trim angle of attack. Take the Cm equation and set
%%%it to zero. Then solve for angle of attack
d = (CLalpha+St/S*CLalphat)*xsmbar - VH*CLalphat;
alphatrim = (-Cmac + CLalpha*alpha0*xsmbar)/d;

%%%In your homework you have Cm = Cm0 + Cma * alpha so alphatrim is just
%%% alphatrim = -Cm0/Cma

%%%Ok so what's this?
%%%You're calculating the actual lift coefficients at trim
CLtrim = CLalpha*(alphatrim-alpha0); %%%Wing
CLttrim = CLalphat*alphatrim; %%%tail
CLwbtrim = CLtrim + St/S*CLttrim; %%%wing body
rho = 1.225; %%%Density of air at sea-level
V = 0:0.1:20; %%%meters/sec

Lift = 0.5*rho*V.^2*S*CLwbtrim;

%%%%What's the point of this?
%%%Figure out cruise speed - when Lift = Weight
figure()
plot(V,Lift)
xlabel('Velocity (m/s)')
ylabel('Lift (N)')

%%%%How do we find required Velocity
Vtrim = sqrt(2*W/(rho*S*CLwbtrim)); %%%Solve L = W = 0.5*rho*V^2*S*CLwbtrim

%%%Compute Thrust Generated
CD0 = 0.02; %%%drag coefficient of the entire thing - drag at zero angle of attack
CDa2 = 5.1; %%%drag coefficient as a function of AoA - drag is nonlinear with AoA
CD = CD0 + CDa2*alphatrim.^2; %%%Total drag coefficent
Thrust = 0.5*rho*Vtrim^2*S*CD %%%Thrust required at that trim velocity

TW = Thrust/W

