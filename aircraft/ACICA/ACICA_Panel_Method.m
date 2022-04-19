%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                           %%
%%%        AIRCRAFT DESIGN CODE               %%
%%%	  ORIGINAL EDIT: EMILY LEYLEK             %%
%%%	  CURRENT EDIT: CARLOS MONTALVO           %%
%%%	                                          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Aerodynamic Coefficient and Inertia Calculator for Aircraft(ACICA)
%%THIS CODE IS BROKEN INTO THREE PARTS.
%1. INPUTS. 
%CHOOSE MAIN WING, TAIL WING AND FUSELAGE GEOMETRIC DIMENSIONS
%AS WELL AS THE MATERIAL DENSITY OF THE AIRCRAFT
%1. CALCULATIONS
%THIS WILL CALCULATE THE CG OF THE AIRCRAFT, THE AREA OF THE WING,
%THE MOMENT OF INERTIA OF THE AIRCRAFT AND THE TOTAL WEIGHT
%2. INPUTS
%PICK A DESIRED CRUISE VELOCITY
%2. CALCULATIONS
%THIS WILL CALCULATE CL_REQ AND Re
%3. INPUTS
%CHOOSE AN AIRFOIL FROM THE FOLDER AIRFOIL_DATA
%3. CALCULATIONS
%THIS WILL CALCULATE: ALL STABILITY DERIVATIVES, ALFA_CRUISE,
%ELEVATOR_DEFLECTION_CRUISE, RUDDER_MAX, AND CREATE A THRUST
%VS. VELOCITY PLOT.
%4. INPUTS
%SET CHORD LENGTHS OF AILERONS,RUDDER,AND ELEVATOR
%4. CALCULATIONS
%CALCULATE CONTROL DERIVATIVES

clear
clc
close all

%%%%%%%%%%%%%%%%%%%%%%FLAGS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PLOTS = 1;
OUTPUT = 0;

%%%%%%%%%%%%%%%%%%GEOMETRIC INPUTS FOR WING AND TAIL%%%%%%%%%%%%%%%%%

N = 10000; %%%Number of elements

%%Main Wing(m)
b = (56.49/12)/3.28
c = (16.22/12)/3.28
maint = 0.09; %thickness of mainwing
taper = 1.0;

%%%%This is what Marzieh gave me
AIRFOIL_CENTER = ((6.39/12)/3.28)/c;
%%%% Center = 13.5

%%%FLAT PLATE
AIRFOIL_CENTER = 0.25;
%%% Center = 11.1 - Note there is no way it could be lower than this. From
%%% what it looks like, flat plates have the furthest forward center of
%%% pressure. Even the NACA 0012 has a 25.75% aerodynamic center. So
%%% putting the aero center at 9" was insanely stable. 10.5 would have
%%% definitely been better. Thing is I wanted a static margin of at least
%%% 2" but I guess a static margin of 0.5" is ok. Unless of course the aero
%%% center is actually 12.25. If that's the case Saturday had a static
%%% margin of -0.75" and the CG at 9" is a static margin of 3.25" which is
%%% way more than I actually wanted. Thing is, I wasn't planning on having
%%% the aircraft break every single test flight. The plan was to try and
%%% fly it. Then move the cg. Fly it again. Then move the cg again. Not
%%% move the cg, break it, etc. It's very wasteful and poorly done in my
%%% opinion. 

%%%NACA 2412
AIRFOIL_CENTER = 0.32;
%%% Center = 12.25

%%%%Supposedly the CG on Saturday was 13"
%%%%The minimum would be 11 if it was a flat plate
%%%%So put it at 10.5" (9 was definitely WAY too conservative)
%%%%Then instead of trying 9.75 and 10.5 we would have tried 11.25 and 12
%%%%Which still would have been less than 13"

%use positive value to indicate aft of nose
xlocation = 0.0; %location of leading edge of wing w.r.t to nose of fuselage
sweep = 30*pi/180;
%%HTail(m)
bt = 0;
ct = 0;
tapert = 1;
tailHt = 0.00; %thickness of horizontaltailwing
alpha_i_t = 0*pi/180; %incidence angle of horizontal tail [rad]
%%VTail(m)
bv = 0.1;
cv = c/1.5;
tailVt = 0.09; %thickness of verticaltailwing
taperv = 1;
%%Fuselage(m)
df = 0;
L = 0.1; %%%This also dictates the position of the horizontal and vertical stabilizers
%Material Density(kg*m^3)
density = 3.2; %%styrofoam

%MASSES OF OTHER PARTS
battmass = 0.0034;
propmass = 0.0152; 
servomass = 0.0093*2; %Two Servos
speedctrlmass = 0.0046;
gpsmass = 0.0;
cameramass = 0.00;
apmodemmass = 0.00; %0.0446; Use GinaMote instead
motormass = 0.0045;
 
%%%%%%%%%%%%%%%%%AERODYNAMIC INPUTS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Vcruise = 8.0; %(m/s)

%%%%%%%%%%%%%%%%%AIRFOIL%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%Directory
rdir = 'Airfoil_Data/';

%%Main Wing
%clcdfile = 'NACA6409clcdRe100k.txt';
%clcmfile = 'NACA6409clcmRe100k.txt';
clcdfile = 'NACA64AclcdRe100k.txt';
clcmfile = 'NACA64AclcmRe100k.txt';

%%Tail Wing
tclcdfile = 'NACA64AclcdRe100k.txt';
tclcmfile = 'NACA64AclcmRe100k.txt';

%%%%%%%%%%%%%%%%%%%%%CONTROL GEOMETRY%%%%%%%%%%%%%%%%%%%%%%%%

elevc = 0.0175;  %%must be > 10% of root chord
aileronc = 0.0175;
rudc = 0.0195;
sr = 0; %%starting length for rudder
er = bv; %%ending length of rudder along elevator
sa = 0;
ea = b/2;
se = 0;
ee = bt/2;
stepe = 0.001; %%section width for elevator
stepr = 0.001; %%rudder
stepa = 0.001; %%aileron

%%%%%%%%%%%%%%%%%%%%STEP 1 CALCULATIONS%%%%%%%%%%%%%%%%%%%%

%%Compute Mass of All parts

%%Wing(assume rectangular)
mainwingthickness = maint*c;
ctip = c*taper;
S = (1/2)*b*(c+taper*c);
Vwing = S*mainwingthickness;
mwing = density*Vwing;

%%Fuselage
AreaF = pi*(df/2)^2;
VolumeF = L*AreaF;
mfuselage = density*VolumeF;

%Horizontal Tail
tailwingthickness = tailHt*ct;
SH = (1/2)*bt*(ct+taper*ct);
cttip = ct*tapert;
Vol_tail = (1/2)*bt*(ct+tapert*ct)*tailwingthickness;
mtail = density*Vol_tail;

%Vertical Tail
vtailwingthickness = tailVt*cv;
SV = (1/2)*bv*(cv+taperv*cv);
cvtip = taperv*cv;
Vol_Vtail = (1/2)*bv*(cv+taperv*cv)*vtailwingthickness;
mVtail = density*Vol_Vtail;

%%Consolidate masses
mass_vec = [mwing;mfuselage;mtail;mVtail;battmass;propmass;servomass;speedctrlmass;gpsmass;cameramass;apmodemmass;motormass];

%%Compute Total Mass of Aircraft
mass = sum(mass_vec);

Weight = 9.81*mass;

%%%%%Compute CG Locations of all Sections(measured from body fixed
%frame at nose of aircraft)

%main wing %%%THIS NEEDS TO GET FIXED FOR SWEPT WINGS
ycoord = linspace(0,b/2,N);
dy = ycoord(2)-ycoord(1);
xwingcg = 0;
for idx = 1:N-1
  ci = c + (2/b)*c*(taper-1)*ycoord(idx);
  xwingcg = xwingcg + xlocation + sin(sweep)*ycoord(idx) + ci/2;
end
maincg = -[xwingcg;0;0];

%Horizontal tail
Htailcg = -[L-ct/2;0;0];

%%Vertical tail
Vtailcg = -[L+cv/2;0;-bv/2];

%Fuselage
fcg = -[L/2;0;0];

%%Propellor
propcg = [0.005;0;0];
%%Motor
motorlength = 0.025; %length of motor [m]
motorcg = -[(0.01+0.5*motorlength);0;0];
%%Battery
battlength = 0.06604; %length of battery [m]
battcg = -[(0.01+motorlength+0.003+0.5*battlength);0;0]; %assuming battery is directly behind motor
%%Speed Control
spdctrllength = 0.034925; %length of speed controller [m]
spdctrlcg = -[(0.01+motorlength+0.003+0.5*spdctrllength);0;-(0.01778/2 + 0.003175/2)]; %speed controller positioned on top of battery
%Place GPS and AP/Modem on top of each other
apmodemcg = maincg - [0.0127;0;-0.0127/2];
gpscg = maincg - [0.0428625/2;0;0.0127/2;];
camcg = gpscg + [- 0.0428625/2 - 0.001 - 0.0127/2;0;0];
servocg = camcg + [- 0.0127/2 - 0.005 - 0.0206375/2;0;0]; %0.0206375 is the length of the servo in [m]

%%Consolidate Cgs
cgs = [maincg,fcg,Htailcg,Vtailcg,battcg,propcg,servocg,spdctrlcg,gpscg,camcg,apmodemcg,motorcg];

%%Now we can compute the cg location of the entire aircraft
ac_cg = [0;0;0];
for ii = 1:length(mass_vec)
  ac_cg = ac_cg + mass_vec(ii)*cgs(:,ii);
end
ac_cg = ac_cg/mass;

%%%DRAW THE AIRCRAFT TO VISUALIZE IT AND MAKE SURE IT LOOKS NORMAL

if PLOTS
  plottool(1,'Aircraft',12,'x','y','z','',[56 18]);
  %%Main Wing
  yleft = 0 - b/2;
  yright = 0 + b/2;
  xmidtop = 0;
  xmidbottom = c;
  xtip_top = xlocation + sin(sweep)*b/2;
  xtip_bottom = xtip_top + c*taper;
  X = -[0 xtip_top xtip_bottom xmidbottom xtip_bottom xtip_top 0];
  Y = [0 yright yright 0 yleft yleft 0];
  patch(X,Y,'b');
  %%Tail Wing
  CubeDraw(ct,bt,tailwingthickness,Htailcg(1),Htailcg(2),Htailcg(3),0,0,0,'blue');
  %%Horizontal Tail
  CubeDraw(cv,bv,vtailwingthickness,Vtailcg(1),Vtailcg(2),Vtailcg(3),pi/2,0,0,'blue');
  %%Fuselage
  [x,y,z] = cylinder; %radius equal to 1
  x = (df/2).*x;y = (df/2).*y;z = L.*z;
  surf(z+2*fcg(1),x+fcg(2),y+fcg(3));
  %%CG of aircraft
  plot3(ac_cg(1),0,df/2,'r*','MarkerSize',20);
  plength = max([L b])/2;
  axis([-2*plength 0 -plength plength -plength plength])
  axis equal
end

%%INERTIAS

%Main Wing
IxxW = (1/12)*mwing*(b^2 + mainwingthickness^2);
IyyW = (1/12)*mwing*(c^2 + mainwingthickness^2);
IzzW = (1/12)*mwing*(b^2 + c^2);

%%Fuselage
Ixxf = (1/2)*mfuselage*(df/2)^2;
Iyyf = (1/12)*mfuselage*(L^2 + 3*(df/2)^2);
Izzf = Iyyf;

%%Horizontal Tail
IxxH = (1/12)*mtail*(bt^2 + tailwingthickness^2);
IyyH = (1/12)*mtail*(ct^2 + tailwingthickness^2);
IzzH = (1/12)*mtail*(bt^2 + ct^2);

%%Vertical Tail
IxxV = (1/12)*mVtail*(bv^2 + vtailwingthickness^2);
IyyV = (1/12)*mVtail*(cv^2 + vtailwingthickness^2);
IzzV = (1/12)*mVtail*(bv^2 + cv^2);

%Inertias of Other Parts (including converting from g-in^2 to kg-m^2)
Ixxprop = 57.475*(1/1000)*0.0254^2; Iyyprop = 28.7375*(1/1000)*0.0254^2; Izzprop = 28.7375*(1/1000)*0.0254^2;
Ixxservo = 0.146069336*(1/1000)*0.0254^2; Iyyservo = 0.236889648*(1/1000)*0.0254^2; Izzservo = 0.164990234*(1/1000)*0.0254^2;
Ixxspdctrl = 0.458447266*(1/1000)*0.0254^2; Iyyspdctrl = 0.976953125*(1/1000)*0.0254^2; Izzspdctrl = 1.419384766*(1/1000)*0.0254^2;
Ixxgps = 0.303385417*(1/1000)*0.0254^2; Iyygps = 1.032552083*(1/1000)*0.0254^2; Izzgps = 1.169270833*(1/1000)*0.0254^2;
Ixxapmodem = 1.161458333*(1/1000)*0.0254^2; Iyyapmodem = 1.161458333*(1/1000)*0.0254^2; Izzapmodem = 1.858333333*(1/1000)*0.0254^2;
Ixxmotor = 5.5625*(1/1000)*0.0254^2; Iyymotor = 5.620442708*(1/1000)*0.0254^2; Izzmotor = 5.620442708*(1/1000)*0.0254^2;
Ixxcam = 0.073958333*(1/1000)*0.0254^2; Iyycam = 0.073958333*(1/1000)*0.0254^2; Izzcam = 0.073958333*(1/1000)*0.0254^2;
Ixxshock = 0;Iyyshock = 0;Izzshock = 0;

battl = battlength; battw = 0.03429; batth = 0.01778;
Ixxbatt = battmass*(battw^2 + batth^2)/12;
Iyybatt = battmass*(battl^2 + batth^2)/12;
Izzbatt = battmass*(battl^2 + battw^2)/12;

%%Sum up inertias including the cg offset
Ixx_vec = [IxxW,Ixxf,IxxH,IxxV,Ixxbatt,Ixxprop,Ixxservo,Ixxspdctrl,Ixxgps,Ixxcam,Ixxapmodem,Ixxmotor];
Iyy_vec = [IyyW,Iyyf,IyyH,IyyV,Iyybatt,Iyyprop,Iyyservo,Iyyspdctrl,Iyygps,Iyycam,Iyyapmodem,Iyymotor];
Izz_vec = [IzzW,Izzf,IzzH,IzzV,Izzbatt,Izzprop,Izzservo,Izzspdctrl,Izzgps,Izzcam,Izzapmodem,Izzmotor];
II = zeros(3,3);
for ii = 1:length(mass_vec)
  IP = [Ixx_vec(ii) 0 0;0 Iyy_vec(ii) 0;0 0 Izz_vec(ii)];
  rcp = cgs(:,ii) - ac_cg;
  rcpM = skew(rcp);
  m = mass_vec(ii);
  IC = IP + m*rcpM*rcpM';
  II = II + IC;
end
Ixx = II(1,1);
Iyy = II(2,2);
Izz = II(3,3);

%%%%%%%%%%%%%%%%%%%%%%%STEP 2 CALCULATIONS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%Using L = W = 1/2*rho*S*Vcr^2*CL we can compute CL_req in cruise

rho = 1.225; %kg*m^3 assume at sea-level
qinf = 1/2*rho*S*(Vcruise^2);
CL_req = Weight/qinf;

%%Compute Reynolds Number
format long g
mu = 1.4607*(10^-5);
Re = (rho*Vcruise*c)/mu
format

%%%%%%%%%%%%%%%%%%%%%%%STEP 3 CALCULATIONS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%Compute Aerodynamic Center

%Main Wing
AR = (b^2)/S; %aspect ratio
cbar = 0;
%mean aerodynamic chord
aero_center = 0;
for idx = 1:N-1
    xi = xlocation + sin(sweep)*ycoord(idx);
    ci = c + (2/b)*c*(taper-1)*ycoord(idx);
    aeroi = xi + ci*AIRFOIL_CENTER;
    aero_center = aero_center + aeroi*dy;
    cbar = cbar + ci*dy;
end
aero_center = aero_center*2/b
aero_center_in = aero_center*3.28*12
aero_over_chord = aero_center/c
cbar = cbar*2/b

%%%Read in Airfoil Data and Extract lift,drag and moment
aclcd = dlmread([rdir,clcdfile]);
aclcm = dlmread([rdir,clcmfile]);

alpha1 = aclcm(:,1)*pi/180;
clift_wing = aclcm(:,2);
cm_wing = aclcm(:,3);

alpha2 = aclcd(:,1)*pi/180;
cl2 = aclcd(:,2);
cd_wing = aclcd(:,3);

%%Use Least Squares to fit CM = CM0 + CM1*alpha
Cm_alpha_wing = polyfit(alpha1, cm_wing, 1);
CMA_wing = Cm_alpha_wing(1);
CMAC_wing = Cm_alpha_wing(2);

%%Use Least Squares to fit CL = CL0 + CL1*alpha
clift_alpha = polyfit(alpha1(1:17), clift_wing(1:17), 1);
a0 = clift_alpha(1);
CL0_wing = clift_alpha(2);
if CL0_wing < 0
  CL0_wing = 0;
end

%%Use least squares to fit drag polar to CD = CD0 + CD2*alpha^2
loc = find(alpha2>0,1);
X = cd_wing(loc:end);
H = [ones(length(X),1) alpha2(loc:end).^2];
thetastar = inv(H'*H)*H'*X;
CD0_wing = thetastar(1);
CD2_wing = thetastar(2);

%Horizontal and Vertical Tails
taclcd = dlmread([rdir,tclcdfile]);
taclcm = dlmread([rdir,tclcmfile]);

alpha1_tail = taclcm(:,1)*pi/180;
clift_tail = taclcm(:,2);
cm_tail = taclcm(:,3);
CMAC_tail = mean(cm_tail);

alpha2_tail = taclcd(:,1)*pi/180;
clift2_tail = taclcd(:,2);
cd_tail = taclcd(:,3);

%%Use least squares to fit CL = CL0 + CL1*alpha
Cl_alpha_tail = polyfit(alpha1_tail(1:17), clift_tail(1:17), 1);
a0_tail = Cl_alpha_tail(1);
CL0_tail = Cl_alpha_tail(2);
if CL0_tail < 0
  CL0_tail = 0;
end

%%Use least Squares to fit CD = CD0 + CD2*alpha^2
loc = find(alpha2_tail>0,1);
X = cd_tail(loc:end);
H = [ones(length(X),1) alpha2_tail(loc:end).^2];
thetastar = inv(H'*H)*H'*X;
CD0_tail = thetastar(1);
CD2_tail = thetastar(2);

%%Fuselage Skin Fristion Drag via Form Factor

f = L/(2*(df/2)); 
FF = 1+60/(f^3) + f/400; %Fuselage form factor
Re_f = (rho*Vcruise*L)/(mu); %Reynolds # of fuselage
Cf = 1.328/sqrt(Re_f); %Skin friction coefficient, laminar
Swet_f = 2*pi*((df/2))*L;
CD0_f = FF*Cf*Swet_f/(2*(S/2));

%%Add in efficiency of Wing to create 3D Lift

e = 1.78*(1-0.045*AR^0.68)-0.64; %From Raymer %0.9899; %oswald's efficiency factor [nd]
CL_alpha_wing = a0/(1+a0/(pi*e*AR));
a_wing = CL_alpha_wing;
ARht = (bt^2)/SH; %aspect ratio
if ARht ~= 0
  CL_alpha_tail = a0_tail/(1+a0_tail/(pi*e*ARht));
else
  CL_alpha_tail = 0;
end
a_tail = CL_alpha_tail;

%%Combine Effects of Tail and Wing to create Total Lift,and Drag

CL0_total = CL0_wing + CL0_tail*(SH/S);
CLA_total = CL_alpha_wing + CL_alpha_tail*(SH/S);
CD0_total = CD0_wing + CD0_tail*(SH/S) + CD0_f;
CD2_total = CD2_wing + CD2_tail*(SH/S);

%%Compute Cm0
lwingcg = xh_nw-xwingcg; %vector pointing from vehicle c.g. to wing a.c.
lhtailcg = xh_nwht-ac_cg(1); %vector pointing from vehicle c.g. to horizontal tail a.c.
lbart = abs(L-ct-xbarht-xlocation-xbar);
Vhbar = lbart*(SH)/(cbar*S);
Vh = (abs(lhtailcg)/cbar)*(SH/S);
eps0 = 0; %%assume no downwash component
deda = 0; %%assume no downwash component

%%%CMO = CMAC_Wing + CMAC_Tail + Lift_Tail*moment_arm
Cm0_total = CMAC_wing + CMAC_tail + a_tail*Vhbar*(eps0-alpha_i_t)*(1-(a_tail/a_wing)*(SH/S)*(1-deda));

%%Compute Cm_alpha_total
%%%Assume M = (1/2)*rho*V^2*S*cbar*(Cm_alpha_total*alpha)
M = 0;
alpha_test = 5*pi/180;
for idx = 1:N-1
  xi = xlocation + sin(sweep)*ycoord(idx);
  ci = c + (2/b)*c*(taper-1)*ycoord(idx);
  aeroi = xi + ci/4;
  deli = xwingcg - aeroi;
  Si = ci*dy;
  M = M + Si*ci*deli*CL_alpha_wing*alpha;
end
Cm_alpha_total = M/(alpha*S*cbar);

if PLOTS
  %plot3(xh_nw,0,-df/2,'r*','MarkerSize',10)
  %plot3(xh_nwht,0,-df/2,'r*','MarkerSize',10)
  %plot3(xh_nwvt,0,bv,'r*','MarkerSize',10)
  plot3(xhn,0,df/2,'g*','MarkerSize',20)
end
if Cm_alpha_total > 0
  disp('Unstable Configuration')
  break;
end

if PLOTS
  %%Plot Main Wing Lift,Drag, Moment
  adeg = alpha1.*180/pi;
  plottool(1,'Wing Lift Coefficient',12,'AoA(deg)','Wing Lift Coefficient')
  plot(adeg,clift_wing,'b-')
  plot(adeg,CL0_wing+a0.*alpha1,'r--')
  plottool(1,'Wing Moment Coefficient',12,'AoA(deg)','Wing Moment Coefficient')
  plot(adeg,cm_wing)
  plot(adeg,CMAC_wing+CMA_wing*alpha1,'r--')
  adeg = alpha2.*180/pi;
  plottool(1,'Wing Drag Coefficient',12,'AoA(deg)','Wing_Drag Coefficient')
  plot(adeg,cd_wing)
  plot(adeg,CD0_wing+CD2_wing.*alpha2.^2,'r--')
  %%%Plot Tail Lift,Drag,Moment
  adeg = alpha1_tail.*180/pi;
  plottool(1,'Tail Lift Coefficient',12,'AoA(deg)','Tail Lift Coefficient')
  plot(adeg,clift_tail,'b-')
  plot(adeg,CL0_tail+a0_tail.*alpha1_tail,'r--')
  plottool(1,'Tail Moment Coefficient',12,'AoA(deg)','Tail Moment Coefficient')
  plot(adeg,cm_tail)
  adeg = alpha2_tail.*180/pi;
  plottool(1,'Tail Drag Coefficient',12,'AoA(deg)','Tail_Drag Coefficient')
  plot(adeg,cd_tail)
  plot(adeg,CD0_tail+CD2_tail.*alpha2_tail.^2,'r--')
end

%Pitch Stiffness Coefficients
CL_q = 2*a_tail*Vh;
Cm_q = -2*a_tail*Vh*abs(lhtailcg)/cbar;

%% Lateral Stability Derivatives
a_vtail = a0_tail/(1+a0_tail/(pi*e*ARvt)); %%3D effects
CL_alpha_vtail = a_vtail;
lf = abs(-L+cv-xbarvt-ac_cg(1));
Vv = SV*lf/(S*b);
Cn_beta_tail = Vv*a_vtail;
Cnbeta_fuse = -1.3*(VolumeF/(S*b)*(df+0.0127+df)/2/(9.81*(mfuselage)));
Cn_beta_total = Cn_beta_tail + Cnbeta_fuse;
Cy_beta = -CL_alpha_vtail*SV/S;
Cl_beta = -CL_alpha_vtail*SV*ybarvt/(S*b);
Cy_p = -CL_alpha_vtail*SV*2*ybarvt/(S*b);
Cn_p = CL_alpha_vtail*Vv*2*ybarvt/(b);
Cl_p = -CL_alpha_vtail*SV*ybarvt^2*2/(S*(b)^2) - CL_alpha_tail*SH*ybarht^2*2/(S*(b)^2);
Cy_r = CL_alpha_vtail*SV*2*lf/(S*b);
Cl_r = CL_alpha_vtail*SV*ybarvt*2*lf/(S*(b)^2);
Cn_r = -CL_alpha_vtail*Vv*2*lf/(b);

%%%%%PERFORMANCE CHARACTERISTICS%%%%%%%%

%%Stall Requirements
CL_max = 0.9*(max(clift_wing) + (SH/S)*max(clift_tail));
Vstall = sqrt(Weight/(S*0.5*rho*CL_max));

%%Plot L/D
if PLOTS
  V = [Vstall:0.1:50];
  CL_req = V.*0;
  aoa = CL_req;
  CD = aoa;
  Induced = aoa;
  Parasite = aoa;
  TReq = aoa;
  for ii = 1:length(V)
    qinf = 1/2*rho*V(ii)^2*S;
    CL(ii) = Weight/qinf;
    aoa(ii) = (CL(ii)-CL0_total)/(CLA_total);
    Cdinduced = (CL(ii)^2)/(pi*e*AR);
    Cdparasite = CD0_total + CD2_total*aoa(ii)^2;
    CD(ii) = Cdparasite + Cdinduced;
    Induced(ii) = Cdinduced*qinf;
    Parasite(ii) = Cdparasite*qinf;
    TReq(ii) = CD(ii)*qinf;
  end
  LDratio = CL./CD;

  plottool(1,'Thrust',12,'Velocity(m/s)','Thrust Required(N)')
  plot(V,TReq,'b-','LineWidth',2)
  plot(V,Induced,'r--')
  plot(V,Parasite,'g--')
  plottool(1,'L/D',12,'Velocity(m/s)','L/D')
  plot(V, LDratio)
  plottool(1,'AoA',12,'Velocity(m/s)','AoA required')
  plot(V,aoa*180/pi)
  plottool(1,'Drag',12,'Velocity(m/s)','CD')
  plot(V,CD)
end

%%%%%%%%%%%%%STEP 4 CALCULATIONS%%%%%%%%%%%%%%%%%

%%From Graph 1
Kfmat = [0.72; 0.7; 0.675; 0.65; 0.61; 0.575; 0.55; 0.525;0.5]; %for elev. deflection of 30 deg.
cKf = [0.10; 0.15; 0.2; 0.25; 0.3; 0.4; 0.5; 0.6; 0.7];
%%From Graph 2
dclddfmat = [1.5; 2.5; 3.2; 3.6; 4; 4.4; 4.75; 5; 5.25; 5.5; 6.0; 6.5]; %[1/rad]
cdcldd = [0.05; 0.1; 0.15; 0.2; 0.25; 0.3; 0.35; 0.4; 0.45; 0.5; 0.6; 0.7];

%ELEVATOR
lift = 0;
momenti = 0;
if elevc/ct < 0.1
  elevc = ct*0.1;
  sprintf('%s \n',['Elevator too small, elevc is now: ',num2str(elevc)])
end
for ii = se:stepe:ee
  chord = ((ct-ct*tapert)/(0-(bt/2)))*ii+ct;
  elevratio = elevc/chord;
  Kf = interp1(cKf,Kfmat,elevratio);
  dclddf = interp1(cdcldd,dclddfmat,elevratio);
  Si = stepe*chord;
  lift = lift + 0.85*Kf*dclddf*Si;
  hi = (-L+ct-ii*(ct-cttip)/(bt/2)-0.25*chord)-ac_cg(1);
  momenti = momenti + 0.85*Kf*dclddf*hi*Si;
end
Clift_delev = 2*lift/S; 
Cd_delev = 2*(Clift_delev/2)^2/(pi*e*ARht);
Cm_delev = 2*momenti/(S*cbar);

%%AILERONS
if aileronc/c < 0.1
  aileronc = c*0.1;
  sprintf('%s \n',['Aileron too small, aileronc is now: ',num2str(aileronc)])
end
lift = 0;
momenti = 0;
roll = 0;
for ii = sa:stepa:ea
  chord = ((c-c*taper)/(0-(b/2)))*ii+c;
  ailrat = aileronc/chord;
  Kf = interp1(cKf,Kfmat,ailrat);
  dclddf = interp1(cdcldd,dclddfmat,ailrat);
  Si = stepa*chord;
  roll = roll + 0.85*Kf*dclddf*ii*Si;
  lift = lift + 0.85*Kf*dclddf*Si;
  hi = (-ii*(ct-cttip)/(bt/2)-0.25*chord)-ac_cg(1);
  momenti = momenti + 0.85*Kf*dclddf*hi*Si;
end
C_l_da = -2*roll/(S*b); 
C_lift_da = 2*lift/S; 
C_d_da = 2*(C_lift_da/2)^2/(pi*e*AR);
C_n_da = -C_d_da*ybar/(2*b);
C_m_da = 2*momenti/(S*cbar);

%%%RUDDER
sumy = 0;
suml = 0;
sumn = 0;
if rudc/cv < 0.1
  rudcc = cv*0.1;
  sprintf('%s \n',['Rudder too small, rudc is now: ',num2str(rudc)])
end
for ii = sr:stepr:er
    chord = ((cv-cv*taperv)/(0-bv))*ii+cv;
    rudrat = rudc/chord;
    Kf = interp1(cKf,Kfmat,rudrat);
    dclddf = interp1(cdcldd,dclddfmat,rudrat);
    Si = stepr*chord;
    vxi = (-L+cv-ii*(cv-cvtip)/bv-0.25*chord)-ac_cg(1);
    sumy = sumy + 0.85*Kf*dclddf*Si;
    suml = suml + 0.85*Kf*dclddf*ii*Si;
    sumn = sumn + 0.85*Kf*dclddf*vxi*Si;
end
Cy_drud = sumy/S;
Cl_drud = suml/(S*b);
Cn_drud = sumn/(S*b);

%%Consolidate Data
MASSFILE = {[num2str(Weight),'	!Weight_(N)_(assume_m1=m2)'];['0	!SLCG(m)'];['0	!BLCG(m)'];['0	!WLCG(m)'];[num2str(Ixx),'        !Ixx(kg*m^2)'];[num2str(Iyy),'	!Iyy'];[num2str(Izz),'	!Izz'];[num2str(II(1,2)),'	!Ixy'];[num2str(II(1,3)),'	!Ixz'];[num2str(II(2,3)),'	!Iyz']}
AEROFILE = {[num2str(CL0_total),'		!C_L_0'];
	    [num2str(CD0_total),' 		!C_D_0'];
	    [num2str(Cm0_total),' 		!C_m_0'];
	    ['0.00000 		!C_D_u'];
	    [num2str(CLA_total),'	        !C_L_alpha'];
	    [num2str(CD2_total),'     	!C_D_alpha2'];
	    [num2str(Cm_alpha_total),' 	!C_m_alpha'];
	    ['0.00000 		!C_m_alpha_dot'];
	    ['0.00000 		!C_m_u'];
	    [num2str(CL_q),' 		!C_L_q'];
	    [num2str(Cm_q),'	 	!C_m_q'];
	    [num2str(Clift_delev),' 		!C_L_de'];
	    [num2str(Cm_delev),'     	!C_m_de'];
	    [num2str(5),' 		!C_x_delThrust'];
	    [num2str(Cy_beta),' 	!C_y_beta'];
	    [num2str(Cl_beta),'      	!C_l_beta'];
	    [num2str(Cn_beta_total),'      	!C_n_beta'];
	    [num2str(Cl_p),'         !C_l_p'];
	    [num2str(Cn_p),'      	!C_n_p'];
	    [num2str(Cl_r),'		!C_l_r'];
	    [num2str(Cn_r),'		!C_n_r'];
	    [num2str(C_l_da),'      	!C_l_da'];
	    [num2str(C_n_da),'		!C_n_da'];
	    [num2str(Cy_drud),'		!C_y_dr'];
	    [num2str(Cl_drud),'      	!C_l_dr'];
	    [num2str(Cn_drud),'		!C_n_dr'];
	    ['0.0057		!C_y_p'];
	    [num2str(S),'		!Reference_Area(m^2)'];
	    [num2str(b),'		!Wingspan(m)'];
	    [num2str(cbar),'		!Mean_chord(m)'];
	    [num2str(Vcruise),'	!Trim_Velocity(m/s)']}

if OUTPUT
  %%OUTPUT EVERYTHING TO A .MASS and .AERO FILE
  fid = fopen('UAV.MASS','wb');
  for ii = 1:length(MASSFILE)
    fprintf(fid,'%s \n',MASSFILE{ii,:});
  end
  fid = fopen('UAV.AERO','wb');
  for ii = 1:length(AEROFILE)
    fprintf(fid,'%s \n',AEROFILE{ii,:});
  end
end

