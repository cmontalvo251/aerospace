
%% Weight Estimation
%Main Wing Airfoil: NACA 6409
%Tail Airfoils: NACA 64A010
%Battery: 7.4V 1250 mAh, 3.1 oz, 3.429 cm x 6.604 cm x 1.778 cm
%Material: ABS 
density = 1040; %material density [kg/m^3]

%Half-Wing Volume Calculations
b = 0.4; %half-span [m]
root = 0.10; %root chord [m]
taper = 0.6;

tip = taper*root; 

area = b*(root+tip)/2; %half-wing area
maxthickness = 0.09*root; %assume thickness is same across span (overestimate)
minthickness = 0.09*tip;
avgthick = (maxthickness+minthickness)/2;
wingvol = area*avgthick; %half wing volume

wingmass = 2*wingvol*density; %full wing mass [kg]

%Fuselage (model as cylinder) & Tail boom

thick = 0.003; 
totallength = 0.41;
tailboomlength = 0.145;
tailboomvol =  (pi*(0.0127+2*thick)^2/4-pi*0.0127^2/4)*tailboomlength;
tailboommass = tailboomvol*density;
din = 0.05; %inner diameter
fuselength = totallength-tailboomlength-0.055;
fusevol = (pi*(din+2*thick)^2/4-pi*din^2/4)*(fuselength) + 1.780576176238355e-005;
fusemass = fusevol*density; %fuselage mass [kg]

%Horizontal Tail:
    %Airfoil: NACA 64A010
htailb = 0.10; %half-span [m]
htailroot = 0.06; %tail root chord [m]
htaper = 0.5;
htailtip = htaper*htailroot;
harea = htailb*(htailroot+htailtip)/2; %half-wing area [m^2]
htailmaxthick = 0.10*htailroot;
htailminthick = 0.10*htailtip;
htailavgthick = (htailmaxthick+htailminthick)/2;
htailvol = htailavgthick*harea; %half-span volume [m^3]
htailmass = 2*htailvol*density; %htailmass [kg]

%Vertical Tail
    %Airfoil: NACA 64A010
vtailb = 0.09;
vtailroot = 0.065;
vtaper = 0.6;
vtailtip = vtaper*vtailroot;
varea = vtailb*(vtailroot+vtailtip)/2;
vtailmaxthick = 0.10*vtailroot;
vtailminthick = 0.10*vtailtip;
vtailavgthick = (vtailmaxthick+vtailminthick)/2;
vtailvol = vtailavgthick*varea;
vtailmass = vtailvol*density;

%Total Structure Weight:
structmass = 0.4352; %(wingmass+fusemass+tailboommass+htailmass+vtailmass); %[N]
xcgstruct = -0.138316941176471;

%Other Parts:
battmass = 0.0878835217;
%From GenMAV:
propmass = 0.0152; 
servomass = 0.0093*2; %Two Servos
speedctrlmass = 0.0246;
gpsmass = 0.016;
cameramass = 0.0071;
apmodemmass = 0.005;%0.0446; Use GinaMote instead
motormass = 0.0445;
shocksmass = 0.091*2; %shockswt

%Total Other Parts Weight [N]:
partsmass = battmass+propmass+servomass+speedctrlmass + gpsmass+cameramass+apmodemmass+motormass+shocksmass;
%Total Estimated Mass [kg]:
TotMass = structmass+partsmass;
%Total Estimated Weight [N]:
Weight = 9.81*TotMass;


%% Trapezoidal Geometry w/ uniform lift coeff. across span, 0 sweep angle
%Main Wing
AR = (2*b)^2/(2*area); %aspect ratio
cbar = 2*root*(1+taper + taper^2)/(3*(1+taper)); %mean aerodynamic chord
ybar = 2*b*(1 + 2*taper)/(6*(1+taper));
dum = (root-tip)/2;  %x-distance from root leading edge to tip leading edge. zero sweep.
xle = 2*b*(1+2*taper)*(dum/b)/(6*(1+taper));
sweep_quart = atan((dum+0.25*tip-0.25*root)/b); %Sweep angle at the quarter chord [rad]
xbar = 0.25*root + ybar*tan(sweep_quart);
h_nw = (xbar-xle)/cbar;

%Horizontal Tail
ARht = (2*htailb)^2/(2*harea); %aspect ratio
cbarht = 2*htailroot*(1+htaper + htaper^2)/(3*(1+htaper)); %mean aerodynamic chord
ybarht = 2*htailb*(1 + 2*htaper)/(6*(1+htaper));
dumht = (htailroot-htailtip);  %x-distance from root leading edge to tip leading edge. zero sweep.
xleht = 2*htailb*(1+2*htaper)*(dumht/htailb)/(6*(1+htaper));
sweep_quartht = atan((dumht+0.25*htailtip-0.25*htailroot)/htailb); %Sweep angle at the quarter chord [rad]
xbarht = 0.25*htailroot + ybarht*tan(sweep_quartht);
h_nwht = (xbarht-xleht)/cbarht;

%Vertical Tail
ARvt = (vtailb)^2/(varea); %aspect ratio
cbarvt = 2*vtailroot*(1+vtaper + vtaper^2)/(3*(1+vtaper)); %mean aerodynamic chord
ybarvt = vtailb*(1 + 2*vtaper)/(6*(1+vtaper));
dumvt = (vtailroot-vtailtip);  %x-distance from root leading edge to tip leading edge. zero sweep.
xlevt = vtailb*(1+2*vtaper)*(dumvt/vtailb)/(6*(1+vtaper));
sweep_quartvt = atan((dumvt+0.25*vtailtip-0.25*vtailroot)/vtailb); %Sweep angle at the quarter chord [rad]
xbarvt = 0.25*vtailroot + ybarvt*tan(sweep_quartvt);
h_nwvt = (xbarvt-xlevt)/cbarvt;

%% CG Estimation -- nose of aircraft is origin
ycg = 0.0; %assume symmetry

%x cg location calculations
xfusecg = -fuselength/2; %x-axis location of the fuselage c.g.
xtailboomcg = -(totallength-tailboomlength/2);
x_rootle_loc = -0.055; %x location of the leading edge of the root chord from the nose[m]
xwingcg = (x_rootle_loc-0.5*root); 
xpropcg = 0.005;
motorlength = 0.025; %length of motor [m]
xmotorcg = -(0.01+0.5*motorlength);
battlength = 0.06604; %length of battery [m]
xbattcg = -(0.01+motorlength+0.003+0.5*battlength); %assuming battery is directly behind motor
xshockcg = x_rootle_loc-0.25*root;
spdctrllength = 0.034925; %length of speed controller [m]
xspdctrlcg = -(0.01+motorlength+0.003+0.5*spdctrllength); %speed controller positioned on top of battery
%Place GPS and AP/Modem on top of each other, behind shocks [Replace with Gina mote?]
xapmodemcg = xwingcg - 0.0127;
xgpscg = xwingcg - 0.0428625/2;
htailside = sqrt((htailroot-htailtip)^2+(htailb)^2);
xhtailcg = -totallength+htailroot - (htailroot/2+(2*htailtip + htailroot)*(htailside^2-htailb^2)/(6*(htailroot^2-htailtip^2)));
vtailside = sqrt((vtailroot-vtailtip)^2+(vtailb)^2);
xvtailcg = -totallength + vtailroot - (vtailroot/2+(2*vtailtip + vtailroot)*(vtailside^2-vtailb^2)/(6*(vtailroot^2-vtailtip^2)));
xcamcg = xgpscg - 0.0428625/2 - 0.001 - 0.0127/2;
xservocg = xcamcg - 0.0127/2 - 0.005 - 0.0206375/2; 

xcg = -0.107991660722541; %(xcgstruct*structmass + xpropcg*propmass + xmotorcg*motormass +...
%     xbattcg*battmass + xshockcg*shocksmass + xspdctrlcg*speedctrlmass + xapmodemcg*apmodemmass +...
%     xgpscg*gpsmass + xcamcg*cameramass+xservocg*servomass)/TotMass;

%z cg location calculations, nose of a/c is origin of z-axis, oriented pos.
%down
%better to have c.g. above wing a.c. for long. stability
zfusecg = 0.0;
zwingcg = 0; %mid wing design, 0 dihedral
zpropcg = 0; 
zmotorcg = 0;
zbattcg = 0;
zshockcg = 0; %?? 
zspdctrlcg = -(0.01778/2 + 0.003175/2);
zapmodemcg = -0.0127/2;
zgpscg = 0.0127/2;
zhtailcg = 0;
zcamcg = 0;
zservocg = 0;
zvtailcg = -(0.0127/2+thick + vtailb*(vtailroot+2*vtailtip)/(3*(vtailtip+vtailroot)));
zcg = -0.002044932453260; %(zfusecg*fusemass + zwingcg*wingmass + zpropcg*propmass + zmotorcg*motormass +...
%     zbattcg*battmass + zshockcg*shocksmass + zspdctrlcg*speedctrlmass + zapmodemcg*apmodemmass +...
%     zgpscg*gpsmass + zhtailcg*htailmass + zvtailcg*vtailmass + zcamcg*cameramass+...
%     zservocg*servomass)/TotMass;

%% Airfoil Data
    %Main Wing
    clcdfile = 'NACA6409clcdRe100k.txt';
    clcmfile = 'NACA6409clcmRe100k.txt';

    aclcd = dlmread(clcdfile);
    aclcm = dlmread(clcmfile);
    
    alpha1 = aclcm(:,1)*pi/180;
    cl1 = aclcm(:,2);
    cm = aclcm(:,3);

    alpha2 = aclcd(:,1)*pi/180;
    cl2 = aclcd(:,2);
    cd = aclcd(:,3);
    ltod = cl2./cd;
    powerrat = cl2.^(3/2)./cd;
    
    Cl_alpha = polyfit(alpha1(1:17), cl1(1:17), 1);
    a0 = Cl_alpha(1);
    CL0 = Cl_alpha(2);
    
    Cd_alpha = polyfit(alpha2, cd, 2);
    CD0 = Cd_alpha(3);
    CDA = Cd_alpha(2);
    CDA2 = Cd_alpha(1);
    CMAC = mean(cm);
    
    %Horizontal and Vertical Tails
    tclcdfile = 'NACA64AclcdRe100k.txt';
    tclcmfile = 'NACA64AclcmRe100k.txt';
    taclcd = dlmread(tclcdfile);
    taclcm = dlmread(tclcmfile);
    
    talpha1 = taclcm(:,1)*pi/180;
    tcl1 = taclcm(:,2);
    tcm = taclcm(:,3);

    talpha2 = taclcd(:,1)*pi/180;
    tcl2 = taclcd(:,2);
    tcd = taclcd(:,3);
    tltod = tcl2./tcd;
    tpowerrat = tcl2.^(3/2)./tcd;
    
    Cl_alpha_tail = polyfit(talpha1(1:17), tcl1(1:17), 1);
    a0_tail = Cl_alpha_tail(1);
    CL0_tail = Cl_alpha_tail(2);

    
    Cd_alpha_tail = polyfit(talpha2, tcd, 2);
    CD0_tail = Cd_alpha_tail(3);
    CDA_tail = Cd_alpha_tail(2);
    CDA2_tail = Cd_alpha_tail(1);
   
    CMAC_tail = mean(tcm);
%% Longitudinal Stability Analysis
rho = 1.225; %air density at sea level, [kg/m^3]
alpha_i_t = -2*pi/180; %incidence angle of horizontal tail [rad]
Vcruise = 16.75;
e = 1.78*(1-0.045*AR^0.68)-0.64; %From Raymer %0.9899; %oswald's efficiency factor [nd]
CLA = a0/(1+a0/(pi*e*AR));
CLA_ht = a0_tail/(1+a0_tail/(pi*e*ARht));
%Cm0
lwingcg = ((x_rootle_loc-xbar)-xcg); %vector pointing from vehicle c.g. to wing a.c.
lhtailcg = -totallength+htailroot-xbarht-xcg; %vector pointing from vehicle c.g. to horizontal tail a.c.
M0_fv = (1/2)*rho*(Vcruise)^2*(2*area)*cbar*CMAC + (1/2)*rho*(Vcruise)^2*(2*area)*CL0*lwingcg +...
    (1/2)*rho*(Vcruise)^2*(2*harea)*(CL0_tail + CLA_ht*alpha_i_t)*lhtailcg;
Cm0_fv = M0_fv/(1/2*rho*Vcruise^2*2*area*cbar);

lbart = abs(-totallength+htailroot-xbarht-(x_rootle_loc-xbar));
Vhbar = lbart*(2*harea)/(cbar*2*area);
Vh = abs(lhtailcg)*2*harea/(cbar*2*area);

hn = h_nw - (1/CLA)*(-Vhbar*CLA_ht);
h = abs(xcg-(x_rootle_loc-xle))/cbar;
Cm_alpha_fv = CLA*(h-hn);

f = fuselength/(din+2*thick); 
FF = 1+60/(f^3) + f/400; %Fuselage form factor
Re_f = (rho*Vcruise*fuselength)/(1.79*10^(-5)); %Reynolds # of fuselage
Cf = 1.328/sqrt(Re_f); %Skin friction coefficient, laminar
% Cf = 0.455/((log10(Re_f))^2.58*(1)^0.65);
Swet_f = 2*pi*(din/2+thick)*fuselength;
CD0_f = FF*Cf*Swet_f/(2*area);

alphatrim = 0.1055;
CLtrim = 0.8252;

%Tail/Fuselage Cx0, Cx_alpha, Cz0, Cz_alpha
M0 = (1/2)*rho*(Vcruise)^2*(2*harea)*(CL0_tail + CLA_ht*alpha_i_t)*lhtailcg;
Cm0 = M0/(1/2*rho*Vcruise^2*2*area*cbar);
Cm_alpha = -Vhbar*CLA_ht;
CL_q = 2*CLA_ht*Vh;
Cm_q = -2*CLA_ht*Vh*abs(lhtailcg)/cbar;

%% Lateral Stability
%Weathercock
CLA_vt = a0_tail/(1+a0_tail/(pi*e*ARvt));
lf = abs(-totallength+vtailroot-xbarvt-xcg);
Vv = varea*lf/(2*area*2*b);
Cn_beta = Vv*CLA_vt;
Cnbeta_fuse = -1.3*(fusevol/(2*area*2*b))*((din+2*thick+0.0127+2*thick)/2/(9.81*(fusemass+tailboommass)));

Cy_beta = -CLA_vt*varea/(2*area);
Cl_beta = -CLA_vt*varea*ybarvt/(2*area*2*b);
Cy_p = -CLA_vt*varea*2*ybarvt/(2*area*2*b);
Cn_p = CLA_vt*Vv*2*ybarvt/(2*b);
Cl_p = -CLA_vt*varea*ybarvt^2*2/(2*area*(2*b)^2) - CLA_ht*2*harea*ybarht^2*2/(2*area*(2*b)^2);
Cy_r = CLA_vt*varea*2*lf/(2*area*2*b);
Cl_r = CLA_vt*varea*ybarvt*2*lf/(2*area*(2*b)^2);
Cn_r = -CLA_vt*Vv*2*lf/(2*b);


%% Control Sizing


%Horizontal Tail:
%Use strip method to calculate roll moment from elevator deflection
    Kfmat = [0.72; 0.7; 0.675; 0.65; 0.61; 0.575; 0.55; 0.525;0.5]; %for elev. deflection of 30 deg.
    crat = [0.10; 0.15; 0.2; 0.25; 0.3; 0.4; 0.5; 0.6; 0.7];
    dclddfmat = [1.5; 2.5; 3.2; 3.6; 4; 4.4; 4.75; 5; 5.25; 5.5; 6.0; 6.5]; %[1/rad]
    crat2 = [0.05; 0.1; 0.15; 0.2; 0.25; 0.3; 0.35; 0.4; 0.45; 0.5; 0.6; 0.7];
    elevc = 0.0175;
    sum = 0;
    summ = 0;
    sumlift = 0;
for i = 0:0.001:htailb
    chord = ((htailroot-htailroot*htaper)/(0-htailb))*i+htailroot;
    elevrat = elevc/chord;
    Kf = interp1(crat,Kfmat,elevrat);
    dclddf = interp1(crat2,dclddfmat,elevrat);
    Si = 0.001*chord;
    sum = sum + 0.85*Kf*dclddf*i*Si;
    sumlift = sumlift + 0.85*Kf*dclddf*Si;
    hxi = (-totallength+htailroot-i*(htailroot-htailtip)/htailb-0.25*chord)-xcg;
    summ = summ + 0.85*Kf*dclddf*hxi*Si;
end
Cl_delev = 2*sum/(2*area*2*b);
Clift_delev = 2*sumlift/(2*area); 
CD_delev = (Clift_delev/2)^2/(pi*e*ARht); %for each half
CN_delev = CD_delev*ybarht/(2*b);
Cm_delev = 2*summ/(2*area*cbar);
sumy = 0;
suml = 0;
sumn = 0;
rudc = 0.0195;
for i = 0:0.001:vtailb
    chord = ((vtailroot-vtailroot*vtaper)/(0-vtailb))*i+vtailroot;
    rudrat = rudc/chord;
    Kf = interp1(crat,Kfmat,rudrat);
    dclddf = interp1(crat2,dclddfmat,rudrat);
    Si = 0.001*chord;
    vxi = (-totallength+vtailroot-i*(vtailroot-vtailtip)/vtailb-0.25*chord)-xcg;
    sumy = sumy + 0.85*Kf*dclddf*Si;
    suml = suml + 0.85*Kf*dclddf*i*Si;
    sumn = sumn + 0.85*Kf*dclddf*vxi*Si;
end
Cy_drud = sumy/(2*area);
Cl_drud = suml/(2*area*2*b);
Cn_drud = sumn/(2*area*2*b);

