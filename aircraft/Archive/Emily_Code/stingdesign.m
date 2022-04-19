l1 = .03625 + 0.1087;  %length of sting inserted into aircraft [m]
l2 = 0.10; %.28; % length from tail end of aircraft to load cell
rho_al = 2700; %kg/m^3

msting = rho_al*(l1+l2)*pi*(0.0127^2)/4;
Wsting = msting*9.81;
ms = 0.187; %[kg]
W = 9.81*ms; %weight [N]

mnc = 0.0106; %mass of nosecone [kg]
xcg_nc = -0.014447; ycgnc = 0.0; zcgnc = 0.0; %cg of nosecone [m], origin is at rear of part.
lnc = 0.035; %total length of nosecone is 35 mm
xcgnc = lnc + xcg_nc;

mcb = 0.051; %mass of cylinder body
%measured from front of cylinder body, pos x towards rear, pos-z up.[m]
xcgcb = 0.055963; ycgcb = 0.0; zcgcb = -0.000086; %cg of cylinder body, [m]
dcb = 0.042; %outer diameter of cylinder body [m]
lcb = 0.1225; %length of cylinder body
hcbz = 0.001505; %z height of corner of zero angle semi-circular hole on fuselage
hcbx = 0.016191; %x length from front of cylinder of front corner of 0 angle semi-circ hole
xcgcb = lnc + xcgcb;

mrc = 0.0106; %mass of rear cone [kg]
%measured from front of cone, x towards the rear, y to the right, z up
xcgrc = 0.012621; ycgrc = 0.0; zcgrc = -0.000002;
lrc = 0.04125; %length of rear cone [m]
xcgrc = lnc + lcb + xcgrc;

mt = 0.0167; %mass of tail boom [kg] %x starts at 
xcgt = 0.240465; ycgt = 0.000033; zcgt = 0.009577; %measured from front of tail boom which is at x = 163.75;
%z positive up, x positive towards rear
xcgt = lnc + lcb + lrc + xcgt-0.16375; 
dt = 0.014025; %outer diameter of tail boom
lt = 0.10875; %length of tail boom [m]

totallength = lnc+lcb+lrc+lt;

mell = 0.0012; %mass of left elevator [kg]
%x starts at 259.375 mm, y starts at -7.013 mm. x pointed towards rear, pos
%z up
xcgell = 0.264126; ycgell = -0.044632; zcgell = 0.00305; 
lell = 0.01311; %chord length of elevator [m]
xcgell = totallength - lell + (xcgell - 0.259375);

melr = 0.0012; 
%x starts at 259.375 mm, y starts at -7.013 mm. x pointed towards rear, pos
%z up
xcgelr = 0.264126; ycgelr = 0.030607; zcgelr = 0.00305;
lelr = 0.01311;
xcgelr = totallength - lelr + (xcgelr - 0.259375);
ycgelr = ycgelr + dt;

mrud = 0.0013;
%x starts at 257.875 mm, Z starts at 7.013 mm 
xcgrud = 0.263076; ycgrud = 0.0; zcgrud = 0.041011; 
lrud = 0.0146460; %chord length of rudder [m]
xcgrud = totallength - lrud + (xcgrud-0.257875);

mrw = 0.0455;
%x = 0, y = 0, z =0 at leading edge wingtip
%x points forward, z up
xcgrw = -0.017219; ycgrw = 0.174568; zcgrw = 0.002911;
xcgrw = lnc + hcbx - 0.01 +0.015 - xcgrw; %measured from nose
ycgrw = dcb/2 + 0.3-ycgrw; %0.3 is the span of the half-wing
zcgrw = hcbz-0.0015 + zcgrw;

mlw = 0.0455;
xcglw = xcgrw; ycglw = -ycgrw; zcglw = zcgrw ;


m = mnc + mcb + mrc + mt + mell + melr + mrud + mrw + mlw; %total mass

%Total CG 
xcg_struct = -(mnc*xcgnc + mcb*xcgcb + mrc*xcgrc + mt*xcgt + mell*xcgell + melr*xcgelr + mrud*xcgrud + mrw*xcgrw + mlw*xcglw)/m;
ycg_struct = (mnc*ycgnc + mcb*ycgcb + mrc*ycgrc + mt*ycgt + mell*ycgell + melr*ycgelr + mrud*ycgrud + mrw*ycgrw + mlw*ycglw)/m;
zcg_struct = -(mnc*zcgnc + mcb*zcgcb + mrc*zcgrc + mt*zcgt + mell*zcgell + melr*zcgelr + mrud*zcgrud + mrw*zcgrw + mlw*zcglw)/m; 

%Full scale
xcg_fs = xcg_struct*4/3; ycg_fs = ycg_struct*4/3; zcg_fs = zcg_struct*4/3;

%Tail/Fuselage Structure CG, Full Scale, standard body axes
mtf = mnc + mcb + mrc + mt + mell + melr + mrud;
xcg_tf = -(4/3)*(mnc*xcgnc + mcb*xcgcb + mrc*xcgrc + mt*xcgt + mell*xcgell + melr*xcgelr + mrud*xcgrud)/mtf;
ycg_tf = (4/3)*(mnc*ycgnc + mcb*ycgcb + mrc*ycgrc + mt*ycgt + mell*ycgell + melr*ycgelr + mrud*ycgrud)/mtf;
zcg_tf = -(4/3)*(mnc*zcgnc + mcb*zcgcb + mrc*zcgrc + mt*zcgt + mell*zcgell + melr*zcgelr + mrud*zcgrud)/mtf;

%Right Wing CG, Full Scale, in wing body axes (origin at leading edge of
%root chord)
b = 0.30; %span of half-wing [m]
xcg_rw = (4/3)*(-0.017219 - 0.015);
ycg_rw = (4/3)*(b - 0.174568);
zcg_rw = (4/3)*(-0.002911);

%Left Wing CG, Full Scale, in wing body axes (origin at leading edge of
%root chord)
xcg_lw = xcg_rw;
ycg_lw = -ycg_rw;
zcg_lw = zcg_rw;

totallength = 0.41*0.75;
Torque = W*(totallength+xcg_struct+l2) + Wsting*(l1+l2)/2; % [N-m]
Torque_E = Torque*39.3700787*0.224808943; %convert from [N-m] to [lbf-in]  

CL0 = 0.4589;
CLA = 3.3234;
CL0_tail = -0.0154;
CLA_ht = 3.7754;
alpha_i_t = -0.0349;
harea = 0.0045;
area = 0.0320;
modarea = 2*area*0.75^2; %area of wind tunnel model which is 75% scale of actual aircraft
V = 17.0;
alpha = 6.0*pi/180;
CLift = CL0 + CL0_tail*harea/area + CLA*alpha + CLA_ht*(alpha+alpha_i_t)*harea/area;
Lift = 0.5*1.225*V^2*modarea*CLift;
T_x = Lift*(totallength+xcg_struct+l2)*39.3700787*0.224808943;
PitchM = -0.035088816211573*39.3700787*0.224808943;
Torque_E + T_x + PitchM