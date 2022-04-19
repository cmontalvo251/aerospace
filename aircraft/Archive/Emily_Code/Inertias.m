clear all
close all
clc

%Nose Cone Inertia, 75% scale. %Inertia about body cg
%Coord. Sys: x pointed rear, y starboard, z up
Ixx_nc = 2.383; Iyy_nc = 2.049; Izz_nc = 2.049; %[kg-mm^2]
Ixy_nc = 0.0; Ixz_nc = 0.0; Iyz_nc = 0.0;

%Cylinder Body Inertia, 75% scale, %Inertia about body cg
%Coord. Sys: x pointed rear, y starboard, z up
Ixx_cb = 13.808; Iyy_cb = 67.409; Izz_cb = 67.255; %[kg-mm^2]
Ixy_cb = 0.0; Ixz_cb = -0.095; Iyz_cb = 0.0;

%Tail Boom/HT/VT combined inertia, 75% scale, inertia about body cg
%Coord. Sys: x pointed rear, y starboard, z up
Ixx_t = 18.231; Iyy_t = 12.727; Izz_t = 19.614; %[kg-mm^2]
Ixy_t = 0.0; Ixz_t = -1.044; Iyz_t = 0.0;

%Rear cone Inertia, 75% scale, inertia about body cg [kg-mm^2]
%Coord. Sys: x pointed rear, y starboard, z up
Ixx_rc = 1.684; Iyy_rc = 1.836; Izz_rc = 1.836;
Ixy_rc = 0.00; Ixz_rc = 0.001; Iyz_rc = 0.00;

%Left Elevator Inertia, 75% scale, inertia about body cg [kg-mm^2]
%Coord. Sys: x pointed rear, y starboard, z up
Ixx_ell = 0.589; Iyy_ell = 0.012; Izz_ell = 0.60;
Ixy_ell = -0.001; Ixz_ell = 0.0; Iyz_ell = -0.003;

%Right Elevator Inertia, 75% scale, inertia about body cg [kg-mm^2]
%Coord. Sys: x pointed rear, y starboard, z up
Ixx_elr = 0.589; Iyy_elr = 0.012; Izz_elr = 0.60;
Ixy_elr = 0.001; Ixz_elr = 0.00; Iyz_elr = 0.003;

%Rudder Inertia, 75% scale, inertia about body cg [kg-mm^2]
%Coord. Sys: x pointed rear, y starboard, z up
Ixx_rud = 0.545; Iyy_rud = 0.561; Izz_rud = 0.017; 
Ixy_rud = 0.0; Ixz_rud = 0.0; Iyz_rud = 0.0; 

%Right Wing Inertia, 75% scale, inertia about body cg
%x points forward, z up
Ixx_rw = 320.030; Iyy_rw = 9.921; Izz_rw = 329.696;
Ixy_rw = -2.579; Ixz_rw = 0.063; Iyz_rw = -1.488;

%Left Wing Inertia, 75% scale, inertia about body cg
%x points forward, z up
Ixx_lw = 320.030; Iyy_lw = 9.921; Izz_lw = 329.696;
Ixy_lw = 2.579; Ixz_lw = 0.063; Iyz_lw = 1.488;

%Form Inertia Matrices, convert to [kg-m^2], convert to full scale

I_nc = (4/3)^5*(1/(1000^2))*[Ixx_nc Ixy_nc Ixz_nc; Ixy_nc Iyy_nc Iyz_nc; Ixz_nc Iyz_nc Izz_nc];
I_cb = (4/3)^5*(1/(1000^2))*[Ixx_cb Ixy_cb Ixz_cb; Ixy_cb Iyy_cb Iyz_cb; Ixz_cb Iyz_cb Izz_cb];
I_t = (4/3)^5*(1/(1000^2))*[Ixx_t Ixy_t Ixz_t; Ixy_t Iyy_t Iyz_t; Ixz_t Iyz_t Izz_t];
I_rc = (4/3)^5*(1/(1000^2))*[Ixx_rc Ixy_rc Ixz_rc; Ixy_rc Iyy_rc Iyz_rc; Ixz_rc Iyz_rc Izz_rc];
I_ell = (4/3)^5*(1/(1000^2))*[Ixx_ell Ixy_ell Ixz_ell; Ixy_ell Iyy_ell Iyz_ell; Ixz_ell Iyz_ell Izz_ell];
I_elr = (4/3)^5*(1/(1000^2))*[Ixx_elr Ixy_elr Ixz_elr; Ixy_elr Iyy_elr Iyz_elr; Ixz_elr Iyz_elr Izz_elr];
I_rud = (4/3)^5*(1/(1000^2))*[Ixx_rud Ixy_rud Ixz_rud; Ixy_rud Iyy_rud Iyz_rud; Ixz_rud Iyz_rud Izz_rud];
I_rw = (4/3)^5*(1/(1000^2))*[Ixx_rw Ixy_rw Ixz_rw; Ixy_rw Iyy_rw Iyz_rw; Ixz_rw Iyz_rw Izz_rw];
I_lw = (4/3)^5*(1/(1000^2))*[Ixx_lw Ixy_lw Ixz_lw; Ixy_lw Iyy_lw Iyz_lw; Ixz_lw Iyz_lw Izz_lw];

%Transform matrices to standard body frame (x point forwards, y to
%starboard, z down

T_nc = [cos(pi) 0 -sin(pi); 0 1 0; sin(pi) 0 cos(pi)];
T_cb = [cos(pi) 0 -sin(pi); 0 1 0; sin(pi) 0 cos(pi)];
T_t = [cos(pi) 0 -sin(pi); 0 1 0; sin(pi) 0 cos(pi)];
T_rc = [cos(pi) 0 -sin(pi); 0 1 0; sin(pi) 0 cos(pi)];
T_ell = [cos(pi) 0 -sin(pi); 0 1 0; sin(pi) 0 cos(pi)];
T_elr = [cos(pi) 0 -sin(pi); 0 1 0; sin(pi) 0 cos(pi)];
T_rud = [cos(pi) 0 -sin(pi); 0 1 0; sin(pi) 0 cos(pi)];
T_rw = [1 0 0; 0 cos(pi) sin(pi); 0 sin(pi) cos(pi)];
T_lw = [1 0 0; 0 cos(pi) sin(pi); 0 sin(pi) cos(pi)];

Inc = T_nc*I_nc*T_nc'; Icb = T_cb*I_cb*T_cb'; Irc = T_rc*I_rc*T_rc';
It = T_t*I_t*T_t'; Iell = T_ell*I_ell*T_ell'; Ielr = T_elr*I_elr*T_elr'; 
Irud = T_rud*I_rud*T_rud'; Irw = T_rw*I_rw*T_rw'; Ilw = T_lw*I_lw*T_lw';

%Inertias of Other Parts (including converting from g-in^2 to kg-m^2)
Ixxprop = 57.475*(1/1000)*0.0254^2; Iyyprop = 28.7375*(1/1000)*0.0254^2; Izzprop = 28.7375*(1/1000)*0.0254^2;
Ixxservo = 0.146069336*(1/1000)*0.0254^2; Iyyservo = 0.236889648*(1/1000)*0.0254^2; Izzservo = 0.164990234*(1/1000)*0.0254^2;
Ixxspdctrl = 0.458447266*(1/1000)*0.0254^2; Iyyspdctrl = 0.976953125*(1/1000)*0.0254^2; Izzspdctrl = 1.419384766*(1/1000)*0.0254^2;
Ixxgps = 0.303385417*(1/1000)*0.0254^2; Iyygps = 1.032552083*(1/1000)*0.0254^2; Izzgps = 1.169270833*(1/1000)*0.0254^2;
Ixxapmodem = 1.161458333*(1/1000)*0.0254^2; Iyyapmodem = 1.161458333*(1/1000)*0.0254^2; Izzapmodem = 1.858333333*(1/1000)*0.0254^2;
Ixxmotor = 5.5625*(1/1000)*0.0254^2; Iyymotor = 5.620442708*(1/1000)*0.0254^2; Izzmotor = 5.620442708*(1/1000)*0.0254^2;
Ixxcam = 0.073958333*(1/1000)*0.0254^2; Iyycam = 0.073958333*(1/1000)*0.0254^2; Izzcam = 0.073958333*(1/1000)*0.0254^2;
Ixxshock = 81.046875*(1/1000)*0.0254^2; Iyyshock = 81.046875*(1/1000)*0.0254^2; Izzshock = 25.59375*(1/1000)*0.0254^2;

 

%Individual CGs
%Structure
xcg_nc = -0.0274040; ycg_nc = 0.0; zcg_nc = 0.0;
xcg_cb = -0.121284; ycg_cb = 0.0; zcg_cb = 0.0001146666666666667;
xcg_rc = -0.226828; ycg_rc = 0.0; zcg_rc = 0.000002666666666666666;
xcg_t = -0.367286666666667; ycg_t =  0.0; zcg_t = -0.012769333333333;
xcg_ell = -0.398854666666667; ycg_ell = -0.059509333333333; zcg_ell = -0.004066666666667;
xcg_elr = -0.398854666666667; ycg_elr = 0.059509333333333; zcg_elr = -0.004066666666667;
xcg_rud = -0.397406666666667; ycg_rud = 0.0; zcg_rud = -0.054681333333333;
xcg_rw = -0.097880;
ycg_rw = 0.195242666666667;
zcg_rw = -0.0038880;
xcg_lw = xcg_rw; ycg_lw = -ycg_rw; zcg_lw = zcg_rw;
grw = 0.0;
glw = 0.0;
ycg_lw = ycg_lw*cos(glw); zcg_lw = zcg_lw + zcg_lw*sin(glw);
ycg_rw = ycg_rw*cos(grw); zcg_rw = zcg_rw - zcg_rw*sin(grw);

%Inner Parts
xpropcg = 0.005;
motorlength = 0.025; %length of motor [m]
xmotorcg = -(0.01+0.5*motorlength);
battlength = 0.06604; %length of battery [m]
xbattcg = -(0.01+motorlength+0.003+0.5*battlength); %assuming battery is directly behind motor

x_rootle_loc = -0.055; root = 0.10; 
xshockcg = x_rootle_loc-0.25*root;
spdctrllength = 0.034925; %length of speed controller [m]
xspdctrlcg = -(0.01+motorlength+0.003+0.5*spdctrllength); %speed controller positioned on top of battery
%Place GPS and AP/Modem on top of each other, behind shocks [Replace with Gina mote?]
xwingcg = (x_rootle_loc-0.5*root);
xapmodemcg = xwingcg - 0.0127;
xgpscg = xwingcg - 0.0428625/2;
xcamcg = xgpscg - 0.0428625/2 - 0.001 - 0.0127/2;
xservocg = xcamcg - 0.0127/2 - 0.005 - 0.0206375/2; %0.0206375 is the length of the servo in [m]

ypropcg = 0.0; ymotorcg = 0.0; ybattcg = 0.0; 
yspdctrlcg = 0.0; yapmodemcg = 0.0; ygpscg = 0.0; ycamcg = 0.0; 
yservorcg = 0.0111125/2+0.0025; yservolcg = -yservorcg;
yshockrcg = 0.025; yshocklcg = -yshockrcg; 

zpropcg = 0.0; zmotorcg = 0.0; zbattcg = 0.0; zshockcg = 0.0; %?? 
zspdctrlcg = -(0.01778/2 + 0.003175/2); zapmodemcg = -0.0127/2;
zgpscg = 0.0127/2; zcamcg = 0.0; zservocg = 0.0;

%Masses
%Structure
mnc = 0.0106*(4/3)^3; %mass of nosecone [kg]
mcb = 0.051*(4/3)^3; mrc = 0.0106*(4/3)^3; mt = 0.0167*(4/3)^3;
mell = 0.0012*(4/3)^3; melr = 0.0012*(4/3)^3; mrud = 0.0013*(4/3)^3;
mrw = 0.0455*(4/3)^3; mlw = 0.0455*(4/3)^3;
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

battl = battlength; battw = 0.03429; batth = 0.01778;
Ixxbatt = battmass*(battw^2 + batth^2)/12;
Iyybatt = battmass*(battl^2 + batth^2)/12;
Izzbatt = battmass*(battl^2 + battw^2)/12;

%Full Body CG
Mtot = mnc+mcb+mrc+mt+mell+melr+mrud+mrw+mlw+battmass+propmass+servomass+speedctrlmass+...
    gpsmass+cameramass+apmodemmass+motormass+shocksmass; 
xcg = (mnc*xcg_nc+mcb*xcg_cb+mrc*xcg_rc+mt*xcg_t+mell*xcg_ell+melr*xcg_elr+mrud*xcg_rud+mrw*xcg_rw+mlw*xcg_lw+battmass*xbattcg+...
    propmass*xpropcg+servomass*xservocg+speedctrlmass*xspdctrlcg+gpsmass*xgpscg+cameramass*xcamcg+apmodemmass*xapmodemcg+...
    motormass*xmotorcg+shocksmass*xshockcg)/Mtot; 
ycg = (mnc*ycg_nc+mcb*ycg_cb+mrc*ycg_rc+mt*ycg_t+mell*ycg_ell+melr*ycg_elr+mrud*ycg_rud+mrw*ycg_rw+mlw*ycg_lw+battmass*ybattcg+...
    propmass*ypropcg+(servomass/2)*yservolcg+(servomass/2)*yservorcg+speedctrlmass*yspdctrlcg+gpsmass*ygpscg+cameramass*ycamcg+apmodemmass*yapmodemcg+...
    motormass*ymotorcg+(shocksmass/2)*yshockrcg+(shocksmass/2)*yshocklcg)/Mtot; 
zcg = (mnc*zcg_nc+mcb*zcg_cb+mrc*zcg_rc+mt*zcg_t+mell*zcg_ell+melr*zcg_elr+mrud*zcg_rud+mrw*zcg_rw+mlw*zcg_lw+battmass*zbattcg+...
    propmass*zpropcg+servomass*zservocg+speedctrlmass*zspdctrlcg+gpsmass*zgpscg+cameramass*zcamcg+apmodemmass*zapmodemcg+...
    motormass*zmotorcg+shocksmass*zshockcg)/Mtot;

%Tail/Fuselage 
Mtot_tf = mnc+mcb+mrc+mt+mell+melr+mrud+battmass+propmass+servomass+speedctrlmass+...
    gpsmass+cameramass+apmodemmass+motormass+shocksmass; 
xcg_tf = (mnc*xcg_nc+mcb*xcg_cb+mrc*xcg_rc+mt*xcg_t+mell*xcg_ell+melr*xcg_elr+mrud*xcg_rud+battmass*xbattcg+...
    propmass*xpropcg+servomass*xservocg+speedctrlmass*xspdctrlcg+gpsmass*xgpscg+cameramass*xcamcg+apmodemmass*xapmodemcg+...
    motormass*xmotorcg+shocksmass*xshockcg)/Mtot_tf; 
ycg_tf = (mnc*ycg_nc+mcb*ycg_cb+mrc*ycg_rc+mt*ycg_t+mell*ycg_ell+melr*ycg_elr+mrud*ycg_rud+battmass*ybattcg+...
    propmass*ypropcg+(servomass/2)*yservolcg+(servomass/2)*yservorcg+speedctrlmass*yspdctrlcg+gpsmass*ygpscg+cameramass*ycamcg+apmodemmass*yapmodemcg+...
    motormass*ymotorcg+(shocksmass/2)*yshockrcg+(shocksmass/2)*yshocklcg)/Mtot_tf; 
zcg_tf = (mnc*zcg_nc+mcb*zcg_cb+mrc*zcg_rc+mt*zcg_t+mell*zcg_ell+melr*zcg_elr+mrud*zcg_rud+battmass*zbattcg+...
    propmass*zpropcg+servomass*zservocg+speedctrlmass*zspdctrlcg+gpsmass*zgpscg+cameramass*zcamcg+apmodemmass*zapmodemcg+...
    motormass*zmotorcg+shocksmass*zshockcg)/Mtot_tf; 

%Vectors from tf CG to each body CG for Tail/Fuse Inertia
rcgx_nc = xcg_nc - xcg_tf; rcgy_nc = ycg_nc - ycg_tf; rcgz_nc = zcg_nc - zcg_tf;    
rcgx_cb = xcg_cb - xcg_tf; rcgy_cb = ycg_cb - ycg_tf; rcgz_cb = zcg_cb - zcg_tf;
rcgx_rc = xcg_rc - xcg_tf; rcgy_rc = ycg_rc - ycg_tf; rcgz_rc = zcg_rc - zcg_tf;
rcgx_t = xcg_t - xcg_tf; rcgy_t = ycg_t - ycg_tf; rcgz_t = zcg_t - zcg_tf;
rcgx_ell = xcg_ell - xcg_tf; rcgy_ell = ycg_ell - ycg_tf; rcgz_ell = zcg_ell - zcg_tf;
rcgx_elr = xcg_elr - xcg_tf; rcgy_elr = ycg_elr - ycg_tf; rcgz_elr = zcg_elr - zcg_tf;
rcgx_rud = xcg_rud - xcg_tf; rcgy_rud = ycg_rud - ycg_tf; rcgz_rud = zcg_rud - zcg_tf;
rcgx_batt = xbattcg - xcg_tf; rcgy_batt = ybattcg - ycg_tf; rcgz_batt = zbattcg - zcg_tf;
rcgx_prop = xpropcg - xcg_tf; rcgy_prop = ypropcg - ycg_tf; rcgz_prop = zpropcg - zcg_tf;
rcgx_servo = xservocg - xcg_tf; rcgy_servol = yservolcg - ycg_tf; rcgy_servor = yservorcg - ycg_tf; rcgz_servo = zservocg - zcg_tf;
rcgx_spdctrl = xspdctrlcg - xcg_tf; rcgy_spdctrl = yspdctrlcg - ycg_tf; rcgz_spdctrl = zspdctrlcg - zcg_tf;
rcgx_gps = xgpscg - xcg_tf; rcgy_gps = ygpscg - ycg_tf; rcgz_gps = zgpscg - zcg_tf;
rcgx_cam = xcamcg - xcg_tf; rcgy_cam = ycamcg - ycg_tf; rcgz_cam = zcamcg - zcg_tf;
rcgx_apmod = xapmodemcg - xcg_tf; rcgy_apmod = yapmodemcg - ycg_tf; rcgz_apmod = zapmodemcg - zcg_tf;
rcgx_motr = xmotorcg - xcg_tf; rcgy_motr = ymotorcg - ycg_tf; rcgz_motr = zmotorcg - zcg_tf;
rcgx_shock = xshockcg - xcg_tf; rcgy_shockl = yshocklcg - ycg_tf; rcgy_shockr = yshockrcg - ycg_tf; rcgz_shock = zshockcg - zcg_tf;

Ixx_tf = Inc(1,1)+mnc*(rcgy_nc^2+rcgz_nc^2)+Icb(1,1)+mcb*(rcgy_cb^2+rcgz_cb^2)+Irc(1,1)+mrc*(rcgy_rc^2+rcgz_rc^2)+...
    It(1,1)+mt*(rcgy_t^2+rcgz_t^2)+Iell(1,1)+mell*(rcgy_ell^2+rcgz_ell^2)+Ielr(1,1)+melr*(rcgy_elr^2+rcgz_elr^2)+...
    Irud(1,1)+mrud*(rcgy_rud^2+rcgz_rud^2)+Ixxbatt+battmass*(rcgy_batt^2+rcgz_batt^2)+Ixxprop+propmass*(rcgy_prop^2+rcgz_prop^2)+...
    Ixxservo+(servomass/2)*(rcgy_servol^2+rcgz_servo^2)+(servomass/2)*(rcgy_servor^2+rcgz_servo^2)+Ixxspdctrl+...
    speedctrlmass*(rcgy_spdctrl^2+rcgz_spdctrl^2)+Ixxgps+gpsmass*(rcgy_gps^2+rcgz_gps^2)+Ixxcam+cameramass*(rcgy_cam^2+rcgz_cam^2)+...
    Ixxapmodem+apmodemmass*(rcgy_apmod^2+rcgz_apmod^2)+Ixxmotor+motormass*(rcgy_motr^2+rcgz_motr^2)+...
    Ixxshock+(shocksmass/2)*(rcgy_shockl^2+rcgz_shock^2)+(shocksmass/2)*(rcgy_shockr^2+rcgz_shock^2);

Iyy_tf = Inc(2,2)+mnc*(rcgx_nc^2+rcgz_nc^2)+Icb(2,2)+mcb*(rcgx_cb^2+rcgz_cb^2)+Irc(2,2)+mrc*(rcgx_rc^2+rcgz_rc^2)+...
    It(2,2)+mt*(rcgx_t^2+rcgz_t^2)+Iell(2,2)+mell*(rcgx_ell^2+rcgz_ell^2)+Ielr(2,2)+melr*(rcgx_elr^2+rcgz_elr^2)+...
    Irud(2,2)+mrud*(rcgx_rud^2+rcgz_rud^2)+Iyybatt+battmass*(rcgx_batt^2+rcgz_batt^2)+Iyyprop+propmass*(rcgx_prop^2+rcgz_prop^2)+...
    Iyyservo+(servomass)*(rcgx_servo^2+rcgz_servo^2)+Iyyspdctrl+...
    speedctrlmass*(rcgx_spdctrl^2+rcgz_spdctrl^2)+Iyygps+gpsmass*(rcgx_gps^2+rcgz_gps^2)+Iyycam+cameramass*(rcgx_cam^2+rcgz_cam^2)+...
    Iyyapmodem+apmodemmass*(rcgx_apmod^2+rcgz_apmod^2)+Iyymotor+motormass*(rcgx_motr^2+rcgz_motr^2)+...
    Ixxshock+(shocksmass)*(rcgx_shock^2+rcgz_shock^2);

Izz_tf = Inc(3,3)+mnc*(rcgy_nc^2+rcgx_nc^2)+Icb(3,3)+mcb*(rcgy_cb^2+rcgx_cb^2)+Irc(3,3)+mrc*(rcgy_rc^2+rcgx_rc^2)+...
    It(3,3)+mt*(rcgy_t^2+rcgx_t^2)+Iell(3,3)+mell*(rcgy_ell^2+rcgx_ell^2)+Ielr(3,3)+melr*(rcgy_elr^2+rcgx_elr^2)+...
    Irud(3,3)+mrud*(rcgy_rud^2+rcgx_rud^2)+Izzbatt+battmass*(rcgy_batt^2+rcgx_batt^2)+Izzprop+propmass*(rcgy_prop^2+rcgx_prop^2)+...
    Izzservo+(servomass/2)*(rcgy_servol^2+rcgx_servo^2)+(servomass/2)*(rcgy_servor^2+rcgx_servo^2)+Izzspdctrl+...
    speedctrlmass*(rcgy_spdctrl^2+rcgx_spdctrl^2)+Izzgps+gpsmass*(rcgy_gps^2+rcgx_gps^2)+Izzcam+cameramass*(rcgy_cam^2+rcgx_cam^2)+...
    Izzapmodem+apmodemmass*(rcgy_apmod^2+rcgx_apmod^2)+Izzmotor+motormass*(rcgy_motr^2+rcgx_motr^2)+...
    Izzshock+(shocksmass/2)*(rcgy_shockl^2+rcgx_shock^2)+(shocksmass/2)*(rcgy_shockr^2+rcgx_shock^2);

Ixy_tf = Inc(1,2)-mnc*(rcgx_nc*rcgy_nc)+Icb(1,2)-mcb*(rcgx_cb*rcgy_cb)+Irc(1,2)-mrc*(rcgy_rc*rcgx_rc)+...
    It(1,2)-mt*(rcgy_t*rcgx_t)+Iell(1,2)-mell*(rcgy_ell*rcgx_ell)+Ielr(1,2)-melr*(rcgy_elr*rcgx_elr)+...
    Irud(1,2)-mrud*(rcgy_rud*rcgx_rud)-battmass*(rcgy_batt*rcgx_batt)-propmass*(rcgy_prop*rcgx_prop)+...
    -(servomass/2)*(rcgy_servol*rcgx_servo)-(servomass/2)*(rcgy_servor*rcgx_servo)+...
    -speedctrlmass*(rcgy_spdctrl*rcgx_spdctrl)-gpsmass*(rcgy_gps*rcgx_gps)-cameramass*(rcgy_cam*rcgx_cam)+...
    -apmodemmass*(rcgy_apmod*rcgx_apmod)-motormass*(rcgy_motr*rcgx_motr)+...
    -(shocksmass/2)*(rcgy_shockl*rcgx_shock)-(shocksmass/2)*(rcgy_shockr*rcgx_shock);

Iyz_tf = Inc(2,3)-mnc*(rcgz_nc*rcgy_nc)+Icb(2,3)-mcb*(rcgz_cb*rcgy_cb)+Irc(2,3)-mrc*(rcgz_rc*rcgx_rc)+...
    It(2,3)-mt*(rcgy_t*rcgz_t)+Iell(2,3)-mell*(rcgy_ell*rcgz_ell)+Ielr(2,3)-melr*(rcgy_elr*rcgz_elr)+...
    Irud(2,3)-mrud*(rcgy_rud*rcgz_rud)-battmass*(rcgy_batt*rcgz_batt)-propmass*(rcgy_prop*rcgz_prop)+...
    -(servomass/2)*(rcgy_servol*rcgz_servo)-(servomass/2)*(rcgy_servor*rcgz_servo)+...
    -speedctrlmass*(rcgy_spdctrl*rcgz_spdctrl)-gpsmass*(rcgy_gps*rcgz_gps)-cameramass*(rcgy_cam*rcgz_cam)+...
    -apmodemmass*(rcgy_apmod*rcgz_apmod)-motormass*(rcgy_motr*rcgz_motr)+...
    -(shocksmass/2)*(rcgy_shockl*rcgz_shock)-(shocksmass/2)*(rcgy_shockr*rcgz_shock);

Ixz_tf = Inc(1,3)-mnc*(rcgx_nc*rcgz_nc)+Icb(1,3)-mcb*(rcgx_cb*rcgz_cb)+Irc(1,3)-mrc*(rcgx_rc*rcgz_rc)+...
    It(1,3)-mt*(rcgx_t*rcgz_t)+Iell(1,3)-mell*(rcgx_ell*rcgz_ell)+Ielr(1,3)-melr*(rcgx_elr*rcgz_elr)+...
    Irud(1,3)-mrud*(rcgx_rud*rcgz_rud)-battmass*(rcgx_batt*rcgz_batt)-propmass*(rcgx_prop*rcgz_prop)+...
    -(servomass)*(rcgx_servo*rcgz_servo)+...
    -speedctrlmass*(rcgx_spdctrl*rcgz_spdctrl)-gpsmass*(rcgx_gps*rcgz_gps)-cameramass*(rcgx_cam*rcgz_cam)+...
    -apmodemmass*(rcgx_apmod*rcgz_apmod)-motormass*(rcgx_motr*rcgz_motr)+...
    -(shocksmass)*(rcgx_shock*rcgz_shock);

%Full Vehicle Inertias
%Vectors from CG to each body CG
rcgx_nc = xcg_nc - xcg; rcgy_nc = ycg_nc - ycg; rcgz_nc = zcg_nc - zcg;    
rcgx_cb = xcg_cb - xcg; rcgy_cb = ycg_cb - ycg; rcgz_cb = zcg_cb - zcg;
rcgx_rc = xcg_rc - xcg; rcgy_rc = ycg_rc - ycg; rcgz_rc = zcg_rc - zcg;
rcgx_t = xcg_t - xcg; rcgy_t = ycg_t - ycg; rcgz_t = zcg_t - zcg;
rcgx_ell = xcg_ell - xcg; rcgy_ell = ycg_ell - ycg; rcgz_ell = zcg_ell - zcg;
rcgx_elr = xcg_elr - xcg; rcgy_elr = ycg_elr - ycg; rcgz_elr = zcg_elr - zcg;
rcgx_rud = xcg_rud - xcg; rcgy_rud = ycg_rud - ycg; rcgz_rud = zcg_rud - zcg;
rcgx_batt = xbattcg - xcg; rcgy_batt = ybattcg - ycg; rcgz_batt = zbattcg - zcg;
rcgx_prop = xpropcg - xcg; rcgy_prop = ypropcg - ycg; rcgz_prop = zpropcg - zcg;
rcgx_servo = xservocg - xcg; rcgy_servol = yservolcg - ycg; rcgy_servor = yservorcg - ycg; rcgz_servo = zservocg - zcg;
rcgx_spdctrl = xspdctrlcg - xcg; rcgy_spdctrl = yspdctrlcg - ycg; rcgz_spdctrl = zspdctrlcg - zcg;
rcgx_gps = xgpscg - xcg; rcgy_gps = ygpscg - ycg; rcgz_gps = zgpscg - zcg;
rcgx_cam = xcamcg - xcg; rcgy_cam = ycamcg - ycg; rcgz_cam = zcamcg - zcg;
rcgx_apmod = xapmodemcg - xcg; rcgy_apmod = yapmodemcg - ycg; rcgz_apmod = zapmodemcg - zcg;
rcgx_motr = xmotorcg - xcg; rcgy_motr = ymotorcg - ycg; rcgz_motr = zmotorcg - zcg;
rcgx_shock = xshockcg - xcg; rcgy_shockl = yshocklcg - ycg; rcgy_shockr = yshockrcg - ycg; rcgz_shock = zshockcg - zcg;

rcgx_lw = xcg_lw - xcg; rcgy_lw = ycg_lw - ycg; rcgz_lw = zcg_lw - zcg;    
rcgx_rw = xcg_rw - xcg; rcgy_rw = ycg_rw - ycg; rcgz_rw = zcg_rw - zcg;    

Ixx = Inc(1,1)+mnc*(rcgy_nc^2+rcgz_nc^2)+Icb(1,1)+mcb*(rcgy_cb^2+rcgz_cb^2)+Irc(1,1)+mrc*(rcgy_rc^2+rcgz_rc^2)+...
    It(1,1)+mt*(rcgy_t^2+rcgz_t^2)+Iell(1,1)+mell*(rcgy_ell^2+rcgz_ell^2)+Ielr(1,1)+melr*(rcgy_elr^2+rcgz_elr^2)+...
    Irud(1,1)+mrud*(rcgy_rud^2+rcgz_rud^2)+Ixxbatt+battmass*(rcgy_batt^2+rcgz_batt^2)+Ixxprop+propmass*(rcgy_prop^2+rcgz_prop^2)+...
    Ixxservo+(servomass/2)*(rcgy_servol^2+rcgz_servo^2)+(servomass/2)*(rcgy_servor^2+rcgz_servo^2)+Ixxspdctrl+...
    speedctrlmass*(rcgy_spdctrl^2+rcgz_spdctrl^2)+Ixxgps+gpsmass*(rcgy_gps^2+rcgz_gps^2)+Ixxcam+cameramass*(rcgy_cam^2+rcgz_cam^2)+...
    Ixxapmodem+apmodemmass*(rcgy_apmod^2+rcgz_apmod^2)+Ixxmotor+motormass*(rcgy_motr^2+rcgz_motr^2)+...
    Ixxshock+(shocksmass/2)*(rcgy_shockl^2+rcgz_shock^2)+(shocksmass/2)*(rcgy_shockr^2+rcgz_shock^2)+...
    Ilw(1,1)+mlw*(rcgy_lw^2+rcgz_lw^2)+Irw(1,1)+mrw*(rcgy_rw^2+rcgz_rw^2);

Iyy = Inc(2,2)+mnc*(rcgx_nc^2+rcgz_nc^2)+Icb(2,2)+mcb*(rcgx_cb^2+rcgz_cb^2)+Irc(2,2)+mrc*(rcgx_rc^2+rcgz_rc^2)+...
    It(2,2)+mt*(rcgx_t^2+rcgz_t^2)+Iell(2,2)+mell*(rcgx_ell^2+rcgz_ell^2)+Ielr(2,2)+melr*(rcgx_elr^2+rcgz_elr^2)+...
    Irud(2,2)+mrud*(rcgx_rud^2+rcgz_rud^2)+Iyybatt+battmass*(rcgx_batt^2+rcgz_batt^2)+Iyyprop+propmass*(rcgx_prop^2+rcgz_prop^2)+...
    Iyyservo+(servomass)*(rcgx_servo^2+rcgz_servo^2)+Iyyspdctrl+...
    speedctrlmass*(rcgx_spdctrl^2+rcgz_spdctrl^2)+Iyygps+gpsmass*(rcgx_gps^2+rcgz_gps^2)+Iyycam+cameramass*(rcgx_cam^2+rcgz_cam^2)+...
    Iyyapmodem+apmodemmass*(rcgx_apmod^2+rcgz_apmod^2)+Iyymotor+motormass*(rcgx_motr^2+rcgz_motr^2)+...
    Ixxshock+(shocksmass)*(rcgx_shock^2+rcgz_shock^2)+Ilw(2,2)+mlw*(rcgx_lw^2+rcgz_lw^2)+...
    Irw(2,2)+mrw*(rcgx_rw^2+rcgz_rw^2);

Izz = Inc(3,3)+mnc*(rcgy_nc^2+rcgx_nc^2)+Icb(3,3)+mcb*(rcgy_cb^2+rcgx_cb^2)+Irc(3,3)+mrc*(rcgy_rc^2+rcgx_rc^2)+...
    It(3,3)+mt*(rcgy_t^2+rcgx_t^2)+Iell(3,3)+mell*(rcgy_ell^2+rcgx_ell^2)+Ielr(3,3)+melr*(rcgy_elr^2+rcgx_elr^2)+...
    Irud(3,3)+mrud*(rcgy_rud^2+rcgx_rud^2)+Izzbatt+battmass*(rcgy_batt^2+rcgx_batt^2)+Izzprop+propmass*(rcgy_prop^2+rcgx_prop^2)+...
    Izzservo+(servomass/2)*(rcgy_servol^2+rcgx_servo^2)+(servomass/2)*(rcgy_servor^2+rcgx_servo^2)+Izzspdctrl+...
    speedctrlmass*(rcgy_spdctrl^2+rcgx_spdctrl^2)+Izzgps+gpsmass*(rcgy_gps^2+rcgx_gps^2)+Izzcam+cameramass*(rcgy_cam^2+rcgx_cam^2)+...
    Izzapmodem+apmodemmass*(rcgy_apmod^2+rcgx_apmod^2)+Izzmotor+motormass*(rcgy_motr^2+rcgx_motr^2)+...
    Izzshock+(shocksmass/2)*(rcgy_shockl^2+rcgx_shock^2)+(shocksmass/2)*(rcgy_shockr^2+rcgx_shock^2)+...
    Ilw(3,3)+mlw*(rcgy_lw^2+rcgx_lw^2)+Irw(3,3)+mrw*(rcgy_rw^2+rcgx_rw^2);

Ixy = Inc(1,2)-mnc*(rcgx_nc*rcgy_nc)+Icb(1,2)-mcb*(rcgx_cb*rcgy_cb)+Irc(1,2)-mrc*(rcgy_rc*rcgx_rc)+...
    It(1,2)-mt*(rcgy_t*rcgx_t)+Iell(1,2)-mell*(rcgy_ell*rcgx_ell)+Ielr(1,2)-melr*(rcgy_elr*rcgx_elr)+...
    Irud(1,2)-mrud*(rcgy_rud*rcgx_rud)-battmass*(rcgy_batt*rcgx_batt)-propmass*(rcgy_prop*rcgx_prop)+...
    -(servomass/2)*(rcgy_servol*rcgx_servo)-(servomass/2)*(rcgy_servor*rcgx_servo)+...
    -speedctrlmass*(rcgy_spdctrl*rcgx_spdctrl)-gpsmass*(rcgy_gps*rcgx_gps)-cameramass*(rcgy_cam*rcgx_cam)+...
    -apmodemmass*(rcgy_apmod*rcgx_apmod)-motormass*(rcgy_motr*rcgx_motr)+...
    -(shocksmass/2)*(rcgy_shockl*rcgx_shock)-(shocksmass/2)*(rcgy_shockr*rcgx_shock)+...
    Ilw(1,2)-mlw*(rcgx_lw*rcgy_lw)+Irw(1,2)-mrw*(rcgx_rw*rcgy_rw);

Iyz = Inc(2,3)-mnc*(rcgz_nc*rcgy_nc)+Icb(2,3)-mcb*(rcgz_cb*rcgy_cb)+Irc(2,3)-mrc*(rcgz_rc*rcgx_rc)+...
    It(2,3)-mt*(rcgy_t*rcgz_t)+Iell(2,3)-mell*(rcgy_ell*rcgz_ell)+Ielr(2,3)-melr*(rcgy_elr*rcgz_elr)+...
    Irud(2,3)-mrud*(rcgy_rud*rcgz_rud)-battmass*(rcgy_batt*rcgz_batt)-propmass*(rcgy_prop*rcgz_prop)+...
    -(servomass/2)*(rcgy_servol*rcgz_servo)-(servomass/2)*(rcgy_servor*rcgz_servo)+...
    -speedctrlmass*(rcgy_spdctrl*rcgz_spdctrl)-gpsmass*(rcgy_gps*rcgz_gps)-cameramass*(rcgy_cam*rcgz_cam)+...
    -apmodemmass*(rcgy_apmod*rcgz_apmod)-motormass*(rcgy_motr*rcgz_motr)+...
    -(shocksmass/2)*(rcgy_shockl*rcgz_shock)-(shocksmass/2)*(rcgy_shockr*rcgz_shock)+...
    Ilw(2,3)-mlw*(rcgz_lw*rcgy_lw)+Irw(2,3)-mrw*(rcgz_rw*rcgy_rw);

Ixz = Inc(1,3)-mnc*(rcgx_nc*rcgz_nc)+Icb(1,3)-mcb*(rcgx_cb*rcgz_cb)+Irc(1,3)-mrc*(rcgx_rc*rcgz_rc)+...
    It(1,3)-mt*(rcgx_t*rcgz_t)+Iell(1,3)-mell*(rcgx_ell*rcgz_ell)+Ielr(1,3)-melr*(rcgx_elr*rcgz_elr)+...
    Irud(1,3)-mrud*(rcgx_rud*rcgz_rud)-battmass*(rcgx_batt*rcgz_batt)-propmass*(rcgx_prop*rcgz_prop)+...
    -(servomass)*(rcgx_servo*rcgz_servo)+...
    -speedctrlmass*(rcgx_spdctrl*rcgz_spdctrl)-gpsmass*(rcgx_gps*rcgz_gps)-cameramass*(rcgx_cam*rcgz_cam)+...
    -apmodemmass*(rcgx_apmod*rcgz_apmod)-motormass*(rcgx_motr*rcgz_motr)+...
    -(shocksmass)*(rcgx_shock*rcgz_shock)+Ilw(1,3)-mlw*(rcgz_lw*rcgx_lw)+Irw(1,3)-mrw*(rcgz_rw*rcgx_rw);