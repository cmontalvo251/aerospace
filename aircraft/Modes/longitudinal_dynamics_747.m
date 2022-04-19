clear
clc
close all

%%%For this 747
rho = 0.00238;
S = 5600;
cbar = 48;
m = 564032/32.2;
Iyy = 32.3*10^6;
CDA = 0.66;
CLA = 5.7;
CMQ = -20.8;
u = 0.25*1116;
alfa = 5.7*pi/180;
CMA = -1.26;
CM0 = 0-CMA*alfa
CD0 = 0.0102-CDA*alfa^2
CL0 = 1.11-CLA*alfa
w = u*tan(alfa);
theta = alfa;
CLQ = 5.4;
g = 32.2;

Asp_DRc = [rho*S*u/(2*m)*(CD0-CLA) -rho*S*u*cbar/(4*m)*CLQ+u;
    rho*S*cbar*u/(2*Iyy)*CMA       rho*S*cbar^2*u/(4*Iyy)*CMQ];

Alon_DRc = [0 0 0 1 sin(theta) 0;
    0 0 -u -sin(theta) 1 0;
    0 0 0 0 0 1;
    0 0 -g*cos(theta) -rho*S*u*CD0/m rho*S*u*CL0/(2*m) -w;
    0 0 -g*sin(theta) -rho*S*u*CL0/m rho*S*u*(CD0-CLA)/(2*m) -rho*S*u*cbar/(4*m)*CLQ+u;
    0 0 0 rho*S*cbar*u*CM0/Iyy rho*S*cbar*u*CMA/(2*Iyy) rho*S*cbar^2*u*CMQ/(4*Iyy)]

eig(Alon_DRc)