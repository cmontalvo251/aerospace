function dxdt = Derivatives490Quad(xstate,t)

waypointcoordinates = [100 100 -20]';

x = xstate(1);
y = xstate(2);
z = xstate(3);
%%%%Kinematics
uvw = xstate(7:9);
phi = xstate(4);
theta = xstate(5);
psi = xstate(6);
T3 = [cos(psi) -sin(psi) 0;sin(psi) cos(psi) 0;0 0 1];
T2 = [cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)];
T1 = [1 0 0;0 cos(phi) -sin(phi); 0 sin(phi) cos(phi)];
%%Eq 1
T321 = T3*T2*T1;
xyzdot = T321*uvw;

%%%Rot Kinematics
pqr = xstate(10:12);
H = [1 sin(phi)*tan(theta) cos(phi)*tan(theta);0 cos(phi)...
    -sin(phi);0 sin(phi)/cos(theta) cos(phi)/cos(theta)];
ptpdot = H*pqr;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%THIS IS A FUNCTION OF YOUR AIRCRAFT

%%%%MASS PROPS HERE
%%% Mass Moment of Inertia Matrix in Body Axis (slug.ft^2)
Ixx = 0.9546;
Iyy = 0.9546;
Izz = 1.3766;
g = 9.81;
m = 1.25;


%%z Geometry and Inertias
I = [Ixx 0 0;0 Iyy 0;0 0 Izz];
Iinv = [1/Ixx 0 0;0 1/Iyy 0;0 0 1/Izz];

% Total airspeed
u = uvw(1);
v = uvw(2);
w = uvw(3);
uaero = u;
vaero = v;
waero = w;
Vaero = norm([uaero vaero waero]);

% Dynamic Pressure
rho = 1.22566578494891;


%Initialize Controls
de = 0;
da = 0;
dT = 0;
dr = 0;



%%% Calculation of body forces
p = pqr(1);
q = pqr(2);
r = pqr(3);


kpphi = 0.006;
kdphi = 0.001;
kptheta = kpphi;
kdtheta = kdphi;
kpz = 15;%
kdz = 18;%
kpwpsi = 100;%
kdwpsi = 300;%
kproll = 150;
kdroll = 230;
kppitch = 150;
kdpitch = 230;


xc = waypointcoordinates(1);
yc = waypointcoordinates(2);
zc = waypointcoordinates(3);
maxangle = 30;
xdotc = 0;
ydotc = 0;
zdotc = 0;

phic = -kpphi*(y-yc)-kdphi*(v-ydotc);
phidotc = 0;
thetac = kptheta*(x-xc)+kdtheta*(u-xdotc);
thetadotc = 0;
psic = 0;
psidotc = 0;

if phic > maxangle*pi/180
    phic = maxangle*pi/180;
end
if phic < -maxangle*pi/180
    phic = -maxangle*pi/180;
end
if thetac > maxangle*pi/180
    thetac = maxangle*pi/180;
end
if thetac < -maxangle*pi/180
    thetac = -maxangle*pi/180;
end
wz = kpz*(z-zc)+kdz*(w-zdotc);
wpsi = -kpwpsi*(psi-psic)-kdwpsi*(r-psidotc);
wroll = kproll*(phi-phic)+kdroll*(p-phidotc);
wpitch = kppitch*(theta-thetac)+kdpitch*(q-thetadotc);
wo = 4000*2*pi/60;


w1 = wo + wz + wpsi - wpitch;
w2 = wo + wz - wpsi - wroll;
w3 = wo + wz + wpsi + wpitch;
w4 = wo + wz - wpsi + wroll;
% w1 = wo;
% w2 = w1;
% w3 = w1;
% w4 = w1;
wbar = (w1+w2+w3+w4)/4;
ct = m*g/(4*wo^2);
cq = ct^(3/2)/sqrt(2);
Alc = .01;
Als = .01;
dx = 0.0426;
dy = 0.0426;

d = 0.254;%%meters
l = 0.254;%%meters 
rrotor = d/2;
Ttot = ct*(w1^2+w2^2+w3^2+w4^2);
if wbar > 0    
Dx = -Ttot*(((Alc/(wbar*rrotor))+dx)*u-(Als/(wbar*rrotor))*v);
Dy = -Ttot*((Als/(wbar*rrotor))*u+((Alc/(wbar*rrotor))+dy)*v);
Dz = 0;
Drag = [Dx Dy Dz]';
else
    Drag = [0 0 0]';
end
Lift = Ttot;
sigma1 = 1;
sigma2 = -1;
sigma3 = 1;
sigma4 = -1;
PHI1 = 0;
PHI2 = PHI1 + (pi/2);
PHI3 = PHI2 + (pi/2);
PHI4 = PHI3 + (pi/2);
Xa = Drag(1);
Ya = Drag(2);
Za = -Lift + Drag(3);
La = ct*((l*sin(PHI1)*w1^2)+(l*sin(PHI2)*w2^2)...
    +(l*sin(PHI3)*w3^2)+(l*sin(PHI4)*w4^2));
Ma = ct*((l*cos(PHI1)*w1^2)+(l*cos(PHI2)*w2^2)...
    +(l*cos(PHI3)*w3^2)+(l*cos(PHI4)*w4^2));
Na = cq*(sigma1*w1^2+sigma2*w2^2+sigma3*w3^2+sigma4*w4^2);

XYZa = [Xa;Ya;Za];
LMNa = [La;Ma;Na];


%%%FORCES AND MOMENTS HERE
g = 9.81;
FgravI = [0;0;m*g];
FgravB = T321'*FgravI;
XYZ = XYZa + FgravB;
LMN = LMNa;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%Translational Dynamics
pqrskew = [0 -r q;r 0 -p;-q p 0];
uvwdot = (1/m)*XYZ - pqrskew*uvw;

%%%%Rotational Dynamics
pqrdot = Iinv*(LMN-pqrskew*I*pqr);


dxdt = [xyzdot;ptpdot;uvwdot;pqrdot];