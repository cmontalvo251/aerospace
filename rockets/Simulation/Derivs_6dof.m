function dxdt = Derivs_6dof(xin,time)
global m Ixx Iyy T

g = 32.2; %ft/s^2
%m = 0.10648/g; %lb
%mfuel = 0.640699223*2.2/32.2; %slugs
mfuel = ((81-67.58)/16)/32.2;
totaltburn = 2.56;%s
mdot = mfuel/totaltburn;
if time > totaltburn
    tmax = totaltburn;
else
    tmax = time;
end

m = (81/16)/32.2-tmax*mdot;
%Ixx=5.46*10^-7;
%Iyy=3.572*10^-5;
%Izz=3.572*10^-5;

%d = .0432; %ft
d = .25;
rad = d/2; %ft
%Len = 3.833; %ft
Len = 3.7976;

rho = 0.00238; %slugs/ft^3
Cx0=0.5;        %zero yaw drag 
%Cx0 = 0.23327635880357467; 
%Cx2 = 5.4691147059580896;
Cx2=4.472;    %yaw drag
%Cna = 5.9431217210808551;
Cna = 2.27253;  %normal force 
Cnpa = 2.04031321285915584*10^-2;   %magnus force of the mass center
Cypa = -4.86262022555776757*10^-3;  %magnus force
Clp = -2.6505781136747761;
%Clp = 0;    %roll damping
%Clp = .048528822396295676; Ballistic
Cmq = .0117454;  %pitch damping

%Cna = 11.303; %%%Estimated from open rocket
Cmq = 0.2008; %%%Estimate from Open rocket

Cg_rear = 17.5920/12 ; %ft
%Cg_rear=0.22486;  
COP_rear = 10.2450/12; %ft
%COP_rear=0.1387;
MCOP_rear = (Cnpa/Cypa)+Cg_rear; %ft
slcop=COP_rear-Cg_rear;
mslcop=MCOP_rear-Cg_rear;

x = xin(1);
y = xin(2);
z = xin(3);
phi = xin(4);
theta =  xin(5);
psi = xin(6);
u = xin(7);
v = xin(8);
w = xin(9);
p = xin(10);
q =  xin(11);
r =  xin(12);

cphi = cos(phi);
ctheta = cos(theta);
cpsi = cos(psi);
sphi = sin(phi);
stheta = sin(theta);
spsi = sin(psi); 
ttheta = tan(theta);

% tmat;

T = [ctheta*cpsi, (sphi*stheta*cpsi)-(cphi*spsi), (cphi*stheta*cpsi)+(sphi*spsi);
    ctheta*spsi, (sphi*stheta*spsi)+(cphi*cpsi), (cphi*stheta*spsi)-(sphi*cpsi);
    -stheta, sphi*ctheta, cphi*ctheta];

%Kinematic Equations
xyzdots = T*[u;v;w];

%%%Dynamic Equations
V = sqrt(u.^2+v.^2+w.^2);
%cd = 0;
XYZg = m*g*[-stheta; sphi*ctheta; cphi*ctheta];


%6dof

%I changed his H mat to R here
R = [1, sphi*ttheta, cphi*ttheta; 0, cphi, -sphi; 0, sphi/ctheta, cphi/ctheta];

phithetapsidots = R*[p, q, r]';

%%Moments of Inertia (slugs-ft^2)
Ixx = (m/12)*(rad^2+rad^2+0.10); %%%0.1 added to include fins
Iyy = (m/12)*(Len^2+d^2);
Izz = Iyy;

Ic = [Ixx, 0, 0; 0, Iyy, 0; 0, 0, Izz;];

Q = (pi/8)*rho*V^2*d^2;
Scp = [0, 0, 0; 0, 0, -slcop; 0, slcop, 0];
Smp = [0, 0, 0; 0, 0, -mslcop; 0, mslcop, 0];


if abs(V) > 1e-4    
    XYZa = Q*[-Cx0-Cx2*((v^2+w^2)/(V^2)); -Cna*(v/V)+Cypa*(w/V)*((p*d)/(2*V)); -Cna*(w/V)+Cypa*(v/V)*((p*d)/(2*V))];
    LMNS = Scp*Q*[-Cx0-Cx2*(v^2+w^2)/V^2; -Cna*(v/V); -Cna*(w/V)];
    LMNU=Smp*Q*[0; Cypa*(w/V)*((p*d)/(2*V)); Cypa*(v/V)*((p*d)/(2*V))];
    LMNP = Q*d*[Clp*((p*d)/(2*V)); Cmq*((q*d)/(2*V)); Cmq* ((r*d)/(2*V))];
    LMN=LMNS+0*LMNU+LMNP;
else
    XYZa = [0;0;0];
    LMN = [0;0;0];
end
%%%%Taken by plotting thrust in open rocket!!
tburn =[0 ,   0.0386,    0.4242,    0.8869,    1.1183,    1.8123,    2.5835,    4.4344];
Tvec = [77.8022,  77.8022,   72.2344,   77.5092,   71.3553,   68.4249,   0,0];
if time > totaltburn
    thrust = 0;
else
    thrust = interp1(tburn,Tvec,time);
end

XYZa(1)=XYZa(1)+thrust;
T = thrust;

L = LMN(1,1);
M = LMN(2,1);
N = LMN(3,1);



pqrdots = inv(Ic)*(LMN-[0, -r, q; r, 0, -p; -q, p, 0]*Ic*[p; q; r]);
uvwdots = (1/m)*(XYZg+XYZa)-[0, -r, q; r, 0, -p; -q, p, 0]*[u; v; w];

dxdt = [xyzdots; phithetapsidots; uvwdots; pqrdots]';




