function dxdt = Derivatives(xstate,t,ustate)
%%%%Parameters of Boat
m = 0.697; %%kg
rho = 1.225; %%kg/m^3
rad = (4/12)/3.28; %%radius of prop in meters
A = pi*rad^2; %%of the prop in m^2
lr = (10/12)/3.28; %%%distance from rudder to cg in meters
%%%These are aero coefficients coming from a NACA 0012 (Online)
Cd0 = 0.006; %%nd - non-dimensional
Cda = 0.3;  %%nd
Cla = 5.7; %%nd
Sr = (2.25)*5/(144*3.28^2); %%%area of the rudder in m^2
a = (7.75/12)/3.28; %%%width of the boat in meters
b = (18/12)/3.28; %%%length of the boat in meters
Izz = (1/12)*m*(a^2+b^2); %%%polar mass moment of inertia of a rectangular prism in m^2

%%%We will figure these out later
crotate = 0.0075; 
cfront = 0.2;
csideways = 0.26;

%%%%Extract all the states and controls from the input variables
x = xstate(1);
y = xstate(2);
psi = xstate(3);
u = xstate(4);
v = xstate(5);
r = xstate(6);
dt = ustate(1); %%%0-100
dr = ustate(2); %%%angle

%%%%Compute xdot and ydot
%%%%This equation relates the body frame components
%%%Of velocity to the inertial frame components
xdot = cos(psi)*u - sin(psi)*v;
ydot = sin(psi)*u + cos(psi)*v;

%%%%psidot 
psidot = r;

%%%Compute the Total forces on the boat

%%%Compute vi which is induced flow from the prop
vi = dt*(6.5)/100; %Performed an experiment and measured this parameter

%%%%Compute Thrust - 1-D Momentum Theory
T = 2*rho*A*vi^2;

%%%Compute the local velocity at the rudder
%%%Velocity at the rudder includes contributions from the boat velocity,
%%%rotation of the boat and the fact that the rudder is offset from the cg
%%%and finally propwash (vi).
ur = u + vi;
vr = v - r*lr;

%%%Compute Angle of attack without rudder deflection - this just comes from geometry
alfa = atan2(vr,ur);

%%%Compute the rudder angle of attack
alfar = alfa + dr;

%%%Total Velocity at the rudder
Vr2 = ur^2+vr^2;

%%%Compute Lift and Drag
L = 0.5*rho*Vr2*Sr*(Cla*(alfar));
D = 0.5*rho*Vr2*Sr*(Cd0 + Cda*alfar^2);

%%%Compute the Rudder Forces - Rotating L and D to the body frame
Xr = -cos(alfa)*D + sin(alfa)*L;
Yr = -sin(alfa)*D - cos(alfa)*L;

%%%Forces include contributions from water drag which is different
%%%depending on front or sideways because the cross section of the boat is
%%%different. Thrust is included as well and rudder forces
X = -cfront*u + T + Xr;
Y = -csideways*v + Yr;

%%%%Dynamics - Newton's 2nd Law - (where do the rv and ru terms come from?
%%%%- The derivative transport theorem - This is a consequence of the body
%%%%frame rotating.
udot = X/m + r*v;
vdot = Y/m - r*u;

%%%Total Moments on the Boat
N = -crotate*r - lr*Yr;

%%%%Rotational Dynamics
rdot = N/Izz;

dxdt = [xdot ydot psidot udot vdot rdot]';










