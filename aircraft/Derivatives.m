function dxdt = Derivatives(xstate,t,ustate)
%%%%Extract all the states and controls from the input variables
x = xstate(1);
z = xstate(2);
theta = xstate(3);
u = xstate(4);
w = xstate(5);
q = xstate(6);
dt = ustate(1); %%%thrust
de = ustate(2); %%%elevator

%%%Compute the Total thrust

%%%Compute vi which is induced flow from the prop
vi = dt*(6.5)/100; %Performed an experiment and measured this parameter

%%%%Compute Thrust - 1-D Momentum Theory
A = (8/12)^2*pi;
rho = 1.225;
T = 2*rho*A*vi^2;

%%%Useful vars
ctheta = cos(theta);
stheta = sin(theta);

%%%%Aerodynamics?
alfa = atan2(w,u);

V = sqrt(u^2+w^2);
b = 2.04;
S = 0.6558;
cbar = 0.3215;
if V > 0
    qhat = q*b/(2*V);
else
    qhat = 0;
end

CL0 = 0.062;
CLA = 5.195;
CD0 = 0.028;
CDA = 1.3537;
CLDE = 0.2167;
CMDE = -0.8551;
CM0 = 0.0598;
CMQ = -5.263;
CLQ = 4.589;
CMALFA = -0.9317;

CL = CL0 + CLA*alfa + CLDE*de + CLQ*qhat;
CD = CD0 + CDA*alfa^2;
CM = CM0 + CMALFA*alfa + CMDE*de + CMQ*qhat;

Lift = 0.5*rho*V^2*S*CL;
Drag = 0.5*rho*V^2*S*CD;
Moment = 0.5*rho*V^2*S*cbar*CM;

%%%%%Vertical and Translational Dynamics
m = 5.6; %%%kg
g = 9.81;

X = Lift*sin(alfa) - Drag*cos(alfa) - m*g*stheta + T;
Z = -Lift*cos(alfa) - Drag*sin(alfa) + m*g*ctheta;

%%% udot + [0  0  q] u = udot + q*w
%%% vdot + [0  0  0] 0 = vdot 
%%% wdot + [-q 0  0] w = wdot - q*u

udot = X/m - q*w;
wdot = Z/m + q*u;

%%%Rotational Dynamics
Iyy = 0.5111;
qdot = Moment/Iyy;

if z >= 0
    if wdot > 0 
        w = 0;
        wdot = 0;
    end
    z = 0;
end

%%%%Kinematics
thetadot = q;
xdot = u*ctheta + w*stheta;
zdot = -u*stheta + w*ctheta;

dxdt = [xdot zdot thetadot udot wdot qdot]';










