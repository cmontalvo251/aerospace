function dxdt = Derivatives(xstate,t,ustate)
%%%%Extract all the states and controls from the input variables
x = xstate(1);
z = xstate(2);
phi = xstate(3);
xdot = xstate(4);
zdot = xstate(5);
p = xstate(6);
r1 = ustate(1); %%%rpm
r2 = ustate(2); %%%rpm

%%%%phidot
phidot = p;

%%%Compute the Total forces on the boat

%%%Compute vi which is induced flow from the prop
vi1 = r1*(2.5)/100; %Performed an experiment and measured this parameter
vi2 = r2*(2.5)/100; 

%%%%Compute Thrust - 1-D Momentum Theory
A = (8/12)^2*pi;
rho = 1.225;
T1 = 2*rho*A*vi1^2;
T2 = 2*rho*A*vi2^2;

%%%Useful vars
cphi = cos(phi);
sphi = sin(phi);

%%%%Aerodynamics?
V = sqrt(xdot^2+zdot^2);
CD = 1.0;
H = (2/12)/3.28;
L = 1/3.28;
FAz = -0.5*rho*V^2*2*A*CD*sign(zdot);
FAx = -0.5*rho*V^2*H*L*CD*sign(xdot);

%%%%%Vertical and Translational Dynamics
m = 1.2; %%%kg
g = 9.81;
zddot = (m*g - (T1+T2)*cphi)/m + FAz;
xddot = (T1+T2)*sphi/m + FAx;

%%%Rotational Dynamics
lr = (6/12)/3.28;
Ixx = m/12*(H^2+L^2);
pdot = (T1*(L/2)-T2*(L/2))/Ixx;

if z >= 0
    if zddot > 0 
        zdot = 0;
        zddot = 0;
    end
    z = 0;
end

dxdt = [xdot zdot phidot xddot zddot pdot]';










