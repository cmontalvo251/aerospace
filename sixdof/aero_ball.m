function [XYZ,LMN] = aero_ball(time,state,inertia)
%state is a 12x1, and inertia is a cell array with inertia{1} =
%mass and inertia{2} = I(moment of inertia matrix)
%LMN and XYZ are column vectors of body moments and forces respectively

x = state(1);
y = state(2);
z = state(3);
phi = state(4);
theta = state(5);
psi = state(6);
u = state(7);
v = state(8);
w = state(9);
p = state(10);
q = state(11);
r = state(12);

rho = 1.225;
V = sqrt(u^2+ v^2 + w^2);
CD = 1;
d = 1;
S = pi*d^2/4;

R = R123(phi,theta,psi);
mass = inertia{1};
g = 9.81;

FxP = 0;
FyP = 0;
zc = 10;
kp = 1000;
FzP = (kp*(zc-z));
Props = [FxP;FyP;FzP];

XYZ = -0.5*rho*V*S*[u;v;w]*CD + mass*g*R'*[0;0;1] + Props;

if abs(V) > 0
    phat = p*d/(2*V);
    qhat = q*d/(2*V);
    rhat = r*d/(2*V);
else
    phat  = 0;
    qhat = 0;
    rhat = 0;
end

C_rotate_damping = -3;

LMN = -0.5*rho*V^2*S*d*C_rotate_damping*[phat;qhat;rhat];