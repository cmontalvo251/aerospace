function dstatedt = dynamics(time,state)
global I invI

%%%Unwrap states
phi = state(1);
theta = state(2);
psi = state(3);
p = state(4);
q = state(5);
r = state(6);

%Compute sines and cosines
ctheta = cos(theta);
stheta = sin(theta);
ttheta = tan(theta);
sphi = sin(phi);
cphi = cos(phi);
spsi = sin(psi);
cpsi = cos(psi);

%%Derivatives
Kptp = [1 sphi*ttheta cphi*ttheta;0 cphi -sphi;0 sphi/ctheta cphi/ctheta];
pqr = [p q r]';
ptpdot = Kptp*pqr;  %%Euler Angle derivatives

%%Get Moments
LMN = thrusters(state,time);  %%Applied moments

%%Dynamics
Kuvw = cross_product(pqr);

Ipqr = I*[p q r]';
wIw = Kuvw*Ipqr;
pqrdot =invI*LMN-invI*wIw;  %%Angular acceleration equations. Does not include probe deployment

%%%Derivatives of state vector
dstatedt = [ptpdot;pqrdot];

% Copyright - Carlos Montalvo 2015
% You may freely distribute this file but please keep my name in here
% as the original owner
