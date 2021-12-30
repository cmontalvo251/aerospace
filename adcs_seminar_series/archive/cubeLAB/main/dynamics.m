function dstatedt = dynamics(time,state)

%%%Unwrap states
phi = state(1);
theta = state(2);
psi = state(3);
p = state(4);
q = state(5);
r = state(6);
pmag0 = state(7); %%Covariance states
pmag1 = state(8);
pmag2 = state(9);

%Compute sines and cosines
ctheta = cos(theta);
stheta = sin(theta);
ttheta = tan(theta);
sphi = sin(phi);
cphi = cos(phi);
spsi = sin(psi);
cpsi = cos(psi);

%%Mass and inertia
m = 5;
Ixx = 1;
Iyy = 2;
Izz = 3;
I = [Ixx 0 0;0 Iyy 0;0 0 Izz]; %%This doesn't have the cross products of inertia

%%Derivatives
Kptp = [1 sphi*ttheta cphi*ttheta;0 cphi -sphi;0 sphi/ctheta cphi/ctheta];
pqr = [p q r]';
ptpdot = Kptp*pqr;  %%Euler Angle derivatives

%%Get Moments
LMN = [0;0;0];  %%Applied moments

%%Dynamics
Kuvw = cross_product(pqr);
invI = inv(I);
pqrdot =invI*LMN-invI*(Kuvw*(I*[p q r]'));  %%Angular acceleration equations. Does not include probe deployment

%%%Covariance Dynamics of Attitude Dynamics
%%%THIS DOES NOT WORK YET
%P = [pmag0 pmag1;pmag1 pmag2];
%%%Q is estimate of model noise
%Q = 0.01*eye(2);
%Pdot = A*P + P*A' + Q;
pmagdot = [0;0;0];

%%%Derivatives of state vector
dstatedt = [ptpdot;pqrdot;pmagdot];

% Copyright - Carlos Montalvo 2015
% You may freely distribute this file but please keep my name in here
% as the original owner
