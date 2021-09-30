function M = KPTP(phi,theta,psi)
%%this will compute the matrix such that
%%ptpdot = M*[p q r];

ctheta = cos(theta);
stheta = sin(theta);
ttheta = tan(theta);
sphi = sin(phi);
cphi = cos(phi);
spsi = sin(psi);
cpsi = cos(psi);

M = [1 sphi*ttheta cphi*ttheta;0 cphi -sphi;0 sphi/ctheta cphi/ctheta];
% Copyright - Carlos Montalvo 2015
% You may freely distribute this file but please keep my name in here
% as the original owner
