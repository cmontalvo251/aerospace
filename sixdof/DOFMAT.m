function [R,M] = DOFMAT(phi,theta,psi)

ctheta = cos(theta);
stheta = sin(theta);
ttheta = tan(theta);
sphi = sin(phi);
cphi = cos(phi);
spsi = sin(psi);
cpsi = cos(psi);

R = [ctheta*cpsi sphi*stheta*cpsi-cphi*spsi cphi*stheta*cpsi+sphi*spsi;ctheta*spsi sphi*stheta*spsi+cphi*cpsi cphi*stheta*spsi-sphi*cpsi;-stheta sphi*ctheta cphi*ctheta];
M = [1 sphi*ttheta cphi*ttheta;0 cphi -sphi;0 sphi/ctheta cphi/ctheta];
% Copyright - Carlos Montalvo 2015
% You may freely distribute this file but please keep my name in here
% as the original owner
