function [phi,theta,psi] = extract_Euler(T)
%%%Assuming R is a 3x3 matrix extract phi,theta,psi Euler angles
%%%assuming a 3-2-1 transformation sequence
%%%Let R be defined such that v(body) = T v(inertial)
%%In Mark Costello's notation, T would TBI

theta = -asin(T(1,3));
sphi  = T(2,3)/cos(theta);
cphi = T(3,3)/cos(theta);
phi = atan2(sphi,cphi);
spsi = T(1,2)/cos(theta);
cpsi = T(1,1)/cos(theta);
psi = atan2(spsi,cpsi);


% Copyright - Carlos Montalvo 2015
% You may freely distribute this file but please keep my name in here
% as the original owner
