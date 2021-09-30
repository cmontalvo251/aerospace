function R = R123(phi,theta,psi)
%R = R123(phi,theta,psi) or R = R123(ptp) where prp is a 3x1 vector
%compute R such that v(inertial) = R v(body)
%In Mark Costello's notation R would be TIB
%IF YOU WANT TO DO THIS WITH QUATERNIONS LOOK FOR RQUAT

if length(phi) == 3
  arg = phi;
  phi = arg(1);
  theta = arg(2);
  psi = arg(3);
end

%Compute sines and cosines
ctheta = cos(theta);
stheta = sin(theta);
sphi = sin(phi);
cphi = cos(phi);
spsi = sin(psi);
cpsi = cos(psi);

%Kinematics
R = [ctheta*cpsi sphi*stheta*cpsi-cphi*spsi cphi*stheta*cpsi+sphi*spsi;ctheta*spsi sphi*stheta*spsi+cphi*cpsi cphi*stheta*spsi-sphi*cpsi;-stheta sphi*ctheta cphi*ctheta];
% Copyright - Carlos Montalvo 2015
% You may freely distribute this file but please keep my name in here
% as the original owner
