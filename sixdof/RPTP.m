function ptp = RPTP(R)
%%this function takes in the R matrix of a 3-2-1 rotation and
%extracts the euler angles phi theta and psi

theta = -asin(R(1,3));
sinpsi = R(1,2)/cos(theta);
cospsi = R(1,1)/cos(theta);
psi = atan2(sinpsi,cospsi);
sinphi = R(2,3)/cos(theta);
cosphi = R(3,3)/cos(theta);
phi = atan2(sinphi,cosphi);

ptp = [phi;theta;psi];

% Copyright - Carlos Montalvo 2015
% You may freely distribute this file but please keep my name in here
% as the original owner
