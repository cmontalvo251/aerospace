function R = ROTMAT(angle,aor)
%%give me an angle that brings the frame from I to B and I will give
%you the rotation matrix R such that v(body) = R'*v(inertial)
%angle must be in radians

%%Create cosines and sines
ca = cos(angle);
sa = sin(angle);

switch aor
 case 1
  %About I1
  R = [1 0 0;0 ca sa;0 -sa ca];
 case 2
  %About I2
  R = [ca 0 -sa;0 1 0;sa 0 ca];
 case 3
  %About I3
  R = [ca sa 0;-sa ca 0;0 0 1];
end

% Copyright - Carlos Montalvo 2015
% You may freely distribute this file but please keep my name in here
% as the original owner
