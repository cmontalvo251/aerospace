function ACCELout = ACCELsim(UVWdotin,PQRin,UVWin,EULERin,PQRdotin,RCGtoA,NOISETYPE,GRAVITY)
%%ACCELout = ACCELsim(UVWdotin,PQRin,UVWin,EULERin,PQRdotin,RCGtoA,NOISETYPE,GRAVITY)
%%%This function will add in ACCELEROMETER errors
%%UVWdotin is the body acceleration vector udot,vdot,wdot
%%PQRin is the body angular velocity p,q,r
%%UVWin is the body velocity u,v,w
%%EULERin is the attiude of the projectile phi,theta,psi
%%PQRdotin is the angular acceleration pdot,qdot,rdot

%%All of the above inputs must be a 3 row matrix with N columns

%%RCGtoA is a 3 column vector denoting the 'actual' displacement
%from the CG

ACCELout = UVWdotin; %%Placeholder

Pin = PQRin(1,:);
Qin = PQRin(2,:);
Rin = PQRin(3,:);
Pdotin = PQRdotin(1,:);
Qdotin = PQRdotin(2,:);
Rdotin = PQRdotin(3,:);
PHI = EULERin(1,:);
THETA = EULERin(2,:);

[r,c] = size(PQRin);
if r ~= 3
  disp('PQRin does not have 3 rows')
  return
end

%%%NOISE PARAMETERS(ALL STANDARD DEVIATIONS)

if NOISETYPE == 2 %%Full Up Noise

  sigmanoise = ones(3,1).*(GRAVITY/1000); %%1 milliG
  sigmaturnon = ones(3,1).*(GRAVITY/100); %%10 milliG
  sigmadrift = ones(3,1).*(GRAVITY/1000); %%1 milliG
  ScaleFactor = ones(1,3) + randn(1,3).*(0.002/100);
  XYcross = randn/1000;  %%cross axis sensitivities;
  YZcross = randn/1000;
  XZcross = randn/1000;
  misalign1 = randn*(1/1000);
  misalign2 = randn*(1/1000);
  misalign3 = randn*(1/1000);

  %%%%Orthogonality and Scale Factor

  Orthogonality = [0 XYcross XZcross;XYcross 0 YZcross;XZcross YZcross 0];
  Ortho_ScaleFactor = diag(ScaleFactor) + Orthogonality;

  %%Transformation (body to gyro - misalignment)

  T = [cos(misalign2)*cos(misalign3) sin(misalign1)*sin(misalign2)*cos(misalign3)-cos(misalign1)*sin(misalign3) cos(misalign1)*sin(misalign2)*cos(misalign3)+sin(misalign1)*sin(misalign3);cos(misalign2)*sin(misalign3) sin(misalign1)*sin(misalign2)*sin(misalign3)+cos(misalign1)*cos(misalign3) cos(misalign1)*sin(misalign2)*sin(misalign3)-sin(misalign1)*cos(misalign3);-sin(misalign2) sin(misalign1)*cos(misalign2) cos(misalign1)*cos(misalign2)];

elseif NOISETYPE == 1 %%White Noise
  sigmadrift = [0;0;0];
  sigmaturnon = [0;0;0];
  sigmanoise = ones(3,1).*(GRAVITY/1000); %%1 mG
  Ortho_ScaleFactor = eye(3);
  T = eye(3);
elseif NOISETYPE == 0 %%No Noise
  sigmadrift = [0;0;0];
  sigmaturnon = [0;0;0];
  sigmanoise = [0;0;0];
  Ortho_ScaleFactor = eye(3);
  T = eye(3);
end

%%%%%%%ACCEL NOISE%%%%%%%%

accelturnon = randn(3,1).*sigmaturnon;
acceldrift = [0;0;0];
for i = 1:c
  acceldrift = acceldrift + sigmadrift.*randn(3,1);
  accelbias = accelturnon + acceldrift;
  accelnoise = randn(3,1).*sigmanoise;
  G = GRAVITY*[-sin(THETA(i));sin(PHI(i))*cos(THETA(i));cos(PHI(i))*cos(THETA(i))];
  omegacross = [0 -Rin(i) Qin(i);Rin(i) 0 -Pin(i);-Qin(i) Pin(i) 0];
  omegadotcross = [0 -Rdotin(i) Qdotin(i);Rdotin(i) 0 -Pdotin(i);-Qdotin(i) Pdotin(i) 0];
  ACCELout(1:3,i) = Ortho_ScaleFactor*T*(UVWdotin(1:3,i) + omegacross*UVWin(1:3,i) + G + omegadotcross*RCGtoA + omegacross*omegacross*RCGtoA + accelbias + accelnoise);
end

%%%%%%%%%%%%%%%%%%%


