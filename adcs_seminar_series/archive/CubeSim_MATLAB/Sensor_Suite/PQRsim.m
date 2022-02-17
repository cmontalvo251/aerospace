function PQRout = PQRsim(PQRin,NOISETYPE)
global tupdate
%%%This function will add in PQR errors
%%The input mus be a vector with 3 rows and N columns
%The rows are p,q,r respectively
PQRout = PQRin;

[r,c] = size(PQRin);
if r ~= 3
  disp('PQRin does not have 3 rows')
  return
end
%%%NOISE PARAMETERS(ALL STANDARD DEVIATIONS)

sigmanoise = ones(3,1).*(1*pi/180); %%deg/s to rad/s
sigmaturnon = ones(3,1).*(100*pi/180*(1/3600)); %%deg/hr to rad/s
sigmadrift = ones(3,1).*(10*pi/180*(1/3600));
ScaleFactor = ones(1,3) + randn(1,3).*(0.001);
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

%%%%%%%PQR NOISE%%%%%%%%

if NOISETYPE == 2 %%All errors
   pqrturnon = randn(3,1).*sigmaturnon;
   pqrdrift = [0;0;0];
   for i = 1:c
      pqrdrift = pqrdrift + sigmadrift.*randn(3,1);
      pqrbias = pqrturnon + pqrdrift;
      pqrnoise = randn(3,1).*sigmanoise;
      PQRout(1:3,i) = Ortho_ScaleFactor*T*(PQRin(1:3,i) + pqrbias + pqrnoise);
   end
end

if NOISETYPE == 1 %%White Noise(Uniform)
  %%Position Noise
  for i = 1:c
    PQRout(1:3,i) = PQRin(1:3,i) - sigmaNoise + (2.*sigmaNoise).*rand(3,1);
  end
end
if NOISETYPE == 0 %%off
  PQRout = PQRin;
end

%%%%%%%%%%%%%%%%%%%



