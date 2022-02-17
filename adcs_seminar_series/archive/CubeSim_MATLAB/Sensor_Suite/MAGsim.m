function [loc MAGout MAGpure MAGdotpure tmeasuredMAG] = MAGsim(phiactual,thetaactual,psiactual,NOISETYPE,MAGupdate,tactual,pactual,qactual,ractual)
global magfield
loc = [];
%%%This function will add in MAG errors
%%The input must be a vector with phi,theta, and psi euler angles

%% Three Axis Magnetometer Modeling Element

%%%Set Magnetometer Frame

xmagfield = [1];
ymagfield = [1];
zmagfield = [1];

magfield = [xmagfield;ymagfield;zmagfield];

magfield = magfield./norm(magfield); %normalize mag vector

%%Sensor Misalignment

misalign1 = randn*(1/1000);
misalign2 = randn*(1/1000);
misalign3 = randn*(1/1000);

%%Transformation (body to gyro - misalignment)

Tbm = [cos(misalign2)*cos(misalign3) sin(misalign1)*sin(misalign2)*cos(misalign3)-cos(misalign1)*sin(misalign3) cos(misalign1)*sin(misalign2)*cos(misalign3)+sin(misalign1)*sin(misalign3);cos(misalign2)*sin(misalign3) sin(misalign1)*sin(misalign2)*sin(misalign3)+cos(misalign1)*cos(misalign3) cos(misalign1)*sin(misalign2)*sin(misalign3)-sin(misalign1)*cos(misalign3);-sin(misalign2) sin(misalign1)*cos(misalign2) cos(misalign1)*cos(misalign2)];

%%%NOISE PARAMETERS(ALL STANDARD DEVIATIONS)

sigmanoise = ones(3,1).*0.01;
sigmabias = ones(3,1).*0.01;
magbias = randn(3,1).*sigmabias;

%%Cross axis sensitivities;

XYcross = randn/1000;  
YZcross = randn/1000;
XZcross = randn/1000;

%%%%Orthogonality and Scale Factor

ScaleFactor = ones(1,3) + randn(1,3).*(0.001/100);
Orthogonality = [0 XYcross XZcross;XYcross 0 YZcross;XZcross YZcross 0];
Ortho_ScaleFactor = diag(ScaleFactor) + Orthogonality;

%%%for loop around time history
MAGout = zeros(3,floor(tactual(end)/MAGupdate+1)); %%<- place holder
MAGpure = phiactual;
MAGdotpure = MAGpure;
tmeasuredMAG = MAGout(1,:);
counter = 1;
MAGupdatetime = 0;

for i = 1:length(phiactual)

  psi = psiactual(i);
  theta = thetaactual(i);
  phi = phiactual(i);
  
  %%%Transformation from body to inertial frame
  
  tbi(1,1) = cos(theta)*cos(psi);
  tbi(1,2) = cos(theta)*sin(psi);
  tbi(1,3) = -sin(theta);
  tbi(2,1) = sin(phi)*sin(theta)*cos(psi) - cos(phi)*sin(psi);
  tbi(2,2) = sin(phi)*sin(theta)*sin(psi) + cos(phi)*cos(psi);
  tbi(2,3) = sin(phi)*cos(theta);
  tbi(3,1) = cos(phi)*sin(theta)*cos(psi) + sin(phi)*sin(psi);
  tbi(3,2) = cos(phi)*sin(theta)*sin(psi) - sin(phi)*cos(psi);
  tbi(3,3) = cos(phi)*cos(theta);
  
  %%Tranform Mag vector from inertial to body frame
  
  magbody = tbi*magfield;
  
  %%Tranform Mag vector from body to sensor
  
  magsensor = Tbm*magbody;
  
  %%%Insert Noise%%%%%
  if (tactual(i) >= MAGupdatetime - 1e-8)
    loc = [loc i];
    if NOISETYPE == 2 %%All errors
      magnoise = randn(3,1).*sigmanoise;
      MAGout(1:3,counter) = Ortho_ScaleFactor*(magsensor + magbias + magnoise);
    end
  
    if NOISETYPE == 1 %%White Noise(Uniform)
	%%Position Noise
	MAGout(1:3,counter) = magsensor - sigmanoise + (2.*sigmanoise).*rand(3,1);
    end
    if NOISETYPE == 0 %%off
      MAGout(1:3,counter) = magsensor;
    end
    MAGupdatetime = MAGupdatetime + MAGupdate;
    tmeasuredMAG(counter) = tactual(i);
    counter = counter + 1;
  end
  %%%%%%%%%%%%%%%%%%
  MAGpure(1:3,i) = magsensor;
  p = pactual(i);
  q = qactual(i);
  r = ractual(i);
  wcross = [0 r -q;-r 0 p;q -p 0];
  MAGdotpure(1:3,i) = wcross*magsensor;
end  










