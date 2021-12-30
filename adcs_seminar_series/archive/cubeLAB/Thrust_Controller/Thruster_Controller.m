purge

fid = fopen('Input_Files/Satellite.THR');
ln = fgetl(fid); %%kp
ln = fgetl(fid); %%kp
ln = fgetl(fid); %%kp
ln = fgetl(fid); %%kd
ln = fgetl(fid); %%kd
ln = fgetl(fid); %%kd
ln = fgetl(fid); %%delay
ln = fgetl(fid); %%duration
ln = fgetl(fid); %%thrust
ln = fgetl(fid); %%Isp
ln = fgetl(fid); %%total mass
ln = fgetl(fid); %%deployment speed
ln = fgetl(fid); %%number of thrusters
NUMTHRUSTERS = 	str2num(ln(1:(find(ln=='!')-1)));

RN = [];
FN = [];

for idx = 1:NUMTHRUSTERS
    ln = fgetl(fid); %%x
    x = 	str2num(ln(1:(find(ln=='!')-1)));
    ln = fgetl(fid); %%y
    y = 	str2num(ln(1:(find(ln=='!')-1)));
    ln = fgetl(fid); %%z
    z = 	str2num(ln(1:(find(ln=='!')-1)));
    rCG2THR = [x;y;z]/1000;
    ln = fgetl(fid); %%yaw
    psi = 	str2num(ln(1:(find(ln=='!')-1)))*pi/180;
    ln = fgetl(fid); %%pitch
    theta = 	str2num(ln(1:(find(ln=='!')-1)))*pi/180;
    ln = fgetl(fid); %%sign
    sgn = 	str2num(ln(1:(find(ln=='!')-1)));
    UnitForceVector = sgn*R123(0,theta,psi)*[1;0;0];
    RN(:,idx) = cross(rCG2THR,UnitForceVector);
    FN(:,idx) = UnitForceVector;
end

Mdesired = [0;0;0.0];

RNGAUSS = RN'*inv(RN*RN')
FNGAUSS = FN'*inv(FN*FN')

break

%%%%Minimize thr'*thr subject to constraint that Md = RN*thrusters
%%%%In standard form we have maximize -thr'*thr with RN*thrusters-Md=0
%%%%thus L = -thr'*thr + lambda'*(RN*thrusters-Md)
%%%%Derivatives dL/dthr = 0
%%%% 0 = -thr + lambda*RN --> thr = lambda'*RN
%%%% dL/lambda = 0
%%%% 0 = RN*thr-Md --> Md = RN*thr --> Md = RN*lambda'*RN
issat = 1;
M0 = Mdesired;
RN0 = RN;
%%%First computation of control
thrusters = RN'*inv(RN*RN')*Mdesired; %%%This gives both positive and negative values of thrust which is not possible
num_sat = 0;
sat_controls = [];
while issat
  %%%First what we do is find all the thrusters bigger than 1.
  sat_controls = [sat_controls;find(thrusters < 0);find(thrusters>1)];
  not_sat = 1:length(thrusters);
  not_sat(sat_controls) = [];
  if (length(sat_controls) == num_sat || isempty(sat_controls))
    issat = 0;
  else
    %%%set saturated control to 1
    thrusters(thrusters>1) = 1;
    thrusters(thrusters<0) = 0;
    %%%Then hit each one of these with RN and reduce the order of
    %RN
    Mdesired = M0 - RN0(:,sat_controls)*thrusters(sat_controls);
    RN = RN0;
    RN(:,sat_controls) = [];
    RNRN = RN*RN';
    %%%Check to make sure RN*RN' is not undefined
    if det(RN*RN') < 1e-8
      issat = 0;
    else
      %%%Re-compute controls
      reduced = RN'*inv(RN*RN')*Mdesired;
      %%%Repopulate controls
      thrusters(not_sat) = reduced;
      %%%Save number of saturated control 
      num_sat = length(sat_controls);
    end
  end
  thrusters
end    

%%%Check computation
RN = RN0;
Mdesired = M0
Mget = RN*thrusters
