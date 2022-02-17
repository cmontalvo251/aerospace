%%%%6DOF Satellite Model
%%%%written by Carlos Montalvo - 10/22/2015
purge

tic

addpath Sensor_Suite

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%MODEL AND SATELLITE SPECIFIC PARAMETERS%%%

parameters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%TIME VECTOR (units are seconds)%%
t0 = 0;
t_vec = t0:(timestep*tskip):tfinal;
tGPSUpdate = 0;
tControlUpdate = 0;
tSaveUpdate = 0;
savectr = 1;

%%%%%%%INITIAL CONDITIONS OF SATELLITE%%%%%%

x0 = 6600*(1000); %%meters
xd0 = 0; %%%meters/second
y0 = 0; %meters
Vsat = sqrt(M*G/x0); %%meters/second
inclination_angle = 52.6*pi/180;
yd0 = Vsat*sin(inclination_angle); %%meters/second
z0 = 0; %%meters
zd0 = Vsat*cos(inclination_angle); %%meters/second
phi0 = 0;
theta0 = 0;
psi0 = 0;
q0123 = euler2quat([phi0 theta0 psi0]);
%p0 = 20*(rand()-0.5); %%rad/s
%q0 = 20*(rand()-0.5);
%r0 = 20*(rand()-0.5);
p0 = 0.1;
q0 = -0.1;
r0 = 0.2;
state_XYZ = [x0 y0 z0 xd0 yd0 zd0]';
state_PTP = [q0123 p0 q0 r0]';
state_XYZ_vec = zeros(6,length(t_vec));
state_PTP_vec = zeros(7,length(t_vec));

%%%%Sensor Setup
if timestep > gpsUpdateRate
  lenGPS_vec = length(t_vec);
else
  lenGPS_vec = tfinal*gpsUpdateRate;
end
GPS_vec = zeros(6,lenGPS_vec);

%%%%%%MODEL VECTORS for Kalman Filter
state_XYZ_tilde = state_XYZ;
state_PTP_tilde = state_PTP;
state_XYZ_tilde_vec = zeros(6,length(t_vec));
state_PTP_tilde_vec = zeros(7,length(t_vec));

%%%%%%MEASURE STATE TO INITIALIZE MODEL%%%%%
t_GPS_vec = zeros(1,lenGPS_vec);
state_XYZ_tilde = GPSsim(state_XYZ,gpsUpdateRate,noiseTypeGPS);
GPS_vec(:,1) = state_XYZ_tilde;
t_GPS_vec(1) = 0;
gpsCounter = 2;
tGPSUpdate = tGPSUpdate + 1/(gpsUpdateRate);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%INTEGRATION LOOP%%%%%%%%%%%%%%%

for time = t0:timestep:tfinal
  
  %%%%%%%%%%%%%%%%%%%MEASURE STATE%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%USING ALL SENSORS%%%%%%%%%%%%%%%%%%%%%%%
  
  %%%Sample GPS
  if (time > tGPSUpdate - 1e-8) && (SIMXYZ)
    GPS_vec(:,gpsCounter) = GPSsim(state_XYZ,gpsUpdateRate,noiseTypeGPS);
    t_GPS_vec(gpsCounter) = time;
    tGPSUpdate = tGPSUpdate + 1/(gpsUpdateRate);
    gpsCounter = gpsCounter + 1;
  end
    
  %%%%%%%%%%%%SAVE STATE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
  if time >= tSaveUpdate
    if (SIMXYZ)
      state_XYZ_vec(:,savectr) = state_XYZ;
      state_XYZ_tilde_vec(:,savectr) = state_XYZ_tilde;
    end
    if (SIMPTP)
      state_PTP_vec(:,savectr) = state_PTP;
      state_PTP_tilde_vec(:,savectr) = state_PTP_tilde;
    end
      savectr = savectr + 1;
      tSaveUpdate = timestep*tskip + tSaveUpdate;
      %%%%%%%%%%%NOTIFY USER OF PROGRESS%%%%%%%%%%%%%%%%
      disp(['Simulation ',num2str(time/tfinal*100),' Percent Complete'])
  end
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %%%%%CONTROL ROUTINE (DICTATED BY SPEED OF PROCESSOR)
  
  if (time > tControlUpdate - 1e-8) && (SIMPTP)
    tControlUpdate = tControlUpdate + 1/(controlUpdateRate);
    LMNcontrol = computeControl(state_PTP_tilde); 
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
  %%%%%%%%%%%%%KALMAN FILTER UPDATE%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%UPDATE IS PERFORMED BY SPECIFIC SENSORS%%%%%%%
  
  if (time > tGPSUpdate - 1e-8) && (SIMXYZ)
    %state_tilde = computeKalmanGPSUpdate(state_tilde);
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %%%%%%%%%%INTEGRATION LOOP (RUNS EVERY TIMESTEP)%%%%%%%%

  if (SIMXYZ)
    state_XYZ = computeRK4Step(@computeXYZDerivs,state_XYZ);
    state_XYZ_tilde = computeRK4Step(@computeXYZDerivs,state_XYZ_tilde);
  end
  if (SIMPTP)
    state_PTP = computeRK4Step(@computePTPDerivs,state_PTP,LMNcontrol);
    state_PTP_tilde = computeRK4Step(@computePTPDerivs,state_PTP_tilde,LMNcontrol);
  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
    %%%%%%%%%NORMALIZE QUATERNIONS%%%%%%%%%%%%%%%%%%%%%%%%%%

    q0123 = state_PTP(1:4);
    q0123 = q0123./norm(q0123);
    state_PTP(1:4) = q0123;
    q0123 = state_PTP_tilde(1:4);
    q0123 = q0123./norm(q0123);
    state_PTP_tilde(1:4) = q0123;
    
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
end

state_PTP_vec(:,end) = state_PTP;
state_PTP_tilde_vec(:,end) = state_PTP_tilde;

%%%%%%%%%%PLOT EVERYTHING%%%%%%%%%

BodyName = '';
dim1 = 'km';
maincolor = 'k-';
colors = {'k-.','k--','k:'};

if (SIMXYZ)
  state_XYZ_vec(:,end) = state_XYZ;
  state_XYZ_tilde_vec(:,end) = state_XYZ_tilde;

  plottool(1,'xy',18,'X Orbit (km)','Y Orbit (km)','Z Orbit (km)');
  plot3(state_XYZ_vec(1,:)./1000,state_XYZ_vec(2,:)./1000,state_XYZ_vec(3,:)./1000)
  plot3(state_XYZ_tilde_vec(1,:)./1000,state_XYZ_tilde_vec(2,:)./1000,state_XYZ_tilde_vec(3,:)./1000,'r--')
  plot3(GPS_vec(1,:)./1000,GPS_vec(2,:)./1000,GPS_vec(3,:)./1000,'ks')
  legend('Numerical','Estimate','GPS')

  Names_XYZ = {['X ',BodyName],['Y ',BodyName],['Z ',BodyName],['U ',BodyName],['V ',BodyName],['W ',BodyName]};

  ylabels_XYZ = {['x (',dim1,') ',BodyName],['y (',dim1,') ',BodyName],['z (',dim1,') ',BodyName],['u (',dim1,'/s) ',BodyName],['v (',dim1,'/s) ',BodyName],['w (',dim1,'/s) ',BodyName]};

  cppdata = dlmread('../../Output_Files/State.OUT');  
  tpp = cppdata(:,1);
  
  %x,y,z,q0,q1,q2,q3,u,v,w ,p ,q , r
  %1,2,3,4 ,5 ,6 ,7 ,8,9,10,11,12,13
  for idx = 1:6
    plottool(1,Names_XYZ{idx},18,'Time (min)',ylabels_XYZ{idx});
    plot(t_vec./60,state_XYZ_vec(idx,:)./1000)
    plot(t_vec./60,state_XYZ_tilde_vec(idx,:)./1000,'r--')
    plot(t_GPS_vec./60,GPS_vec(idx,:)./1000,'ks')
    plot(tpp./60,cppdata(:,idx+1)./1000,'g-')
    legend('Numerical','Estimate','Sensor')
  end
  
end

if (SIMPTP)
  
  if CPPDATA
      cppdata = dlmread('../../Output_Files/State.OUT');
      tpp = cppdata(:,1);
  end
  
    
  ylabels_PTP = {['q0 (rad) ',BodyName],['q1 (rad) ',BodyName],['q2 (rad) ',BodyName],['q3 (rad) ',BodyName],['p (rad/s) ',BodyName],['q (rad/s) ',BodyName],['r (rad/s) ',BodyName]};
  Names_PTP = {['q0 ',BodyName],['q1 ',BodyName],['q2 ',BodyName],['q3 ',BodyName],['P ',BodyName],['Q ',BodyName],['R ',BodyName]};
  
  for idx = 1:7
    plottool(1,Names_PTP{idx},18,'Time (min)',ylabels_PTP{idx});
    plot(t_vec./60,state_PTP_vec(idx,:))
    plot(t_vec./60,state_PTP_tilde_vec(idx,:),'r--')
    if CPPDATA
        plot(tpp./60,cppdata(:,idx+7),'g-')
        legend('Numerical','Estimate','CPP')
    end
  end

  q0123_vec = state_PTP_vec(1:4,:);
  ptp_vec = quat2euler(q0123_vec');
  q0123_tilde_vec = state_PTP_tilde_vec(1:4,:);
  ptp_tilde_vec = quat2euler(q0123_tilde_vec');
  if CPPDATA
      q0123cpp_vec = cppdata(:,8:11);
      ptpcpp_vec = quat2euler(q0123cpp_vec);
  end
  
  eangle = {['\phi (deg) ',BodyName],['\theta (deg) ',BodyName],['\psi (deg) ',BodyName]};

  for idx = 1:3
    plottool(1,eangle{idx},18,'Time (nd)',eangle{idx});
    plot(t_vec./60,ptp_vec(:,idx)*180/pi)
    plot(t_vec./60,ptp_tilde_vec(:,idx)*180/pi,'r--')
    if CPPDATA
        plot(tpp./60,ptpcpp_vec(:,idx)*180/pi,'g-')
    end
  end
end

toc