clear
clc
close all

%%%SIMULATE FOR 0 to TMAX seconds
TMAX = 10;
tspan = [0 TMAX];

%%%Covariance of Attitude Model - THis does not work yet
pmag0 = 0;pmag1 = 0;pmag2 = 0;

%%%INITIAL CONDITIONS
%phi,theta,psi are the angles
%p,q,r are the angular rates
%xinitial = [phi0,theta0,psi0,p0,q0,r0,covariance states]
xinitial = [0;0;0;0.5;0.2;0.1;pmag0;pmag1;pmag2];
%%%Timestep of RK4
timestep = 0.01;


%%%%SETUP RK4
next = 10;
[N,flag] = size(xinitial);
tout = tspan(1):timestep:tspan(end);
integrationsteps = length(tout);
xout = zeros(N,integrationsteps);
x = xinitial;
threshold = 0;

%%%SETUP Rate gyro
gyro_update_rate = (1/10); %%%10 hz
tgyro_next = 0;
gyro_ctr = 1;
rate_gyro = zeros(3,ceil(tspan(2)/gyro_update_rate));
trate_gyro = zeros(1,ceil(tspan(2)/gyro_update_rate));

%%SETUP Magnetometer
mag_update_rate = (1/1); %%%1 hz
tmag_next = 0;
mag_ctr = 1;
mag = zeros(3,ceil(tspan(2)/mag_update_rate));
tmag = zeros(1,ceil(tspan(2)/mag_update_rate));

%%%RK4 SIMULATION
for ii = 1:length(tout)
    
  %%%NOTIFY USER
  time = tout(ii);
  if next
    percent = (floor(100*time/tout(end)));
    if percent >= threshold
      disp(['Simulation ',num2str(percent),'% Complete'])
      threshold = threshold + next;
    end
  end
  
  %%Save States
  xout(:,ii) = x;
  
  %%Integrate 4 times like RK4
  xdot1 = dynamics(time, x);
  xdot2 = dynamics(time + (.5*timestep), x + (xdot1*.5*timestep));
  xdot3 = dynamics(time + (.5*timestep), x + (xdot2*.5*timestep));
  xdot4 = dynamics(time + timestep, x + (xdot3*timestep));
  xdotRK4 = (1/6) * (xdot1 + (2*xdot2) + (2*xdot3) + xdot4);
  nextstate = x + (timestep * xdotRK4);
  x = nextstate;
  
  %%%Check for Sensor Update
  pqr_truth = x(4:6);
  ptp_truth = x(1:3);
  if time > tmag_next
    tmag_next = tmag_next + mag_update_rate;
    tmag(mag_ctr) = time+timestep;
    mag(:,mag_ctr) = ptp_truth + 0.01*(0.5-rand(3,1));
    
    %%%%Kalman Update -- THIS DOES NOT WORK YET
    %K = eye(3); %%%THIS IS WRONG
    %And so is this below so I'm commenting it out
    %x(1:3) = x(1:3) + K*(mag(:,mag_ctr)-x(1:3));
    
    mag_ctr = mag_ctr + 1;
  end
  
  if time > tgyro_next
    tgyro_next = tgyro_next + gyro_update_rate;
    trate_gyro(gyro_ctr) = time+timestep;
    rate_gyro(:,gyro_ctr) = pqr_truth + 0.01*(0.5-rand(3,1));
    
    %%%Kalman Update
    %K = eye(3);
    %x(4:6) = x(4:6) + K*(rate_gyro(:,gyro_ctr)-x(4:6));
    
    gyro_ctr = gyro_ctr + 1;
  end
end

%%%GET READY FOR PLOTTING
mag = mag';
rate_gyro = rate_gyro';
xout = xout'; %%This changes everything to a more file friendly format
tout = tout';

%%%PLOT EVERYTHING
BodyName = 'Satellite';
dim1 = 'deg';
Names_XYZ = {['\phi ',BodyName],['\theta ',BodyName],['\psi ',BodyName],['P ',BodyName],['Q ',BodyName],['R ',BodyName]};
ylabels_XYZ = {['\phi (',dim1,') ',BodyName],['\theta (',dim1,') ',BodyName],['\psi (',dim1,') ',BodyName],['p (',dim1,'/s) ',BodyName],['q (',dim1,'/s) ',BodyName],['r (',dim1,'/s) ',BodyName]};

for idx = 1:6
  fancy_plotting(1,Names_XYZ{idx},18,'Time (min)',ylabels_XYZ{idx});
  plot(tout,180/pi*xout(:,idx),'b-')
  hold on
  if idx <= 3
    plot(tmag,180/pi*mag(:,idx),'rs')
  else
    plot(trate_gyro,180/pi*rate_gyro(:,idx-3),'rs')
  end
  legend('Truth Signal','Sensor Measurement')
end
