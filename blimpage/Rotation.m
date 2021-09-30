function Rotation()
%%% Calibrate Rotation %%%
close all

% Thrust on
yawvideoON = [0 pi/2 pi 3*pi/2];
tvideoON = [8.43 11.17 13.17 14.9];
tvideoON(end)-tvideoON(1)

% Thrust off
yawvideoOFF = [0 pi/2 pi 3*pi/2 2*pi]+yawvideoON(end)+pi/2;
tvideoOFF = [16.37 18 20 22.7 25.73];

% Total Timehistory
yawvideo = [yawvideoON yawvideoOFF];
tvideo = [tvideoON tvideoOFF];

tvideo = tvideo - tvideo(1);

plottool(1,'Rotational Estimation',14,'Time (sec)','\psi (deg)');
plot(tvideo,yawvideo*180/pi,'b*','MarkerSize',10)

% Yaw vector
yaw0 = yawvideo(1);
yawdot0 = 0;
state = [yaw0;yawdot0];

% Time vector
t0 = tvideo(1);
tfinal = tvideo(end);
timestep = 0.01;
time = t0:timestep:tfinal;

StateOUT = zeros(2,length(time));

for idx = 1:length(time)
  StateOUT(:,idx) = state;
  statedot = Derivatives(state,time(idx));
  state = state + statedot*timestep;
end

yawEuler = StateOUT(1,:);

hold on
plot(time,yawEuler*180/pi,'g-','LineWidth',2)
%legend('Experimental Data','Simulated Data') 

yawEulerInterp = interp1(time,yawEuler,tvideo);

cost = sum((yawvideo-yawEulerInterp).^2)


function statedot = Derivatives(state,t)
yaw = state(1); 
yawdot = state(2);

b = (38/12)/3.28; %%%diameter of blimp
Vol = (4/3)*pi*(b/2)^3; %%%volume
mstructure = 32/1000; %%%%mass of structure
mgondola = 68/1000; %%%mass of gondola and batteries
density_helium = 0.18; %%%%kg/m^3 of helium
m = density_helium*Vol + mstructure + mgondola; %%%total mass
I =  (2/5)*m*(b/2)^2; %%%Inertia of a sphere
d = (5/12)/3.28; %%%Moment arm
k = (8/3)*(b/2)^3;
k1 = k;
k2 = k;
k3 = 0;
rho = 1.225;

if t > 16.37-8.43
    deltat = 0;         % Thruster is on (0-1)
else
    deltat = 1;
end

% sig = .0655;
% Cnr = .02;
sig = .045;
Cnr = .023;

Na = -0.5*rho*(Vol)^(5/3)*abs(yawdot)*Cnr*yawdot;

N = (sig*deltat*d) + Na;    % The total moment

yawdbldot = N/(k3 + I);    % EOMs

statedot = [yawdot;yawdbldot];
