function Drag()
%%% Calibrate Drag %%%
close all
clc

% Drag (x in feet/3.28 = meters)
xvideo = [0 .25 .5 .75 1 1.25 1.5 1.75 2 2.25 2.5 2.75 3 3.25 3.5 3.75 4 4.25 4.5 4.75 5]/3.28;
%%%Tvideo in seconds
tvideo = [4.37 4.43 4.53 4.63 4.73 4.83 4.93 5.03 5.13 5.23 5.37 5.47 5.57 5.73 5.87 6.0 6.1 6.17 6.23 6.37 6.53];

tvideo = tvideo - tvideo(1);        % Initialize at 0 
xvideo = xvideo - xvideo(1); 

[h1,ax1] = plottool(1,'Drag Estimation',18,'Time (sec)','x (m)');
plot(tvideo,xvideo,'b*','MarkerSize',10)

%%%Numerical Differentiation of Video
xdotvideo = xvideo;
for idx = 1:length(xvideo)-1
    xdotvideo(idx) = (xvideo(idx+1)-xvideo(idx))/(tvideo(idx+1)-tvideo(idx));
end

%%%%%Polyfit position data
coeff = polyfit(tvideo,xvideo,2);
t0 = tvideo(1);
tfinal = tvideo(end);
timestep = 0.01;
time = t0:timestep:tfinal;

%%%%Plotting derivative
[h2,ax2] = plottool(1,'Drag Estimation',18,'Time (sec)','xdot (m/s)');
plot(tvideo,xdotvideo,'b*','MarkerSize',10)

%%%%Differentiating polyfit data to get velocity
a1 = coeff(2);
a2 = coeff(1);
xdotest = a1 + 2*a2*time;

%%%%Getting initial conditions
x0 = xvideo(1);
xdot0 = xdotest(1);

state = [x0;xdot0];
StateOUT = zeros(2,length(time));       % Initialize StateOUT vector

%%%%Simulation using Eulers Method
for idx = 1:length(time)
  StateOUT(:,idx) = state;
  statedot = Derivatives(state,time(idx));
  state = state + statedot*timestep;
end

%%%%Extract State information
xsimulation = StateOUT(1,:);
xdotsimulation = StateOUT(2,:);

hold on
plot(ax1,time,xsimulation,'g-','LineWidth',2)
plot(ax2,time,xdotsimulation,'g-','LineWidth',2)
legend(ax1,'Experimental Data','Simulated Data')
legend(ax2,'Experimental Data','Simulated Data')
%legend(ax1,'Experimental Data','PolyFit Data','Simulated Data')
%legend(ax2,'Experimental Data','PolyFit Data','Simulated Data')

%%%Get cost estimate to tune C_D
xInterp = interp1(time,xsimulation,tvideo); 
cost = sum((xvideo-xInterp).^2)*3.28*12*length(xvideo)

function statedot = Derivatives(state,t)

x = state(1);
xdot = state(2);

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

C_D = 0.125;
rho = 1.225;
Drag = -0.5*rho*abs(xdot)*Vol^(2/3)*C_D*xdot;

X = Drag;

udot = X/(m+k1);

xdbldot = udot;

statedot = [xdot;xdbldot];
