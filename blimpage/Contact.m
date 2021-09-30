function cost = Contact(make_plots)
global K C del_max
%%% Calibrate Drag %%%
close all

del_max = 0;

% Drag (x in feet/3.28)
xvideoforward = [0 .5 1 1.5 2 2.5 3 3.5]/3.28;
tvideoforward = [35.3 35.37 35.43 35.5 35.57 35.63 35.73 35.8];

%%%feet/3.28
xvideoreturn = [3.5 3 2.5 2 1.5 1 .5 0]/3.28;
tvideoreturn = [35.83 36.03 36.3 36.63 36.83 37.1 37.43 37.83];

xvideo = [xvideoforward xvideoreturn];
tvideo = [tvideoforward tvideoreturn];
tvideo = tvideo - tvideo(1);           % Initialize at 0 
xvideo = xvideo - xvideo(1); 

if make_plots
    [h1,ax1] = plottool(1,'Drag Estimation',14,'Time (sec)','x (m)');
    plot(tvideo,xvideo,'b*','MarkerSize',10)
end

xdotvideo = xvideo;
for idx = 1:length(xvideo)-1
    xdotvideo(idx) = (xvideo(idx+1)-xvideo(idx))/(tvideo(idx+1)-tvideo(idx));
end

N = 10;
coeff = polyfit(tvideo,xvideo,N);
t0 = tvideo(1);
tfinal = tvideo(end);
timestep = 0.01;
time = t0:timestep:tfinal;

xest = polyval(coeff,time);
% plot(time,xest,'r--','LineWidth',2)

if make_plots
    [h2,ax2] = plottool(1,'Drag Estimation',14,'Time (sec)','xdot (m/s)');
    plot(tvideo,xdotvideo,'b*','MarkerSize',10)
end

xdotest = (xest(2:end)-xest(1:end-1))./(time(2:end)-time(1:end-1));
xdotest(end+1) = xdotest(end);

% plot(time,xdotest,'r--','LineWidth',2)

x0 = xvideo(1);
xdot0 = mean(xdotvideo(1:4));

state = [x0;xdot0];
StateOUT = zeros(2,length(time));       % Initialize StateOUT vector

for idx = 1:length(time)
  StateOUT(:,idx) = state;
  statedot = Derivatives(state,time(idx));
  state = state + statedot*timestep;
end

xsimulation = StateOUT(1,:);
xdotsimulation = StateOUT(2,:);

if make_plots
    hold on
    plot(ax1,time,xsimulation,'g-','LineWidth',2)
    plot(ax2,time,xdotsimulation,'g-','LineWidth',2)
    legend(ax1,'Experimental Data','Simulated Data')
    legend(ax2,'Experimental Data','Simulated Data')
    del_max*3.28*12
end

xInterp = interp1(time,xsimulation,tvideo); 

cost = sum((xvideo-xInterp).^2)*3.28*12*length(xvideo);

function statedot = Derivatives(state,t)
global K C del_max
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

xwall = (50/12)/3.28; %%%Estimated from the video
% K = 20;
% C=10;
% C = 0.0541;
if x+b/2 > xwall
    delx = (x+b/2-xwall);
    if abs(delx) > del_max
        del_max = delx;
    end
    if xdot > 0
        delxdot = xdot;
    else
        delxdot = 0;
    end
    Spring = -K*delx-C*delxdot;
else
    Spring = 0;
end
   
C_D = 0.125;
rho = 1.225;
X = -0.5*rho*abs(xdot)*Vol^(2/3)*C_D*xdot + Spring;

udot = X/(m+k1);

xdbldot = udot;

statedot = [xdot;xdbldot];
