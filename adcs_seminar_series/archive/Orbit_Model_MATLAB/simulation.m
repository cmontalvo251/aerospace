function [Got_99_problems_and]= simulation()
close all
clear
clc

disp('Orbital Rendezvous Program')

% Constants and Initial Variables
L = 42164*(1000);
T = 86400;
M = 5.972*10^24 ;
G = 6.6742*10^-11;
mu = G*M;
gaia = 6371*1000;
anakin= gaia/L;

x0 = 42164*(1000);
x0hat = x0/L;
xd0 = 0;             
xd0hat = xd0*T/L;

y0 = 0;
y0hat = y0/L;
yd0 = sqrt(M*G/x0);
yd0hat = yd0*T/L;

tinitial = 0;
tfinal = 2;
timestep = 0.01;
time = tinitial:timestep:tfinal;

t1 = 2*pi*sqrt(x0.^3/mu);
phi = pi/2;
tphi = phi*t1/(2*pi);
tmeet = t1 + tphi;
semi = (mu*(tmeet/(2*pi))^2)^(1/3);
ve = sqrt(mu*((2/x0)-(1/semi)));
vehat = ve*T/L;

% Super cool book-keeping method:
satellites(1).mass = 1.33;
satellites(2).mass = 1;
satellites(1).statehat = [x0hat y0hat xd0hat yd0hat]';
satellites(2).statehat = [y0hat x0hat -vehat xd0hat]';

for n = 1:length(satellites)
    satellites(n).muhat = -G*(satellites(n).mass+M)*T^2/L^3;
    satellites(n).statedhat = [xd0hat yd0hat 0 0]';
    satellites(n).stateouthat = zeros(4,length(time));
end

u = zeros(length(time), 1);
rendezvous = 1;

%%% Integration
for idx = 1:length(time)
    for n = 1:length(satellites)
        satellites(n).stateouthat(:,idx)=satellites(n).statehat;   
          u1= Control(satellites(n).statehat,time(idx));
          k1 = Derivatives(satellites(n).statehat,satellites(n).muhat);
          k2 = Derivatives(satellites(n).statehat + k1*timestep/2, satellites(n).muhat);
          k3 = Derivatives(satellites(n).statehat + k2*timestep/2, satellites(n).muhat);
          k4 = Derivatives( satellites(n).statehat + k3*timestep, satellites(n).muhat);
          phi = (1/6)*(k1 + 2*k2 + 2*k3 + k4);
          satellites(n).statehat = satellites(n).statehat + phi*timestep;
          u (idx ) = u1;   
    end
    x1 = satellites(1).statehat(1);
    y1 = satellites(1).statehat(2);
    x2 = satellites(2).statehat(1);
    y2 = satellites(2).statehat(2);
    kitty = atan2(y2,x2);
    if abs(kitty -pi/2) < (1/43) && rendezvous && time(idx) > 0 %change velocity of satellite. Resent velocity xd0hat. %1/43 -pi/2
        satellites(2).statehat(3:4) = satellites(1).statehat(3:4);% [-yd0hat 0];
        rendezvous = 0;
    end
end


%%% Re-Dimensionalize stuff
% Position
x1hat = satellites(1).stateouthat(1,:);
y1hat = satellites(1).stateouthat(2,:);
x1 = x1hat *L;
y1 = y1hat *L;
x2hat = satellites(2).stateouthat(1,:);
y2hat = satellites(2).stateouthat(2,:);
x2 = x2hat*L;
y2 = y2hat*L;
% Velocity
xd1hat = satellites(1).stateouthat(3,:);
yd1hat = satellites(1).stateouthat(4,:);
xd2hat = satellites(2).stateouthat(3,:);
yd2hat = satellites(2).stateouthat(4,:);
xd1 = xd1hat*L/T;
yd1 = yd1hat*L/T;
v1hat = sqrt(xd1hat.^2+yd1hat.^2);
v1 = sqrt(xd1.^2+yd1.^2);
xd2 = xd2hat*L/T;
yd2 = yd2hat*L/T;
v2hat = sqrt(xd2hat.^2+yd2hat.^2);
v2 = sqrt(xd2.^2+yd2.^2);
% Time
time= time*24;


%%% PLOTS

% Phase angle
phi1 =unwrap(atan2(y1hat,x1hat));
phi2 =unwrap(atan2(y2hat,x2hat));
phidel = phi1-phi2;
fig = figure();
set(fig,'color','white');
plot(time, phidel, 'b-','Linewidth',2)
set(gca,'Fontsize',12)
xlabel( 'Time(hours)')
ylabel( 'Difference in Angle')
title('Phase Angle Offset between two Satellites')
xlim([0 60])
grid on

% Velocity
fig = figure();
set(fig, 'color','white')
plot(time, v1hat,'g-','linewidth',2)
set(gca,'Fontsize',18)
xlabel('Time (hours)')
ylabel('Velocity of Satellites (nd)')
title('Velocity')
grid on
xlim([0 60])
hold on
plot(time, v2hat,'r-','linewidth',2)
legend('Velocity of Satellite 1','Velocity of Satellite 2')

% Simulation
fig = figure();
A = imread('earth.png');
set(fig,'color','white');
sky = imread('world-map.jpg');
[X,Y,Z] = ellipsoid(0,0,0,anakin,anakin,anakin);
for n = 1:length(x1hat)
    cla 
    hold on 
    surf(X,Y,Z);
    h = findobj('Type','surface');
    set(h,'CData',flipdim(sky,1),'FaceColor','texturemap','edgecolor','none','FaceLighting','Gouraud','Clipping','off')
    colors = {'k-','m--'};
    for m = length(satellites):-1:1
        xhat = satellites(m).stateouthat(1,:);
        yhat = satellites(m).stateouthat(2,:);
        plot(xhat, yhat, colors{m},'Linewidth',1)
        CubeDraw(.15,.15,.15,xhat(n), yhat(n),0,0,0,0,[1 0 0])
    end
    view(-180,90)
    axis equal 
    axis ([-2. 2. -2. 2.])
    set(gca,'Fontsize',14)
    set(gca,'ZTickLabel',{'-0.1' '' '0.1'})
    xlabel( '$\hat{X}$(nd)','interpreter','LaTeX')
    ylabel( '$\hat{Y}$(nd)','interpreter','LaTeX')
    title('Unit Position')
    grid on
    pause(.1)
    f = getfilename(n,5);
    saveas(gcf,['Frames/Frame_',f,'.jpg'])
end 
system('cd Frames/; makemovie 30 ~/Desktop/Movie.avi jpg');