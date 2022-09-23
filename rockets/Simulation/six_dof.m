%6 dof simulation
%3/3/2015 - Created by Andrew Tindell
%8/18/2015 - Edited by Kent Lino
%1/2/2016 - Edited By Carlos Montalvo
clc,clear
close all

global m Iyy Ixx T

saveplots = 0;

%givens
x = 0; 
y = 0;
z = 0;
u = 0; 
v = 0;
w = 0;

phi = 0*pi/180;
theta = 45*pi/180; 
psi = 0*pi/180;

p = 0;
q = 0;
r = 0;
%Hr = 0;5

x0 = [x, y, z, phi, theta, psi, u, v, w, p, q, r];

dt = 5*10^-3; %1*10^-5
t = 0;
tfinal = 25.8; %0.023
nsteps = round((tfinal-t)/dt);

%begin
k = 1;
out = zeros(nsteps,length(x0));
out(k,:) = x0;
zdot= zeros (nsteps,1);
xdot0 = Derivs_6dof(x0,0);
zdot(1)= xdot0(3);
latvel=zeros (nsteps,1);
latvel(1)=xdot0(1);
t = linspace(0,tfinal,nsteps);
mass_usa(1) = m;
Iyy_vec(1) = Iyy;
Ixx_vec(1) = Ixx;
thrust_curve(1) = T;
for k = 1:nsteps-1
    x = out(k,:);
    xdot1 = Derivs_6dof(x,t(k));
    xdot2 = Derivs_6dof(x+0.5*xdot1*dt,t(k));
    xdot3 = Derivs_6dof(x+0.5*xdot2*dt,t(k));
    xdot4 = Derivs_6dof(x+xdot3*dt,t(k));
    xdotRK4 = (1/6)*(xdot1 + 2*xdot2 + 2*xdot3 + xdot4);
    out(k+1,:) = out(k,:) + xdotRK4*dt;
    zdot(k+1)= xdotRK4(3); %vertical velocity
    latvel(k+1) = xdotRK4(1); %lateral velocity...?
    if abs(mod(t(k),1))<dt
      % disp(['T = ',num2str(t(k)),' out of ',num2str(t(end))])
    end
    mass_usa(k+1) = m;
    Iyy_vec(k+1) = Iyy;
    Ixx_vec(k+1) = Ixx;
    thrust_curve(k+1) = T;
end

C=180/pi;

x = out(:,2);
y = out(:,1);
z = out(:,3);
phi = out(:,4)*C;
theta = out(:,5)*C;
psi = out(:,6)*C;
u=out(:,7);
v = out(:,8);
w = out(:,9);
p = out(:,10)*C;
q = out(:,11)*C;
r = out(:,12)*C;

%%%%Process data
filename = '../Flight_Data/cmontalvo1';
%system(['./process_flight_data.py ',filename]);
data = csvread([filename,'_processed.csv']);

time = data(:,1); 
position_east_of_launch = data(:,2);
lateral_distance = data(:,4);

latitude = data(:,5); 
longitude = data(:,6);
[xlat,ylat] = convertLATLON(latitude,longitude,[latitude(1),longitude(1)]); %% Convert to cartesian
xlat = xlat.*3.28;
ylat = ylat.*3.28;

figure
hold on
plot(t,y, 'LineWidth', 2)       %%%%lateral distance in openrocket
xlabel('Time (sec)','FontSize',18), ylabel('East Coordinate (ft)','FontSize',18)
plot(time,position_east_of_launch, 'r--','LineWidth',2)
%plot(time,ylat,'g--','LineWidth',2)
%plot(time,lateral_distance,'g--','LineWidth',2)
legend ('FASTrocket','OpenRocket')
grid on

position_north_of_launch = data(:,3);

figure
hold on
plot(t,x, 'LineWidth', 2)           %%%should stay zero or close to it
xlabel('Time (sec)','FontSize',18), ylabel('North Coordinate (ft)','FontSize',18)
plot(time,position_north_of_launch, 'r--','LineWidth',2)
%plot(time,xlat,'g--','LineWidth',2)
legend ('FASTrocket','OpenRocket')
grid on

altitude = data(:,7); %altitude

figure
hold on                         %altitude
plot(t,-z,'b-', 'LineWidth', 2)
xlabel('Time (sec)','FontSize',18), ylabel('Altitude (ft)','FontSize',18)
plot(time,altitude,'r--','LineWidth',2)
legend ('FASTrocket', 'OpenRocket')
grid on

vertical_velocity = data(:,8);  %vertical velocity

figure
hold on                   %vertical velocity
plot(t,-zdot,'b-', 'LineWidth', 2)
xlabel('Time (sec)','FontSize',18), ylabel('Vertical Velocity (ft/s)','FontSize',18)
plot(time,vertical_velocity,'r--','LineWidth',2)
legend ('FASTrocket', 'OpenRocket')
grid on

lateral_velocity = data(:,9);

figure
hold on                   %lateral velocity
plot(t,latvel,'b-', 'LineWidth', 2)
xlabel('Time (sec)','FontSize',18), ylabel('Translational Velocity (ft/s)','FontSize',18)
plot(time,lateral_velocity,'r--','LineWidth',2)
legend ('FASTrocket', 'OpenRocket')
grid on

total_velocity = data(:,10);

figure
hold on
plot(t,sqrt(u.^2+v.^2+w.^2),'b-', 'LineWidth', 2)
xlabel('Time (sec)','FontSize',18), ylabel('Total Velocity (ft/s)','FontSize',18)
plot(time,total_velocity,'r--','LineWidth',2)
legend ('FASTrocket', 'OpenRocket')
grid on

pitch_angle = data(:,11);

figure
hold on
plot(t, theta, 'b-','LineWidth', 2)
xlabel('Time (sec)','FontSize',18), ylabel('Pitch Angle (deg)','FontSize',18)
plot(time,pitch_angle,'r--','LineWidth',2)
legend ('FASTrocket','OpenRocket')
grid on

yaw_angle = data(:,12);

figure
hold on
plot(t, psi, 'LineWidth', 2)
xlabel('Time (sec)','FontSize',18), ylabel('Yaw Angle (deg)','FontSize',18)
plot(time,yaw_angle,'r--','LineWidth',2)
legend ('FASTrocket','OpenRocket')
grid on

roll_rate = data(:,13)*C;

figure
hold on
plot(t,p,'b-', 'LineWidth', 2)
xlabel('Time (sec)','FontSize',18), ylabel('p (deg/s)','FontSize',18)
plot(time,roll_rate,'r--','LineWidth',2)
legend ('FASTrocket', 'OpenRocket')
grid on

pitch_rate = data(:,14)*C;

figure
hold on
plot(t,q,'b-', 'LineWidth', 2)
xlabel('Time (sec)','FontSize',18), ylabel('q (deg/s)','FontSize',18)
plot(time,pitch_rate,'r--','LineWidth',2)
legend ('FASTrocket', 'OpenRocket')
grid on

yaw_rate = data(:,15)*C;

figure
hold on
plot(t,r, 'LineWidth', 2)
xlabel('Time (sec)','FontSize',18), ylabel('r(deg/s)','FontSize',18)
plot(time,yaw_rate,'r--','LineWidth',2)
legend ('FASTrocket','OpenRocket')
grid on

angle_of_attack = data(:,17);
AOA = atan2(sqrt(v.^2+w.^2),u);

figure
hold on
plot(t,AOA, 'LineWidth', 2)
xlabel('Time (sec)','FontSize',18), ylabel('Angle of Attack (deg)','FontSize',18)
plot(time,angle_of_attack,'r--','LineWidth',2)
legend ('FASTrocket','OpenRocket')
grid on

mass = data(:,18);

figure()
hold on
plot(t,mass_usa*32.2*16,'LineWidth',2)
xlabel('Time (sec)','FontSize',18), ylabel('Mass of Rocket (oz)','FontSize',18)
plot(time,mass,'r--','LineWidth',2)
legend('FASTrocket','OpenRocket')

longitudinal_moment_inertia = data(:,19);

figure()
hold on
plot(t,Iyy_vec*32.2,'LineWidth',2)
xlabel('Time (sec)','FontSize',18), ylabel('Longitudinal MOI (lbf-ft^2)','FontSize',18)
plot(time,longitudinal_moment_inertia,'r--','LineWidth',2)
legend('FASTrocket','OpenRocket')

rotational_moment_inertia = data(:,20);

figure()
hold on
plot(t,Ixx_vec*32.2,'LineWidth',2)
xlabel('Time (sec)','FontSize',18), ylabel('Rotational MOI (lbf-ft^2)','FontSize',18)
plot(time,rotational_moment_inertia,'r--','LineWidth',2)
legend('FASTrocket','OpenRocket')

thrust = data(:,21);

figure()
hold on
plot(t,thrust_curve,'LineWidth',2)
xlabel('Time (sec)','FontSize',18), ylabel('Thrust (lbf)','FontSize',18)
plot(time,thrust/4.44,'r--','LineWidth',2)
legend('FASTrocket','OpenRocket')

figure()
hold on
plot(t,thrust_curve*4.44,'LineWidth',2)
xlabel('Time (sec)','FontSize',18), ylabel('Thrust (N)','FontSize',18)
plot(time,thrust,'r--','LineWidth',2)
legend('FASTRocket','OpenRocket')

if saveplots
   system('rm Frames/*.jpg');
   system('rm Frames/*.pdf');
   system('rm Frames/*.png');
   ii = 50;
   while ~isempty(get(0,'CurrentFigure') )
       f = getfilename(ii,5);
       saveas(gcf,['Frames/Frame_',f,'.eps'],'epsc')
       close(gcf)
       ii = ii - 1;
   end
   disp('Converting All Pdfs to State.pdf')
   system('convert Frames/*.eps State.pdf')
end

% figure
% plot3(x,-y,-z, 'LineWidth', 2)
% title('x, y, z Coordinates')
% xlabel('x (ft)','FontSize',18), ylabel('y (ft)','FontSize',18), zlabel('z (ft)','FontSize',12)
% grid on

% figure
% hold on
% plot(t, phi, 'LineWidth', 2)
% title('Roll angle','FontSize',18)
% xlabel('Time (sec)','FontSize',12), ylabel('phi(deg)','FontSize',12)
% grid on