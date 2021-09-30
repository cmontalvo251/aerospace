close all
clear
clc

global uk WPctr

WPctr = 1;

%%%Time vector
dt = 0.03;
time = 0:dt:200;

%%%Initial Conditions
xstate = zeros(6,length(time)); %x,z,theta,u,w,q
ustate = zeros(2,length(time)); %dt, de

x0 = [0 0 0 0 0 0]'; %%On the ground chilling

ANIMATE = 1;
AMODE = 'AUTO'; %%%choices are MANUAL

if ANIMATE
    animation_setup();
    pause
end

xstate(:,1) = x0;

%%%Runge_Kutta
for k = 1:length(time)-1
    xk = xstate(:,k);
    tk = time(k);

    %%%%COMPUTE CONTROL
    uk = Control(xk,tk,AMODE,uk);

    %%%%SEND CURRENT STATE TO ANIMATION
    if ANIMATE
        aircraft_animation(xk,uk,tk);
        %pause(0.1)
    end

    %%%%%DERIVATIVES CALLS
    k1 = Derivatives(xk,tk,uk);
    k2 = Derivatives(xk+k1*dt/2,tk+dt/2,uk);
    k3 = Derivatives(xk+k2*dt/2,tk+dt/2,uk);
    k4 = Derivatives(xk+k3*dt,tk+dt,uk);
    phi = (1/6)*(k1 + 2*k2 + 2*k3 + k4);
    
    %%%%SAVE XSTATE AND USTATE
    xstate(:,k+1) = xstate(:,k) + phi*dt;
    ustate(:,k+1) = uk;
end

%%%EXTRACT STATES
x = xstate(1,:);
z = xstate(2,:);
theta = xstate(3,:);
u = xstate(4,:);
w = xstate(5,:);
q = xstate(6,:);
dt = ustate(1,:);
de = ustate(2,:);

%%%PLOT STATES
fig = figure;
plot(time,x,'b-','LineWidth',2)
hold on
xlabel('Time (sec)')
ylabel('X (m)')
%plot(tuvideo,xvideo,'rs')



figure()
%plot(tvideo,yvideo,'rs')
hold on
plot(time,z)
xlabel('Time (sec)')
ylabel('Z (m)')

figure()
%plot(tpvideo,psivideo*180/pi,'rs')
hold on
plot(time,theta*180/pi)
xlabel('Time (sec)')
ylabel('\theta (deg)')

fig = figure;
set(gca,'FontSize',18)
plot(time,u,'b-','LineWidth',2)
hold on
xlabel('Time (sec)')
ylabel('u (m/s)')

fig = figure;
set(gca,'FontSize',18)
plot(time,w,'b-','LineWidth',2)
hold on
xlabel('Time (sec)')
ylabel('w (m/s)')

fig = figure;
set(gca,'FontSize',18)
plot(time,q,'b-','LineWidth',2)
hold on
xlabel('Time (sec)')
ylabel('Q (rad/s)')

fig = figure;
plot(x,z)
xlabel('X (m)')
ylabel('Z (m)')


fig = figure;
set(gca,'FontSize',18)
plot(time,dt,'b-','LineWidth',2)
hold on
xlabel('Time (sec)')
ylabel('Thrust Percent (0-100%)')


fig = figure;
set(gca,'FontSize',18)
plot(time,de*180/pi,'b-','LineWidth',2)
hold on
xlabel('Time (sec)')
ylabel('Elevator Deflection Angle (deg)')

%%%Output Data to text file
dlmwrite('Aircrat_Data.txt',[time' xstate' ustate'])
