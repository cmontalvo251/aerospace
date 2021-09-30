close all
clear
clc

global uk WPctr

WPctr = 1;

%%%Time vector
dt = 0.05;
time = 0:dt:200;

%%%Initial Conditions
xstate = zeros(6,length(time)); %x,z,phi,xdot,zdot,p
ustate = zeros(4,length(time)); %r1 and r2 (0-100)

x0 = [0 0 0 0 0 0]'; %%On the ground chilling

ANIMATE = 1;
QMODE = 'AUTO'; %%%choices are RATE_HOLD,RATE,STD,AUTO

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
    uk = Control(xk,tk,QMODE,uk);

    %%%%SEND CURRENT STATE TO ANIMATION
    if ANIMATE
        quad_animation(xk,uk,tk);
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
phi = xstate(3,:);
xdot = xstate(4,:);
zdot = xstate(5,:);
p = xstate(6,:);
r1 = ustate(1,:);
r2 = ustate(2,:);

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
plot(time,phi*180/pi)
xlabel('Time (sec)')
ylabel('\phi (deg)')

fig = figure;
set(gca,'FontSize',18)
plot(time,xdot,'b-','LineWidth',2)
hold on
xlabel('Time (sec)')
ylabel('Xdot (m/s)')

fig = figure;
set(gca,'FontSize',18)
plot(time,zdot,'b-','LineWidth',2)
hold on
xlabel('Time (sec)')
ylabel('Zdot (m/s)')

fig = figure;
set(gca,'FontSize',18)
plot(time,p,'b-','LineWidth',2)
hold on
xlabel('Time (sec)')
ylabel('P (rad/s)')

fig = figure;
plot(x,z)
xlabel('X (m)')
ylabel('Z (m)')


fig = figure;
set(gca,'FontSize',18)
plot(time,r1,'b-','LineWidth',2)
hold on
xlabel('Time (sec)')
ylabel('Left Percent (0-100%)')


fig = figure;
set(gca,'FontSize',18)
plot(time,r2,'b-','LineWidth',2)
hold on
xlabel('Time (sec)')
ylabel('Right Percent (0-100%)')

%%%Output Data to text file
dlmwrite('Quad_Data.txt',[time' xstate' ustate'])
