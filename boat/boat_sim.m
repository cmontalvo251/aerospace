close all
clear
clc

global uk WPctr

WPctr = 1;

%%%%USE THIS FOR C_V
%tvideo = [28/30 1+29/30 3+27/30 5+12/30]-28/30;
%yvideo = ([1 2 3 3.5]-1)./3.28;
%ydot0 = 1.2*(yvideo(2)-yvideo(1))/(tvideo(2)-tvideo(1))

%%%%USE THIS FOR C_R
%tpvideo = [1+22/30 2+7/30 3 4 5+4/30]-(1+22/30);
%psivideo = [0  45 90 135 180]*pi/180;
%psidot0 = 1.2*(psivideo(2)-psivideo(1))/(tpvideo(2)-tpvideo(1))

%%%%USE THIS FOR C_U
%tuvideo = [5+6/30 5+24/30 6+12/30 7+6/30 8+4/30]-(5+6/30);
%xvideo = [0 1 2 3 4]./3.28;
%xdot0 = 1.2*(xvideo(2)-xvideo(1))/(tuvideo(2)-tuvideo(1))

%%%Time vector
dt = 0.03;
time = 0:dt:200;

%%%Initial Conditions
xstate = zeros(6,length(time)); %x,y,psi,u,v,r
ustate = zeros(2,length(time));

x0 = [0 0 0 0 0 0]';
%xeq = [0 0 0 4.278 0 0]';

ANIMATE = 1;
AUTONOMOUS = 0;

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
    if AUTONOMOUS
        uk = Control(xk,tk);
    end

    %%%%SEND CURRENT STATE TO ANIMATION
    if ANIMATE
        boat_animation(xk,uk,tk);
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
y = xstate(2,:);
psi = xstate(3,:);
u = xstate(4,:);
v = xstate(5,:);
r = xstate(6,:);
dt = ustate(1,:);
dr = ustate(2,:);

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
plot(time,y)
xlabel('Time (sec)')
ylabel('Y (m)')

figure()
%plot(tpvideo,psivideo*180/pi,'rs')
hold on
plot(time,psi*180/pi)
xlabel('Time (sec)')
ylabel('\psi (deg)')

fig = figure;
set(gca,'FontSize',18)
plot(time,u,'b-','LineWidth',2)
hold on
xlabel('Time (sec)')
ylabel('U (m/s)')

fig = figure;
set(gca,'FontSize',18)
plot(time,v,'b-','LineWidth',2)
hold on
xlabel('Time (sec)')
ylabel('V (m/s)')

fig = figure;
set(gca,'FontSize',18)
plot(time,r,'b-','LineWidth',2)
hold on
xlabel('Time (sec)')
ylabel('R (rad/s)')

fig = figure;
plot(x,y)
xlabel('X (m)')
ylabel('Y (m)')


fig = figure;
set(gca,'FontSize',18)
plot(time,dt,'b-','LineWidth',2)
hold on
xlabel('Time (sec)')
ylabel('Thrust (0-100%)')


fig = figure;
set(gca,'FontSize',18)
plot(time,dr*180/pi,'b-','LineWidth',2)
hold on
xlabel('Time (sec)')
ylabel('\delta_R (deg)')

%%%Output Data to text file
dlmwrite('Boat_Data.txt',[time' xstate' ustate'])
