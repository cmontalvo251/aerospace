%%%Initialize
clear
clc
close all

tic

%%%Globals
global BB m I invI lastMagUpdate nextMagUpdate lastSensorUpdate 
global nextSensorUpdate BfieldMeasured pqrMeasured BfieldNav pqrNav
global BfieldNavPrev pqrNavPrev

BfieldNavPrev = [0;0;0];
pqrNavPrev = [0;0;0];

%%%%Simulation of a Low Earth Satellite
disp('Simulation Started')

%%%Setup the IGRF Model
disp(['You must download the igrf model from MathWorks in order to' ...
      ' use this software'])
disp('https://www.mathworks.com/matlabcentral/fileexchange/34388-international-geomagnetic-reference-field-igrf-model')
addpath 'igrf/'

%%%Get Planet Parameters
planet

%%%Get mass and inertia properties
inertia

%%%Initial Conditions Position and Velocity
altitude = 600*1000; %%meters
x0 = R + altitude;
y0 = 0;
z0 = 0; 
xdot0 = 0;
inclination = 56*pi/180;
semi_major = norm([x0;y0;z0]);
vcircular = sqrt(mu/semi_major);
ydot0 = vcircular*cos(inclination);
zdot0 = vcircular*sin(inclination);

%%%Intitial Conditions for Attitude and Angular Velocity
phi0 = 0;
theta0 = 0;
psi0 = 0;
ptp0 = [phi0;theta0;psi0];
q0123_0 = EulerAngles2Quaternions(ptp0);
p0 = 0.08;
q0 = -0.02;
r0 = 0.03;

state = [x0;y0;z0;xdot0;ydot0;zdot0;q0123_0;p0;q0;r0];

%%%Need time window
period = 2*pi/sqrt(mu)*semi_major^(3/2);
number_of_orbits = 12;
tfinal = period*number_of_orbits;
timestep = 1;
tout = 0:timestep:tfinal;
stateout = zeros(length(tout),length(state));
%%%This is where we integrate the equations of motion
%[tout,stateout] = ode45(@Satellite,tspan,stateinitial);

%%%Loop through time to integrate
BxBout = 0*stateout(:,1);
ByBout = BxBout;
BzBout = BxBout;
BxBm = 0*stateout(:,1);
ByBm = BxBout;
BzBm = BxBout;
pqrm = zeros(length(tout),3);

BxBN = 0*stateout(:,1);
ByBN = BxBout;
BzBN = BxBout;
pqrN = zeros(length(tout),3);

nextMagUpdate = 10;
lastMagUpdate = 0;

%%%Sensor Parameters
lastSensorUpdate = 0;
sensor_params

%%%Print Next
next = 100;
lastPrint = 0;
for idx = 1:length(tout)
    %%%Save the current state
    stateout(idx,:) = state';
    
    if tout(idx) > lastPrint
        disp(['Time = ',num2str(tout(idx)),' out of ',num2str(tfinal)])
        lastPrint = lastPrint + next;
    end
    
    %%%%Then we make our 4 function calls for the RK4
    k1 = Satellite(tout(idx),state);
    k2 = Satellite(tout(idx)+timestep/2,state+k1*timestep/2);
    k3 = Satellite(tout(idx)+timestep/2,state+k2*timestep/2);
    k4 = Satellite(tout(idx)+timestep,state+k3*timestep);
    k = (1/6)*(k1 + 2*k2 + 2*k3 + k4);
    state = state + k*timestep;
    
    %%%Save the magnetic field
    BxBout(idx) = BB(1);
    ByBout(idx) = BB(2);
    BzBout(idx) = BB(3);
    BxBm(idx) = BfieldMeasured(1);
    ByBm(idx) = BfieldMeasured(2);
    BzBm(idx) = BfieldMeasured(3);
    BxBN(idx) = BfieldNav(1);
    ByBN(idx) = BfieldNav(2);
    BzBN(idx) = BfieldNav(3);
    %%%Save the polluted pqr signal
    pqrm(idx,:) = pqrMeasured';
    %%%Save the filtered pqr signal
    pqrN(idx,:) = pqrNav';
end

disp('Simulation Complete')


%%%Convert state to kilometers
stateout(:,1:6) = stateout(:,1:6)/1000;

%%%Extract the state vector
xout = stateout(:,1);
yout = stateout(:,2);
zout = stateout(:,3);
q0123out = stateout(:,7:10);
ptpout = Quaternions2EulerAngles(q0123out);
pqrout = stateout(:,11:13);

%%%Make an Earth
[X,Y,Z] = sphere(100);
X = X*R/1000;
Y = Y*R/1000;
Z = Z*R/1000;

%%%Plot X,Y,Z as a function of time
fig0 = figure();
set(fig0,'color','white')
plot(tout,xout,'b-','LineWidth',2)
hold on
grid on
plot(tout,yout,'r-','LineWidth',2)
plot(tout,zout,'g-','LineWidth',2)
xlabel('Time (sec)')
ylabel('Position (m)')
legend('X','Y','Z')

%%%Plot 3D orbit
fig = figure();
set(fig,'color','white')
plot3(xout,yout,zout,'b-','LineWidth',4)
xlabel('X')
ylabel('Y')
zlabel('Z')
grid on
hold on
surf(X,Y,Z,'EdgeColor','none')
axis equal

%%%Plot Magnetic Field
fig2 = figure();
set(fig2,'color','white')
plot(tout,BxBout,'b-','LineWidth',2)
hold on
grid on
plot(tout,ByBout,'y-','LineWidth',2)
plot(tout,BzBout,'g-','LineWidth',2)
plot(tout,BxBm,'b--','LineWidth',2)
plot(tout,ByBm,'y--','LineWidth',2)
plot(tout,BzBm,'g--','LineWidth',2)
plot(tout,BxBN,'r-','LineWidth',2)
plot(tout,ByBN,'m-','LineWidth',2)
plot(tout,BzBN,'k-','LineWidth',2)
xlabel('Time (sec)')
ylabel('Mag Field (T)')
legend('X','Y','Z')

%%%And Norm
Bnorm = sqrt(BxBout.^2 + ByBout.^2 + BzBout.^2);
fig3 = figure();
set(fig3,'color','white')
plot(tout,Bnorm,'LineWidth',2)
grid on

%%%plot Euler Angles
fig4 = figure();
set(fig4,'color','white')
plot(tout,ptpout,'LineWidth',2)
grid on
xlabel('Time (sec)')
ylabel('Angles (rad)')

%%%Plot Angular Velocity
fig5 = figure();
set(fig5,'color','white')
plot(tout,pqrout,'LineWidth',2)
hold on
plot(tout,pqrm,'--','LineWidth',2)
plot(tout,pqrN,'LineWidth',2)
grid on
xlabel('Time (sec)')
ylabel('Angular Velocity (rad/s)')

toc


