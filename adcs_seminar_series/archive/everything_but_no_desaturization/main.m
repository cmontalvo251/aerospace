%%Initialize
clear
clc
close all

tic

%%%Globals
global BB m I Is invI mu lastMagUpdate nextMagUpdate lastSensorUpdate 
global nextSensorUpdate BfieldMeasured pqrMeasured ptpMeasured BfieldNav pqrNav ptpNav
global BfieldNavPrev BfieldCtlPrev pqrNavPrev ptpNavPrev current Ir1Bcg Ir2Bcg Ir3Bcg n1 n2 n3
global maxSpeed maxAlpha Ir1B Ir2B Ir3B rwalphas Bdot
global fsensor MagFieldBias AngFieldBias EulerBias R Amax lmax CD
global MagFieldNoise AngFieldNoise EulerNoise IrR Jinv

%%%%Simulation of a Low Earth Satellite
disp('Simulation Started')

%%%Setup the IGRF Model
disp(['You must download the igrf model from MathWorks in order to' ...
      ' use this software'])
disp('https://www.mathworks.com/matlabcentral/fileexchange/34388-international-geomagnetic-reference-field-igrf-model')
addpath '../igrf/'
nextMagUpdate = 1;
lastMagUpdate = 0;

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
p0 = 0.8;
q0 = -0.2;
r0 = 0.3;
%%%Initial conditions of my reaction wheels
w10 = 0;
w20 = 0;
w30 = 0;

state = [x0;y0;z0;xdot0;ydot0;zdot0;q0123_0;p0;q0;r0;w10;w20;w30];

%%%Need time window
period = 2*pi/sqrt(mu)*semi_major^(3/2);
number_of_orbits = 1;
tfinal = period*number_of_orbits;
%tfinal = 100;
next = 10;
timestep = 0.1;
tout = 0:timestep:tfinal;
stateout = zeros(length(tout),length(state));

%%%Create empty array for simulation
BxBout = 0*stateout(:,1);
ByBout = BxBout;
BzBout = BxBout;
BxBm = 0*stateout(:,1);
ByBm = BxBout;
BzBm = BxBout;
pqrm = zeros(length(tout),3);
ptpm = zeros(length(tout),3);
ptpN = 0*ptpm;
BxBN = 0*stateout(:,1);
ByBN = BxBout;
BzBN = BxBout;
pqrN = zeros(length(tout),3);
ix = 0*stateout(:,1);
iy = ix;
iz = ix;
rwa = 0*ptpm;
BfieldNavPrev = [-99;0;0];
pqrNavPrev = [0;0;0];
ptpNavPrev = [0;0;0];
Bdot = [0;0;0];

%%%Sensor Parameters
lastSensorUpdate = 0;
sensor_params

%%%%Call the Derivatives Routine to initialize variables
k1 = Satellite(tout(1),state);

%%%Print Next
lastPrint = 0;

%%%%Control Parameters
lastControl = -timestep;
nextControl = 0.1;

%%%Loop through time to integrate
for idx = 1:length(tout)
    %%%Save the current state
    stateout(idx,:) = state';
    
    %%%Save the Current
    ix(idx) = current(1);
    iy(idx) = current(2);
    iz(idx) = current(3);
    
    %%%%Save reaction wheel acceleration
    rwa(idx,:) = rwalphas';
    
    %%%Save the magnetic field
    BxBout(idx) = BB(1);
    ByBout(idx) = BB(2);
    BzBout(idx) = BB(3);
    %%%Save the Measured Magnetic field
    BxBm(idx) = BfieldMeasured(1);
    ByBm(idx) = BfieldMeasured(2);
    BzBm(idx) = BfieldMeasured(3);
    %%%Save the Nav magnetic field
    BxBN(idx) = BfieldNav(1);
    ByBN(idx) = BfieldNav(2);
    BzBN(idx) = BfieldNav(3);
    %%%THe actual pqr truth signal is embedded in The state vector. Save
    %%%the measured pqr signal
    pqrm(idx,:) = pqrMeasured';
    %%%Save the Nav pqr signal
    pqrN(idx,:) = pqrNav';
    %%%Save ptp
    ptpm(idx,:) = ptpMeasured';
    ptpN(idx,:) = ptpNav';
    
    if tout(idx) > lastPrint
        disp(['Time = ',num2str(tout(idx)),' out of ',num2str(tfinal)])
        lastPrint = lastPrint + next;
    end
    
    %%%CONTROL BLOCK
    if tout(idx) > lastControl
        [current,rwalphas] = Control(BfieldNav,pqrNav,ptpNav);
        lastControl = lastControl + nextControl;
    end
    
    %%%This is where we integrate the equations of motion and
    %%%%make our 4 function calls for the RK4
    k1 = Satellite(tout(idx),state);
    k2 = Satellite(tout(idx)+timestep/2,state+k1*timestep/2);
    k3 = Satellite(tout(idx)+timestep/2,state+k2*timestep/2);
    k4 = Satellite(tout(idx)+timestep,state+k3*timestep);
    k = (1/6)*(k1 + 2*k2 + 2*k3 + k4);
    state = state + k*timestep;
    
end

%%%Save original State
stateout_original = stateout;

%%%
disp('Simulation Complete')
toc
tic
disp('Initiation plotting routine....')

%%Convert state to kilometers
stateout(:,1:6) = stateout_original(:,1:6)/1000;

%%%Extract the state vector
xout = stateout(:,1);
yout = stateout(:,2);
zout = stateout(:,3);
q0123out = stateout(:,7:10);
ptpout = Quaternions2EulerAngles(q0123out);
pqrout = stateout(:,11:13);
w123 = stateout(:,14:16);

%%%Make an Earth
[X,Y,Z] = sphere(100);
X = X*R/1000;
Y = Y*R/1000;
Z = Z*R/1000;

%%%Plot X,Y,Z as a function of time
fig0 = figure();
set(fig0,'color','white')
plot(tout,sqrt(xout.^2+yout.^2+zout.^2)-R/1000,'b-','LineWidth',2)
%hold on
grid on
%plot(tout,yout,'r-','LineWidth',2)
%plot(tout,zout,'g-','LineWidth',2)
xlabel('Time (sec)')
ylabel('Altitude (km)')
%legend('X','Y','Z')

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
p1 = plot(tout,BxBout,'b-','LineWidth',2);
hold on
grid on
p2 = plot(tout,ByBout,'r-','LineWidth',2);
p3 = plot(tout,BzBout,'g-','LineWidth',2);
p1m = plot(tout,BxBm,'b-s','LineWidth',2);
p2m = plot(tout,ByBm,'r-s','LineWidth',2);
p3m = plot(tout,BzBm,'g-s','LineWidth',2);
p1N = plot(tout,BxBN,'b--','LineWidth',2);
p2N = plot(tout,ByBN,'r--','LineWidth',2);
p3N = plot(tout,BzBN,'g--','LineWidth',2);
legend([p1,p2,p3,p1m,p2m,p3m,p1N,p2N,p3N],'Bx','By','Bz','Bx Measured','By Measured','Bz Measured','Bx Nav','By Nav','Bz Nav')
xlabel('Time (sec)')
ylabel('Mag Field (T)')

%%%And Norm
Bnorm = sqrt(BxBout.^2 + ByBout.^2 + BzBout.^2);
fig3 = figure();
set(fig3,'color','white')
plot(tout,Bnorm,'LineWidth',2)
xlabel('Time (sec)')
ylabel('Norm of Magnetic Field (T)')
grid on

%%%plot Euler Angles
fig4 = figure();
set(fig4,'color','white')
p1 = plot(tout,ptpout*180/pi,'-','LineWidth',2);
hold on
p2 = plot(tout,ptpm*180/pi,'-s','LineWidth',2);
p3 = plot(tout,ptpN*180/pi,'--','LineWidth',2);
grid on
xlabel('Time (sec)')
ylabel('Euler Angles (deg)')
legend('Phi','Theta','Psi')
legend([p1(1),p2(1),p3(1)],'Actual','Measured','Nav')

%%%Plot Angular Velocity
fig5 = figure();
set(fig5,'color','white')
p1=plot(tout,pqrout,'-','LineWidth',2);
hold on
p2=plot(tout,pqrm,'-s','LineWidth',2);
p3=plot(tout,pqrN,'--','LineWidth',2);
grid on
xlabel('Time (sec)')
ylabel('Angular Velocity (rad/s)')
legend([p1(1),p2(1),p3(1)],'Actual','Measured','Nav')

%%%Plot the current in the magnetorquers
fig6 = figure();
set(fig6,'color','white')
plot(tout,ix*1000,'LineWidth',2)
hold on
plot(tout,iy*1000,'LineWidth',2)
plot(tout,iz*1000,'LineWidth',2)
grid on
xlabel('Time (sec)')
ylabel('Current Magnetorquers (mA)')
legend('X','Y','Z')

%%%Plot the Total current in magnetorquers
mag_current_total = abs(ix)+abs(iy)+abs(iz);
fig6 = figure();
set(fig6,'color','white')
plot(tout,(mag_current_total)*1000,'LineWidth',2)
grid on
xlabel('Time (sec)')
ylabel('Total Current Magnetorquers (mA)')

%%%Plot the acceleration in reaction wheels
fig10 = figure();
set(fig10,'color','white')
plot(tout,rwa,'LineWidth',2)
grid on
xlabel('Time (sec)')
ylabel('Angular Acceleration RWs (rad/s^2)')
legend('X','Y','Z')

%%%Plot the angular velocity of the RWs
fig7 = figure();
set(fig7,'color','white')
plot(tout,w123,'LineWidth',2)
grid on
xlabel('Time (sec)')
ylabel('Angular Velocity of RWs (rad/s)')
legend('X','Y','Z')

%%%Plot the current in the reaction wheels
current_rwa = Amps2Alpha * rwa; 
fig11 = figure();
set(fig11,'color','white')
plot(tout,current_rwa*1000,'LineWidth',2)
grid on
xlabel('Time (sec)')
ylabel('Current RWs (mA)')
legend('X','Y','Z')

%%%Plot the current in the reaction wheels
rw_current_total = abs(current_rwa(:,1)) + abs(current_rwa(:,2)) + abs(current_rwa(:,3));
fig14 = figure();
set(fig14,'color','white')
plot(tout,rw_current_total*1000,'LineWidth',2)
grid on
xlabel('Time (sec)')
ylabel('Total Current RWs (mA)')

%%%Plot the total current in the entire system 
fig12 = figure();
set(fig12,'color','white')
plot(tout,(mag_current_total + rw_current_total)*1000,'LineWidth',2)
grid on
xlabel('Time (sec)')
ylabel('Total Current Magnetorquers+RWs (mA)')

toc


