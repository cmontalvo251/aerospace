%%%%Simulation of a Low Earth Satellite
disp('Simulation Started')

%%%Get Planet Parameters
planet

%%%Initial Conditions
altitude = 10*254*1.6*1000; %%meters
x0 = R + altitude;
y0 = 0;
z0 = 0; 
xdot0 = 0;
inclination = 51.6*pi/180;
semi_major = norm([x0;y0;z0]);
vcircular = sqrt(mu/semi_major);
ydot0 = vcircular*cos(inclination);
zdot0 = -vcircular*sin(inclination);
stateinitial = [x0;y0;z0;xdot0;ydot0;zdot0];

%%%Need time window
period = 2*pi/sqrt(mu)*semi_major^(3/2);
number_of_orbits = 1;
tspan = [0 period*number_of_orbits];

%%%This is where we integrate the equations of motion
[tout,stateout] = ode45(@Satellite,tspan,stateinitial);

%%%Convert state to kilometers
stateout = stateout/1000;

%%%Extract the state vector
xout = stateout(:,1);
yout = stateout(:,2);
zout = stateout(:,3);

%%%Make an Earth
[X,Y,Z] = sphere(100);
X = X*R/1000;
Y = Y*R/1000;
Z = Z*R/1000;

%%%Plot 3D orbit
fig = figure();
set(fig,'color','white')
plot3(xout,yout,zout,'b-','LineWidth',4)
grid on
hold on
surf(X,Y,Z,'EdgeColor','none')
axis equal