purge
load Flight_Data/altimeter.mat

timestep = 0.1;
time = 0:timestep:(timestep*(length(Pressure)-1));

%%%%Launch Time
pressure = Pressure;
tlaunch = 3360;
tland = 3500;
s = find(time > tlaunch,1);
e = find(time > tland,1);
time = time(s:e);
pressure = Pressure(s:e);

figure()
plot(time,pressure)
xlabel('Time (sec)')
ylabel('Pressure (mb)')

%%%Convert to height in feet
height = 3.28*44330.0*(1.0 - (pressure/pressure(1)).^0.1903);

figure()
plot(time,height)
xlabel('Time (sec)')
ylabel('Height (ft)')

%%%First Order Derivative to get Velocity
velocity = (height(2:end)-height(1:end-1))./timestep;

%%%Throw it through a complimentary filter
velocity_filtered = Complimentary(velocity,0.95);

figure
plot(time(1:end-1),velocity)
hold on
plot(time(1:end-1),velocity_filtered,'r-')
xlabel('Time (sec)')
ylabel('Velocity (ft/s)')
legend('First Order Differentiation','Complimentary Filter')
