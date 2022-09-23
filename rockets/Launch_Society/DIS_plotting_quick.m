purge

disp('Loading Data')
load All_the_datas
disp('Loaded')

%%%%Since we did this experiment indoors we won't have any GPS data. Thus I
%%%%will be using LastPrint as the time vector although remind me to give
%%%%you a MATLAB code that will fuse ArduinoTime with GPS time. 
%%%%time_vec = LastPrint-LastPrint(1); 
Time_FS = linspace(0,98,length(Pressure_FS))';
    
skipFS = 10;
skipTM = 100;

%%%PLOTTING

%%%Fix Telemega time
%Time_TM = Time_TM - Time_TM(1);

%%%Filter Magnetometer Data
sig = 0.05;
X_Mag_FS = filter_signal(X_Mag_FS,sig);
Y_Mag_FS = filter_signal(Y_Mag_FS,sig);
Z_Mag_FS = filter_signal(Z_Mag_FS,sig);
MAG_FS = [X_Mag_FS,Y_Mag_FS,Z_Mag_FS];
X_Mag_Shifted = filter_signal(X_Mag_Shifted,sig);
Y_Mag_Shifted = filter_signal(Y_Mag_Shifted,sig);
Z_Mag_Shifted = filter_signal(Z_Mag_Shifted,sig);
MAG_Shifted = [X_Mag_Shifted,Y_Mag_Shifted,Z_Mag_Shifted];
X_Mag_TM = filter_signal(X_Mag_TM,sig);
Y_Mag_TM = filter_signal(Y_Mag_TM,sig);
Z_Mag_TM = filter_signal(Z_Mag_TM,sig);
MAG_TM = [X_Mag_TM,Y_Mag_TM,Z_Mag_TM];

%%%Plot Magnetometer Data
plottool(1,'X-Mag',18,'Time (sec)','X-Magnetometer (uT)');
plot(Time_FS,X_Mag_Shifted,'k--','linewidth',2)
plot(Time_TM,X_Mag_TM*100,'k-','linewidth',2)
xlim([0 100])
legend('FASTsensor','TeleMega')

plottool(1,'Y-Mag',18,'Time (sec)','Y-Magnetometer (uT)');
plot(Time_FS,Y_Mag_Shifted,'k--','linewidth',2)
plot(Time_TM,Y_Mag_TM*100,'k-','linewidth',2)
xlim([0 100])
legend('FASTsensor','TeleMega')

plottool(1,'Z-Mag',18,'Time (sec)','Z-Magnetometer (uT)');
plot(Time_FS,Z_Mag_Shifted,'k--','linewidth',2)
plot(Time_TM,Z_Mag_TM*100,'k-','linewidth',2)
xlim([0 100])
legend('FASTsensor','TeleMega')

%%%Filter Accelerometer Data
sig = 0.05;
X_Accel_FS = filter_signal(X_Accel_FS,sig);
Y_Accel_FS = filter_signal(Y_Accel_FS,sig);
Z_Accel_FS = filter_signal(Z_Accel_FS,sig);
ACCEL_FS = [X_Accel_FS,Y_Accel_FS,Z_Accel_FS];
X_Accel_TM = filter_signal(X_Accel_TM,sig);
Y_Accel_TM = filter_signal(Y_Accel_TM,sig);
Z_Accel_TM = filter_signal(Z_Accel_TM,sig);
ACCEL_TM = [X_Accel_TM,Y_Accel_TM,Z_Accel_TM];

%%%Plot Accelerometer Data
plottool(1,'X-Accel',18,'Time (sec)','X-Accelerometer (ft/s^2)');
plot(Time_FS,X_Accel_FS*3.28,'k--','linewidth',2)
plot(Time_TM,X_Accel_TM*3.28,'k-','linewidth',2)
xlim([0 100])
legend('FASTsensor','TeleMega')

plottool(1,'Y-Accel',18,'Time (sec)','Y-Accelerometer (ft/s^2)');
plot(Time_FS,Y_Accel_FS*3.28,'k--','linewidth',2)
plot(Time_TM,Y_Accel_TM*3.28,'k-','linewidth',2)
xlim([0 100])
legend('FASTsensor','TeleMega')

plottool(1,'Z-Accel',18,'Time (sec)','Z-Accelerometer (ft/s^2)');
plot(Time_FS,Z_Accel_FS*3.28,'k--','linewidth',2)
plot(Time_TM,Z_Accel_TM*3.28,'k-','linewidth',2)
xlim([0 100])
legend('FASTsensor','TeleMega')

%%%filter Rate data
sig_gyro = 0.1;
Roll_Rate_FS = filter_signal(Roll_Rate_FS,sig_gyro);
Pitch_Rate_FS = filter_signal(Pitch_Rate_FS,sig_gyro);
Yaw_Rate_FS = filter_signal(Yaw_Rate_FS,sig_gyro);
RATE_FS = [Roll_Rate_FS,Pitch_Rate_FS,Yaw_Rate_FS];
Roll_Rate_TM = filter_signal(Roll_Rate_TM,sig_gyro);
Pitch_Rate_TM = filter_signal(Pitch_Rate_TM,sig_gyro);
Yaw_Rate_TM = filter_signal(Yaw_Rate_TM,sig_gyro);
RATE_TM = [Roll_Rate_TM,Pitch_Rate_TM,Yaw_Rate_TM];

%%%Plot Angular Rate Data
[fig_roll_rate,ax_roll_rate]=plottool(1,'Roll Rate',18,'Time (sec)','Roll Rate (deg/s)');
plot(Time_FS,-Roll_Rate_FS*(180/pi),'k--','linewidth',2)
% plot(Time_OR,Roll_Rate_OR*(180/pi),'k-','linewidth',2)
plot(Time_TM,Roll_Rate_TM,'k-','linewidth',2)
% plot(Time_FR,Roll_Rate_FR,'g-','linewidth',2)
xlim([0 100])
legend('FASTsensor','TeleMega')
%legend('FASTsensor','OpenRocket','TeleMega','FASTrocket')

plottool(1,'Pitch Rate',18,'Time (sec)','Pitch Rate (deg/s)');
plot(Time_FS,Pitch_Rate_FS*(180/pi),'k--','linewidth',2)
% plot(Time_OR,Pitch_Rate_OR*(180/pi),'k-','linewidth',2)
plot(Time_TM,Pitch_Rate_TM,'k-','linewidth',2)
% plot(Time_FR,Pitch_Rate_FR,'g-','linewidth',2)
xlim([0 100])
legend('FASTsensor','TeleMega')
%legend('FASTsensor','OpenRocket','TeleMega','FASTrocket')

plottool(1,'Yaw Rate',18,'Time (sec)','Yaw Rate (deg/s)');
plot(Time_FS,-Yaw_Rate_FS*(180/pi),'k--','linewidth',2)
% plot(Time_OR,Yaw_Rate_OR*(180/pi),'k-','linewidth',2)
plot(Time_TM,Yaw_Rate_TM,'k-','linewidth',2)
% plot(Time_FR,Yaw_Rate_FR,'g-','linewidth',2)
xlim([0 100])
legend('FASTsensor','TeleMega')
%legend('FASTsensor','OpenRocket','TeleMega','FASTrocket')

%%%%Compute Phi,Theta and Psi using Multiple Methods.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Phi, theta and Psi can be computed using the rate gyro by integrating
%%%the signals. This will be very noise and will most likely require a
%%%filter. I suggest a complimentary filter but we shall see how that goes.
%%%We're gonna need initial conditions. Let's also use trapezoidal rule for
%%%integration to reduce noise
%%%Phi and Theta can be computed using the accelerometer using the
%%%relationship of tranforming z_inertial into the body frame. 
%%%phi_accel = atan2(Y_Accel,Z_Accel);
Phi_FS = atan(Y_Accel_FS./Z_Accel_FS);
Phi_TM = atan(Y_Accel_TM./Z_Accel_TM);
Theta_FS = atan2(-X_Accel_FS,Y_Accel_FS.*sin(Phi_FS) + Z_Accel_FS.*cos(Phi_FS));
Theta_TM = atan2(-X_Accel_TM,Y_Accel_TM.*sin(Phi_TM) + Z_Accel_TM.*cos(Phi_TM));

Theta_FS2 = 180/pi*atan(-X_Accel_FS./(Y_Accel_FS.*sin(Phi_FS) + Z_Accel_FS.*cos(Phi_FS)));
Theta_TM2 = 180/pi*atan(-X_Accel_TM./(Y_Accel_TM.*sin(Phi_TM) + Z_Accel_TM.*cos(Phi_TM)));

loc_FS = (Time_FS<70 & Theta_FS2<-30);
Theta_FS2(loc_FS) = -Theta_FS2(loc_FS);

loc_TM = (Time_TM<70 & Theta_TM2<-30);
Theta_TM2(loc_TM) = -Theta_TM2(loc_TM);

%%%%After running this code a few times I realized that when the pitch
%%%%angle approaches 90 degrees we run into a singularity problem. As such
%%%%I multiplied phi by cos(theta_accel)^2
Phi_FS = Phi_FS.*cos(Theta_FS).^2;
Phi_TM = Phi_TM.*cos(Theta_TM).^2;

%%%Get Roll rate from Derivative Filter
[tfilter,y] = DerivativeFilter(Phi_FS*180/pi,Time_FS,1,0,200);
plot(ax_roll_rate,Time_FS,y,'b-','LineWidth',2)

[tfilter,y] = DerivativeFilter(Phi_TM*180/pi,Time_TM,1,0,200);
plot(ax_roll_rate,Time_TM,y,'r-','LineWidth',2)


%%%It is not possible to compute heading with the accelerometer

%%%%So the magnetometer can be used to obtain psi using
%%%%the dot product relationship between the initial value of the mag field
%%%%and the mag field obtained during flight. Turns out Adafruit only gets
%%%%psi rather than getting phi or theta. I guess that makes sense since a
%%%%dot product only returns one angle

%%%One thing I noticed with the magnetometer is that the sensor took about
%%%0.2 seconds to calibrate. So I'm going to calibrate the magnetometer by
%%%averaging the first 5 seconds of data and then setting the first 5
%%%seconds to that value
loc_FS = find(Time_FS > 5,1);
Z_Mag_AVG_FS = mean(Z_Mag_FS(1:loc_FS));
Y_Mag_AVG_FS = mean(Y_Mag_FS(1:loc_FS));
X_Mag_AVG_FS = mean(X_Mag_FS(1:loc_FS));
X_Mag_FS(1:loc_FS) = X_Mag_AVG_FS;
Y_Mag_FS(1:loc_FS) = Y_Mag_AVG_FS;
Z_Mag_FS(1:loc_FS) = Z_Mag_AVG_FS;

loc_TM = find(Time_TM > 5,1);
Z_Mag_AVG_TM = mean(Z_Mag_TM(1:loc_TM));
Y_Mag_AVG_TM = mean(Y_Mag_TM(1:loc_TM));
X_Mag_AVG_TM = mean(X_Mag_TM(1:loc_TM));
X_Mag_FS(1:loc_TM) = X_Mag_AVG_TM;
Y_Mag_FS(1:loc_TM) = Y_Mag_AVG_TM;
Z_Mag_FS(1:loc_TM) = Z_Mag_AVG_TM;

Psi_FS = atan2(Z_Mag_Shifted .* sin(Phi_FS) - Y_Mag_Shifted .* cos(Phi_FS), X_Mag_Shifted .* cos(Theta_FS) + Y_Mag_Shifted .* sin(Theta_FS) .* sin(Phi_FS) + Z_Mag_Shifted .* sin(Theta_FS) .* cos(Phi_FS));
Psi_TM = atan2(Z_Mag_TM .* sin(Phi_TM) - Y_Mag_TM .* cos(Phi_TM), X_Mag_TM .* cos(Theta_TM) + Y_Mag_TM .* sin(Theta_TM) .* sin(Phi_TM) + Z_Mag_TM .* sin(Theta_TM) .* cos(Phi_TM));

%%%Finally I'm pretty sure this equation assumes a certain magnetic field reference
%%%point so I'm going to calibrate the reading to zero initially
Psi_FS = Psi_FS - Psi_FS(1);
Psi_TM = Psi_TM - Psi_TM(1);

%%%The singularity in pitch causes the yaw angle to go kind of wonky as
%%%well so I've multiplied yaw_mag by cos(theta)^2
Psi_FS = Psi_FS.*cos(Theta_FS).^2;
Psi_TM = Psi_TM.*cos(Theta_TM).^2;

%%%Compute Phi Theta and Psi from integrating Gyro Data
phi_gyro_FS = zeros(length(X_Gyro_FS)-1,1);
theta_gyro_FS = phi_gyro_FS;
psi_gyro_FS = theta_gyro_FS;

phi_gyro_FS(1) = Phi_FS(1);
theta_gyro_FS(1) = Theta_FS(1);
psi_gyro_FS(1) = Psi_FS(1);

for idx = 1:length(Roll_Rate_FS)-1
    phi_gyro_FS(idx+1) = phi_gyro_FS(idx) + 0.5*(Roll_Rate_FS(idx) + Roll_Rate_FS(idx+1))*(Time_FS(idx+1)-Time_FS(idx));
    theta_gyro_FS(idx+1) = theta_gyro_FS(idx) + 0.5*(Pitch_Rate_FS(idx) + Pitch_Rate_FS(idx+1))*(Time_FS(idx+1)-Time_FS(idx));
    psi_gyro_FS(idx+1) = psi_gyro_FS(idx) + 0.5*(Yaw_Rate_FS(idx) + Yaw_Rate_FS(idx+1))*(Time_FS(idx+1)-Time_FS(idx));
end

%%%Plot Phi
%%%FASTsensor Phi is from accelerometer
plottool(1,'Phi',18,'Time (sec)','\phi (deg)');
plot(Time_FS,Phi_FS*180/pi,'k--','LineWidth',2);
plot(Time_TM,Phi_TM*180/pi,'k-','LineWidth',2);
% plot(Time_FR,Phi_FR,'g-','LineWidth',2);
xlim([0 100])
%legend('FASTsensor','TeleMega','FASTrocket')
legend('FASTsensor','TeleMega')

%%%Plot Theta
%%%FASTsensor Theta is from accelerometer
plottool(1,'Theta',18,'Time (sec)','\theta (deg)');
plot(Time_FS,Theta_FS2,'k--','LineWidth',2);
% plot(Time_OR,Theta_OR,'k-','LineWidth',2);
plot(Time_TM,Theta_TM2,'k-','LineWidth',2);
% plot(Time_FR,Theta_FR,'g-','LineWidth',2);
xlim([0 100])
%legend('FASTsensor','OpenRocket','TeleMega','FASTrocket')
legend('FASTsensor','TeleMega')

%%%Plot Psi
%%%FASTsensor Psi is from magnetometer
loc_FS = (Time_FS<20 & Psi_FS<-200*pi/180);
Psi_FS(loc_FS) = 0;

loc_TM = (Time_TM<20 & Psi_TM<-200*pi/180);
Psi_TM(loc_TM) = 0;

plottool(1,'Psi',18,'Time (sec)','\psi (deg)');
plot(Time_FS,Psi_FS*180/pi,'k--','LineWidth',2);
%plot(Time_OR,Psi_OR,'k-','LineWidth',2);
plot(Time_TM,Psi_TM*180/pi,'k-','LineWidth',2);
%plot(Time_FR,Psi_FR,'g-','LineWidth',2);
xlim([0 100])
%legend('FASTsensor','OpenRocket','TeleMega','FASTrocket')
legend('FASTsensor','TeleMega')

%%%%%It should then be possible to stitch the data together using a
%%%%%weighted average technique. Unfortunately though you can see when you
%%%%%plot this that the Z_Gyro signal for some reason had a ton of bias and
%%%%%unfortunately we integrated that bias so it just got worse. This is
%%%%%why magnetometers are used to get rid of that bias. It looks like roll
%%%%%and pitch are in good agreement but even they have some integration
%%%%%bias. I would say for now just use the phi and theta computed by accel
%%%%%and yaw from the mag. With more time I can figure out how to create a
%%%%%filtering algorithm that attempts to measure the bias in the
%%%%%magnetometer and correct for it.

%%% Calculating Height given Barometer Pressure
Po = Ground_Pressure; %Pressure at Ground
R = 8.31447; %Gas Constant
L = 0.0065; %Temperature Lapse Rate
g = 9.81; %Acceleration due to gravity
M = 0.0289644; %Molar Mass of Dry Air
To = Temperature_FS(1)+273.13; %Temperature at Ground

Altitude_FS = zeros(length(Pressure_FS),1);

for iter = 1:length(Pressure_FS)
    Altitude_FS(iter) = (1 - (Pressure_FS(iter)/Po)^((R*L/(g*M))))*(To/L);
end

%%%Calculate vertical velocity profile through differentiating pressure 
%%%curve
Vertical_Velocity_FS = zeros(length(Pressure_FS),1);

for iter = 1:length(Pressure_FS)-1
    Vertical_Velocity_FS(iter+1) = (Altitude_FS(iter+1)-Altitude_FS(iter))/(Time_FS(iter+1)-Time_FS(iter));
end

Peak = find(abs(Vertical_Velocity_FS) > 1000);
Vertical_Velocity_FS(Peak) = Vertical_Velocity_FS(Peak(1)-1);

%%%Filter Vertical Velocity
sig = 0.12;
Vertical_Velocity_FS = filter_signal(Vertical_Velocity_FS,sig);

%%%Plot Altitude
plottool(1,'Altitude',18,'Time (sec)','Altitude (ft)');
plot(Time_FS,Altitude_FS*3.28,'k-','linewidth',2)
plot(Time_FS(1:skipFS:end),Altitude_FS(1:skipFS:end)*3.28,'ks','linewidth',2)
p1 = plot(101+Time_FS(1:skipFS:end),Altitude_FS(1:skipFS:end)*3.28,'k-s','linewidth',2);

p2 = plot(Time_OR,Altitude_OR,'k--','linewidth',2);

plot(Time_TM,Altitude_TM*3.28,'k--','linewidth',2)
plot(Time_TM(1:skipTM:end),Altitude_TM(1:skipTM:end)*3.28,'k*','linewidth',2)
p3 = plot(Time_TM(1:skipTM:end),Altitude_TM(1:skipTM:end)*3.28,'k--*','linewidth',2);

p4 = plot(Time_FR,Altitude_FR*3.28,'k-','linewidth',2);

xlim([0 100])
%grid mi
legend([p1,p2,p3,p4],'FASTsensor','OpenRocket','TeleMega','FASTrocket')

%%%Plot Vertical Velocity
plottool(1,'Vertical Velocity',18,'Time (sec)','Vertical Velocity (ft/s)');

plot(Time_FS,Vertical_Velocity_FS*3.28,'k-','linewidth',2);
plot(Time_FS(1:skipFS:end),Vertical_Velocity_FS(1:skipFS:end)*3.28,'ks','linewidth',2);
p1 = plot(101+Time_FS(1:skipFS:end),Vertical_Velocity_FS(1:skipFS:end)*3.28,'k-s','linewidth',2);

p2 = plot(Time_OR,Vertical_Velocity_OR,'k--','linewidth',2);

plot(Time_TM,Vertical_Velocity_TM*3.28,'k--','linewidth',2);
plot(Time_TM(1:skipTM:end),Vertical_Velocity_TM(1:skipTM:end)*3.28,'k*','linewidth',2);
p3 = plot(101+Time_TM(1:skipTM:end),Vertical_Velocity_TM(1:skipTM:end)*3.28,'k--*','linewidth',2);

p4 = plot(Time_FR,Vertical_Velocity_FR*3.28,'k-','linewidth',2);

xlim([0 100])
legend([p1,p2,p3,p4],'FASTsensor','OpenRocket','TeleMega','FASTrocket')



