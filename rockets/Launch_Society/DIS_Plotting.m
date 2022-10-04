clear
clc
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%FASTSensor Data
%filename = 'FAST_01.TXT';
%delimiter = ' ';
disp('Opening FastSensor Data')
filename = 'Flight_Data/Faust_Wiesneth_Brown_Data/FASTSensor_Data_NO_FIRST_LINE.TXT';
delimiter = ',';
formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'EmptyValue' ,NaN, 'ReturnOnError', false);
fclose(fileID);
disp('Processing Data')
Hour = dataArray{:, 1};
Minute = dataArray{:, 2};
Second = dataArray{:, 3};
LastPrint = dataArray{:, 4};
Fix = dataArray{:, 5};
Latitude_FS = dataArray{:, 6};
Longitude_FS = dataArray{:, 7};
Speed = dataArray{:, 8};
Angle = dataArray{:, 9};
Altitude_GPS = dataArray{:, 10};
X_Accel_FS = -dataArray{:, 11};
Y_Accel_FS = -dataArray{:, 12};
Z_Accel_FS = -dataArray{:, 13};
Z_Mag_FS = dataArray{:, 14};
Z_Mag_Shifted = dataArray{:, 14}-15;
X_Mag_FS = dataArray{:, 15};
X_Mag_Shifted = dataArray{:, 15}+38;
Y_Mag_FS = -dataArray{:, 16};
Y_Mag_Shifted = -dataArray{:, 16}-95;
X_Gyro_FS = dataArray{:, 17};
Y_Gyro_FS = -dataArray{:, 18};
Z_Gyro_FS = dataArray{:, 19};
Temperature_FS = dataArray{:, 20};
Pressure_FS = dataArray{:, 21};
Ground_Pressure = mean(Pressure_FS(1:100,:));

X_Accel_FS = X_Accel_FS(3.76e4:3.857e4,:);
Y_Accel_FS = Y_Accel_FS(3.76e4:3.857e4,:);
Z_Accel_FS = Z_Accel_FS(3.76e4:3.857e4,:);
Pressure_FS = Pressure_FS(3.76e4:3.857e4,:);
X_Mag_FS = X_Mag_FS(3.76e4:3.857e4,:);
X_Mag_Shifted = X_Mag_Shifted(3.76e4:3.857e4,:);
Y_Mag_FS = Y_Mag_FS(3.76e4:3.857e4,:);
Y_Mag_Shifted = Y_Mag_Shifted(3.76e4:3.857e4,:);
Z_Mag_FS = Z_Mag_FS(3.76e4:3.857e4,:);
Z_Mag_Shifted = Z_Mag_Shifted(3.76e4:3.857e4,:);
X_Gyro_FS = X_Gyro_FS(3.76e4:3.857e4,:);
Roll_Rate_FS = X_Gyro_FS;
Y_Gyro_FS = Y_Gyro_FS(3.76e4:3.857e4,:);
Pitch_Rate_FS = Y_Gyro_FS;
Z_Gyro_FS = Z_Gyro_FS(3.76e4:3.857e4,:);
Yaw_Rate_FS = Z_Gyro_FS;
Temperature_FS = Temperature_FS(3.76e4:3.857e4,:);
Latitude_FS = Latitude_FS(3.76e4:3.857e4,:);
Longitude_FS = Longitude_FS(3.76e4:3.857e4,:);
Speed = Speed(3.76e4:3.857e4,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%OpenRocket Data
%filename = 'Open_Rocket.csv';
filename = 'Flight_Data/Faust_Wiesneth_Brown_Data/Open_Rocket_Data.csv';
delimiter = ',';
formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%s%s%f%f%[^\n\r]';
fileID = fopen(filename,'r');
disp('Opening OpenRocket Data')
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN, 'ReturnOnError', false);
fclose(fileID);
disp('Processing data')
%%% Allocate imported array to column variable names
Time_OR = dataArray{:, 1};
Altitude_OR = dataArray{:, 2};
Vertical_Velocity_OR = dataArray{:, 3};
Vertical_Acceleration_OR = dataArray{:, 4};
Total_Velocity_OR = dataArray{:, 5};
Total_Acceleration_OR = dataArray{:, 6};
Position_East_of_Launch_OR = dataArray{:, 7};
Position_North_of_Launch_OR = dataArray{:, 8};
Lateral_Distance_OR = dataArray{:, 9};
Lateral_Velocity_OR = dataArray{:, 10};
Lateral_Acceleration_OR = dataArray{:, 11};
Latitude_OR = dataArray{:, 12};
Longitude_OR = dataArray{:, 13};
Angle_of_Attack_OR = dataArray{:, 14};
Roll_Rate_OR = dataArray{:, 15};
Pitch_Rate_OR = dataArray{:, 16};
Yaw_Rate_OR = dataArray{:, 17};
Mass_OR = dataArray{:, 18};
Longitudinal_Moment_of_Inertia_OR = dataArray{:, 19};
Rotational_Moment_of_Inertia_OR = dataArray{:, 20};
Thrust_OR = dataArray{:, 21};
Pitch_Moment_Coefficient_OR = dataArray{:, 22};
Yaw_Moment_Coefficient_OR = dataArray{:, 23};
Theta_OR = dataArray{:, 24};
Psi_OR = dataArray{:, 25};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%TeleMega Data
%filename = 'Jan_Flight.csv'; %%% Set correct file path
filename = 'Flight_Data/Faust_Wiesneth_Brown_Data/TeleMega_Data.csv';
startRow = 51; %%% Change this to exclude rows that have numbers that can't be data i.e. 2.15E+09.
endRow = 2693; %%% Change based on how many rows are in the file
formatSpec = '%2s%6s%3s%8s%10s%10s%12s%7s%9s%9s%11s%9s%9s%9s%9s%6s%6s%6s%6s%8s%8s%8s%8s%8s%8s%8s%8s%8s%3s%3s%4s%13s%13s%9s%6s%4s%4s%4s%4s%4s%10s%10s%5s%5s%7s%7s%7s%4s%4s%4s%4s%4s%4s%4s%4s%4s%4s%4s%4s%4s%4s%4s%4s%4s%4s%4s%4s%4s%4s%4s%4s%4s%4s%4s%4s%4s%4s%4s%s%[^\n\r]';
fileID = fopen(filename,'r');
disp('Opening TeleMega Data')
dataArray = textscan(fileID, formatSpec, endRow-startRow+1, 'Delimiter', '', 'WhiteSpace', '', 'HeaderLines', startRow-1, 'ReturnOnError', false);
fclose(fileID);


raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
disp('Initial Post-Process')
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = dataArray{col};
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));
disp('Secondary Post-Process....this may take a while')
for col=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79]
    rawData = dataArray{col};
    for row=1:size(rawData, 1);
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData{row}, regexstr, 'names');
            numbers = result.numbers;
            invalidThousandsSeparator = false;
            if any(numbers==',');
                thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(thousandsRegExp, ',', 'once'));
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            if ~invalidThousandsSeparator;
                numbers = textscan(strrep(numbers, ',', ''), '%f');
                numericData(row, col) = numbers{1};
                raw{row, col} = numbers{1};
            end
        catch me
        end
    end
end
disp('Final Post-Process')

R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw);
raw(R) = {NaN};

VarName1 = cell2mat(raw(:, 1));
Config = cell2mat(raw(:, 2));
ve = cell2mat(raw(:, 3));
rsion1 = cell2mat(raw(:, 4));
Time_TM = cell2mat(raw(:, 5));
VarName6 = cell2mat(raw(:, 6));
VarName7 = cell2mat(raw(:, 7));
VarName8 = cell2mat(raw(:, 8));
VarName9 = cell2mat(raw(:, 9));
acceleration_TM = cell2mat(raw(:, 10));
pressure_TM = cell2mat(raw(:, 11));
sea_altitude_TM = cell2mat(raw(:, 12));
Altitude_TM = cell2mat(raw(:, 13));
accel_speed_TM = cell2mat(raw(:, 14));
baro_speed_TM = cell2mat(raw(:, 15));
Vertical_Velocity_TM = baro_speed_TM;
Temperature_TM = cell2mat(raw(:, 16));
drogue_voltage_TM = cell2mat(raw(:, 17));
main_voltage_TM = cell2mat(raw(:, 18));
battery_voltage_TM = cell2mat(raw(:, 19));
X_Accel_TM = cell2mat(raw(:, 20));
Y_Accel_TM = cell2mat(raw(:, 21));
Z_Accel_TM = cell2mat(raw(:, 22));
X_Gyro_TM = cell2mat(raw(:, 23));
Roll_Rate_TM = X_Gyro_TM;
Y_Gyro_TM = cell2mat(raw(:, 24));
Pitch_Rate_TM = Y_Gyro_TM;
Z_Gyro_TM = cell2mat(raw(:, 25));
Yaw_Rate_TM = Z_Gyro_TM;
X_Mag_TM = cell2mat(raw(:, 26));
Y_Mag_TM = cell2mat(raw(:, 27));
Z_Mag_TM = cell2mat(raw(:, 28));
VarName29 = cell2mat(raw(:, 29));
VarName30 = cell2mat(raw(:, 30));
VarName31 = cell2mat(raw(:, 31));
Latitude_TM = cell2mat(raw(:, 32));
Longitude_TM = cell2mat(raw(:, 33));
gps_altitude_TM = cell2mat(raw(:, 34));
VarName35 = cell2mat(raw(:, 35));
VarName36 = cell2mat(raw(:, 36));
VarName37 = cell2mat(raw(:, 37));
VarName38 = cell2mat(raw(:, 38));
VarName39 = cell2mat(raw(:, 39));
VarName40 = cell2mat(raw(:, 40));
VarName41 = cell2mat(raw(:, 41));
VarName42 = cell2mat(raw(:, 42));
VarName43 = cell2mat(raw(:, 43));
VarName44 = cell2mat(raw(:, 44));
pdop_TM = cell2mat(raw(:, 45));
hdop_TM = cell2mat(raw(:, 46));
vdop_TM = cell2mat(raw(:, 47));
VarName48 = cell2mat(raw(:, 48));
VarName49 = cell2mat(raw(:, 49));
VarName50 = cell2mat(raw(:, 50));
VarName51 = cell2mat(raw(:, 51));
VarName52 = cell2mat(raw(:, 52));
VarName53 = cell2mat(raw(:, 53));
VarName54 = cell2mat(raw(:, 54));
VarName55 = cell2mat(raw(:, 55));
VarName56 = cell2mat(raw(:, 56));
VarName57 = cell2mat(raw(:, 57));
VarName58 = cell2mat(raw(:, 58));
VarName59 = cell2mat(raw(:, 59));
VarName60 = cell2mat(raw(:, 60));
VarName61 = cell2mat(raw(:, 61));
VarName62 = cell2mat(raw(:, 62));
VarName63 = cell2mat(raw(:, 63));
VarName64 = cell2mat(raw(:, 64));
VarName65 = cell2mat(raw(:, 65));
VarName66 = cell2mat(raw(:, 66));
VarName67 = cell2mat(raw(:, 67));
VarName68 = cell2mat(raw(:, 68));
VarName69 = cell2mat(raw(:, 69));
VarName70 = cell2mat(raw(:, 70));
VarName71 = cell2mat(raw(:, 71));
VarName72 = cell2mat(raw(:, 72));
VarName73 = cell2mat(raw(:, 73));
VarName74 = cell2mat(raw(:, 74));
VarName75 = cell2mat(raw(:, 75));
VarName76 = cell2mat(raw(:, 76));
VarName77 = cell2mat(raw(:, 77));
VarName78 = cell2mat(raw(:, 78));
VarName79 = cell2mat(raw(:, 79));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%FASTRocket Data
%filename = 'FASTRocket_Data.csv';
filename = 'Flight_Data/Faust_Wiesneth_Brown_Data/FASTRocket_Data.csv';
delimiter = ',';
startRow = 2;
formatSpec = '%s%s%s%s%s%s%s%s%s%[^\n\r]';
fileID = fopen(filename,'r');
disp('Opening FASTRocket Data')
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false);
fclose(fileID);
disp('Initial Post-Process')
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = dataArray{col};
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));
disp('Secondary Post-Process.....this may take a while')
for col=[1,2,3,4,5,6,7,8,9]
    rawData = dataArray{col};
    for row=1:size(rawData, 1);
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData{row}, regexstr, 'names');
            numbers = result.numbers;
            invalidThousandsSeparator = false;
            if any(numbers==',');
                thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(numbers, thousandsRegExp, 'once'));
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            if ~invalidThousandsSeparator;
                numbers = textscan(strrep(numbers, ',', ''), '%f');
                numericData(row, col) = numbers{1};
                raw{row, col} = numbers{1};
            end
        catch me
        end
    end
end
disp('Final post-process')
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
raw(R) = {NaN}; % Replace non-numeric cells

Time_FR = cell2mat(raw(:, 1));
Altitude_FR = -cell2mat(raw(:, 2))/3.23;
Roll_Rate_FR = cell2mat(raw(:, 3));
Pitch_Rate_FR = cell2mat(raw(:, 4));
Yaw_Rate_FR = cell2mat(raw(:, 5));
Phi_FR = cell2mat(raw(:, 6));
Theta_FR = cell2mat(raw(:, 7));
Psi_FR = cell2mat(raw(:, 8));
Vertical_Velocity_FR = -cell2mat(raw(:, 9))/3.23;

clearvars filename delimiter startRow formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp me R;

disp('Plotting')

save All_the_datas

%%%%Since we did this experiment indoors we won't have any GPS data. Thus I
%%%%will be using LastPrint as the time vector although remind me to give
%%%%you a MATLAB code that will fuse ArduinoTime with GPS time. 
%%%%time_vec = LastPrint-LastPrint(1); 
Time_FS = linspace(0,98,length(Pressure_FS))';
    
%%%PLOTTING

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
figure1 = figure();
set(figure1,'color','white')
plot(Time_FS,X_Mag_FS,'r-*','linewidth',2)
hold on
plot(Time_TM,X_Mag_TM*100,'b--','linewidth',2)
grid on
title('X Axis Magnetometer','fontsize',20,'fontweight','bold')
xlabel('Time(s)','fontsize',16,'fontweight','bold')
ylabel('Magnetic Field(uT)','fontsize',16,'fontweight','bold')
legend('FASTSensor','TeleMega')

figure2 = figure();
set(figure2,'color','white')
plot(Time_FS,Y_Mag_FS,'r-*','linewidth',2)
hold on
plot(Time_TM,Y_Mag_TM*100,'b--','linewidth',2)
grid on
title('Y Axis Magnetometer','fontsize',20,'fontweight','bold')
xlabel('Time(s)','fontsize',16,'fontweight','bold')
ylabel('Magnetic Field(uT)','fontsize',16,'fontweight','bold')
legend('FASTSensor','TeleMega')

figure3 = figure();
set(figure3,'color','white')
plot(Time_FS,Z_Mag_FS,'r-*','linewidth',2)
hold on 
plot(Time_TM,Z_Mag_TM*100,'b--','linewidth',2)
grid on
title('Z Axis Magnetometer','fontsize',20,'fontweight','bold')
xlabel('Time(s)','fontsize',16,'fontweight','bold')
ylabel('Magnetic Field(uT)','fontsize',16,'fontweight','bold')
legend('FASTSensor','TeleMega')

% plot(time_vec,MAG)
% xlabel('Time (sec)')
% ylabel('Mag Field (nF)')
% legend('X','Y','Z')

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
figure4 = figure();
set(figure4,'color','white')
plot(Time_FS,X_Accel_FS,'r-*','linewidth',2)
hold on
plot(Time_TM,X_Accel_TM,'b--','linewidth',2)
grid on
title('X Axis Accelerometer','fontsize',20,'fontweight','bold')
xlabel('Time(s)','fontsize',16,'fontweight','bold')
ylabel('Acceleration(m/s)','fontsize',16,'fontweight','bold')
legend('FASTSensor','TeleMega')

figure5 = figure();
set(figure5,'color','white')
plot(Time_FS,Y_Accel_FS,'r-*','linewidth',2)
hold on
plot(Time_TM,Y_Accel_TM,'b--','linewidth',2)
grid on
title('Y Axis Accelerometer','fontsize',20,'fontweight','bold')
xlabel('Time(s)','fontsize',16,'fontweight','bold')
ylabel('Acceleration(m/s)','fontsize',16,'fontweight','bold')
legend('FASTSensor','TeleMega')

figure6 = figure();
set(figure6,'color','white')
plot(Time_FS,Z_Accel_FS,'r-*','linewidth',2)
hold on
plot(Time_TM,Z_Accel_TM,'b--','linewidth',2)
grid on
title('Z Axis Accelerometer','fontsize',20,'fontweight','bold')
xlabel('Time(s)','fontsize',16,'fontweight','bold')
ylabel('Acceleration(m/s)','fontsize',16,'fontweight','bold')
legend('FASTSensor','TeleMega')

%%%filter Rate data
sig_gyro = 0.05;
Roll_Rate_FS = filter_signal(Roll_Rate_FS,sig_gyro);
Pitch_Rate_FS = filter_signal(Pitch_Rate_FS,sig_gyro);
Yaw_Rate_FS = filter_signal(Yaw_Rate_FS,sig_gyro);
RATE_FS = [Roll_Rate_FS,Pitch_Rate_FS,Yaw_Rate_FS];
Roll_Rate_TM = filter_signal(Roll_Rate_TM,sig_gyro);
Pitch_Rate_TM = filter_signal(Pitch_Rate_TM,sig_gyro);
Yaw_Rate_TM = filter_signal(Yaw_Rate_TM,sig_gyro);
RATE_TM = [Roll_Rate_TM,Pitch_Rate_TM,Yaw_Rate_TM];

%%%Plot Angular Rate Data
figure7 = figure();
set(figure7,'color','white')
plot(Time_FS,Roll_Rate_FS*(180/pi),'r-*','linewidth',2)
hold on
plot(Time_OR,Roll_Rate_OR*(180/pi),'k-.','linewidth',4)
plot(Time_TM,Roll_Rate_TM,'b--','linewidth',2)
plot(Time_FR,Roll_Rate_FR,'g-','linewidth',2)
grid on
title('Roll Rate','fontsize',20,'fontweight','bold')
xlabel('Time(s)','fontsize',16,'fontweight','bold')
ylabel('Angular Velocity(degrees/sec)','fontsize',16,'fontweight','bold')
legend('FASTSensor','OpenRocket','TeleMega','FASTRocket')

figure8 = figure();
set(figure8,'color','white')
plot(Time_FS,Pitch_Rate_FS*(180/pi),'r-*','linewidth',2)
hold on
plot(Time_OR,Pitch_Rate_OR*(180/pi),'k-.','linewidth',4)
plot(Time_TM,Pitch_Rate_TM,'b--','linewidth',2)
plot(Time_FR,Pitch_Rate_FR,'g-','linewidth',2)
grid on
title('Pitch Rate','fontsize',20,'fontweight','bold')
xlabel('Time(s)','fontsize',16,'fontweight','bold')
ylabel('Angular Velocity(degrees/sec)','fontsize',16,'fontweight','bold')
legend('FASTSensor','OpenRocket','TeleMega','FASTRocket')

figure9 = figure();
set(figure9,'color','white')
plot(Time_FS,Yaw_Rate_FS*(180/pi),'r-*','linewidth',2)
hold on
plot(Time_OR,Yaw_Rate_OR*(180/pi),'k-.','linewidth',4)
plot(Time_TM,Yaw_Rate_TM,'b--','linewidth',2)
plot(Time_FR,Yaw_Rate_FR,'g-','linewidth',2)
grid on
title('Yaw Rate','fontsize',20,'fontweight','bold')
xlabel('Time(s)','fontsize',16,'fontweight','bold')
ylabel('Angular Velocity(degrees/sec)','fontsize',16,'fontweight','bold')
legend('FASTSensor','OpenRocket','TeleMega','FASTRocket')

% figure()
% plot(time_vec,RATE)
% xlabel('Time (sec)')
% ylabel('Rate Gyro (rad/s)')
% legend('X','Y','Z')

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

%%%%After running this code a few times I realized that when the pitch
%%%%angle approaches 90 degrees we run into a singularity problem. As such
%%%%I multiplied phi by cos(theta_accel)^2
Phi_FS = Phi_FS.*cos(Theta_FS).^2;
Phi_TM = Phi_TM.*cos(Theta_TM).^2;

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

Psi_FS = atan2(Z_Mag_FS .* sin(Phi_FS) - Y_Mag_FS .* cos(Phi_FS), X_Mag_FS .* cos(Theta_FS) + Y_Mag_FS .* sin(Theta_FS) .* sin(Phi_FS) + Z_Mag_FS .* sin(Theta_FS) .* cos(Phi_FS));
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

% Phi_TM = zeros(length(X_Gyro_TM)-1,1);
% Theta_TM = Phi_TM;
% Psi_TM = Theta_TM;
% 
% Phi_TM(1) = Phi_FS(1);
% Theta_TM(1) = Theta_FS(1);
% Psi_TM(1) = Psi_FS(1);

% for idx = 1:length(Roll_Rate_TM)-1
%     Phi_TM(idx+1) = Phi_TM(idx) + 0.5*(Roll_Rate_TM(idx) + X_Gyro_TM(idx+1))*(Time_TM(idx+1)-Time_TM(idx));
%     Theta_TM(idx+1) = Theta_TM(idx) + 0.5*(Pitch_Rate_TM(idx) + Pitch_Rate_TM(idx+1))*(Time_TM(idx+1)-Time_TM(idx));
%     Psi_TM(idx+1) = Psi_TM(idx) + 0.5*(Yaw_Rate_TM(idx) + Yaw_Rate_TM(idx+1))*(Time_TM(idx+1)-Time_TM(idx));
% end

%%%Plot Phi
%%%FASTSensor Phi is from accelerometer
figure10 = figure();
set(figure10,'color','white')
plot(Time_FS,Phi_FS*180/pi,'r-*','LineWidth',2);
hold on
plot(Time_TM,Phi_TM*180/pi,'b--','LineWidth',2);
plot(Time_FR,Phi_FR,'g-','LineWidth',2);
grid on
title('Phi','fontsize',20,'fontweight','bold')
xlabel('Time(s)','fontsize',16,'fontweight','bold')
ylabel('Phi(degrees)','fontsize',16,'fontweight','bold')
legend('FASTSensor','TeleMega','FASTRocket')

%%%Plot Theta
%%%FASTSensor Theta is from accelerometer
figure11 = figure();
set(figure11,'color','white')
plot(Time_FS,Theta_FS*180/pi,'r-*','LineWidth',2);
hold on
plot(Time_OR,Theta_OR,'k-.','LineWidth',4);
plot(Time_TM,Theta_TM*180/pi,'b--','LineWidth',2);
plot(Time_FR,Theta_FR,'g-','LineWidth',2);
grid on
title('Theta','fontsize',20,'fontweight','bold')
xlabel('Time(s)','fontsize',16,'fontweight','bold')
ylabel('Theta(degrees)','fontsize',16,'fontweight','bold')
legend('FASTSensor','OpenRocket','TeleMega','FASTRocket')

%%%Plot Psi
%%%FASTSensor Psi is from magnetometer
figure12 = figure();
set(figure12,'color','white')
plot(Time_FS,Psi_FS*180/pi,'r-*','LineWidth',2);
hold on
plot(Time_OR,Psi_OR,'k-.','LineWidth',4);
plot(Time_TM,Psi_TM*180/pi,'b--','LineWidth',2);
plot(Time_FR,Psi_FR,'g-','LineWidth',2);
grid on
title('Psi','fontsize',20,'fontweight','bold')
xlabel('Time(s)','fontsize',16,'fontweight','bold')
ylabel('Psi(degrees)','fontsize',16,'fontweight','bold')
legend('FASTSensor','OpenRocket','TeleMega','FASTRocket')

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
%%%Plot Altitude
figure13 = figure();
set(figure13,'color','white')
plot(Time_FS,Altitude_FS,'r-*','linewidth',2)
hold on
plot(Time_OR,Altitude_OR/3.23,'k-.','linewidth',4)
plot(Time_TM,Altitude_TM,'b--','linewidth',2)
plot(Time_FR,Altitude_FR,'g-','linewidth',2)
grid on
title('Altitude','fontsize',20,'fontweight','bold')
xlabel('Time(s)','fontsize',16,'fontweight','bold')
ylabel('Altitude(m)','fontsize',16,'fontweight','bold')
legend('FASTSensor','OpenRocket','TeleMega','FASTRocket')

%%%Getting Xlat and Ylat data
bad_data_FS = find(abs(Latitude_FS) < 28);
Latitude_FS(bad_data_FS)= [];
Longitude_FS(bad_data_FS) = [];
Altitude_FS(bad_data_FS) = [];
bad_data_TM = find(abs(Latitude_TM) < 28);
Latitude_TM(bad_data_TM)= [];
Longitude_TM(bad_data_TM) = [];
Altitude_TM(bad_data_TM) = [];

[xlat_FS,ylat_FS] = convertLATLON(Latitude_FS,Longitude_FS,[Latitude_FS(1),Longitude_FS(1)]); %% Convert to cartesian
% xlat_FS = xlat_FS.*3.28;
% ylat_FS = ylat_FS.*3.28;
[xlat_OR,ylat_OR] = convertLATLON(Latitude_OR,Longitude_OR,[Latitude_OR(1),Longitude_OR(1)]); %% Convert to cartesian
% xlat_FS = xlat_FS.*3.28;
% ylat_FS = ylat_FS.*3.28;
[xlat_TM,ylat_TM] = convertLATLON(Latitude_TM,Longitude_TM,[Latitude_TM(1),Longitude_TM(1)]); %% Convert to cartesian
% xlat_FS = xlat_FS.*3.28;
% ylat_FS = ylat_FS.*3.28;


%%%Plot Xlat and Ylat
figure14 = figure();
set(figure14,'color','white')
plot(xlat_FS,ylat_FS,'r-*','linewidth',2)
hold on
plot(xlat_OR,ylat_OR,'k-.','linewidth',3)
plot(xlat_TM,ylat_TM,'b--','linewidth',2)
grid on
title('Ylat vs Xlat','fontsize',20,'fontweight','bold')
xlabel('Xlat(m)','fontsize',16,'fontweight','bold')
ylabel('Ylat(m)','fontsize',16,'fontweight','bold')
legend('FASTSensor','OpenRocket','TeleMega')


%%%Plot Vertical Velocity
figure15 = figure();
set(figure15,'color','white')
plot(Time_FS,Vertical_Velocity_FS,'r-*','linewidth',2)
hold on
plot(Time_OR,Vertical_Velocity_OR/3.23,'k-.','linewidth',4)
plot(Time_TM,Vertical_Velocity_TM,'b--','linewidth',2)
plot(Time_FR,Vertical_Velocity_FR,'g-','linewidth',2)
grid on
title('Vertical Velocity','fontsize',20,'fontweight','bold')
xlabel('Time(s)','fontsize',16,'fontweight','bold')
ylabel('Velocity(m/s)','fontsize',16,'fontweight','bold')
legend('FASTSensor','OpenRocket','TeleMega','FASTRocket')


%%%Plot Temperature
figure16 = figure();
set(figure16,'color','white')
plot(Time_FS,Temperature_FS,'r-*','linewidth',2)
hold on
plot(Time_TM,Temperature_TM,'b--','linewidth',2)
grid on
title('Temperature','fontsize',20,'fontweight','bold')
xlabel('Time(s)','fontsize',16,'fontweight','bold')
ylabel('Temperature(�C)','fontsize',16,'fontweight','bold')
legend('FASTSensor','TeleMega')


%%%Plot Shifted Mag 

figure17 = figure();
set(figure17,'color','white')
plot(Time_FS,X_Mag_Shifted,'r-*','linewidth',2)
hold on
plot(Time_TM,X_Mag_TM*100,'b--','linewidth',2)
grid on
title('X Axis Magnetometer','fontsize',20,'fontweight','bold')
xlabel('Time(s)','fontsize',16,'fontweight','bold')
ylabel('Magnetic Field(uT)','fontsize',16,'fontweight','bold')
legend('FASTSensor','TeleMega')

figure18= figure();
set(figure18,'color','white')
plot(Time_FS,Y_Mag_Shifted,'r-*','linewidth',2)
hold on
plot(Time_TM,Y_Mag_TM*100,'b--','linewidth',2)
grid on
title('Y Axis Magnetometer','fontsize',20,'fontweight','bold')
xlabel('Time(s)','fontsize',16,'fontweight','bold')
ylabel('Magnetic Field(uT)','fontsize',16,'fontweight','bold')
legend('FASTSensor','TeleMega')

figure19 = figure();
set(figure19,'color','white')
plot(Time_FS,Z_Mag_Shifted,'r-*','linewidth',2)
hold on 
plot(Time_TM,Z_Mag_TM*100,'b--','linewidth',2)
grid on
title('Z Axis Magnetometer','fontsize',20,'fontweight','bold')
xlabel('Time(s)','fontsize',16,'fontweight','bold')
ylabel('Magnetic Field(uT)','fontsize',16,'fontweight','bold')
legend('FASTSensor','TeleMega')


