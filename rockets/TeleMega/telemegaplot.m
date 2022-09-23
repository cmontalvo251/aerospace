clear
close all
clc

%filename = 'C:\Users\wdbro\Desktop\FASTSensor\flight_data.csv'; %%% Set correct file path
filename = '../Flight_Data/2014-05-17-serial-1254-flight-0002.csv';
startRow = 26; %%% Change this to exclude rows that have numbers that can't be data i.e. 2.15E+09.
endRow = 2034; %%% Change based on how many rows are in the file

formatSpec = '%2s%6s%3s%8s%10s%10s%12s%7s%9s%9s%11s%9s%9s%9s%9s%6s%6s%6s%6s%8s%8s%8s%8s%8s%8s%8s%8s%8s%3s%3s%4s%13s%13s%9s%6s%4s%4s%4s%4s%4s%10s%10s%5s%5s%7s%7s%7s%4s%4s%4s%4s%4s%4s%4s%4s%4s%4s%4s%4s%4s%4s%4s%4s%4s%4s%4s%4s%4s%4s%4s%4s%4s%4s%4s%4s%4s%4s%4s%s%[^\n\r]';

fileID = fopen(filename,'r');

dataArray = textscan(fileID, formatSpec, endRow-startRow+1, 'Delimiter', '', 'WhiteSpace', '', 'HeaderLines', startRow-1, 'ReturnOnError', false);

fclose(fileID);


raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = dataArray{col};
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));

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

R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw);
raw(R) = {NaN};

VarName1 = cell2mat(raw(:, 1));
Config = cell2mat(raw(:, 2));
ve = cell2mat(raw(:, 3));
rsion1 = cell2mat(raw(:, 4));
time = cell2mat(raw(:, 5));
VarName6 = cell2mat(raw(:, 6));
VarName7 = cell2mat(raw(:, 7));
VarName8 = cell2mat(raw(:, 8));
VarName9 = cell2mat(raw(:, 9));
acceleration = cell2mat(raw(:, 10));
pressure = cell2mat(raw(:, 11));
altitude = cell2mat(raw(:, 12));
height = cell2mat(raw(:, 13));
accel_speed = cell2mat(raw(:, 14));
baro_speed = cell2mat(raw(:, 15));
temperature = cell2mat(raw(:, 16));
drogue_voltage = cell2mat(raw(:, 17));
main_voltage = cell2mat(raw(:, 18));
battery_voltage = cell2mat(raw(:, 19));
accel_x = cell2mat(raw(:, 20));
accel_y = cell2mat(raw(:, 21));
accel_z = cell2mat(raw(:, 22));
gyro_x = cell2mat(raw(:, 23));
gyro_y = cell2mat(raw(:, 24));
gyro_z = cell2mat(raw(:, 25));
mag_x = cell2mat(raw(:, 26));
mag_y = cell2mat(raw(:, 27));
mag_z = cell2mat(raw(:, 28));
VarName29 = cell2mat(raw(:, 29));
VarName30 = cell2mat(raw(:, 30));
VarName31 = cell2mat(raw(:, 31));
latitude = cell2mat(raw(:, 32));
longitude = cell2mat(raw(:, 33));
gps_altitude = cell2mat(raw(:, 34));
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
pdop = cell2mat(raw(:, 45));
hdop = cell2mat(raw(:, 46));
vdop = cell2mat(raw(:, 47));
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

clearvars filename startRow endRow formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp me R;


fig1 = figure()
plot(time,acceleration)
xlabel('Time (s)')
ylabel('Acceleration (m/s^2)')
title('Acceleration vs Time')
grid on

fig2 = figure()
plot(time,pressure)
xlabel('Time (s)')
ylabel('Pressure (Pa)')
title('Pressure vs Time')
grid on

fig3 = figure()
plot(time,altitude)
xlabel('Time (s)')
ylabel('Altitude (m)')
title('Altitude vs Time')
grid on

fig4 = figure()
plot(time,height)
xlabel('Time (s)')
ylabel('Height (m)')
title('Height vs Time')
grid on

fig5 = figure()
plot(time,accel_speed)
xlabel('Time (s)')
ylabel('Speed (m/s)')
title('Speed vs Time')
grid on

fig6 = figure()
plot(time,baro_speed)
xlabel('Time (s)')
ylabel('Speed (m/s)')
title('Speed vs Time')
grid on

fig7 = figure()
plot(time,temperature)
xlabel('Time (s)')
ylabel('Temperature (�C)')
title('Temperature vs Time')
grid on

fig8 = figure()
plot(time,drogue_voltage)
xlabel('Time (s)')
ylabel('Drogue Voltage (V)')
title('Drogue Voltage vs Time')
grid on

fig9 = figure()
plot(time,main_voltage)
xlabel('Time (s)')
ylabel('Main Voltage (V)')
title('Main Voltage vs Time')
grid on

fig10 = figure()
plot(time,battery_voltage)
xlabel('Time (s)')
ylabel('Battery Voltage (V)')
title('Battery Voltage vs Time')
grid on

fig11 = figure()
plot(time,accel_x)
xlabel('Time (s)')
ylabel('Acceleration (m/s^2)')
title('X-Axis Acceleration')
grid on


fig12 = figure()
plot(time,accel_y)
xlabel('Time (s)')
ylabel('Acceleration (m/s^2)')
title('Y-Axis Acceleration')
grid on

fig13 = figure()
plot(time,accel_z)
xlabel('Time (s)')
ylabel('Acceleration (m/s^2)')
title('Z-Axis Acceleration')
grid on

fig14 = figure()
plot(time,gyro_x)
xlabel('Time (s)')
ylabel('Rotation Rate (�/s)')
title('X-Axis Gyro')
grid on

fig15 = figure()
plot(time,gyro_y)
xlabel('Time (s)')
ylabel('Rotation Rate (�/s)')
title('Y-Axis Gyro')
grid on

fig16 = figure()
plot(time,gyro_z)
xlabel('Time (s)')
ylabel('Rotation Rate (�/s)')
title('Z-Axis Gyro')
grid on

fig17 = figure()
plot(time,mag_x)
xlabel('Time (s)')
ylabel('Acceleration (m/s^2)')
title('X-Axis Magnetometer')
grid on

fig18 = figure()
plot(time,mag_y)
xlabel('Time (s)')
ylabel('Acceleration (m/s^2)')
title('Y-Axis Magnetometer')
grid on

fig19 = figure()
plot(time,mag_z)
xlabel('Time (s)')
ylabel('Acceleration (m/s^2)')
title('Z-Axis Magnetometer')
grid on

fig20 = figure()
plot(time,gps_altitude)
xlabel('Time (s)')
ylabel('Altitude (m)')
title('GPS Altitude')
grid on

fig21 = figure()
plot(time,pdop)
xlabel('Time (s)')
ylabel('Dilution of Precision')
title('Position DOP')
grid on

fig22 = figure()
plot(time,hdop)
xlabel('Time (s)')
ylabel('Dilution of Precision')
title('Horizontal DOP')
grid on

fig23 = figure()
plot(time,vdop)
xlabel('Time (s)')
ylabel('Dilution of Precision')
title('Vertical DOP')
grid on

fig24 = figure()
plot3(latitude,longitude,altitude)
xlabel('Latitude')
ylabel('Longitude')
zlabel('Altitude (m)')
title('Latitude and Longitude vs Altitude')
