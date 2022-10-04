close all

data = dlmread('Velocity_Data.csv');

x = data(:,1);
y = data(:,2);

xfit = x(x>1.6);
yfit = y(x>1.6);

xfit = xfit(xfit < 9.7);
yfit = yfit(xfit < 9.7);

time = xfit;
velocity = yfit;

t0 = time(1);
v0 = velocity(1);

dt = 0.1;

teuler = t0:dt:time(end);
veuler = 0*teuler;
for idx = 1:length(teuler)-1
    %vprime = .....
    veuler(idx+1) = veuler(idx) + vprime*dt;
end
hold on
plot(teuler,veuler,'r-')





