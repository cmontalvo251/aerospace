%%%Compute Center of Pressure
clear
clc
close all

in = importdata('naca0012_cpx.txt');

data = in.data;

x = data(:,1);
cp = data(:,2);

%%%Split into lower and upper surface
found = 0;
idx = 1;
while ~found
    if x(idx) < x(idx+1)
       found = 1;
    end
   idx = idx + 1;
end

xupper = x(1:idx-2);
xlower = x(idx-1:end);
cpupper = cp(1:idx-2);
cplower = cp(idx-1:end);

%%%Flip Xupper around
xupper = xupper(end:-1:1);
cpupper = cpupper(end:-1:1);

plot(xupper,-cpupper,'r-')
hold on
plot(xlower,-cplower,'g-')

%%%Integrate Pressure Curve
N = 1000;
xcoord = linspace(0,1,N);
dx = xcoord(2)-xcoord(1);
cm = 0;
cl = 0;
for idx = 1:N
    if xcoord(idx) <= 0
        cpu = 0;
        cpl = 0;
    else
        cpu = interp1(xupper,cpupper,xcoord(idx));
        cpl = interp1(xlower,cplower,xcoord(idx));
    end
    delcp = cpu-cpl;
    cm = cm + delcp*xcoord(idx)*dx;
    cl = cl + delcp*dx;
end
center = cm/cl
%%%NACA 0012 = 0.2575
%%%NACA 2212 = 0.3005
%%%NACA 2412 = 0.3201

