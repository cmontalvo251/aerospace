clear
clc
close all

xeq = [0 0 0 4.278 0 0 ]';
ueq = [50.52 0]';


feq = Derivatives(xeq,0,ueq)

dx1 = 1e-6;
xeqp1 = xeq + [dx1 0 0 0 0 0]';
xeqm1 = xeq - [dx1 0 0 0 0 0]';
A1stcolumn = (Derivatives(xeqp1,0,ueq)-Derivatives(xeqm1,0,ueq))/(2*dx1)

dx2 = 1e-6;
xeqp2 = xeq + [0 dx2 0 0 0 0]';
xeqm2 = xeq - [0 dx2 0 0 0 0]';
A2ndcolumn = (Derivatives(xeqp2,0,ueq)-Derivatives(xeqm2,0,ueq))/(2*dx2)

dx3 = 1e-6;
xeqp3 = xeq + [0 0 dx3 0 0 0]';
xeqm3 = xeq - [0 0 dx3 0 0 0]';
A3rdcolumn = (Derivatives(xeqp3,0,ueq)-Derivatives(xeqm3,0,ueq))/(2*dx3)

dx4 = 1e-6;
xeqp4 = xeq + [0 0 0 dx4 0 0]';
xeqm4 = xeq - [0 0 0 dx4 0 0]';
A4thcolumn = (Derivatives(xeqp4,0,ueq)-Derivatives(xeqm4,0,ueq))/(2*dx4)

dx5 = 1e-6;
xeqp5 = xeq + [0 0 0 0 dx5 0]';
xeqm5 = xeq - [0 0 0 0 dx5 0]';
A5thcolumn = (Derivatives(xeqp5,0,ueq)-Derivatives(xeqm5,0,ueq))/(2*dx5)

dx6 = 1e-6;
xeqp6 = xeq + [0 0 0 0 0 dx6]';
xeqm6 = xeq - [0 0 0 0 0 dx6]';
A6thcolumn = (Derivatives(xeqp6,0,ueq)-Derivatives(xeqm6,0,ueq))/(2*dx6)


A = [A1stcolumn A2ndcolumn A3rdcolumn A4thcolumn A5thcolumn A6thcolumn]