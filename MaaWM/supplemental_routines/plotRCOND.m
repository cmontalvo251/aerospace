clear
clc
close all

A = importdata('Output_Files/Data_Point_Studies/RCOND.txt');

RCOND = A.data(:,end);
dx = 1000./(30:40);


semilogy(dx,RCOND,'LineWidth',2)
plottool(0,'RCOND',12,'Distance between Data Points (m)','RCOND');

dRCONDdx = RCOND(2:end)-RCOND(1:end-1);




