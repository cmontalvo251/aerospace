%%This function will test the Sensor Sims
clear
clc
close all

%%Quant Sim
time = -10:0.001:10;
BITS = 12;
actualsignal = sin(time);
quantsignal = QUANTsim(actualsignal,2,BITS,[-1 1]);
h1 = figure();
set(h1,'color','white')
plot(time,actualsignal,'b--')
hold on
plot(time,quantsignal,'b-','LineWidth',2)
grid on