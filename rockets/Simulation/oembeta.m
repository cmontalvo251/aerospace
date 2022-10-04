function [teuler, veuler]=oembeta(time,V,Tk)
close all
clc,clear
data1 = dlmread('udata.csv');
data2 = xlsread('vel_with_time.csv');
x = data1(:,1);
y = data1(:,2);

xfit = x(x>0.01);
yfit = y(x>0.01);

xfit = xfit(xfit < 25.8);
yfit = yfit(xfit < 25.8);

time = x;
V = y;



% figure()
% plot(time,V)
% hold on
% xlabel('Time'), ylabel('Velocity')

Tk = 0.2;

[veuler, teuler] = eulerbeta(time,V,Tk);


vinterp = interp1(teuler,veuler,time);

time2=data2(:,1);
vel=data2(:,4);

N = length(V);
err = (1/N)*sum((V-vinterp).^2);
iters = 0;
while iters < 10
    Tk = Tk - errprime(Tk,V,time)/err2prime(Tk,time,V);
    [veuler, teuler] = eulerbeta(time,V,Tk);
    vinterp = interp1(teuler,veuler,time);
    err = (1/N)*sum((V-vinterp));
    iters = iters+1;
    hold on
    plot(teuler,veuler,'b-', 'LineWidth', 2)
    plot(time2,vel,'r--','LineWidth',2)
    title('Axial Velocity','Fontsize',18)
    xlabel('t (s)','FontSize',18), ylabel('velocity (ft/s)','FontSize',18)
    legend('6DOF','OpenRocket')
    grid on
   

end












