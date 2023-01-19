function Expectation()
global m

close all

dt = 0.01;
tauc = 10;
tau = -tauc:dt:tauc;
t = 0:dt:5;
tend = t(end)

%%%Generate X
x = xvar(t);

%plottool(1,'X',18,'Time (sec)','X');
%plot(t,x,'LineWidth',2)

%%%%Expectation Operator
R_exp = 0.*tau;
R_interp = 0.*R_exp;
twindow = 50; %percent
sides = round((twindow/100)*length(t));
s = find(tau >= -t(1),1);
e = find(tau >= t(end)-t(end-sides),1)-1;
tlength = (t(end-sides)-t(1));

for ii = s:e
    for tdx = 1:(length(t)-sides)
        xinterp = interp1(t,x,t(tdx)+tau(ii));
        R_interp(ii) = R_interp(ii) + dt*xvar(t(tdx))*xinterp;
        xactual = xvar(t(tdx)+tau(ii));
        R_exp(ii) = R_exp(ii) + dt*x(tdx)*xactual;
    end
    disp(tau(ii))
end
R_exp = (1/tlength)*R_exp;
R_interp = (1/tlength)*R_interp;

%%%Generate R
R_analytic =( m^2/tlength)*(tlength^3/3 + tau.*tlength^2/2);

plottool(1,'R',18,'Tau (sec)','R');
plot(tau,R_analytic,'b-','LineWidth',2)
plot(tau,R_exp,'r-','LineWidth',2)
plot(tau,R_interp,'g-','LineWidth',2)
xlim([tau(s) tau(e)])
legend('Analytic','Numerical','Interpolation')

keyboard

function out = xvar(t)
global m

m = 2.3;

out = m*t;
