function Example_Problem()
global S0 w1 w2

close all
clc

MAKE_R = 1;
MAKE_S = 0;
MAKE_X = 1;

%%%Parameters
w1 = 2;
w2 = 3;
S0 = 1;
wc = 4;
tauc = 20;
dtau = 0.01;
dw = 0.01;
tau = -tauc:dtau:tauc;
w = -wc:dw:wc;

%%%Now create a time series using S(w) and Shinozuka's formula
if MAKE_X

    plottool(1,'X',18,'Time (sec)','x');
    colors = {'b-','r-','g-'};
    for ll = 1:1
    dt = dtau;
    tend = 10;
    NW = 1000;
    O1L = -wc;
    O1U = wc;
    dO1 = (O1U - O1L)/NW;
    delO1 = 0.1*dO1;
    dO = dO1;
    t = 0:dt:tend;
    x = 0.*t;
    phin = zeros(NW,1);
    
    for ii = 1:length(t)
        for k1 = 1:NW
            if ii == 1
                phin(k1) = 2*pi*rand;
            end
            O1k1 = O1L + (k1-0.5)*dO1;
            Ok = O1k1;
            w1k1 = O1k1 + delO1*rand;
            wk = w1k1;
            %%%In one dimension Pmn = P11, thus THETA = S = P11*P11 = P11^2
            %%%So P11 = sqrt(S)
            P = sqrt(S(Ok));
            x(ii) = x(ii) + abs(P)*sqrt(2*dO)*cos(wk*t(ii) + phin(k1));
        end
    end
    
    plot(t,x,colors{ll},'LineWidth',2)
    drawnow
    end
    
end

if MAKE_R || MAKE_S
    %%%Create R from S
    %Rxx(tau) = int S(w) e^-iw tau dw

    %%%Analytic Expression
    %R_analytic = (4*S0./tau).*cos(0.5*(w1+w2).*tau).*sin(0.5*(w2-w1).*tau);
    R_analytic = 2*S0./tau.*(sin(w2.*tau)-sin(w1.*tau));

    %%%Numerical and Expectation Operator
%     R_numerical = 0*R_analytic;
%     for ii = 1:length(tau)
%          %%%Numerical
%          for wdx = 1:length(w)
%              wi = w(wdx);
%              R_numerical(ii) = R_numerical(ii) + real(dw*S(wi)*exp(-1i*wi*tau(ii)));
%          end
%          disp(tau(ii))
%      end
    
    %%%%Expectation Operator
    R_interp = 0.*tau;
    twindow = 50; %percent
    sides = round((twindow/100)*length(t));
    s = find(tau >= -t(1),1);
    e = find(tau >= t(end)-t(end-sides),1)-1;
     tlength = (t(end-sides)-t(1));
     for ii = s:e
         for tdx = 1:(length(t)-sides)
             xinterp = interp1(t,x,t(tdx)+tau(ii));
             R_interp(ii) = R_interp(ii) + dt*x(tdx)*xinterp;
         end
         disp(tau(ii))
     end
     R0 = 0;
     for ii = 1:length(w)
         R0 = R0 + dw*S(w(ii));
     end
     scale = mean(x(1:(length(t)-sides)).^2);
     R_interp = (R0/tlength)*R_interp/scale;

     plottool(1,'R',18,'\tau','R');
     plot(tau,R_analytic,'b-','LineWidth',2);
%    plot(tau,R_numerical,'r-','LineWidth',2);
     plot(tau,R_interp,'g-','LineWidth',2)
     legend('Analytic','Expectation')
%    legend('Analytic','Numerical')
xlim([tau(s) tau(e)])

end

%%%Recreate S from R
if MAKE_S
    
    S_analytic = 0*w;
    S_numerical = 0*w;
    for ii = 1:length(w)
        S_analytic(ii) = S(w(ii));
         for taudx = 1:length(tau)
             taut = tau(taudx);
             R = R_analytic(taudx);
             if ~isnan(R)
                 S_numerical(ii) = S_numerical(ii) + (1/(2*pi))*real(dtau*R_analytic(taudx)*exp(-1i*w(ii)*taut));
             end
         end
         disp(w(ii))
   end
    
    plottool(1,'S',18,'\Omega_x(rad/s)','\Theta');
    plot(w,S_analytic,'b-','LineWidth',2);
    plot(w,S_numerical,'r-','LineWidth',2);
    legend('Analytic','Numerical')
    
end

keyboard

function val = S(w)
global w1 w2 S0
%%%Assume S is a narrow band process

if (-w2 < w && w < -w1) || (w1 < w && w < w2)
    val = S0;
else
    val = 0;
end
%val = sin(w).^2;
