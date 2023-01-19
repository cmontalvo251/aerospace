%function uvw = Dryden(x,N,wupper,wlower,h)
% function uvw = Dryden(x,N,wupper,wlower,h)
% x is a 3xP vector or points to sample
% wupper is a vector containing the upper limits of the waves
% wlower is a vector containing the lower limits of the waves
% Ni is the number of intervals such that the waves sampled are:
% wi = linspace(wlower,wupper,Ni)
% h is the height above ground
%
% The Dryden model uses the energy function
%
% E(W) = 8*o^2/(pi)*(1/(1+(L*W)^2)^3)
%
% This is placed in the Spectral Density Function
%
% Sij(w) = E(W)/(4*pi*W^4)*(W^2*Kdij-wi*wj)
%
% Note W = norm(w) where w = [w1,....wm] and Kdij is the Kronecker
% delta function
%
% This is then factorized into H(w) where
%
% S(w) = H*conj(H)^T
%
% Finally the system is realized for a given set of random numbers
%
% Using Shinozuka we arrive at the equation for uvw
%
% uvw = [f1(x),f2(x),...,fm(x)];
%
% N=(N1,...,Nm)
%
% fj(x) = sum(m=1,j)sum(l=1,N)norm(Hjm(wl))sqrt(2*dw)*cos(wlp'*x+Ojm(wl)+PHIml)
%
% wlp = wl + delwl
%
% delwl = random number between -dwl/2,dwl/2
%
% wl = wlower + (l-1/2)*dwl
%
% dw = dw1*...*dwm
%
% Ojm(wl) = atan(Im(Hjm(wl))/Re(Hjm(wl)))
%
% PHIml = random number between 0,2pi

%%Example problem
purge
P = 10000;
x = linspace(0,100,P);
h = 1;
N = 100;
wupper = 5;
wlower = -5;

%%1D
[r,P] = size(x);
f = 0.*x;
S = 0;
tic

%%Generate Length Scales and Variance
sigu = 1;
Lu = h;

%%Generate PHI
PHI = zeros(N(1),1);
dwp_rand = zeros(N(1),1);
for ii = 1:N(1)
  PHI(ii) = 2*pi*rand*0;
  dwp_rand(ii) = rand*0;
end
dw = (wupper-wlower)./N;
dwp = 0.01.*dw;

%%Generate f
for idx = 1:P
  for k1 = 1:N(1)
    w1k1 = wlower(1) + (k1-1/2)*dw(1);
    delw1 = -dwp(1)/2 + dwp(1)*dwp_rand(k1,1);
    wp = [w1k1+delw1];
    W = sqrt(wp(1)^2);
    S(1,1) = ((sigu^2)*Lu/(pi))*1/(1+(Lu*W)^2);
    H = Cholesky(S);
    magH = abs(H);
    f(idx) = f(idx) + magH*cos(wp*x(idx) + 0*PHI(k1));
  end
  f(idx) = sqrt(2*dw)*f(idx);
end
plottool(1,'Name',12,'x','f');
plot(x,f)

%%Regenerate R
% R = 0.*f;
% Rregen = 0.*f;
% for idx = 1:P-1
%   zed = x(:,idx)-x(:,1);
%   dzed = x(:,idx+1)-x(:,idx);
%   R(idx) = exp(-zed./Lu);
%   for jdx = 1:P
%     xzed = x(:,jdx) + zed;
%     locxzed = find(xzed <= [x,x+x(end)],1);
%     if locxzed > P
%       locxzed = locxzed - P;
%     end
%     frzed = f(locxzed);
%     fr = f(jdx);
%     Rregen(idx) = Rregen(idx) + fr*frzed;
%   end
% end
% Rregen = Rregen./Rregen(1);
% plottool(1,'Name',12,'x','R');
% plot(x,R)
% plot(x,Rregen,'r-')
% legend('R Actual','R regenerated')

%%Regenerate S
% R(zeta) = exp(-zeta/L);
Svec = 0.*PHI;
Wvec = Svec;
Sregen = Svec;
Sregen2 = Svec;
Rvec = Svec;
zvec = Svec;
delta = x(2)-x(1);
for ii = 1:N(1)
  w1k1 = wlower(1) + (ii-1/2)*dw(1);
  w = [w1k1];
  W = sqrt(w(1)^2);
  Wvec(ii) = w1k1;
  Svec(ii) = ((sigu^2)*Lu/(pi))*1/(1+(Lu*W)^2);
  for idx = 1:P-1
    zed = x(:,idx);
    zvec(idx) = zed;
    dzed = x(:,idx+1)-x(:,idx);
    R = exp(-abs(zed)./Lu);
    Rvec(idx) = R;
    Sregen2(ii) = Sregen2(ii) + (delta/(2*pi))*R*exp(-w*i*zed);
    Sregen(ii) = Sregen(ii) + (delta/(pi))*f(idx)*cos(w*x(idx))/(x(end)-x(1));
  end
  Sregen2(ii) = norm(Sregen2(ii));
end
% Sregen2 = fft(f,N(1));
% Sregen2 = [Sregen2(round(end/2)+1:end),Sregen2(1:round(end/2))];
% Sregen2 = (Sregen2(1:N(1)));
% Sregen2 = ((Sregen2.*conj(Sregen2)))./N(1);
plottool(1,'Name',12,'W','S');
plot(Wvec,Svec)
plot(Wvec,Sregen,'r-')
plot(Wvec,Sregen2,'g-')
legend('Actual','Regeneration Using x','Regeneration Using R')

% plottool(1,'Name',12,'W','S');
% plot(Wvec,abs(Svec-Sregen2))
% ylim([0 max(Svec)])

toc