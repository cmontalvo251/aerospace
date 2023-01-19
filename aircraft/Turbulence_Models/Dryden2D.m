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
P = 10;
xlimit = 10;
x = linspace(0,xlimit,P);
y = linspace(0,xlimit,P);
h = 1;
N = [100,100];
wupper = [5,5];
wlower = [-5,-5];

%%1D
[r,P] = size(x);
[xx,yy] = meshgrid(x,y);
u = zeros(P,P);
v = u;
S = [0 0;0 0];
tic

%%Generate Length Scales and Variance
sigu = 1;
Lu = h;

%%Generate PHI
PHI1 = zeros(N(1),N(2));
PHI2 = zeros(N(1),N(2));
dwp_rand = zeros(N(1),N(2));
for ii = 1:N(1)
  for jj = 1:N(2)
    PHI1(ii,jj) = 2*pi*rand;
    PHI2(ii,jj) = 2*pi*rand;
    dwp_rand(ii,jj) = rand;
  end
end
dw = (wupper-wlower)./N;
dwp = 0.01.*dw;

%%Generate f
alfa = 4*(sigu^2)*(Lu^5)/(pi^2);
for idx = 1:P
  idx
  for jdx = 1:P
    state = [x(idx);y(jdx)];
    for k1 = 1:N(1)
      w1k1 = wlower(1) + (k1-1/2)*dw(1);
      for k2 = 1:N(2)
	w2k2 = wlower(2) + (k2-1/2)*dw(2);
	delw1 = -dwp(1)/2 + dwp(1)*dwp_rand(k1,k2);
	delw2 = -dwp(2)/2 + dwp(2)*dwp_rand(k1,k2);
	wp = [w1k1+delw1;w2k2+delw2];
	W = sqrt(wp(1)^2+wp(2)^2);
	e = 1/(1+(Lu*W)^2)^3;
	%S = alfa*e*[W^2-wp(1)^2 -wp(1)*wp(2);-wp(1)*wp(2) W^2-wp(2)^2];
	%H = Cholesky(S);
	%f1 = cos(wp'*state + PHI1(k1,k2));
	%u(idx,jdx) = u(idx,jdx) + abs(H(1,1))*f1;
      	%v(idx,jdx) = v(idx,jdx) + abs(H(2,1))*f1;% + abs(H(2,2))*cos(wp'*state + PHI2(k1,k2));;
	f = sqrt(e)*cos(wp'*state+PHI1(k1,k2));
	u(idx,jdx) = u(idx,jdx) + abs(wp(2))*f;
	v(idx,jdx) = v(idx,jdx) + abs(wp(1))*f;
      end
    end
    f = sqrt(alfa*2*dw(1)*dw(2));
    u(idx,jdx) = f*u(idx,jdx);
    v(idx,jdx) = f*v(idx,jdx);
  end
end
plottool(1,'Name',12,'x','y','u');
mesh(xx,yy,u)
plottool(1,'Name',12,'x','y','v');
mesh(xx,yy,v)
drawnow

SpectralDensity

%%Regenerate Spectral Density
delta1 = (x(2)-x(1));
delta2 = (y(2)-y(1));
S11regen = 0.*S11;
for ii = 1:N(1)
  ii
  w1k1 = wlower(1) + (ii-1/2)*dw(1);
  for jj = 1:N(2)
    w2k2 = wlower(2) + (jj-1/2)*dw(2);
    w = [w1k1;w2k2];
    for idx = 1:P
      for jdx = 1:P
	state = [x(idx);y(jdx)];
	S11regen(ii,jj) = S11regen(ii,jj) + delta1*delta2/(pi^2)*u(idx,jdx)^2*cos(w'*state);
      end
    end
  end
end
mesh(ww1,ww2,S11regen)


toc