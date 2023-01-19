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
%purge
xlimit = 10; %%%%How big is our grid
res = 1; %%%and at what resolution
wupper = [1,1,1].*20; %%Then where do we want to truncate the frequencies
wlower = [-1,-1,-1].*20;

%%%THIS IS ALL SET%%%%
x = 0:res:xlimit;
y = 0:res:xlimit;
h = 5000/3.28;
z = 200; %%let z be a constant since it does not change that much
N = [20,20,20]; %%you may as well make this big. You only have to
                %do it once

%%1D
[r,P] = size(x);
[xx,yy] = meshgrid(x,y);
u = zeros(P,P);
v = u;
w = u;
S = [0 0 0;0 0 0;0 0 0];

%%Generate Length Scales and Variance
sigu = 1;
Lu = h;

%%Generate PHI
PHI1 = zeros(N(1),N(2),N(3));
PHI2 = zeros(N(1),N(2),N(3));
PHI3 = zeros(N(1),N(2),N(3));
dwp_rand = zeros(N(1),N(2),N(3));
for ii = 1:N(1)
  for jj = 1:N(2)
    for kk = 1:N(3)
      PHI1(ii,jj,kk) = 2*pi*rand;
      PHI2(ii,jj,kk) = 2*pi*rand;
      PHI3(ii,jj,kk) = 2*pi*rand;
      dwp_rand(ii,jj,kk) = rand;
    end
  end
end
dw = (wupper-wlower)./N;
dwp = 0.01.*dw;
disp('Random Numbers Generated')

tic

%%Generate f
alfa = 4*(sigu^2)*(Lu^5)/(pi^2);
factor = sqrt(2*dw(1)*dw(2)*dw(3));
for idx = 1:P
  for jdx = 1:P
    state = [x(idx);y(jdx);z];
    disp(['X = ',num2str(state(1)),' Y = ',num2str(state(2))])
    for k1 = 1:N(1)
      w1k1 = wlower(1) + (k1-1/2)*dw(1);
      for k2 = 1:N(2)
	w2k2 = wlower(2) + (k2-1/2)*dw(2);
	for k3 = 1:N(3)
	  w3k3 = wlower(3) + (k3-1/2)*dw(3);
	  delw1 = -dwp(1)/2 + dwp(1)*dwp_rand(k1,k2,k3);
	  delw2 = -dwp(2)/2 + dwp(2)*dwp_rand(k1,k2,k3);
	  delw3 = -dwp(3)/2 + dwp(3)*dwp_rand(k1,k2,k3);
	  wp = [w1k1+0*delw1;w2k2+0*delw2;w3k3+0*delw3];
	  W = sqrt(wp(1)^2+wp(2)^2+wp(3)^2);
	  e = 1/(1+(Lu*W)^2)^3;
	  f = alfa*e;
	  %S
	  W2 = W^2;
	  S11 = f*(W2-wp(1)^2);
	  S12 = f*-wp(1)*wp(2);
	  S22 = f*(W2-wp(2)^2);
	  S31 = f*-wp(1)*wp(3);
	  S32 = f*-wp(2)*wp(3);
	  S33 = f*(W2-wp(3)^2);
	  %H
	  H11 = sqrt(S11);
	  H21 = S12/H11;
	  H22 = sqrt(S22-H21^2);
	  H31 = S31/H11;
	  H32 = (S32-H31*H21)/H22;
	  H33 = sqrt(S33-H31^2-H32^2);
	  wpstate = wp'*state;
	  f1 = cos(wpstate+PHI1(k1,k2,k3));
	  f2 = cos(wpstate+PHI2(k1,k2,k3));
	  f3 = cos(wpstate+PHI3(k1,k2,k3));
	  u(idx,jdx) = u(idx,jdx) + abs(H11)*f1;
	  v(idx,jdx) = v(idx,jdx) + abs(H21)*f1 + abs(H22)*f2;
	  w(idx,jdx) = w(idx,jdx) + abs(H31)*f1 + abs(H32)*f2 + abs(H33)*f3;
	end
      end
    end
    u(idx,jdx) = factor*u(idx,jdx);
    v(idx,jdx) = factor*v(idx,jdx);
    w(idx,jdx) = factor*w(idx,jdx);
  end
end
toc
plottool(1,'Name',12,'x','y','u','',[-27,30]);
mesh(xx,yy,u)
plottool(1,'Name',12,'x','y','v','',[-27,30]);
mesh(xx,yy,v)
plottool(1,'Name',12,'x','y','w','',[-27,30]);
mesh(xx,yy,w)
drawnow

