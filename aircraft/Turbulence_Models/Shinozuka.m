purge

%%Method 1
%%Inputs
% N = 500;
% x1 = 1;
% x2 = -1;
% x3 = 2;
% x = [x1;x2];

% %%Compute sigma
% o = 0;
% %So = 5/exp(-w1^2 - w2^2);
% w1 = -5:0.1:5;
% w2 = -5:0.1:5;
% dw1 = w1(2)-w1(1);
% dw2 = w2(2)-w2(1);
% [ww1,ww2] = meshgrid(w1,w2);
% So = ww1;
% for ii = 1:length(w1)
%   w1k = w1(ii);
%   for jj = 1:length(w2)
%     w2k = w2(jj);
%     So(ii,jj) = 5/(exp(w1k^2 + w2k^2));
%     o = o + So(ii,jj)*dw1*dw2;
%   end
% end
% mesh(ww1,ww2,So)

% f = 0;
% for k = 1:N
%   phik = 2*pi*rand;
%   w1k = 0; %%How do we pick these random frequencies??
%   w2k = 0;
%   f = f + cos(w1k*x1 + w2k*x2 + phik);
% end
% f = f*o*sqrt(2/N)

%%Method 2
% n = 2;
% %%pick upper and lower bounds
% w1u = 5;
% w2u = 5;
% w1l = -5;
% w2l = -5;
% wu = [w1u,w2u];
% wl = [w1l,w2l];
% %pick N
% N1 = 100;
% N2 = 100;
% N = [N1,N2];
% dw = (wu-wl)./N;
% dwp = 0.01.*dw;
% %Pick x location
% Nt = 1;
% Nx = 100;
% f = zeros(Nt,Nx);
% x = 10*rand(2,Nx);
% tic
% for ii = 1:Nt
%   %%Compute f
%   f(ii,:) = 0;
%   for k1 = 1:N1
%     w1k1 = w1l + (k1-1/2)*dw(1);
%     for k2 = 1:N2
%       w2k2 = w2l + (k2-1/2)*dw(2);
%       delw1 = -dwp(1)/2 + dwp(1)*rand;
%       delw2 = -dwp(2)/2 + dwp(2)*rand;
%       w1p = [w1k1+delw1;w2k2+delw2];
%       So = 5/(exp(w1k1^2 + w2k2^2));
%       phik = 2*pi*rand;
%       for jj = 1:Nx
% 	f(ii,jj) = f(ii,jj) + sqrt(So*dw(1)*dw(2))*cos(w1p'*x(:,jj) + phik);
%       end
%     end
%   end
%   f(ii,:) = f(ii,:)*sqrt(2);
% end
% toc
% plot(f)

%%Method 3 Compute Sqrt(So*dw) ahead of time.
% n = 2;
% %%pick upper and lower bounds
% w1u = 5;
% w2u = 5;
% w1l = -5;
% w2l = -5;
% wu = [w1u,w2u];
% wl = [w1l,w2l];
% %pick N
% N1 = 100;
% N2 = 100;
% N = [N1,N2];
% dw = (wu-wl)./N;
% dwp = 0.01.*dw;

% %%Compute Sodw
% Sodw = 0;
% for k1 = 1:N1
%   w1k1 = w1l + (k1-1/2)*dw(1);
%   for k2 = 1:N2
%     w2k2 = w2l + (k2-1/2)*dw(2);
%     So = 5/(exp(w1k1^2 + w2k2^2));
%     Sodw = Sodw + sqrt(So*dw(1)*dw(2))*dw(1)*dw(2);
%   end
% end

% %Pick x location
% f = zeros(100,3);
% tic
% for ii = 1:100
%   x1 = [-1;1];
%   x2 = [0;1];
%   x3 = [1;1];
%   %%Compute f
%   f(ii,:) = 0;
%   x = [x1;x2];
%   for k1 = 1:N1
%     w1k1 = w1l + (k1-1/2)*dw(1);
%     for k2 = 1:N2
%       w2k2 = w2l + (k2-1/2)*dw(2);
%       delw1 = -dwp(1)/2 + dwp(1)*rand;
%       delw2 = -dwp(2)/2 + dwp(2)*rand;
%       w1p = [w1k1+delw1;w2k2+delw2];
%       phik = 2*pi*rand;
%       f(ii,1) = f(ii,1) + cos(w1p'*x1 + phik);
%       % delw1 = -dwp(1)/2 + dwp(1)*rand;
%       % delw2 = -dwp(2)/2 + dwp(2)*rand;
%       % w1p = [w1k1+delw1;w2k2+delw2];
%       % phik = 2*pi*rand;
%       f(ii,2) = f(ii,2) + cos(w1p'*x2 + phik);
%       % delw1 = -dwp(1)/2 + dwp(1)*rand;
%       % delw2 = -dwp(2)/2 + dwp(2)*rand;
%       % w1p = [w1k1+delw1;w2k2+delw2];
%       % phik = 2*pi*rand;
%       f(ii,3) = f(ii,3) + cos(w1p'*x3 + phik);
%     end
%   end
%   f(ii,1) = f(ii,1)*sqrt(2)*Sodw;
%   f(ii,2) = f(ii,2)*sqrt(2)*Sodw;
%   f(ii,3) = f(ii,3)*sqrt(2)*Sodw;
%   ii
% end
% toc
% hold on
% plot(f)

%%Method 4 - Compute all random numbers ahead of time
n = 2;
%%pick upper and lower bounds
w1u = 5;
w2u = 5;
w1l = -5;
w2l = -5;
wu = [w1u,w2u];
wl = [w1l,w2l];
%pick N
N1 = 100;
N2 = 100;
N = [N1,N2];
dw = (wu-wl)./N;
dwp = 0.01.*dw;
%Pick x location
Nx = 100;
x = 10*rand(2,Nx);
f = zeros(1,Nx);
for k1 = 1:N1
  w1k1(k1) = w1l + (k1-1/2)*dw(1);
  for k2 = 1:N2
    w2k2(k2) = w2l + (k2-1/2)*dw(2);
    delw1(k1,k2) = -dwp(1)/2 + dwp(1)*rand;
    delw2(k1,k2) = -dwp(2)/2 + dwp(2)*rand;
    w1p1(k1,k2) = w1k1(k1)+delw1(k1,k2);
    w2p2(k1,k2) = w2k2(k2)+delw2(k1,k2);
    w1p = [w1p1(k1,k2);w2p2(k1,k2)];
    So(k1,k2) = 5/(exp(w1k1(k1)^2 + w2k2(k2)^2));
    phik(k1,k2) = 2*pi*rand;
  end
end
p = 0.*f;
tic
for k1 = 1:N1
  for k2 = 1:N2
    Sodw = sqrt(So(k1,k2)*dw(1)*dw(2));
    w1p = [w1p1(k1,k2);w2p2(k1,k2)];
    for jj = 1:Nx
      p(:,jj) = p(:,jj) + Sodw*cos(w1p'*x(:,jj) + phik(k1,k2));
    end
  end
end
toc
plot(p)