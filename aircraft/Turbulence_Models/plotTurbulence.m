%Dryden3D
purge

%system('compilec Drydencpp.cpp -w -O3');
%system('./Run.exe');

U = dlmread('Uturb.txt');
V = dlmread('Vturb.txt');
W = dlmread('Wturb.txt');

[r,c] = size(U);

x = 1:c;
y = 1:r;
[xx,yy] = meshgrid(x,y);

plottool(1,'U',12,'x','y','U','',[-27 30]);
mesh(xx,yy,U)

plottool(1,'V',12,'x','y','V','',[-27 30]);
mesh(xx,yy,V)

plottool(1,'W',12,'x','y','W','',[-27 30]);
mesh(xx,yy,W)

%%Sample along a line
% close all

% x = 1:500;
% y = U(:,10)';

% temp1 = y(1);
% temp2 = y(end);

% temp = 0.5*(y(1)+y(end));

% y(1) = temp;
% y(end) = temp;

% %%March in from the left
% filt = 0.9;
% l = length(y);
% m = round(0.1*l);
% for ii = 2:m
%   y(ii) = y(ii-1)*(1-filt) + y(ii)*filt;
% end
% %%March in from right
% for ii = (l-1):-1:(l-m)
%   y(ii) = y(ii+1)*(1-filt) + y(ii)*filt;
% end

% plot([x,x+x(end)],[y,y])

% x = 0:0.1:pi/2;
% y = 0:0.1:pi/2;
% [xx,yy] = meshgrid(x,y);
% z = sin(xx.*yy);

z = U;
temp = 0.5.*(z(:,1)+z(:,end));
z(:,1) = temp;
z(:,end) = temp;

temp = 0.5.*(z(1,:)+z(end,:));
z(1,:) = temp;
z(end,:) = temp;
U = z;

z = V;
temp = 0.5.*(z(:,1)+z(:,end));
z(:,1) = temp;
z(:,end) = temp;

temp = 0.5.*(z(1,:)+z(end,:));
z(1,:) = temp;
z(end,:) = temp;
V = z;

z = W;
temp = 0.5.*(z(:,1)+z(:,end));
z(:,1) = temp;
z(:,end) = temp;

temp = 0.5.*(z(1,:)+z(end,:));
z(1,:) = temp;
z(end,:) = temp;
W = z;

% dlmwrite('UturbSmooth.txt',U,'delimiter',' ');
% dlmwrite('VturbSmooth.txt',V,'delimiter',' ');
% dlmwrite('WturbSmooth.txt',W,'delimiter',' ');

% %%March in from left
% filt = 0.5;
% l = length(y);
% m = round(0.5*l);
% for ii = 2:m
%   z(:,ii) = z(:,ii-1)*(1-filt) + z(:,ii)*filt;
% end
% %%March in from right
% for ii = (l-1):-1:(l-m)
%   z(:,ii) = z(:,ii+1)*(1-filt) + z(:,ii)*filt;
% end
% %%Marhc in from bottom
% for ii = 2:m
%   z(ii,:) = z(ii-1,:)*(1-filt) + z(ii,:)*filt;
% end
% %%March in from top
% for ii = (l-1):-1:(l-m)
%   z(ii,:) = z(ii+1,:)*(1-filt) + z(ii,:)*filt;
% end
% x = 1:(2*c);
% y = 1:(2*r);
% [xx,yy] = meshgrid(x,y);
%mesh(xx,yy,[[z,z];[z,z]])