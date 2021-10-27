clear
clc
close all

%%%Time vector
dt = 0.01;
time = 0:dt:60;

%%%Initial Conditions
xstate = zeros(12,length(time));

x0 = 0;
y0 = 0;
z0 = 1;
phi0 = 0;
theta0 = 0;
psi0 = 0;
u0 = 0;
v0 = 0;
w0 = 0;
p0 = 0;
q0 = 0;
r0 = 0;

A = zeros(12,12);
xvec0 = [x0;y0;z0;phi0;theta0;psi0;u0;v0;w0;p0;q0;r0];

dx = 1e-4;

for idx = 1:12
    xvec0(idx) = xvec0(idx) + dx;
    f2 = Derivatives490Quad(xvec0,0);
    xvec0(idx) = xvec0(idx) - 2*dx;
    f1 = Derivatives490Quad(xvec0,0);
    xvec0(idx) = xvec0(idx) + dx;
    A(:,idx) = (f2-f1)/(2*dx);
end    

A;
[s,v] = eig(A);

% break

xstate(:,1) = xvec0;

% figure();

for k = 1:length(time)-1
    %%%Runge_Kutta
    xk = xstate(:,k);
    tk = time(k);
    k1 = Derivatives490Quad(xk,tk);
    k2 = Derivatives490Quad(xk+k1*dt/2,tk+dt/2);
    k3 = Derivatives490Quad(xk+k2*dt/2,tk+dt/2);
    k4 = Derivatives490Quad(xk+k3*dt,tk+dt);
    phi = (1/6)*(k1 + 2*k2 + 2*k3 + k4);
    xstate(:,k+1) = xstate(:,k) + phi*dt;
%     cla;
%     x = xk(1);
%     y = xk(2);
%     z = xk(3);
%     phi = xk(4);
%     theta = xk(5);
%     psi = xk(6);
%     CubeDraw(2,2,0.2,x,y,z,-phi,-theta,psi,[1 0 0])
%     drawnow
%     grid on
%     xlim([x-10 x+10])
%     ylim([y-10 y+10])
%     zlim([z-10 z+10])
%     view(-28,22)
%     axis equal
end

ylabels = {'x','y','z','\phi','\theta','\psi','u','v','w','p','q','r'}';

for idx = 1:12
    fig = figure();
    plot(time,xstate(idx,:))
    xlabel('Time (sec)')
    ylabel(ylabels{idx})
end

% %%%Angular Momentum = I*w
% Ixxp = 1*xstate(10,:);
% Iyyq = 2*xstate(11,:);
% Izzr = 3*xstate(12,:);
% 
% Hmag = sqrt(Ixxp.^2 + Iyyq.^2 + Izzr.^2);
% 
% figure
% plot(time,Hmag)
% ylim([0 10])
