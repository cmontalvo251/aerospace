clear
close all
clc

t = linspace(1, 200, 1000);

C = tf([1, 226, 11475, 11250],[1, 1200, 470000, 60000000]);
%zeroes = -1, -75, -150
%poles = -300, -400, -500
A = tf([1.316],[1, 1.316]);
G = tf([1883],[1, 0, 0]);
H = tf([1000],[1, 1000]);

figure()
bode(G)
margin(G)

figure()
rlocus(G)

GOL = A*G*H;

figure()
bode(GOL)
margin(GOL)


figure()
rlocus(GOL)

%controlSystemDesigner('bode',GOL)

GCL = C*A*G*H;

figure()
bode(GCL)
margin(GCL)


figure()
rlocus(GCL)

figure()
pzmap(GCL)
grid on


%controlSystemDesigner('bode',GCL)

%State Space


num = [2.478e6, 5.6e8, 2.844e10, 2.788e10];
den = [1, 2201, 1.673e6, 5.322e8, 6.07e10, 7.896e10, 0, 0];

[W,X,Y,Z] = tf2ss(num,den) 


% y = 0.35309*t - 0.00941695*exp(-997.983*t) + 0.131082*exp(-506.407*t) - 0.197895*exp(-393.823*t) + 0.0710376*exp(-301.471*t) - 0.083554*exp(-1.31594*t) + 0.0887465;
% 
% figure()
% plot(t, y)
% grid on
