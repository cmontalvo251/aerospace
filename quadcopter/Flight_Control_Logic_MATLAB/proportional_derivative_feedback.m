%%%%Open Loop Control
clear
clc
close all

%%%Parameters
d = 0.1; %distance from cg to prop in meters
m = 0.6; %mass of quad kg
w = 0.1; %width of quad
h = 0.075; %height of quad
J = (m/12)*(w^2 + h^2); %%Moment of inertia
C = ((7.269-1.422)/(1215-1170))/2; %%%N/us %%Divide by 2 for left and right rotors
alfa = 2*d*C/J;

%%% G = alfa/(s^2) - This transfer is 
%%% PHI/DELTAF => PHI - Radians, DELTAF => Newtons

G = tf([alfa],[1 0 0])
%impulse(G)

%%%Feedback Control - Proportional Feedback
%G_KP = alfa*kp/(s^2 + alfa*kp)
%%%inputs = commanded roll angle in radians
%%%output = actual roll in radians
kp = 233;
G_KP = tf([alfa*kp],[1 0 alfa*kp])

%close all
%step(G_KP)

%%%Proportional Derivative Control
kp = 16/alfa
%for kd = 0:0.1:100
kd = 8/alfa
kp = 10;
kd = 2;
%    kd
    G_KP_KD = tf([alfa*kp],[1 alfa*kd alfa*kp])
    step(G_KP_KD)
    %pause
%end



