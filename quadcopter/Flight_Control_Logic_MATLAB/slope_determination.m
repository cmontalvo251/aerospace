%%%%Figuring out slope to ms impulse
clear
clc
close all

%%%Mass of quad
m = 0.6; %mass of quad kg
Weight = m*9.81

%%%Rotors off 
Nominal_weight = 941*9.81/1000;

Force_On_Scale = [796 670 510 330 200]*9.81/1000;

Throttle_impulse = [1170 1180 1195 1204 1215];

Thrust_from_rotors = -(Force_On_Scale-Nominal_weight);

plot(Throttle_impulse,Thrust_from_rotors)

xlabel('us')

ylabel('N')


