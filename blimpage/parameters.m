%%%%Parameters
clear
clc

b = (38/12)/3.28 %%%diameter of blimp

Vol = (4/3)*pi*(b/2)^3 %%%volume

mstructure = 32/1000; %%%%mass of structure

mgondola = 68/1000; %%%mass of gondola and batteries

density_helium = 0.18; %%%%kg/m^3 of helium

mass = density_helium*Vol + mstructure + mgondola %%%total mass

Inertia =  (2/5)*mass*(b/2)^2 %%%Inertia of a sphere

d = (5/12)/3.28 %%%Moment arm

%%%%Apparent Masses The derivation for apparent mass can be found in
%%%%Tuckerman and Mueller. The airship considered in this paper is a
%%%%special case where the major and minor axes are equal giving an
%%%%eccentricity of zero. Thus the rotational apparent mass is zero (k3=0).
%%%%In addition the translation apparent masses are equal such that k1 = k2
%%%%= k and is given by.
k = (8/3)*(b/2)^3




