%%%inertia and mass
ms = 2.6; %%kilograms

%%%Inertia of satellite in kg-m^2
Is = [0.9 0 0 ;0 0.9 0;0 0 0.3];

%%%Call the reaction wheel params
reaction_wheel_params

%%%Add everything up
m = ms + 3*mr;
I = Is + Ir1Bcg + Ir2Bcg + Ir3Bcg;

%%%Invert the matrix
invI = inv(I);

