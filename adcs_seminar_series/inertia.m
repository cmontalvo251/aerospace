%%%inertia and mass
ms = 2.6; %%kilograms
lx = 10/100; %%meters
ly = 10/100; %%meters
lz = 20/100; %%meters
l = [lx;ly;lz];
lsort = sort(l);
Amax = lsort(2)*lsort(3);
lmax = l(3);
CD = 1.0;
%%%Inertia of satellite in kg-m^2
Is = (ms/12)*[(ly^2+lz^2) 0 0 ;0 (lx^2+lz^2) 0;0 0 (lx^2+ly^2)];

%%%Call the reaction wheel params
reaction_wheel_params

%%%Add everything up
m = ms + 3*mr;
I = Is + Ir1Bcg + Ir2Bcg + Ir3Bcg;

%%%Invert the matrix
invI = inv(I);

