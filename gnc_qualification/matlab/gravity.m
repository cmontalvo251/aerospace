function Mgrav = gravity(r,max_moment_arm,m)
constants

%%%Compute Gravity at the bottom and top
bottom = r - max_moment_arm/2;
top = r + max_moment_arm/2;
accel_bottom = (G*MEarth)./(bottom.^2);
accel_top = (G*MEarth)./(top.^2);
delta_accel = (accel_bottom-accel_top);
force_delta = delta_accel * m;
Mgrav = force_delta * max_moment_arm;
