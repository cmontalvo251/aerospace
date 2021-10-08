function Mrad = solar_radiation(r,max_area,max_moment_arm)
constants

solar_rad_force = solar_rad_pressure*max_area;
Mrad = 0*r + solar_rad_force*max_moment_arm;

