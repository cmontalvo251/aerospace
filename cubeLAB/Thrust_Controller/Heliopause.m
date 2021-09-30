
%%%Calibration Exp
force_tether_length = 90e-9; %nN/meter
length = 16000;
force_tether = force_tether_length*length
force_tether = 0.0144
mass = 24;
tether_accel = force_tether/mass;

%%%% xfinal = 0.5*a*t^2
helio_distance = 0.5*tether_accel*(10*365*24*3600)^2;

%%% Estimate Time for Solar Sail
force = 720e-6;
mass = 14; 
solar_sail_accel = force/mass;
time_solar_sail = sqrt(helio_distance*2/solar_sail_accel)
time_solar_sail_years = (((time_solar_sail/3600)/24)/365)