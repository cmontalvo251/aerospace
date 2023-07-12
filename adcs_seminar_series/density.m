function d = density(altitude)
ds = 1.225; %%kg/m^3
sigma = 0.1354; %%scale height;
d = ds * exp(-sigma*altitude/1000);