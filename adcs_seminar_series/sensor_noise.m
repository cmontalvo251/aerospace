%%%Sensor Noise Parameters
global fsensor

MagscaleNoise = (1e-5)*fsensor; %%T
MagFieldNoise = MagscaleNoise*(2*rand()-1); %%% -1 to 1 (0 to 1) (0 to 2

AngscaleNoise = 0.001*fsensor; %%rad/s
AngFieldNoise = AngscaleNoise*(2*rand()-1); %%% -1 to 1 (0 to 1) (0 to 2

EulerScaleNoise = 1*pi/180*fsensor;
EulerNoise = EulerScaleNoise*(2*rand()-1);