%%%Update Rate
nextSensorUpdate = 1.0;
fsensor=1;
%%%Bias and Noise
MagscaleBias = (4e-7)*fsensor; %%T
MagFieldBias = MagscaleBias*(2*rand()-1); %%% -1 to 1 (0 to 1) (0 to 2

AngscaleBias = 0.01*fsensor; %%rad/s
AngFieldBias = AngscaleBias*(2*rand()-1); %%% -1 to 1 (0 to 1) (0 to 2

EulerBias = 2*pi/180*fsensor; %%rad
EulerBias = EulerBias*(2*rand()-1);