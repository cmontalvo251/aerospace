%%%Update Rate
nextSensorUpdate = 1.0;
f=1;
%%%Bias and Noise
MagscaleBias = (4e-7)*f; %%T
MagFieldBias = MagscaleBias*(2*rand()-1); %%% -1 to 1 (0 to 1) (0 to 2

MagscaleNoise = (1e-5)*f; %%T
MagFieldNoise = MagscaleNoise*(2*rand()-1); %%% -1 to 1 (0 to 1) (0 to 2

AngscaleBias = 0.01*f; %%rad/s
AngFieldBias = AngscaleBias*(2*rand()-1); %%% -1 to 1 (0 to 1) (0 to 2

AngscaleNoise = 0.001*f; %%rad/s
AngFieldNoise = AngscaleNoise*(2*rand()-1); %%% -1 to 1 (0 to 1) (0 to 2

EulerBias = 2*pi/180*f; %%rad
EulerNoise = EulerBias*(2*rand()-1);