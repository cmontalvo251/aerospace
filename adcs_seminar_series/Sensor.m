function [BB,pqr,ptp] = Sensor(BB,pqr,ptp)
global fsensor MagFieldBias AngFieldBias EulerBias
global MagFieldNoise AngFieldNoise EulerNoise

for idx = 1:3
    %%%Get our sensor params
    sensor_noise
    %%%Pollute the data
    BB(idx) = BB(idx) + MagFieldBias + MagFieldNoise;
    pqr(idx) = pqr(idx) + AngFieldBias + AngFieldNoise;
    ptp(idx) = ptp(idx) + EulerBias + EulerNoise;
end


