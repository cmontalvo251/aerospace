function [BB,pqr] = Sensor(BB,pqr)

for idx = 1:3
    %%%Get our sensor params
    sensor_params
    %%%Pollute the data
    BB(idx) = BB(idx) + MagFieldBias + MagFieldNoise;
    pqr(idx) = pqr(idx) + AngFieldBias + AngFieldNoise;
end


