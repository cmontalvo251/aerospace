%%%What does the atan2 function output?

purge


argx = -1000:10:1000;
argy = argx;

[xx,yy] = meshgrid(argx,argy);

zz = atan2(yy,xx);

mesh(xx,yy,zz)

figure()

delpsi = -pi:0.1:pi;

f = (pi-abs(delpsi))/pi;


plot(delpsi,f)





