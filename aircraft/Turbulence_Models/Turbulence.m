purge

%%Plot some spectra

%%Inputs
hm = 533; %meters
O = (0:0.001:(2.5*3.28)); %rad/m
u20m = 1; %speed at 20 feet (m/s)

h = hm*3.28; %convert to feet
u20 = u20m*3.28;

if (h > 1000)
  sigu = 1;
  sigw = 1;
  sigv = 1;
  Lu = 1750;
  Lv = Lu;
  Lw = Lu;
else
  %%Compute Length Scales
  Lw = h;
  Lv = h/(0.177+0.000823*h)^1.2;
  Lu = Lv;
  %%Compute Turbulence Scales
  sigw = 0.1*u20;
  sigu = sigw/(0.177+0.000823*h)^0.4;
  sigv = sigu;
end

O = O./3.28; %%Convert to 1/feet

PHIU = (2*(sigu^2)*Lu./(pi)).*(1./(1+(Lu.*O).^2));

plottool(1,'PHIU',12,'Omega(mrad/ft)','PHIU(Omega)(ft^3/s^2)');
plot(O.*1000,PHIU)

%%Spectral Density of Random Input
N1 = 1;
dx = 10;
PHIO = ((N1^2)./(dx.*pi.*O.^2)).*(1-cos(dx.*O));
PHIAPPROX = (N1^2).*dx./(2*pi).*ones(length(O),1);

plottool(1,'PHIO',12,'Omega(rad/ft)','PHI(Omega)(ft^3/s^2)');
plot(O,PHIO)
plot(O,PHIAPPROX,'r-')

%%Correlation functions
zeta1 = -20:1:20;
zeta2 = -20:1:20;
dy = 10;
R = ones(length(zeta1),length(zeta2));
for ii = 1:length(zeta1)
  for jj = 1:length(zeta2)
    R(ii,jj) = (N1^2)*(1-abs(zeta1(ii)/dx))*(1-abs(zeta2(jj)/dy));
    if abs(zeta1(ii)) >= dx || abs(zeta2(jj)) >= dy
      R(ii,jj) = 0;
    end
  end
end
[zz1,zz2] = meshgrid(zeta1,zeta2);
plottool(1,'R',12)
mesh(zz1,zz2,R)