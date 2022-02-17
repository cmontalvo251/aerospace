%%%Initialize
clear
clc
close all

%%%Setup the IGRF Model
disp(['You must download the igrf model from MathWorks in order to' ...
      ' use this software'])
disp('https://www.mathworks.com/matlabcentral/fileexchange/34388-international-geomagnetic-reference-field-igrf-model')
addpath '../igrf/'

%%%Get Planet Parameters
planet
%% R is the radius of the planet

distance_from_surface = linspace(0,60000,100)*1000; %Times 1000 to get to km

%%%Loop through stateout to extract Magnetic Field
BxIout = 0*distance_from_surface;
ByIout = BxIout;
BzIout = BxIout;
for idx = 1:length(distance_from_surface)
	%%%Call the magnetic field model
	%%%Convert Cartesian x,y,z into Lat,Lon, Alt
	phiE = 0;
	z = 0; %%assume that z and y are zero
	y = 0;
	x = distance_from_surface(idx)+R;
	r = [x;y;z];
	rho = norm(r);
	thetaE = acos(z/rho);
	psiE = atan2(y,x);
	latitude = 90-thetaE*180/pi;
	longitude = psiE*180/pi;
	rhokm = (rho)/1000;
	[BN,BE,BD] = igrf('01-Jan-2020',latitude,longitude,rhokm,'geocentric');
    BxIout(idx) = BN;
    ByIout(idx) = BE;
    BzIout(idx) = BD;
end

%%%Plot Magnetic Field
fig2 = figure();
set(fig2,'color','white')
plot(distance_from_surface/1000,BxIout,'b-','LineWidth',2)
hold on
grid on
plot(distance_from_surface/1000,ByIout,'y-','LineWidth',2)
plot(distance_from_surface/1000,BzIout,'g-','LineWidth',2)
xlabel('Distance From Surface (km)')
ylabel('Mag Field (nT)')
legend('X','Y','Z')

%%%And Norm
Bnorm = sqrt(BxIout.^2 + ByIout.^2 + BzIout.^2);

%%%Now we want to fit an empirical model
%%% Bnorm_fit = A*exp(-sd)
%%% when d = 0 Bnorm_fit(d=0) = A
A = Bnorm(1)
%%% Bnorm/A = exp(-s*d)
%%% ln(Bnorm/A) = -s*d
%%% -ln(Bnorm/A)/d = s
s = mean(-log(Bnorm(2:end)./A)./distance_from_surface(2:end))
Bnorm_fit_exp = A*exp(-s*distance_from_surface);

%%%Ok whoops. It's actually r^2
%%%Bnorm_fit = B/r^2
%%%when r = R , Bnorm(r=R) = Bnorm(1) = B/R^2
B = Bnorm(1)*R^2
Bnorm_fit_r2 = B./(distance_from_surface+R).^2;

%%%Let's try r^3
C = Bnorm(1)*R^3
Bnorm_fit_r3 = C./(distance_from_surface+R).^3;

fig3 = figure();
set(fig3,'color','white')
plot(distance_from_surface/1000,Bnorm,'LineWidth',2)
hold on
plot(distance_from_surface/1000,Bnorm_fit_exp,'r--','LineWidth',2)
plot(distance_from_surface/1000,Bnorm_fit_r2,'k-.','LineWidth',2)
plot(distance_from_surface/1000,Bnorm_fit_r3,'g-','LineWidth',2)
grid on
xlabel('Distance From Surface (m)')
ylabel('Mag Field (nT)')

%%%%Let's assume we are JagSAT
N = 84
I = 0.04
A = (20./100)*(20/100.)
mu = N*I*A
bfield_strength = C/(400*1000+R)^3
bT = bfield_strength*10^(-9)
torque = bT*mu

%%%GTO
bfield_strength = C/(35000*1000+R)^3
bT = bfield_strength*10^(-9)
torque = bT*mu