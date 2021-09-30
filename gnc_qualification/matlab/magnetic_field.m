function [BxI,ByI,BzI] = magnetic_field(x,y,z)
%%%Make dummy variables
BxI = 0*x;
ByI = BxI;
BzI = BxI;

%%%Distance from the center of the Earth
rho = sqrt(x.^2+y.^2+z.^2);

%%%Call the magnetic field model
%%%Convert Cartesian x,y,z into Lat,Lon, Alt
phiE = 0;
thetaE = acos(z./rho);
psiE = atan2(y,x);
latitude = 90-thetaE*180/pi;
longitude = psiE*180/pi;
rhokm = (rho)/1000;
disp(['Computing Magnetic Field for ',num2str(length(x)),' points....'])
for idx = 1:length(x)
	[BN,BE,BD] = igrf('28-Sep-2019',latitude(idx),longitude(idx),rhokm(idx),'geocentric'); 
	%%%Convert NED (North East Down to X,Y,Z in ECI frame)
	%%%First we need to create a rotation matrix from the NED frame to the 
	%%%inertial frame
	BNED = [BN;BE;BD]; 
	BI = TIB(phiE,thetaE(idx)+pi,psiE(idx))*BNED;
	BxI(idx) = BI(1);
	ByI(idx) = BI(2);
	BzI(idx) = BI(3);
end
disp('Done')