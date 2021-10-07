function [BxI,ByI,BzI,B500] = magnetic_field(x,y,z,r)
%%%Make dummy variables
BxI = 0*x;
ByI = BxI;
BzI = BxI;

%%%Call the magnetic field model
%%%Convert Cartesian x,y,z into Lat,Lon, Alt
phiE = 0;
thetaE = acos(z./r);
psiE = atan2(y,x);
latitude = 90-thetaE*180/pi;
longitude = psiE*180/pi;
rkm = (r)/1000;
disp(['Computing Magnetic Field for ',num2str(length(x)),' points....'])
for idx = 1:length(x)
	[BN,BE,BD] = igrf('28-Sep-2019',latitude(idx),longitude(idx),rkm(idx),'geocentric'); 
	%%%Convert NED (North East Down to X,Y,Z in ECI frame)
	%%%First we need to create a rotation matrix from the NED frame to the 
	%%%inertial frame
	BNED = [BN;BE;BD]; 
	BI = TIB(phiE,thetaE(idx)+pi,psiE(idx))*BNED;
	BxI(idx) = BI(1);
	ByI(idx) = BI(2);
	BzI(idx) = BI(3);
end
constants
[BN,BE,BD] = igrf('28-Sep-2019',0,0,500+REarth/1000,'geocentric');
B500 = sqrt(BN^2+BE^2+BD^2);
disp('Done Computing Magnetic Field')