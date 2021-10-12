function maneuver_time(Power,MaxT,maneuver_angle,f,Ixx,Iyy,Izz)

%%Maximum Acceleration per axis
alfax = MaxT/Ixx;
alfay = MaxT/Iyy;
alfaz = MaxT/Izz;

%%%Time for first and last manuever
t1x = sqrt(2*maneuver_angle*pi/180*f/alfax);
t1y = sqrt(2*maneuver_angle*pi/180*f/alfay);
t1z = sqrt(2*maneuver_angle*pi/180*f/alfaz);

%%%Total Time
tx = 2*t1x+maneuver_angle*pi/180*(1-2*f)/(alfax*t1x);
ty = 2*t1y+maneuver_angle*pi/180*(1-2*f)/(alfay*t1y);
tz = 2*t1z+maneuver_angle*pi/180*(1-2*f)/(alfaz*t1z);

disp(['Total Time For Maneuver Through ',num2str(maneuver_angle),['' ...
		    ' degrees (sec) = ',num2str(tx),' ',num2str(ty),['' ...
		    ' ',num2str(tz)]]])

%%Energy for each axis
Ex = tx*Power;
Ey = ty*Power;
Ez = tz*Power;

disp(['Total Energy Expended Per Maneuver Through ',num2str(maneuver_angle),['' ...
		    ' degrees (Wh) = ',num2str(Ex),' ',num2str(Ey),['' ...
		    ' ',num2str(Ez)]]])
