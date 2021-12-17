function [Inertia,max_moment_arm,max_area] = inertia(L,W,D,total_mass,mass_solar_panel,number_of_panels,LWD_sp,xyz_sp)
%%This code will compute the inertia of ta CubeSAT with Solar Panels
% L is the length of the satellite in cm
% W is the width of the satellite in cm
% D is the depth of the satellite in cm
% total_mass is the total mass of the satellite including solar panels
% mass_solar_panel is the mass of solar panels
% number_of_panels is the number of solar panels
% LWD_sp is the size of 1 panel in meters
% xyz_sp is the position of the solar panels from the center of mass of the cuboid in meters

max_moment_arm = max([L,W,D]/200);
disp(['Maximum Moment Arm (m) = ',num2str(max_moment_arm)])
area1 = L*W/(100^2);
area2 = W*D/(100^2);
area3 = L*D/(100^2);
max_area = max([area1,area2,area3]);
disp(['Maximum Area (m^2) = ',num2str(max_area)])

%%%Mass/Inertia of Cuboid
mass_cuboid = total_mass - mass_solar_panel*number_of_panels;
Ixx = (mass_cuboid/12)*(L^2+W^2)/(100^2);
Iyy = (mass_cuboid/12)*(L^2+D^2)/(100^2);
Izz = (mass_cuboid/12)*(W^2+D^2)/(100^2);
Inertia_cuboid = [[Ixx,0,0];[0,Iyy,0];[0,0,Izz]];

%%%Inertia of Solar panel about centroid
Length_sp = LWD_sp(1);
Width_sp = LWD_sp(2);
Depth_sp = LWD_sp(3);
Inertia_sp_centroid = mass_solar_panel/12.0 * [[Length_sp^2 + Width_sp^2,0,0];[0,Length_sp^2+Depth_sp^2,0];[0,0,Width_sp^2+Depth_sp^2]];

%disp(['Inertia of Cuboid (kg-m^2) = '])
%disp(Inertia_cuboid)
%disp(['Inertia of Solar Panel about Centroid (kg-m^2) = '])
%disp(Inertia_sp_centroid)

%##Use the parallel axis theorem and move the sp to the centroid of the cuboid
Inertia = Inertia_cuboid;
for i = 1:number_of_panels
	lx = xyz_sp(i,1);
	ly = xyz_sp(i,2);
	lz = xyz_sp(i,3);
	rsp = [[0,-lz,ly];[lz,0,-lx];[-ly,lx,0]];
	rspT = rsp';
	Inertia_sp = Inertia_sp_centroid + mass_solar_panel*rsp*rspT;
	%disp(['Inertia of Solar Panel about Cuboid Center (kg-m^2) = '])
	%disp(Inertia_sp)
	Inertia = Inertia + Inertia_sp;
end

disp('Total Inertia of Satellite (kg-m^2)')
disp(Inertia)