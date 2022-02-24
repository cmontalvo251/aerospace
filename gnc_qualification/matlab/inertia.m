function [Inertia,max_moment_arm,max_area] = inertia(Lx,Wy,Dz,total_mass,mass_solar_panel,number_of_panels,LWD_sp,xyz_sp)
%%This code will compute the inertia of ta CubeSAT with Solar Panels
% L is the length of the satellite in cm
% W is the width of the satellite in cm
% D is the depth of the satellite in cm
% total_mass is the total mass of the satellite including solar panels
% mass_solar_panel is the mass of solar panels
% number_of_panels is the number of solar panels
% LWD_sp is the size of 1 panel in meters
% xyz_sp is the position of the solar panels from the center of mass of the cuboid in meters

max_moment_arm = max([Lx,Wy,Dz]/2);
disp(['Maximum Moment Arm (m) = ',num2str(max_moment_arm)])
area1 = Lx*Wy; %%%Top face
area2 = Wy*Dz; %%Ram face
area3 = Lx*Dz; %%%Side face
max_area = max([area1,area2,area3]);
disp(['Maximum Area (No Panels) (m^2) = ',num2str(max_area)])

%%%Mass/Inertia of Cuboid
mass_cuboid = total_mass - mass_solar_panel*number_of_panels;
Ixx = (mass_cuboid/12)*(Wy^2+Dz^2);
Iyy = (mass_cuboid/12)*(Lx^2+Dz^2);
Izz = (mass_cuboid/12)*(Lx^2+Wy^2);
Inertia_cuboid = [[Ixx,0,0];[0,Iyy,0];[0,0,Izz]];

%%%Inertia of Solar panel about centroid
Lx_sp = LWD_sp(1,1);
Wy_sp = LWD_sp(2,1);
Dz_sp = LWD_sp(3,1);
Inertia_sp_centroid = mass_solar_panel/12.0 * [[Wy_sp^2 + Dz_sp^2,0,0];[0,Lx_sp^2+Dz_sp^2,0];[0,0,Lx_sp^2+Wy_sp^2]];

%%%Compute Max Area of Solar Panels
area1sp = 0; %%Top face (basically nothing)
area2sp = (Wy_sp * Dz_sp)*number_of_panels;  %%%Ram face
area3sp = 0; %%side face (basically nothing)

%%%Compute Area of Faces with panels
areaT1 = area1 + area1sp;
areaT2 = area2 + area2sp;
areaT3 = area3 + area3sp;

max_areaT = max([areaT1,areaT2,areaT3]);
disp(['Maximum Area (with panels) (m^2) = ',num2str(max_areaT)])

%disp(['Inertia of Cuboid (kg-m^2) = '])
%disp(Inertia_cuboid)
%disp(['Inertia of Solar Panel about Centroid (kg-m^2) = '])
%disp(Inertia_sp_centroid)

%##Use the parallel axis theorem and move the sp to the centroid of the cuboid
Inertia = Inertia_cuboid;
for i = 1:number_of_panels
	lx = xyz_sp(1,i);
	ly = xyz_sp(2,i);
	lz = xyz_sp(3,i);
	rsp = [[0,-lz,ly];[lz,0,-lx];[-ly,lx,0]];
	rspT = rsp';
	Inertia_sp = Inertia_sp_centroid + mass_solar_panel*rsp*rspT;
	%disp(['Inertia of Solar Panel about Cuboid Center (kg-m^2) = '])
	%disp(Inertia_sp)
	Inertia = Inertia + Inertia_sp;
end

disp('Total Inertia of Satellite (kg-m^2)')
disp(Inertia)