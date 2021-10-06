function [Ixx,Iyy,Izz,max_moment_arm] = inertia(L,W,D,m)

max_moment_arm = max([L,W,D]/200);
disp(['Maximum Moment Arm (m) = ',num2str(max_moment_arm)])
area1 = L*W/(100^2);
area2 = W*D/(100^2);
area3 = L*D/(100^2);
max_area = max([area1,area2,area3]);
disp(['Maximum Area (m^2) = ',num2str(max_area)])

Ixx = (m/12)*(L^2+W^2)/(100^2);
Iyy = (m/12)*(L^2+D^2)/(100^2);
Izz = (m/12)*(W^2+D^2)/(100^2);

disp(['Inertia (kg-m^2) = ',num2str(Ixx),' ',num2str(Iyy),' ',num2str(Izz)])