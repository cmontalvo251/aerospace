clear
clc
close all

xyz1 = [0;0;0];
xyz2 = [1;1;1];

r = xyz1-xyz2;

x = r(1);
y = r(2);
z = r(3);

psi = atan2(y,x);
theta = atan2(-z,sqrt(x^2+y^2));
phi = 0;

ptp = [phi;theta;psi];

TIB = R123(ptp(1),ptp(2),ptp(3));

xhatB = TIB(:,1);
yhatB = TIB(:,2);
zhatB = TIB(:,3);

%%%Draw a vector between satellite 1 and 2
plot3([xyz1(1) xyz2(1)],[xyz1(2) xyz2(2)],[xyz1(3) xyz2(3)],'g-','LineWidth',2)
hold on

%Artificially Plot an inertial frame
plot3([0 1],[0 0],[0 0],'b-','LineWidth',2)
plot3([0 0],[0 1],[0 0],'b-','LineWidth',2)
plot3([0 0],[0 0],[0 1],'b-','LineWidth',2)
view(45,45)
grid on


plot3([xyz2(1) xyz2(1)+TIB(1,1)],[xyz2(2) xyz2(2)+TIB(2,1)],[xyz2(3) xyz2(3)+TIB(3,1)],'r-','LineWidth',3)
plot3([xyz2(1) xyz2(1)+TIB(1,2)],[xyz2(2) xyz2(2)+TIB(2,2)],[xyz2(3) xyz2(3)+TIB(3,2)],'r-','LineWidth',3)
plot3([xyz2(1) xyz2(1)+TIB(1,3)],[xyz2(2) xyz2(2)+TIB(2,3)],[xyz2(3) xyz2(3)+TIB(3,3)],'r-','LineWidth',3)

axis equal

