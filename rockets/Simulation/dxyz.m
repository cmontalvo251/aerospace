function uvw = dxyz(xdot, ydot,zdot,T)

xdot=0;
ydot=0;
zdot=0;
T=

xyzdots=[xdot; ydot; zdot];
uvw=T\xyzdots; 
u=uvw(1,:);
v=uvw(2,:);
w=uvw(3,:);
uvw=[u; v; w];