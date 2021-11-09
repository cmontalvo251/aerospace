model = stlread('sat_model_2.stl');

points = model.Points;
cList = model.ConnectivityList;

ax = 0.2;
ay = 0.1;
az = 2.0;

Rx = [1 0 0 ; 0 cos(ax) -sin(ax) ; 0 sin(ax) cos(ax)];
Ry = [cos(ay) 0 sin(ay) ; 0 1 0 ; -sin(ay) 0 cos(ay)];
Rz = [cos(az) -sin(az) 0 ; sin(az) cos(az) 0 ; 0 0 1];


pointR = points*Rx;
pointR = points*Ry;
pointR = points*Rz;


figure
axis equal
h1 = trimesh(model);
h1.FaceColor = 'flat';
rotate(h1,[1 0 0],25);
% xlim([-13.4 -10.4])
% ylim([-28.5 -25.5])
% zlim([-1.5 1.5])