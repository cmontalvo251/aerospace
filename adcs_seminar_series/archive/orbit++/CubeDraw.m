function CubeDraw(dx,dy,dz,xc,yc,zc,phi,theta,psi,color)
%function CubeDraw(dx,dy,dz,xc,yc,zc,phi,theta,psi,color,ax)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xb = [ 1  1 -1 -1  1  1 -1 -1].*dx/2;
yb = [-1  1  1 -1 -1  1  1 -1].*dy/2;
zb = [-1 -1 -1 -1  1  1  1  1].*dz/2;
Xb = xb;
Yb = yb;
Zb = zb;

T = R123(phi,theta,psi);

[n,m] = size(xb);

rc = [xc;yc;zc];

for ii = 1:n;
  for jj = 1:m
    xs = xb(ii,jj);
    ys = yb(ii,jj);
    zs = zb(ii,jj);

    xyz = [xs,ys,zs]';

    xyz_trans = rc + T*xyz;

    Xb(ii,jj) = xyz_trans(1);
    Yb(ii,jj) = xyz_trans(2);
    Zb(ii,jj) = xyz_trans(3);
  end
end

for k = 1:1
    Face1 = patch(Xb(k,1:4),Yb(k,1:4),Zb(k,1:4),Zb(k,1:4),'facecolor',color,'facelighting','gouraud');
    Face2 = patch([Xb(k,1:2) Xb(k,6) Xb(k,5)],[Yb(k,1:2) Yb(k,6) Yb(k,5)],[Zb(k,1:2) Zb(k,6) Zb(k,5)],Zb(k,1:4),'facecolor',color,'facelighting','gouraud');
    Face3 = patch([Xb(k,2:3) Xb(k,7) Xb(k,6)],[Yb(k,2:3) Yb(k,7) Yb(k,6)],[Zb(k,2:3) Zb(k,7) Zb(k,6)],Zb(k,1:4),'facecolor',color,'facelighting','gouraud');
    Face4 = patch([Xb(k,3:4) Xb(k,8) Xb(k,7)],[Yb(k,3:4) Yb(k,8) Yb(k,7)],[Zb(k,3:4) Zb(k,8) Zb(k,7)],Zb(k,1:4),'facecolor',color,'facelighting','gouraud');
    Face5 = patch([Xb(k,1) Xb(k,4) Xb(k,8) Xb(k,5)],[Yb(k,1) Yb(k,4) Yb(k,8) Yb(k,5)],[Zb(k,1) Zb(k,4) Zb(k,8) Zb(k,5)],Zb(k,1:4),'facecolor',color,'facelighting','gouraud');
    Face6 = patch(Xb(k,5:8),Yb(k,5:8),Zb(k,5:8),Zb(k,5:8),'facecolor',color,'facelighting','gouraud');
end







% Copyright - Carlos Montalvo 2015
% You may freely distribute this file but please keep my name in here
% as the original owner
