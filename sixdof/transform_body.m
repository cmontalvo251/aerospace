function [Xn,Yn,Zn] = transform_body(X,Y,Z,pos,att,SF)
SFMAT = [SF(1) 0 0;0 SF(2) 0;0 0 SF(3)];

x = pos(1);
y = pos(2);
z = pos(3);

ph = att(1);
th = att(2);
ps = att(3);

T_IB = R123(ph,th,ps);
   
R = [x,y,z]';
            
[n,m] = size(X);

Xn = zeros(n,m);
Yn = zeros(n,m);
Zn = zeros(n,m);

for i = 1:n;
    for j = 1:m
        xs = X(i,j);
        ys = Y(i,j);
        zs = Z(i,j);
        
        xyz = [xs,ys,zs]';
        
        xyz_trans = SFMAT*T_IB*xyz + R;
        
        Xn(i,j) = xyz_trans(1);
        Yn(i,j) = xyz_trans(2);
        Zn(i,j) = xyz_trans(3);
    end
end


% Copyright - Carlos Montalvo 2015
% You may freely distribute this file but please keep my name in here
% as the original owner
