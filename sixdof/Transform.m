function out = Transform(in,rB)
[r,c] = size(in);
out = in;
%%%First translate xyz by rB
xyzI = in(:,1:3);
ptp = in(:,4:6);
for ii = 1:r
    TIB = R123(ptp(ii,1),ptp(ii,2),ptp(ii,3));
    out(ii,1:3) = (-TIB*rB)' + xyzI(ii,:);
end
%%%Ptp are identical
out(:,4:6) = ptp;

%%%Calculate UVW
uvwB = in(:,7:9);
pqrB = in(:,10:12);
for ii = 1:r
    wcross = [0 -pqrB(ii,3) pqrB(ii,2);pqrB(ii,3) 0 -pqrB(ii,1);-pqrB(ii,2) pqrB(ii,1) 0];
    out(ii,7:9) = uvwB(ii,:) + (wcross*rB)';
end
   
%%%PQR are identical
out(:,10:12) = pqrB;


% Copyright - Carlos Montalvo 2015
% You may freely distribute this file but please keep my name in here
% as the original owner
