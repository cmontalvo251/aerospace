minx = -BOUNDARY + offset; %%Add a safety factor
maxx = BOUNDARY - offset;
miny = -BOUNDARY + offset;
maxy = BOUNDARY - offset;
minz = min(alts);
maxz = max(alts);
xnew = linspace(minx,maxx,100);
ynew = linspace(miny,maxy,100);
znew = linspace(minz,maxz,100);
[xx,zz] = meshgrid(xnew,znew);
uappx = 0.*xx;
vappx = uappx;
wappx = vappx;
umesh = zeros(length(xnew),length(ynew));
vmesh = zeros(length(xnew),length(ynew));
wmesh = zeros(length(xnew),length(ynew));

arc_length = 0*x;

for ii = 1:length(xnew)
    for jj = 1:length(znew)
        xi = xx(ii,jj);
        yi = 0; %%%Change this later!
        zi = zz(ii,jj);

        state = [xi yi zi];
        ux = RBFhack2(state,xknot,alfu_uno,IKIND);
        ux = ux + Vbuild(state,uveeze,K35);
        vx = RBFhack2(state,xknot,alfv_dos,IKIND);
        vx = vx + Vbuild(state,vveeze,K35);
        wx = RBFhack2(state,xknot,alfw_tres,IKIND);
        wx = wx + Vbuild(state,wveeze,K35);
        uvwwrf = uvwout(xi,yi,-zi,0,dataloc,1);
        uappx(ii,jj) = ux;
        vappx(ii,jj) = vx;
        wappx(ii,jj) = wx;
        umesh(ii,jj) = uvwwrf(1);
        vmesh(ii,jj) = uvwwrf(2);
        wmesh(ii,jj) = uvwwrf(3);
    end
end

myplot('U Fit',18,'X (m)','U (m/s)','Z (m)',[-17 20]);
mesh(xx,uappx,zz,uappx)
%surf(xx,umesh,zz,umesh)
reverse(['y','z'])
zlim([-500 -100])

myplot('V Fit',18,'X (m)','V (m/s)','Z (m)',[-17 20]);
mesh(xx,vappx,zz,vappx)
%surf(xx,vmesh,zz,vmesh)
reverse(['y','z'])
zlim([-500 -100])

myplot('W Fit',18,'X (m)','W (m/s)','Z (m)',[-17 20]);
mesh(xx,wappx,zz,wappx)
%surf(xx,wmesh,zz,wmesh)
reverse(['y','z'])
zlim([-500 -100])

Sr = sum(sum((umesh-uappx).^2));
St = sum(sum((umesh-mean(mean(umesh))).^2));
r2u = ((St-Sr)/St);
Sr = sum(sum((vmesh-vappx).^2));
St = sum(sum((vmesh-mean(mean(vmesh))).^2));
r2v = ((St-Sr)/St);
Sr = sum(sum((wmesh-wappx).^2));
St = sum(sum((wmesh-mean(mean(wmesh))).^2));
r2w = ((St-Sr)/St);
rsVertical = [r2u,r2v,r2w]




