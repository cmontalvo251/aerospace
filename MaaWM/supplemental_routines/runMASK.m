close all
%%%%Compute Average Error as a function of number of data points
%%%%%Create a mask
areamask = 2500; %%m^2
dx = xnew(2)-xnew(1);
dy = ynew(2)-ynew(1);
Lmask = sqrt(areamask);
rx = ceil(Lmask/dx);
ry = ceil(Lmask/dy);
if rx == 1
    rx = 2;
end
if ry == 1
    ry = 2;
end
mask = ones(rx,ry);
S = sum(sum(mask));
xmask = rx*dx;
ymask = ry*dy;
%%%%Now move through error function and extract average error
errorV = abs(vmesh-vappx);
meanV = mean(mean(vmesh));
%%%Compute Local Correlation as well
plottool(1,'Actual Error',18,'X (m)','Y (m)','Error (m/s)');
mesh(xx,yy,errorV)
r = 1:(rx-1):(length(xnew)-rx+1);
c = 1:(ry-1):(length(ynew)-ry+1);
avg_error = zeros(length(r),length(c));
rU = 0*avg_error;
rV = 0*avg_error;
rW = 0*avg_error;
num_data_pts = avg_error;
xctr = 0;
yctr = 0;
plottool(1,'Data Points',18,'x (m)','y (m)');
hold on
for idx = 1:(rx-1):(length(xnew)-rx+1)
    xctr = xctr + 1;
    yctr = 0;
    for jdx = 1:(ry-1):(length(ynew)-ry+1)
        yctr = yctr + 1;
        sx = idx:(idx+rx-1);
        sy = jdx:(jdx+ry-1);
        avg_error(xctr,yctr) = sum(sum(mask.*errorV(sx,sy)))/S;
        %%%Compute Correlation
        uxy = umesh(sx,sy);
        uaxy = uappx(sx,sy);
        Sr = sum(sum((uxy-uaxy).^2));
        St = sum(sum((uxy-mean(mean(uxy))).^2));
        rU(xctr,yctr) = (St-Sr)/St;
        
        vxy = vmesh(sx,sy);
        vaxy = vappx(sx,sy);
        Sr = sum(sum((vxy-vaxy).^2));
        St = sum(sum((vxy-mean(mean(vxy))).^2));
        rV(xctr,yctr) = (St-Sr)/St;
        
        wxy = wmesh(sx,sy);
        waxy = wappx(sx,sy);
        Sr = sum(sum((wxy-waxy).^2));
        St = sum(sum((wxy-mean(mean(wxy))).^2));
        rW(xctr,yctr) = (St-Sr)/St;
        %%%We also need to compute the number of data points in this
        %%%mask
        xspace = xx(sx,sy);
        yspace = yy(sx,sy);
        minx = min(min(xspace));
        maxx = max(max(xspace));
        miny = min(min(yspace));
        maxy = max(max(yspace));
        rectangle('Position',[minx miny maxx-minx maxy-miny])
        for ldx = 1:length(x)
            if minx <= x(ldx) && x(ldx) <= maxx && miny <= y(ldx) && y(ldx) <= maxy
                num_data_pts(xctr,yctr) = num_data_pts(xctr,yctr) + 1;
            end
        end
    end
end
plot(x,y,'b*')
plottool(1,'Number of Data Points',18,'x(m)','y(m)','Data Points/2500 m^2');
xspace = linspace(min(x),max(x),length(r));
yspace = linspace(min(y),max(y),length(c));
[xxspace,yyspace] = meshgrid(xspace,yspace);
mesh(xxspace,yyspace,num_data_pts)
plottool(1,'Average Error',18,'x (m)','y (m)','Average Error (m/s)');
mesh(xxspace,yyspace,avg_error)

%%%%Now Convert the two matrices to vectors
num_data_vec = unwrapmatrix(num_data_pts)';
avg_error_vec = unwrapmatrix(avg_error)';
rU_vec = unwrapmatrix(rU)';
rV_vec = unwrapmatrix(rV)';
rW_vec = unwrapmatrix(rW)';
%data_pt_density = 1000*num_data_vec/areamask;
data_pt_density = num_data_vec/areamask;

%%%Convert to Distance between data points
X = data_pt_density;
L = X==0;
X(L) = [];
Y = avg_error_vec;
Y(L) = [];
%X = 1000./sqrt(X*1000);
XD = 50./sqrt(X*2500);

%%%%Plot correlation
plottool(1,'r^2 versus Ndata points',18,'Distance Between Data Points (m)','r^2');
rU_vec(L) = [];
rU_vec(rU_vec < 0) = 0;
rV_vec(L) = [];
rV_vec(rV_vec < 0) = 0;
rW_vec(L) = [];
rW_vec(rW_vec < 0) = 0;
plot(XD,rU_vec,'b*')
%plot(X,rV_vec,'r*')
%plot(X,rW_vec,'g*')

%%%%Create a Fit
%%% Y = a0 + a1*(X-x0)^2
[XS,I] = sort(XD);
a0 = 1;
x0 = min(XD);
H = [(XD-x0).^2];
AU = inv(H'*H)*H'*(rU_vec-a0);
YU = H*AU+a0;
YU = YU(I);
plot(XS,YU,'b-','LineWidth',2)
AV = inv(H'*H)*H'*(rV_vec-a0);
YV = H*AV+a0;
YV = YV(I);
%plot(XS,YV,'r-','LineWidth',2)
AW = inv(H'*H)*H'*(rW_vec-a0);
YW = H*AW+a0;
YW = YW(I);
%plot(XS,YW,'g-','LineWidth',2)
ylim([0 1])

%%%Plot it
plottool(1,'Error versus Number of Data Points',18,'Data Point Density (N/m^2)','Average Error (m/s)');
%plot(data_pt_density,avg_error_vec,'b*','LineWidth',2)
plot(X,Y,'b*','LineWidth',2)

%%%%Compute a regression curve through the points
a0 = max(Y);
%Y = avg_error_vec-a0;
%H = [1./(data_pt_density+1).^2];
%astar = inv(H'*H)*H'*Y;
xsorted = linspace(min(data_pt_density),max(data_pt_density),100)';

%H = [1./(xsorted+1).^2];
%yfit = H*astar+a0;
N = 0;
Nstep = 1;
Sr = realmax;
iter = 1;
while iter < 10000
    N = N + Nstep;
    yfit = a0./(data_pt_density+1).^N;
    Sr_new = sum((yfit-avg_error_vec).^2);
    if Sr_new > Sr
        Nstep = -Nstep/2;
    else
        Sr = Sr_new;
    end
    iter = iter + 1;
end
yfit = a0./(xsorted+1).^N;
plot(xsorted,yfit,'r-','LineWidth',2)

plottool(1,'Error versus Distance Between Points',18,'Distance Between Data Points (m)','Average Error (m/s)');
plot(XD,Y,'b*','LineWidth',2)
%%%Y = H A
x0 = min(XD);
a0 = Y(find(XD==x0));
H = [(XD-x0).^2];
Astar = inv(H'*H)*H'*(Y-a0);
Ytilde = H*Astar+a0;
[X,I] = sort(XD);
Ytilde = Ytilde(I);
plot(X,Ytilde,'b-','LineWidth',2)
%ylim([0 10])


plottool(1,'Percent Error versus Distance Between Points',18,'Distance Between Data Points (m)','Average Percent Error (%)');
plot(XD,(Y./meanV)*100,'b*','LineWidth',2)
plot(X,(Ytilde/meanV)*100,'b-','LineWidth',2)


return