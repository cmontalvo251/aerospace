%%%%Compute 3D Correlation
Xcenfile = 'Output_Files/Xcen.dat';
Ucenfile = 'Output_Files/Ualf.dat';
Vcenfile = 'Output_Files/Valf.dat';
Wcenfile = 'Output_Files/Walf.dat';
inputfile = 'Output_Files/Winds.OUT';
Uveezefile = 'Output_Files/Uveeze.dat';
Vveezefile = 'Output_Files/Vveeze.dat';
Wveezefile = 'Output_Files/Wveeze.dat';

read_some_data;

if ~dataflag
    disp('Error in input files')
    return
end

Adata = dlmread(inputfile);

%Adata(1,:) = []; %%Throw out first row

[r,c] = size(Adata);

x = Adata(:,1);
y = Adata(:,2);
z = Adata(:,3);
%%%Extract Altitudes
alts = zeros(1,NUM_AIRCRAFT);
zsort = sort(z);
zsortdiff = zsort(2:end)-zsort(1:end-1);
splits = [1;find(zsortdiff>20);length(zsort)];
for idx = 1:length(splits)-1
    s = splits(idx);
    e = splits(idx+1);
    alts(idx) = mean(zsort(s:e));
end
%%%%%%%%%%%%
t = Adata(:,4);
u = Adata(:,5);
v = Adata(:,6);
w = Adata(:,7);

dataloc = 'WRF_Wind_Data/';
UVWfrontend
offset = 100;
minx = -BOUNDARY + offset; %%Add a safety factor
maxx = BOUNDARY - offset;
miny = -BOUNDARY + offset;
maxy = BOUNDARY - offset;
minz = min(alts);
maxz = max(alts);
xnew = linspace(minx,maxx,100);
ynew = linspace(miny,maxy,100);
znew = linspace(minz,maxz,100);
uappx = 0.*xx;
vappx = uappx;
wappx = vappx;
umesh = zeros(length(xnew),length(ynew),length(znew));
vmesh = zeros(length(xnew),length(ynew),length(znew));
wmesh = zeros(length(xnew),length(ynew),length(znew));

    
for ii = 1:length(xnew)
    xi = xnew(ii)
    for jj = 1:length(ynew)
        yi = ynew(jj);
        for kk = 1:length(znew)
            zi = znew(kk);
            %disp([num2str(ii),' out of ',num2str(length(x))])
            state = [xi yi zi];
            ux = RBFhack2(state,xknot,alfu_uno,IKIND);
            ux = ux + Vbuild(state,uveeze,K35);
            vx = RBFhack2(state,xknot,alfv_dos,IKIND);
            vx = vx + Vbuild(state,vveeze,K35);
            wx = RBFhack2(state,xknot,alfw_tres,IKIND);
            wx = wx + Vbuild(state,wveeze,K35);
            uvwwrf = uvwout(xi,yi,-zi,0,dataloc,1);
            uappx(ii,jj,kk) = ux;
            vappx(ii,jj,kk) = vx;
            wappx(ii,jj,kk) = wx;
            umesh(ii,jj,kk) = uvwwrf(1);
            vmesh(ii,jj,kk) = uvwwrf(2);
            wmesh(ii,jj,kk) = uvwwrf(3);
        end
    end
end
Sr = sum(sum(sum((umesh-uappx).^2)));
St = sum(sum(sum((umesh-mean(mean(mean(umesh)))).^2)));
r2u = ((St-Sr)/St)

Sr = sum(sum(sum((vmesh-vappx).^2)));
St = sum(sum(sum((vmesh-mean(mean(mean(vmesh)))).^2)));
r2v = ((St-Sr)/St)

Sr = sum(sum(sum((wmesh-wappx).^2)));
St = sum(sum(sum((wmesh-mean(mean(mean(wmesh)))).^2)));
r2w = ((St-Sr)/St)