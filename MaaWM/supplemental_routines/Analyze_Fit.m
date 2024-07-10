%function Analyze_Fit(IMESH,IKIND)

root_directory = 'Output_Files/';
%root_directory = 'Input_Files/Example_Fit/';
[dataflag,xknot,alfu_uno,alfv_dos,alfw_tres,uveeze,vveeze,wveeze] = read_some_data(root_directory,K35);

if ~dataflag
    disp('Error in input files')
    return
end

inputfile = 'Output_Files/Winds.OUT';
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

if length(alts) > 1
    close all
end

if IMESH
  switch GRIDTYPE
   case 'WRF'
    dataloc = getdataloc('Input_Files/MaaWM.ifiles');
    UVWfrontend
   case 'SINE'
    freqy = freqx;
  end
end

for Nac = 1:length(alts)
    disp(['Current Altitude = ',num2str(alts(Nac))])

    if IMESH %%%%%%%%%THIS WILL BREAK HERE
        offset = OFFSET;
        minx = -BOUNDARY + offset; %%Add a safety factor
        maxx = BOUNDARY - offset;
        miny = -BOUNDARY + offset;
        maxy = BOUNDARY - offset;
        xnew = linspace(minx,maxx,100);
        ynew = linspace(miny,maxy,100);
        [xx,yy] = meshgrid(xnew,ynew);
        uappx = 0.*xx;
        umesh = zeros(length(xnew),length(ynew));
        vmesh = zeros(length(xnew),length(ynew));
        wmesh = zeros(length(xnew),length(ynew));
    else
        xnew = x;
        ynew = 0;
        uappx = 0*x;
    end
    vappx = uappx;
    wappx = vappx;
    
    arc_length = 0*x;
    
    for ii = 1:length(xnew)
        for jj = 1:length(ynew)
            %disp([num2str(ii),' out of ',num2str(length(x))])
            if IMESH
                xi = xx(ii,jj);
                yi = yy(ii,jj);
                zi = alts(Nac);
            else
                xi = x(ii);
                yi = y(ii);
                zi = z(1);
            end
            state = [xi yi zi];
            ux = RBFhack2(state,xknot,alfu_uno,IKIND);
            ux = ux + Vbuild(state,uveeze,K35);
            vx = RBFhack2(state,xknot,alfv_dos,IKIND);
            vx = vx + Vbuild(state,vveeze,K35);
            wx = RBFhack2(state,xknot,alfw_tres,IKIND);
            wx = wx + Vbuild(state,wveeze,K35);
            if IMESH
                uvwwrf = 3.28084*uvwout(xi*0.3048,yi*0.3048,-zi*0.3048,0,dataloc,1); %%Convert input position in feet to meters, Then conver output velocity from m/s to ft/s
                uappx(ii,jj) = ux;
                vappx(ii,jj) = vx;
                wappx(ii,jj) = wx;
                umesh(ii,jj) = uvwwrf(1);
                vmesh(ii,jj) = uvwwrf(2);
                wmesh(ii,jj) = uvwwrf(3);
            else
                uappx(ii) = ux;
                vappx(ii) = vx;
                wappx(ii) = wx;
            end
            % try
            %   dx = x(ii)-x(ii-1);
            %   dy = y(ii)-y(ii-1);
            %   dz = z(ii)-z(ii-1);
            %   arc_length(ii) = arc_length(ii-1) + sqrt(dx^2+dy^2+dz^2);
            % end
        end
    end
    
    if IMASK
        runMASK
    end
    
    % format long g
    % arc_length(end)
    % t = arc_length(end)/20
    
    %%%Plot everything
    myplot('U Fit',18,'X (ft)','Y (ft)','U (ft/s)',[-27 30]);
    if IMESH
        %mesh(xx,yy,uappx,'EdgeColor','r','LineWidth',2,'FaceColor','k')
        %mesh(xx,yy,umesh,'EdgeColor','b','LineWidth',2,'FaceColor','k')
        mesh(xx,yy,uappx)
        surf(xx,yy,umesh)
        %legend('RBF Fit','Sampled Data')
        %plot3(x,y,z-mean(z),'b-','LineWidth',2)
    else
        plot3(x,y,u,'b*')
        plot3(x,y,uappx,'r*')
    end
    
    myplot('V Fit',18,'X (ft)','Y (ft)','V (ft/s)',[-27 30]);
    if IMESH
        surf(xx,yy,vmesh)
        mesh(xx,yy,vappx)
    else
        plot3(x,y,v,'b*')
        plot3(x,y,vappx,'r*')
        legend('Sampled Data','RBF Fit')
    end
    
    myplot('W Fit',18,'X (ft)','Y (ft)','W (ft/s)',[-43 50]);
    if IMESH
        surf(xx,yy,wmesh)
        mesh(xx,yy,wappx)
    else
        plot3(x,y,w,'b*','LineWidth',1)
        plot3(x,y,wappx,'r*')
        legend('Sampled Data','RBF Fit')
    end
    
    myplot('U Error',18,'X (ft)','Y (ft)','UError (ft/s)',[-27 30]);
    if IMESH
        mesh(xx,yy,umesh-uappx)
        Sr = sum(sum((umesh-uappx).^2));
        St = sum(sum((umesh-mean(mean(umesh))).^2));
        uappx_lin = interp2(xnew,ynew,uappx,x,y);
        for udx = 1:length(uappx_lin)
            if isnan(uappx_lin(udx))
                uappx_lin(udx) = 0;
                u(udx) = 0;
            end
        end
    else
        plot3(x,y,u-uappx,'b*')
        %%%Compute r^2
        Sr = sum((u-uappx).^2);
        St = sum((u-mean(u)).^2);
    end
    r2u = ((St-Sr)/St);
    
    myplot('V Error',18,'X (ft)','Y (ft)','VError (%)',[-27 30]);
    if IMESH
        mesh(xx,yy,vmesh-vappx)
        Sr = sum(sum((vmesh-vappx).^2));
        St = sum(sum((vmesh-mean(mean(vmesh))).^2));
        vappx_lin = interp2(xnew,ynew,vappx,x,y);
        for vdx = 1:length(vappx_lin)
            if isnan(vappx_lin(vdx))
                vappx_lin(vdx) = 0;
                v(vdx) = 0;
            end
        end
    else
        errorV = 100*abs(v-vappx)./v;
        l = find(abs(v) < 1e-1);
        errorV(l) = 0;
        plot3(x,y,errorV,'b*')
        %%%Compute r^2
        Sr = sum((v-vappx).^2);
        St = sum((v-mean(v)).^2);
    end
    r2v = ((St-Sr)/St);
    
    myplot('W Error',18,'X (ft)','Y (ft)','WError (ft/s)',[-27 30]);
    % errorW = 100*abs(w-wappx)./w;
    % l = find(abs(w) < 1e-1);
    % errorW(l) = 0;
    if IMESH
        mesh(xx,yy,wmesh-wappx)
        Sr = sum(sum((wmesh-wappx).^2));
        St = sum(sum((wmesh-mean(mean(wmesh))).^2));
        wappx_lin = interp2(xnew,ynew,wappx,x,y);
        for wdx = 1:length(wappx_lin)
            if isnan(wappx_lin(wdx))
                wappx_lin(wdx) = 0;
                w(wdx) = 0;
            end
        end
    else
        errorW = abs(w-wappx);
        plot3(x,y,errorW,'b*')
        Sr = sum((w-wappx).^2);
        St = sum((w-mean(w)).^2);
    end
    r2w = ((St-Sr)/St);
    
    rs = [r2u,r2v,r2w]
    
    if IMESH
        NDATA = length(u);
        TEND = t(end);
        avgU = mean(mean(abs(umesh-uappx)));
        maxUE = max(max(abs(umesh-uappx)));
        avgV = mean(mean(abs(vmesh-vappx)));
        maxVE = max(max(abs(vmesh-vappx)));
        avgW = mean(mean(abs(wmesh-wappx)));
        maxWE = max(max(abs(wmesh-wappx)));
        fprintf(montefile,'%f %f %f %f %f %f %f %f %f %f %f %f \n',NUM_AIRCRAFT,r2u,r2v,r2w,avgU,avgV,avgW,maxUE,maxVE,maxWE,NDATA,TEND);
    end
end

if length(alts) > 1 && IMESH
    runVerticalSlice
end
    
