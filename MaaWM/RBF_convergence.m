addpath 'supplemental_routines/'

purge

global NSTATES

dataloc = getdataloc('Input_Files/MaaWM.ifiles');
UVWfrontend

[DYNSKIP,BOUNDARY,IKIND,K35] = read_boundary('Input_Files/MaaWM.SIM');
GRIDTYPE = 'WRF'; %%%THIS IS GRIDTYPE OF WIND MODEL. SET TO WRF to simply have the routine grab the .ATM file dataloc
%system('make'); %%This won't work on Windows. You need to simulate ahead of time
root_directory = 'Output_Files/';
inputfile = 'Output_Files/Winds.OUT';
ncens_vec =[];
data_vec = [];
rs_vec = [];

for ncens = [5:10]%[5:5:60 64] %%Do I have to do 64?

    simfile = slurp('Input_Files/MaaWM.SIM');
    simfile{14} = [num2str(ncens),' !Number of centers'];
    writedata('Input_Files/MaaWM.SIM',simfile);

    rs = [-10,-10,-10];

    data_pts = ncens; %%change to ncen next time
    while sum(rs < 0) > 0 && data_pts < 10*ncens
        data_pts = data_pts + 1
        freqx = grid_space(data_pts,BOUNDARY,GRIDTYPE,K35,IKIND);

        system('./Run.exe Input_Files/MaaWM.ifiles'); %%%This won't work on windows either

        [dataflag,xknot,alfu_uno,alfv_dos,alfw_tres,uveeze,vveeze,wveeze] = read_some_data(root_directory,K35);

        Adata = dlmread(inputfile);
        [r,c] = size(Adata);
    
        x = Adata(:,1);
        y = Adata(:,2);
        z = Adata(:,3);

        offset = 0;
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

        vappx = uappx;
        wappx = vappx;

        arc_length = 0*x;

        for ii = 1:length(xnew)
            for jj = 1:length(ynew)
                xi = xx(ii,jj);
                yi = yy(ii,jj);
                zi = z(1);
                state = [xi yi zi];
                ux = RBFhack2(state,xknot,alfu_uno,IKIND);
                ux = ux + Vbuild(state,uveeze,K35);
                vx = RBFhack2(state,xknot,alfv_dos,IKIND);
                vx = vx + Vbuild(state,vveeze,K35);
                wx = RBFhack2(state,xknot,alfw_tres,IKIND);
                wx = wx + Vbuild(state,wveeze,K35);
            uvwwrf = 3.28*uvwout(xi/3.28,yi/3.28,-zi/3.28,0,dataloc,1); %%Convert input position in feet to meters, Then conver output velocity from m/s to ft/s
            uappx(ii,jj) = ux;
            vappx(ii,jj) = vx;
            wappx(ii,jj) = wx;
            umesh(ii,jj) = uvwwrf(1);
            vmesh(ii,jj) = uvwwrf(2);
            wmesh(ii,jj) = uvwwrf(3);
            end
        end

        %myplot('U Error',18,'X (ft)','Y (ft)','UError (ft/s)',[-27 30]);
        %mesh(xx,yy,umesh-uappx)
        Sr = sum(sum((umesh-uappx).^2));
        St = sum(sum((umesh-mean(mean(umesh))).^2));
        uappx_lin = interp2(xnew,ynew,uappx,x,y);
        for udx = 1:length(uappx_lin)
            if isnan(uappx_lin(udx))
                uappx_lin(udx) = 0;
                u(udx) = 0;
            end
        end
        r2u = ((St-Sr)/St);

        %myplot('V Error',18,'X (ft)','Y (ft)','VError (%)',[-27 30]);
        %mesh(xx,yy,vmesh-vappx)
        Sr = sum(sum((vmesh-vappx).^2));
        St = sum(sum((vmesh-mean(mean(vmesh))).^2));
        vappx_lin = interp2(xnew,ynew,vappx,x,y);
        for vdx = 1:length(vappx_lin)
            if isnan(vappx_lin(vdx))
                vappx_lin(vdx) = 0;
                v(vdx) = 0;
            end
        end
        r2v = ((St-Sr)/St);

        %myplot('W Error',18,'X (ft)','Y (ft)','WError (ft/s)',[-27 30]);
        %mesh(xx,yy,wmesh-wappx)
        Sr = sum(sum((wmesh-wappx).^2));
        St = sum(sum((wmesh-mean(mean(wmesh))).^2));
        wappx_lin = interp2(xnew,ynew,wappx,x,y);
        for wdx = 1:length(wappx_lin)
            if isnan(wappx_lin(wdx))
                wappx_lin(wdx) = 0;
                w(wdx) = 0;
            end
        end

        r2w = ((St-Sr)/St);

        rs = [r2u,r2v,r2w]
    end
    ncens_vec = [ncens_vec;ncens];
    data_vec = [data_vec;data_pts];
    rs_vec = [rs_vec;rs];
end

ncens_vec
data_vec
rs_vec

save wave1_results.mat
