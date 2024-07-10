%%%Check Measurements
purge
dataloc = getdataloc('Input_Files/MaaWM.ifiles');
UVWfrontend

root_directory = 'Output_Files/';
[DYNSKIP,BOUNDARY,IKIND,K35] = read_boundary('Input_Files/MaaWM.SIM');
[dataflag,xknot,alfu_uno,alfv_dos,alfw_tres,uveeze,vveeze,wveeze] = read_some_data(root_directory,K35);

inputfile = 'Output_Files/Winds.OUT';
Adata = dlmread(inputfile);
[r,c] = size(Adata);

x = Adata(:,1);
y = Adata(:,2);
z = Adata(:,3);
t = Adata(:,4);
vx = Adata(:,5);
vy = Adata(:,6);
vz = Adata(:,7);
plot3(x,y,z,'b*')

%break

figure()
plot3(x,y,vx,'r*') %%%Red is straight from Winds.OUT file
hold on
for ii = 1:length(x)
    xi = x(ii);
    yi = y(ii);
    zi = z(ii);
    uvwwrf = 3.28084*uvwout(xi*0.3048,yi*0.3048,-zi*0.3048,0,dataloc,1);
    state = [xi,yi,zi];
    vx = RBFhack2(state,xknot,alfv_dos,IKIND);
    vx = vx + Vbuild(state,vveeze,K35);
    
    ux = RBFhack2(state,xknot,alfu_uno,IKIND);
    ux = ux + Vbuild(state,uveeze,K35);

    %plot3(xi,yi,vx,'g*')
    plot3(xi,yi,uvwwrf(1),'b*') %%%Blue is straight from MATLAB UVWRF Model
end



