%purge

%root_directory = 'Input_Files/Example_Fit/';
root_directory = 'Output_Files/';
K35 = 0;
IKIND = 1;
[dataflag,xknot,alfu_uno,alfv_dos,alfw_tres,uveeze,vveeze,wveeze] = read_some_data(root_directory,K35);

datax = 41;
datay = datax;
boundary = 500;

X = linspace(-boundary,boundary,datax);
Y = linspace(-boundary,boundary,datay);
Z = -200;

%plottool(1,'Grid Space',18,'X (m)','Y (m)');
for idx = 1:length(X)
    for jdx = 1:length(Y)
      state = [X(idx) Y(jdx) Z];
      ux = RBFhack2(state,xknot,alfu_uno,IKIND);
      ux = ux + Vbuild(state,uveeze,K35);
      vx = RBFhack2(state,xknot,alfv_dos,IKIND);
      vx = vx + Vbuild(state,vveeze,K35);
      wx = RBFhack2(state,xknot,alfw_tres,IKIND);
      wx = wx + Vbuild(state,wveeze,K35);
      uvwwrf = [ux,vx,wx];
      u = uvwwrf(1);
      v = uvwwrf(2);
      w = uvwwrf(3);
      plot3(X(idx),Y(jdx),u,'r*')
    end
end

view(-27,50)