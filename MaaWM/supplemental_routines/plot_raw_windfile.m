%%%Plot Raw Data
purge

myplot('Raw Data',18,'X (m)','Y (m)','U (m/s)',[-27 30]);

inputfiles = {'Output_Files/Winds.OUT';'Input_Files/Example_Fit/Winds.OUT'};

colors = {'b*','r*'};

for idx = 1:length(inputfiles)
    
    inputfile = inputfiles{idx};
    
    Adata = dlmread(inputfile);

    %Adata(1,:) = []; %%Throw out first row

    [r,c] = size(Adata);

    x = Adata(:,1);
    y = Adata(:,2);
    z = Adata(:,3);
    t = Adata(:,4);
    u = Adata(:,5);
    v = Adata(:,6);
    w = Adata(:,7);


    plot3(x,y,u,colors{idx})

end

%{'Output_Files/Winds.OUT';'Input_Files/Example_Fit/Winds.OUT'};
legend('Output Files','Example Fit')