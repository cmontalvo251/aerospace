purge

all = importdata('Output_Files/Multi_Aircraft/MonteCarlo.txt');

data = all.data;

[r,c] = size(data);

datax = 1:20;

rU = zeros(20,1);
rV = zeros(20,1);
rW = zeros(20,1);
data_pts = rW;
time = rW;

%%%Average all aircraft
ctr = 1;
for idx = 1:20
    s = (idx-1)*3 + 1;
    e = s + 2;
    rU3 = data(s:e,2);
    rU3(rU3<0) = [];
    rU(ctr) = mean(rU3);
    
    rV3 = data(s:e,3);
    rV3(rV3<0) = [];
    rV(ctr) = mean(rV3);
    
    rW3 = data(s:e,4);
    rW3(rW3<0) = [];
    rW(ctr) = mean(rW3);
    
    data3 = data(s:e,11);
    data_pts(ctr) = mean(data3);
    
    time3 = data(s:e,12);
    time(ctr) = mean(time3);
    
    ctr = ctr + 1;
end

plottool(1,'Error',18,'Number of Aircraft','Mean of Correlation Value');
%plot(datax,rU,'b-','LineWidth',2)
%plot(datax,rV,'r-','LineWidth',2)
%plot(datax,rW,'g-','LineWidth',2)
meanR = (1/3)*(rU+rV+rW);
plot(datax,meanR,'b-','LineWidth',2)
ylim([0 1])

plottool(1,'Data Points',18,'Number of Aircraft','Average Data Points');
plot(datax,data_pts,'b-','LineWidth',2)

plottool(1,'Time',18,'Number of Aircraft','Time to Convergence (min)')
plot(datax,time/60,'b-','LineWidth',2)