purge

%all = importdata('Output_Files/Data_Point_Studies/Data_Point_MonteCarlo.txt');
all = importdata('MonteCarlo.txt');

data = all.data;

datax = data(:,1);
myxlabel = 'Number of Grid Points Along X';

datax = data(:,1).^2;
myxlabel = 'Total Number of Grid Points';

datax = 1000./data(:,1);
myxlabel = 'Distance Between Grid Points (m)';

rU = data(:,2);
rV = data(:,3);
rW = data(:,4);
avgU = data(:,5);
avgV = data(:,6);
avgW = data(:,7);
maxU = data(:,8);
maxV = data(:,9);
maxW = data(:,10);
ndata = data(:,11);

%[1000./sqrt(ndata),rU,rV,rW,ndata]
[rU,rV,rW,ndata]

plottool(1,'Error',18,myxlabel,'Average Error (m/s)');
plot(datax,avgU,'b-','LineWidth',2)
plot(datax,avgV,'r-','LineWidth',2)
plot(datax,avgW,'g-','LineWidth',2)
normUVW = sqrt(avgU.^2+avgV.^2+avgW.^2);
%plot(datax,normUVW,'k-','LineWidth',3)
avgUVW = (1/3)*(avgU+avgV+avgW);
%plot(datax,avgUVW,'k--','LineWidth',3)
legend('U','V','W','Norm','Average')

plottool(1,'Error',18,myxlabel,'Max Error (m/s)');
plot(datax,maxU,'LineWidth',2)
plot(datax,maxV,'r-','LineWidth',2)
plot(datax,maxW,'g-','LineWidth',2)
maxnormUVW = sqrt(maxU.^2+maxV.^2+maxW.^2);
%plot(datax,maxnormUVW,'k-','LineWidth',3)
maxavgUVW = (1/3)*(maxU+maxV+maxW);
%plot(datax,maxavgUVW,'k--','LineWidth',3)
legend('U','V','W','Norm','Average')

plottool(1,'Error',18,myxlabel,'Correlation Function');
plot(datax,rU,'b-','LineWidth',2)
plot(datax,rV,'r-','LineWidth',2)
plot(datax,rW,'g-','LineWidth',2)
normR = sqrt(rU.^2+rV.^2+rW.^2);
%plot(datax,normUVW,'k-','LineWidth',3)
avgR = (1/3)*(rU+rV+rW);
%plot(datax,avgUVW,'k--','LineWidth',3)
legend('U','V','W','Norm','Average')
