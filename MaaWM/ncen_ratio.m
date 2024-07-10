purge

%%100 meter windfield

%%% 100 meter has more humps so it needs more centers to
%converge. whereas the WRF model is only 200
%meters so it can converge with only 5 centers.

%%%100 WF = Minimum number of data points to converge is 29 data points with
%7 centers. 

%%% 29 / 7

%%%WRF WF = Minimum number of data points to converge is 11 data
%points with 5 centers. 

%%% 11 / 5 -- Way less dense

ncens_vec_100 = [1 2 3 4  5 6 7 8 9  10   15   20   25   30   35   40   45   50   55   60   64];
data_vec_100 =  [10 20 30 40 50 60 29 27 23 23   18   23   28   33   38   43   48   53   58   63   67];
rs_vec_100 = [  
    0.0083  -31.0415   -0.0007; %1
    0.0083  -31.0415   -0.0007; %2
    0.0083  -31.0415   -0.0007; %3
    0.0083  -31.0415   -0.0007; %4
    -0.0370   -0.0335   -0.0324; %5
    -0.0006    0.0015    0.0036 %6
    0.0003    0.0001    0.0040 %7 --- first converges at 7
    0.0015    0.0000    0.0046 %8
    0.0017    0.0004    0.0049 %9
    0.0024    0.0012    0.0056; %10 -- 
    0.0182    0.0099    0.0143;
    0.4181    0.3622    0.4220;
    0.6457    0.6497    0.6419;
    0.7807    0.7955    0.7790;
    0.8650    0.8710    0.8598;
    0.9143    0.9182    0.9089;
    0.9423    0.9486    0.9394;
    0.9588    0.9670    0.9583;
    0.9699    0.9789    0.9711;
    0.9769    0.9859    0.9789;
    0.9809    0.9896    0.9833];

%%WRF MODEL

ncens_vec = [1  2  3  4  5  10  15  20  25  30  35  40  45  50  55  60   64];
data_vec =  [10 20 30 40 11 13 18  23  28  33  38  43  48  53  58  63  67];
rs_vec =[   
    0.0083  -31.0415   -0.0007; %1
    0.0272  -20.2120    0.0359; %2
    0.1016   -1.8119    0.0736; %3
    0.1916   -0.5554    0.1712; %4
    0.2713    0.0000    0.2584; %5 -- first converges at 5
    0.7284    0.7040    0.7402;
    0.9003    0.8870    0.8954;
    0.9561    0.9482    0.9510;
    0.9769    0.9744    0.9735;
    0.9857    0.9849    0.9843;
    0.9906    0.9898    0.9896;
    0.9940    0.9946    0.9937;
    0.9956    0.9960    0.9954;
    0.9968    0.9970    0.9967;
    0.9977    0.9978    0.9976;
    0.9981    0.9982    0.9981;
    0.9985    0.9986    0.9984];



%%%Remove all rs_vec
loc = sum((rs_vec<0)')==0;
loc_100 = sum((rs_vec_100<0)')==0;

ncens_vec(~loc) = [];
ncens_vec_100(~loc_100) = [];

data_vec(~loc) = [];
data_vec_100(~loc_100) = [];

rs_vec(~loc,:) = [];
rs_vec_100(~loc_100,:) = [];

%%%PLot number of data points
plottool(1,'Data Required',18,'Number of Centers','Number of Data Points')
plot(ncens_vec,data_vec,'b-')
plot(ncens_vec_100,data_vec_100,'r-')
legend('WRF','Single Wavelength')

%%%Plot ratio
plottool(1,'Ratio',18,'Number of Centers','Number of Data Points/Number of Centers')
plot(ncens_vec,data_vec./ncens_vec,'b-')
plot(ncens_vec_100,data_vec_100./ncens_vec_100,'r-')
legend('WRF','Single Wavelength')

%%%Get mean of rs_vec
plottool(1,'Correlation Function',18,'Number of Centers','Correlation Value')
plot(ncens_vec,mean(rs_vec'),'b-')
plot(ncens_vec_100,mean(rs_vec_100'),'r-')
legend('WRF','Single Wavelength','Location','southeast')

