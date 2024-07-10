function plotstates(arg1,filename)
global NSTATES DELTA NOISE DUAL mass I_b I_b_inv num_ac rcg2ac

%close all

if exist('arg1','var')
    states = arg1;
else
    states = 1:15;
end
if ~exist('filename','var')
    filename = 'source/Out_Files/Meta3_3';
end

%%%PLOTTING ROUTINE
data = dlmread(filename);
xout = data(:,2:end);
xout = xout';
tout = data(:,1);
dim1 = 'ft';
[c,r] = size(xout);
num_ac = round(c/NSTATES);

%%%Skip rate
skip = 1;
xout = xout(:,1:skip:end);
tout = tout(1:skip:end);
units = 1;
Names = {'X','Y','Z','PHI','THETA','PSI','U','V','W','P','Q','R'};
ylabels = {['x(',dim1,')'],['y(',dim1,')'],['z(',dim1,')'],'\phi(deg)','\theta(deg)','\psi(deg)',['u(',dim1,'/s)'],['v(',dim1,'/s)'],['w(',dim1,'/s)'],'p(rad/s)','q(rad/s)','r(rad/s)'};
LineWidth = 2;
%colors = {'b','r','g','m','c','k','y'};
colors = {'k','k'};
linetype = {'-','--','-.','-','--','-.','-','--','-.','-'};
for ii = states
  if ii <= 15
    if ii >= 4 && ii <= 6
      factor = 180/pi;
    else
      factor = units;
    end
    if ii >= 10
      factor = 1;
    end
    h1 = plottool(1,Names{ii},12,'Time(sec)',ylabels{ii});
    p(1) = plot(tout,factor.*xout(ii,:),[colors{1},linetype{1}],'LineWidth',LineWidth);
    for jj = 2:num_ac
      val = (jj-1)*NSTATES+ii;
      jdx = jj;
      while jdx > 7
          jdx = jdx - 7;
      end
      p(jj) = plot(tout,factor.*xout(val,:),[colors{1},linetype{jj}],'LineWidth',LineWidth);
    end
    if exist('LegendNames','var')
      legend(p,LegendNames)
    end
  end
end

