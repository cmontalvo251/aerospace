function plot3dim(filename)
global NSTATES DELTA NOISE DUAL mass I_b I_b_inv num_ac rcg2ac

close all

%%%PLOTTING ROUTINE
data = dlmread(filename);
xout = data(:,2:end);
xout = xout';
tout = data(:,1);
dim1 = 'm';
[c,r] = size(xout);
num_ac = round(c/NSTATES);

%%%Skip rate
skip = 1;
xout = xout(:,1:skip:end);
tout = tout(1:skip:end);
units = 1;
LineWidth = 2;
%colors = {'b','r','g','m','c','k','y'};
colors = {'k','k'};
linetype = {'-','--','-.','-','--','-.','-','--','-.','-'};
h1 = plottool(1,'3D',12,'X(m)','Y(m)',);
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

