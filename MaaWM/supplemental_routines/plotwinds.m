function plotwinds(filename)

windout = dlmread(filename);
windout(1,:) = []; %%Throw out first row
[r,c] = size(windout);
windsUVW = zeros(r,4,1);
state = zeros(r,3,1);
idx = ones(1,10);
ac_max = 1;
for ii = 1:r
    %%%Grab Aircraft Counter
    num_ac = windout(ii,8);
    windsUVW(idx(num_ac),:,num_ac) = windout(ii,[4:7]);
    state(idx(num_ac),:,num_ac) = windout(ii,1:3);
    idx(num_ac) = idx(num_ac) + 1;
    if num_ac > ac_max
        ac_max = num_ac;
    end
end
windsUVW(max(idx):end,:,:) = [];
state(max(idx):end,:,:) = [];

colors = {'b','r','g','m','c','k','y'};
linetype = {'-','--','-.','-','--','-.'};
LineWidth = 2;

Names = {'U','V','W'};
ylabels = Names;
for ii = 1:3
  ctr1 = 0;
  ctr2 = 1;
  h1 = plottool(1,Names{ii},12,'Time(sec)',ylabels{ii});
  for jj = 1:ac_max
      ctr1 = ctr1 + 1;
      if ctr1 > 7
          ctr1 = 1;
          ctr2 = ctr2 + 1;
      end
      plot(windsUVW(:,1,jj),windsUVW(:,ii+1,jj),[colors{ctr1},linetype{ctr2}],'LineWidth',LineWidth);
  end
  if exist('LegendNames','var')
    legend(p,LegendNames,'Location','NorthEastOutside')
  end
end
