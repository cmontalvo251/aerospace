function [fig,ax] = plottool(newfigure,Name,FontSize,xaxislabel,yaxislabel,zaxislabel,titlename,viewparameters,limtype,limits)
%%Makes plots look cool
%plottool(newfigure,Name,FontSize,xaxislabel,yaxislabel,zaxislabel,titlename,viewparameters);
global SAVE2PDF

if isempty(SAVE2PDF)
  SAVE2PDF = 0;
end

if newfigure
  if SAVE2PDF
    saveas(gcf,[Name,'.jpg']);
  end
  fig = figure('Name',Name);
  if nargin >= 3
    set(axes,'FontSize',FontSize)
  end
else
  fig = gcf;
  set(fig,'Name',Name);
  cla;
end
set(fig,'color','white')
hold on
grid on
if nargin >=4
  xlabel(xaxislabel,'FontSize',FontSize)
end
if nargin >=5
  if ~isempty(yaxislabel)
    if iscell(yaxislabel)
      %{'$\dot p \$ (rad/s^2)','latex'} <- this will create pdot
      ylabel(yaxislabel{1},'FontSize',FontSize,'interpreter',yaxislabel{2})
    else
      ylabel(yaxislabel,'FontSize',FontSize)
    end
  end
end
if nargin >=6
  if ~strcmp('',zaxislabel)
    zlabel(zaxislabel)
  end
end
if nargin >=7
  title(titlename)
end
if nargin >=8
  if ~isempty(viewparameters)
    view(viewparameters)
  end
end
if nargin >=10
  if strcmp(limtype,'xlim')
    xlim(limits)
  end
  if strcmp(limtype,'ylim')
    ylim(limits)
  end
end
ax = gca;
% Copyright - Carlos Montalvo 2015
% You may freely distribute this file but please keep my name in here
% as the original owner
