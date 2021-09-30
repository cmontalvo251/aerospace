function theta = thetafunct(t,td)

% if t>=td+0.5
% theta = 90;
% elseif t>=0.5
%   m = 90/(t-0.5);
%  theta = m*t;
% else
%     theta = 0;
% end

theta = (45*(1 + tanh((t-td)/.05)))*pi/180;