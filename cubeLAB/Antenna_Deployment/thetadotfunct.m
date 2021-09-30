function thetadot = thetadotfunct(t, td)

% if t>=td+0.5
% thetadot = 0;
% elseif t>=0.5
% thetadot = 90/(t-0.5);
% else
%     thetadot = 0;
% end


thetadot = ((45/.05)*sech((t-td)/.05).^2)*pi/180;