% Derivative function
% % % Is the u term really necessary - can I just get rid of it?
function StateDhat = Derivatives(Statehat,muhat)

xhat=Statehat(1);
yhat=Statehat(2);
xdhat=Statehat(3);
ydhat= Statehat(4);
rhat = sqrt((xhat^2) +(yhat^2));

xddhat = muhat*xhat/(rhat^3);
yddhat = muhat*yhat/(rhat^3);
StateDhat = [xdhat ydhat xddhat yddhat]';
