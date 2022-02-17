function statedot = computeXYZDerivs(state,LMNcontrol)
%%%DERIVATIVE OF SATELLITE XYZ
parameters

x = state(1);
y = state(2);
z = state(3);
xd = state(4);
yd = state(5);
zd = state(6);

rSat = sqrt((x^2) +(y^2) + (z^2));

xdd = muSat*x/(rSat^3);
ydd = muSat*y/(rSat^3);
zdd = muSat*z/(rSat^3);

statedot = [xd;yd;zd;xdd;ydd;zdd];