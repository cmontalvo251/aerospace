%%%Compute Inertia
purge

%%%E-sail
%m = 12;
%ax = 10/100;
%by = 20/100;
%cz = 30/100;

%%%E-sail joined
m = 24;
ax = 20/100;
by = 20/100;
cz = 30/100;

Ixx = (m/12)*(by^2+cz^2)
Iyy = (m/12)*(ax^2+cz^2)
Izz = (m/12)*(ax^2+by^2)