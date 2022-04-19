%%%%USE THIS CODE TO DETERMINE THE COP FOR A SWEPT WING


clear
clc
close all


%%Main Wing(m)
b = 4
c_root = 1.0
maint = 0.09; %thickness of mainwing (percentage)
taper = 0.7;

%%%%This is the center of pressure of the airfoil cross section
%%This code assumes that the airfoill is constant throughout the wing
%%% For a FLAT PLATE use this
AIRFOIL_CENTER = 0.25;

%%%It seems like a Static Margin of |0.5"| is ok.

%%%NACA 2412
%AIRFOIL_CENTER = 0.32;

%use positive value to indicate aft of nose
xlocation = 0.0; %location of leading edge of wing w.r.t to nose of fuselage
sweep = 30*pi/180; %%This is something you can change
 
%%%%%%%%%%%%%%%%%%%%STEP 1 CALCULATIONS%%%%%%%%%%%%%%%%%%%%

%%Compute Mass of All parts

%%Wing(assume rectangular)
mainwingthickness = maint*c_root;
ctip = c_root*taper;
S = (1/2)*b*(c_root+taper*c_root);

%%%%%Compute CG Locations of all Sections(measured from body fixed
%frame at nose of aircraft)

%main wing %%%THIS NEEDS TO GET FIXED FOR SWEPT WINGS
N = 10000;
ycoord = linspace(0,b/2,N);
dy = ycoord(2)-ycoord(1);

%%%%%%%%%%%%%%%%%%%%%%%STEP 3 CALCULATIONS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%Compute Aerodynamic Center

%Main Wing
AR = (b^2)/S; %aspect ratio
cbar = 0;
%mean aerodynamic chord
aero_center = 0;
for idx = 1:N-1
    xi = xlocation + sin(sweep)*ycoord(idx);
    ci = c_root + (2/b)*c_root*(taper-1)*ycoord(idx);
    aeroi = xi + ci*AIRFOIL_CENTER;
    aero_center = aero_center + aeroi*dy;
    cbar = cbar + ci*dy;
end
aero_center = aero_center*2/b
aero_over_chord = aero_center/c_root
cbar = cbar*2/b

%%%Plot Everything
figure()
%%%Plotting the wing
plot([0,c_root],[0,0])
hold on
plot([0,sin(sweep)*b/2],[0,b/2])
plot([sin(sweep)*b/2,sin(sweep)*b/2+ctip],[b/2,b/2])
plot([c_root,sin(sweep)*b/2+ctip],[0,b/2])
%%%Plotted the A/C line
plot([0+AIRFOIL_CENTER*c_root,sin(sweep)*b/2+ctip*AIRFOIL_CENTER],[0,b/2],'b--')
%%%Plot the otherside
%%%Plotting the wing
plot([0,c_root],[0,0])
hold on
plot([0,sin(sweep)*b/2],-[0,b/2])
plot([sin(sweep)*b/2,sin(sweep)*b/2+ctip],-[b/2,b/2])
plot([c_root,sin(sweep)*b/2+ctip],-[0,b/2])
%%%Plotted the A/C line
plot([0+AIRFOIL_CENTER*c_root,sin(sweep)*b/2+ctip*AIRFOIL_CENTER],-[0,b/2],'b--')
%%%Plot the Aerodynamic Center
plot([aero_center,aero_center],[0,0],'bo','MarkerSize',20)
axis equal