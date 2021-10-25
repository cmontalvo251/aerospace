%%%%In order to tune the autopilot we need a few things
function auto_tune()
global coeff I W0 us_left us_right %%%need some globals for the function at the bottom

%%%%Clear stuff but don't clear variables
clc 
close all

%%%%First we need to get the moment of inertia
%%%%For a mass with finite inertia and hinged at a distance d from the cg
%%%%The period is equal to
%I = (T/(2*pi))*m*g*d - m*d^2;

%%%Need to get Period from an experiment
T = 2*97/240; %%seconds
%%%Weigh the quad
m = 735/1000; %%%kilograms
%%%Need to measure the distance from the cg to the rotation point
d = ((3+3/8)/12)/3.28; %%%inches converted to meters
%%%Assume gravity
g = 9.81;

%%%%This will output I in kg-m^2
I = (T/(2*pi))*m*g*d - m*d^2;

%%%%%The next thing we need to do is measure the thrust of one rotor
%W0 = m*g; %%%kilograms converted to Newtons - this is the weight of the quad with the rotors off

%%%%Ok so we take data and get data from the microsecond pulse and the
%%%%weight from the scale
%W = [1130 1116 1087 1066 1016 975 908 835 740 640 540 380]*g/1000; %%%You stop taking data when the quad almost lifts off
W = [0 80 110 150 205 243 315 380 450 550]*g/1000; %%Use this if the quad is pushing down
us = [1019 1290 1350 1405 1470 1522 1590 1660 1720 1800];

%%%Just for kicks lets plot this
figure()
plot(us,W,'b*')
xlabel('Throttle Position (us)')
ylabel('Weight (N)')

%%%This will give us a polynomial that we can use to evaluate thrust
%W=f(us)
coeff = polyfit(us,W,2);

%%%Let's then plot the trend line to make sure this makes sense.
x_fit = linspace(us(1),us(end),100);
y_fit = polyval(coeff,x_fit);
hold on
plot(x_fit,y_fit,'r-')

%%%ok so now we have a way to get Thrust from the us pulse.
%%Now we can kick off our simulation
%%%First we need initial conditions
%%%Let's assume the quad is at 12 degrees from nominal with zero angular
%%%velocity
x0 = [12*pi/180;0];

%%%Then we use the easy ode45 to integrate our equations of motion
[tout,xout] = ode45(@QuadDerivs,[0 2],x0);

us_R_vec = zeros(length(tout),1);
us_L_vec = zeros(length(tout),1);
for idx = 1:length(tout)
    dxdt = QuadDerivs(tout(idx),xout(idx,:));
    us_R_vec(idx) = us_right;
    us_L_vec(idx) = us_left;
end

%%%Once ode45 finishes we can plot the roll angle as a function of time
phi = xout(:,1);
figure()
plot(tout,phi*180/pi) %%Convert it to degrees
xlabel('Time (sec)')
ylabel('Roll Angle (deg)')

figure()
plot(tout,us_R_vec,'b-')
hold on
plot(tout,us_L_vec,'r-')
xlabel('Time (sec)')
ylabel('us Pulse')
legend('Right Rotors','Left Rotors')

%%%Down here we just need to code the derivatives
function dxdt = QuadDerivs(t,x)
global coeff I W0 us_left us_right %%same globals from above
%%%The input is x here which is the angle and angular velocity
phi = x(1);
phidot = x(2);

%%%We need to measure the distance from the cg to the rotors on the left
%%%and the rotors on the right. We can just take the average of the front
%%%and back rotors
%front_rotor = ;   %%%do this in inches
%rear_rotor = ;
%avg_rotor = (front_rotor+rear_rotor)/2;
%dl = (avg_rotor/12)/3.28; %%%Inches converted to feet
dl = 0.0762;

%%%%Now we need to use our fancy PD controller and compute the us pulse to
%%%%each rotor
kp = 40; %%%We can play with these numbers until we get what we want - us/deg
kd = 300; %%%This too - us/(rad/s) - CHECK THIS!!!! It might be deg/s
%%%However theta is in degrees and thetadot is in rad/s so we need to put
%%%that into our code here
roll_command = kp*phi*180/pi + kd*phidot;

%%%Then we need to send a us pulse to the left and right rotors
nominal_us = 1560; %%%just enough to keep us in the air
us_left = nominal_us - roll_command;
us_right = nominal_us + roll_command;

%%%%We need a saturation block so we don't break anything
us_max = 1700;
us_min = 1400;
if us_left > us_max
    us_left = us_max;
end
if us_left < us_min
    us_left = us_min;
end
if us_right > us_max
    us_right = us_max;
end
if us_right < us_min
    us_right = us_min;
end

%%%THen using our fancy new polynomial fit we can get thrust on the rotors
%TL = -(polyval(coeff,us_left)-W0)/cosd(9.7); %%%The 9.7 is there because the quad was at an angle when we powered up
%TR = -(polyval(coeff,us_right)-W0)/cosd(9.7);
%%%use this if you are pushing down.
TL = polyval(coeff,us_left);
TR = polyval(coeff,us_right);

%%%So all we need to do from here is compute the angular acceleration
phidbldot = 2*(TL - TR)*dl/I; %%%The two is there because their are two rotors on each side

%%%The derivative is then just
dxdt = [phidot;phidbldot];



