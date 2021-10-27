#!/usr/bin/python

import numpy as np
import scipy.integrate as S
import matplotlib.pyplot as plt


#%%%Up here we just need to code the derivatives
def QuadDerivs(x,t):
    global coeff,I,W0,us_left,us_right #%%same globals from above
    #%%%The input is x here which is the angle and angular velocity
    phi = x[0]
    phidot = x[1]

    #%%%We need to measure the distance from the cg to the rotors on the left
    #%%%and the rotors on the right. We can just take the average of the front
    #%%%and back rotors
    #%front_rotor = ;   %%%do this in inches
    #%rear_rotor = ;
    #%avg_rotor = (front_rotor+rear_rotor)/2.;
    #%dl = (avg_rotor/12)/3.28; %%%Inches converted to feet converted meters
    #dl = 0.0762;
    dl = 0.10115

    #%%%%Now we need to use our fancy PD controller and compute the us pulse to
    #%%%%each rotor
    kp = 40; #%%%We can play with these numbers until we get what we want - us/deg
    kd = 300; #%%%This too - us/(rad/s) - CHECK THIS!!!! It might be deg/s
    #%%%However theta is in degrees and thetadot is in rad/s so we need to put
    #%%%that into our code here
    roll_command = kp*phi*180./np.pi + kd*phidot;

    #%%%Then we need to send a us pulse to the left and right rotors
    nominal_us = 1560; #%%%just enough to keep us in the air
    us_left = nominal_us - roll_command;
    us_right = nominal_us + roll_command;

    #%%%%We need a saturation block so we don't break anything
    us_max = 1700;
    us_min = 1400;
    if (us_left > us_max):
        us_left = us_max;
    if (us_left < us_min):
        us_left = us_min;
    if (us_right > us_max):
        us_right = us_max;
    if (us_right < us_min):
        us_right = us_min;
        
    #%%%THen using our fancy new polynomial fit we can get thrust on the rotors
    #%TL = -(polyval(coeff,us_left)-W0)/cosd(9.7); %%%The 9.7 is there because the quad was at an angle when we powered up
    #%TR = -(polyval(coeff,us_right)-W0)/cosd(9.7);
    #%%%use this if you are pushing down.
    TL = np.polyval(coeff,us_left);
    TR = np.polyval(coeff,us_right);

    #%%%So all we need to do from here is compute the angular acceleration
    phidbldot = 2*(TL - TR)*dl/I; #%%%The two is there because their are two rotors on each side
    
    #%%%The derivative is then just
    dxdt = np.asarray([phidot,phidbldot]);

    return dxdt

    #%%%%In order to tune the autopilot we need a few things
    #global coeff I W0 us_left us_right %%%need some globals for the function at the bottom

#MAIN ROUTINE
#%%%%Clear stuff but don't clear variables
#%%%%First we need to get the moment of inertia
#%%%%For a mass with finite inertia and hinged at a distance d from the cg
#%%%%The period is equal to
#%I = (T/(2.*np.pi))*m*g*d - m*d**2;

#%%%Need to get Period from an experiment
T = 0.9
#%%%Weigh the quad
m = 735./1000.; #%%%kilograms
#%%%Need to measure the distance from the cg to the rotation point
d = ((4.)/12.)/3.28; #%%%inches converted to meters
#%%%Assume gravity
g = 9.81;

#%%%%This will output I in kg-m^2
I = (T/(2.*np.pi))*m*g*d - m*d**2;
#dl = 0.0762;
dl = 0.10115

print('Inertia = ',I)
print('Moment Arm = ',dl)

##Period for Pitch Channel was 97 frames at 240 fps
#T = 2*97./240.; #%%seconds
##The distance from the center of rotation to the cg was 3 + 3/8 of an inch
#d = ((3.+3./8.)/12.)/3.28; #%%%inches converted to meters

##Using that data
Iy = 0.07413
Ix = 0.0973
dy = 0.0762
dx = 0.10115
ky = 8.0

kx = (dy/dx)*(Ix/Iy)*ky

print('ky = ',ky,' kx = ',kx)

##Period for the roll channel was 0.9 seconds
#The distance from the center of rotation to the cg was 4 inches

#%%%%%The next thing we need to do is measure the thrust of one rotor
#W0 = m*g; #%%%kilograms converted to Newtons - this is the weight of the quad with the rotors off

#%%%%Ok so we take data and get data from the microsecond pulse and the
#%%%%weight from the scale
#W = [1130 1116 1087 1066 1016 975 908 835 740 640 540 380]*g/1000; %%%You stop taking data when the quad almost lifts off
W = np.asarray([0,80,110,150,205,243,315,380,450,550])*g/1000.; #%%Use this if the quad is pushing down
us = np.asarray([1019,1290,1350,1405,1470,1522,1590,1660,1720,1800]);

#%%%Just for kicks lets plot this
plt.figure()
plt.plot(us,W,'b*')
plt.xlabel('Throttle Position (us)')
plt.ylabel('Weight (N)')

#%%%This will give us a polynomial that we can use to evaluate thrust
#%W=f(us)
coeff = np.polyfit(us,W,2);

#%%%Let's then plot the trend line to make sure this makes sense.
x_fit = np.linspace(us[0],us[-1],100);
y_fit = np.polyval(coeff,x_fit);
plt.plot(x_fit,y_fit,'r-')

#%%%ok so now we have a way to get Thrust from the us pulse.
#%%Now we can kick off our simulation
#%%%First we need initial conditions
#%%%Let's assume the quad is at 12 degrees from nominal with zero angular
#%%%velocity
x0 = np.asarray([12*np.pi/180.0,0]);

#%%%Then we use the easy ode45 to integrate our equations of motion
tout = np.linspace(0,2,1000)
xout = S.odeint(QuadDerivs,x0,tout);

us_R_vec = []
us_L_vec = []
for idx in range(0,len(tout)):
    dxdt = QuadDerivs(xout[idx,:],tout[idx]);
    us_R_vec.append(us_right);
    us_L_vec.append(us_left);

us_R_vec = np.asarray(us_R_vec)
us_L_vec = np.asarray(us_L_vec)

#%%%Once ode45 finishes we can plot the roll angle as a function of time
phi = xout[:,0];
plt.figure()
plt.plot(tout,phi*180/np.pi) #%%Convert it to degrees
plt.xlabel('Time (sec)')
plt.ylabel('Roll Angle (deg)')

plt.figure()
plt.plot(tout,us_R_vec,'b-')
plt.plot(tout,us_L_vec,'r-')
plt.xlabel('Time (sec)')
plt.ylabel('us Pulse')
plt.legend()

plt.show()
