%%%Rotor Test Data
g = 9.81;

%%%%%The next thing we need to do is measure the thrust of one rotor
W0 = (1130)/1000*g; %%%kilograms converted to Newtons - this is the weight of the quad with the rotors off

%%%%Ok so we take data and get data from the microsecond pulse and the
%%%%weight from the scale
W = [1130 1116 1087 1066 1016 975 908 835 740 640 540 380]*g/1000; %%%You stop taking data when the quad almost lifts off
us = [1100 1151 1209 1252 1310 1354 1412 1456 1514 1572 1615 1659];

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


%%%THen using our fancy new polynomial fit we can get thrust on the rotors
us_left = 1100:1900;
Weight_one_leg = polyval(coeff,us_left);
Thrust_TOTAL = -4*(Weight_one_leg - W0)/cosd(9.7);

m = 1.482; %%%kilograms
Weight_of_AC = m*g;

figure()
plot(us_left,Thrust_TOTAL)
hold on
plot(us_left,Weight_of_AC,'r-')
