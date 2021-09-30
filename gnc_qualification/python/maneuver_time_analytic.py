import numpy as np
import matplotlib.pyplot as plt

###Double check all of these numbers and comment the code heavily
power = 6.0
f = np.linspace(0,0.5,100)
max_torque = 0.025
final_angle = 180.0*np.pi/180.0

###Do this for the X-axis first
inertiax = 0.15
angular_accelerationx = max_torque/inertiax
t1x = np.sqrt(2*final_angle*f/angular_accelerationx)
energyx = 2*power*t1x
timex = 2*t1x + final_angle*(1.0-2.0*f)/(angular_accelerationx*t1x)
plt.plot(timex,energyx)

###Do this for Y-Axis
inertiay = 0.0

###Do it again for Z-axis

###add a legend, x and y axis labels(units), a grid
###a title with the manuever angle


plt.show()
