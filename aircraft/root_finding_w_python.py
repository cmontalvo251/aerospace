##Import the optimization toolbox
import scipy.optimize as RF
##Import the numeric python toolbox
import numpy as np
##Import the plotting toolbox
import matplotlib.pyplot as plt

##This is the function we're trying to find the zero of
def f(alpha):
    CD0 = 0.0119903112
    CDA = 1.2468520598
    CL0 = 0.
    CLA = 4.3068025604
    y = CD0 + CDA*alpha**2 - np.tan(alpha)*CL0 - np.tan(alpha)*CLA*alpha
    return y

###Let's plot it first
##So let's use linspace like MATLAB and make a vector from -10 to 10 degrees
alpha = np.linspace(-10,10,100)*np.pi/180.0

###Compute y just by calling the function we defined above
y = f(alpha)

####Close all figures and change the fontsize to 18
plt.close("all")
plt.rcParams.update({'font.size': 18})

###Plot the function so we can arbitrarily determine our bounds by glancing 
###at the function in graph form
plt.plot(alpha*180.0/np.pi,y)
plt.xlabel('Alpha (deg)')
plt.ylabel('Root Finding Function')
plt.grid()

##Minimize function f using the minimization toolbox
##Specifically the bisection methodd.
sol = RF.bisect(f,0.0,0.15)

print 'Alpha = ',sol

##Plot a red x where the solution in
plt.plot(sol*180.0/np.pi,0,'rx')

##Show plot
plt.show()
