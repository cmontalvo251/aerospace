import numpy as np
import matplotlib.pyplot as plt

##From Airfoiltools
Cl = 1.0
##Vary the AR
AR = np.linspace(2,10,100)
##TOtal Lift
CL = Cl/(1+Cl/(np.pi*AR))
plt.plot(AR,CL)

##Based on the plot we set AR
AR = 6.0
##Re compute the lift coefficient
CL = Cl/(1+Cl/(np.pi*AR))
print('CL = ',CL)

###Then we vary W/S
WS = np.linspace(1/8.,2.0,100)
rho = 0.00238
Vs = np.sqrt(2*WS/(rho*CL))
plt.figure()
plt.plot(WS,Vs)

##Based on the plot
WS = 1.5 ## lb / ft^2

##Based on weight estimate
W = 12.0 #lbs
##Compute hte area
S = W/WS
print('Area = ',S)

##Recompute our stall speed
Vs = np.sqrt(2*WS/(rho*CL))
print('Vs = ',Vs)

##Compute Wingspan
b = np.sqrt(AR*S)
print('b = ',b)

##Compute Chord
c = S/b
print('C = ',c)


plt.show()
