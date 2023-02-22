
import numpy as np
import matplotlib.pyplot as plt

a0 = 4.0
a = np.linsapce(-4,12,1000)
Cd0 = 0.025

#Cd = Cd0 + Cda*(a-a0)^2
##Grab a data point
astar = 12
Cdstar = 0.04
Cda = (Cd-Cd0)/(astar-a0)**2
print('Cda = ',Cda)
Cd = Cd0 + Cda*(a-a0)**2

plt.plot(a,Cd)
plt.grid()
plt.show()