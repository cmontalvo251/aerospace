
import numpy as np
import matplotlib.pyplot as plt

a0 = 4.0
a = np.linspace(-4,12,1000)
Cd0 = 0.025

#Cd = Cd0 + Cda*(a-a0)^2
##Grab a data point
astar = 12
Cdstar = 0.04
Cda = (Cdstar-Cd0)/(astar-a0)**2
print('Cda (/deg^2) = ',Cda)
Cda_rad = Cda*(180/np.pi)**2
print('Cda (/rad^2) = ',Cda_rad)
Cd = Cd0 + Cda*(a-a0)**2

plt.plot(a,Cd)
plt.grid()
plt.show()