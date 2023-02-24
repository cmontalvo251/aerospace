import numpy as np
import matplotlib.pyplot as plt

a0 = 0.0 #angle of attack at minimum drag
a = np.linspace(-6,6,1000) #range for your drag polar
Cd0 = 0.01 #drag coefficient at minimum drag

#Cd = Cd0 + Cda*(a-a0)^2
##Grab a data point
astar = 6 
Cdstar = 0.02
Cda = (Cdstar-Cd0)/(astar-a0)**2
print('Cda /deg**2 = ',Cda)
print('Cda /rad**2 = ',Cda*(180./np.pi)**2)
Cd = Cd0 + Cda*(a-a0)**2

plt.plot(a,Cd)
plt.grid()
plt.show()
