
import numpy as np
import matplotlib.pyplot as plt

##Lift vs alfa
a = np.linspace(-10,10,1000)
Cl0 = 0.0
Cla = (1.075-Cl0)/(10-0)
print('Cla (/deg) = ',Cla)
print('Cla (/rad) = ',Cla*(180/np.pi))
print('Cl0 = ',Cl0)
Cl = Cl0 + Cla*a
plt.figure()
plt.plot(a,Cl)
plt.xlabel('Angle of Attack (deg)')
plt.ylabel('Lift Coefficient')
plt.grid()

##Drag Polar
Cd0 = 0.00832
print('Cd0 = ',Cd0)
#Cd = Cd0 + Cda*a^2
##Grab a data point
astar = 13
Cdstar = 0.0167
Cda = (Cdstar-Cd0)/(astar**2)
print('Cda (/deg^2) = ',Cda)
Cda_rad = Cda*(180/np.pi)**2
print('Cda (/rad^2) = ',Cda_rad)
Cd = Cd0 + Cda*a**2
plt.figure()
plt.plot(a,Cd)
plt.grid()
plt.xlabel('Angle of Attack (deg)')
plt.ylabel('Drag Coefficient')

plt.show()
