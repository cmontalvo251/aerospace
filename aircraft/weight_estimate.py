import numpy as np
import matplotlib.pyplot as plt

##THIS IS YOUR DATA
wingspans = np.asarray([48,55.2,44.5,86,58.5,63])
weights = np.asarray([3.37,5.64,2.5,25,6.61,6.28])

##DO NOT TOUCH
coeff = np.polyfit(wingspans,weights,1)
wingspans_fit = np.linspace(np.min(wingspans),np.max(wingspans),1000)
weights_fit = np.polyval(coeff,wingspans_fit)

###NO NEED TO CHANGE EXCEPT FOR THE UNITS on X AND Y AXES
plt.plot(wingspans,weights,'b*')
plt.plot(wingspans_fit,weights_fit,'r--')
plt.xlabel('Wingspan (in)')
plt.ylabel('Weight (lb)')
plt.grid()

plt.show()