import numpy as np
import matplotlib.pyplot as plt

##THIS IS YOUR DATA
wingspans = np.asarray([45.24,42.5,59.1,66.92,98.43,19.89,43.31,45.28,59.06,60.24])
weights = np.asarray([61,39.5,144.62,190.69,672,21.38,49.03,61.09,65.54,83.71])

##DO NOT TOUCH
coeff = np.polyfit(wingspans,weights,1)
wingspans_fit = np.linspace(np.min(wingspans),np.max(wingspans),1000)
weights_fit = np.polyval(coeff,wingspans_fit)

###NO NEED TO CHANGE EXCEPT FOR THE UNITS on X AND Y AXES
plt.plot(wingspans,weights,'b*')
plt.plot(wingspans_fit,weights_fit,'r--',label='Linear Fit')
plt.xlabel('Wingspan (in)')
plt.ylabel('Weight (oz)')
plt.grid()

##Plot weight and wingspan of designed aircraft
wingspan_star = 45.0
weight_star = np.polyval(coeff,wingspan_star)
plt.plot(wingspan_star,weight_star,'go',markersize=10)

###ADD a regression model that is weight = c1 * exp (c2 * wingspan)
## ln (weight) = ln(c1) + c2 * wingspan
##  ln_weights = b1 + c2* wingspan
ln_weights = np.log(weights)
coeff_exp = np.polyfit(wingspans, ln_weights, 1)
c2 = coeff_exp[0]
b1 = coeff_exp[1]
c1 = np.exp(b1)
weights_fit_exp = c1 * np.exp(c2 * wingspans_fit)
plt.plot(wingspans_fit, weights_fit_exp, 'g--', label='Exponential Fit')
weight_star_exp = c1 * np.exp(c2 * wingspan_star)
plt.plot(wingspan_star, weight_star_exp, 'ro', markersize=10)
plt.legend()

plt.title('Weight(Lin,Exp) = ' + str(round(weight_star,2)) + ', ' + str(round(weight_star_exp,2)) + ' : Wingspan =  ' + str(round(wingspan_star,2)))

plt.show()