import numpy as np
import matplotlib.pyplot as plt

voltage = np.asarray([0.05,0.52,1.03,1.5,2,2.56])
length = np.asarray([0,0.5,1,1.5,2,2.5])

coeff = np.polyfit(voltage,length,1)
print(coeff)

voltage_fit = np.linspace(voltage[0],voltage[-1],1000)
length_fit = np.polyval(coeff,voltage_fit)

length_vals = np.polyval(coeff,voltage)

residuals = length - length_vals

plt.plot(voltage,length,'b*')
plt.plot(voltage_fit,length_fit,'r-')
plt.xlabel('Voltage (V)')
plt.ylabel('Length (cm)')
plt.grid()

plt.figure()
plt.plot(residuals,'b*')
plt.grid()
plt.show()