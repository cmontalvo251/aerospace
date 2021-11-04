import numpy as np
import matplotlib.pyplot as plt

pwm = np.array([1502,1600,1700,1800,1900,(1502-992)+1502])
force = np.array([0,1.2,1.33,1.49,1.52,1.60])

dpwm = pwm - pwm[0]

A = 1.7
s = -np.log(1-force/A)/dpwm
s = np.mean(s[1:])

dpwm_fit = np.linspace(dpwm[0],dpwm[-1],100)
force_fit = A*(1-np.exp(-s*dpwm_fit))

plt.plot(dpwm,force,'b*')
plt.plot(dpwm_fit,force_fit,'r-')
plt.show()