import numpy as np
import sys
import matplotlib.pyplot as plt


alpha_deg = np.linspace(0,12.0,1000)
alpha = np.pi/180.0 * alpha_deg

CL0 = 0.
CLA = 4.3068025604

CL = CL0 + CLA*alpha

##This way
CD0 = 0.0119903112
CDA = 1.2468520598
CD = CD0 + CDA*alpha**2

plt.close("all")
plt.plot(alpha,CL/CD)
plt.show()