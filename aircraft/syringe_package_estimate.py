import numpy as np
import matplotlib.pyplot as plt

Vfast = 130.0

distance = 2500*2+2*100 ##feet
tfast = distance/Vfast
Nfast = 1.0
M2fast = Nfast/tfast
print(M2fast)

Vus = np.arange(1,130,1)
Nus = np.arange(1,10,1)
[VV,NN] = np.meshgrid(Vus,Nus)
tus = distance / VV
M2us = NN/tus
fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
ax.plot_surface(VV,NN,M2us-M2fast)
#ax.plot_surface(VV,NN,0.0*NN)
plt.grid()

Vreq = distance * M2fast / Nus
plt.figure()
plt.plot(Nus,Vreq)

plt.show()