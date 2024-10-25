import numpy as np
import matplotlib.pyplot as plt

### Optimize weight / flight time
chord = np.linspace(4,24,100)/12
wingspan = 6
S = wingspan*chord
WS = 16 #oz/ft^3  weight / S^(3/2)
weight = S**(1.5)*WS #ounces
weight /= 16.0 #lb

##Compute aerodynamics
AR = wingspan**2/S
alfa = 5.0*np.pi/180.0
Cl = 0.9227
e = 0.9
Clmax = 1.4
CLmax = Clmax / (1 + Clmax / (np.pi*e*AR))
CL = Cl / (1 + Cl/(np.pi*e*AR))
Cd0 = 0.01016 + 0.6
Cd = 0.01266
CD0 = Cd0 / (1 + Cd0/(np.pi*e*AR))
CD = Cd / (1 + Cd/(np.pi*e*AR))
#K = (CD - CD0)/CL**2
#print(K)

###Now compute max speed
rho = 0.00238 #slugs/ft**3
stall_speed = np.sqrt(2*weight / (rho*CLmax*S))
Tmax = 22.2
Power = Tmax * stall_speed
TW = Tmax/weight
#speed = np.sqrt((TW*WS+WS*np.sqrt(TW**2-4*CD0*K))/(rho*CD0)) ##Eq 5.5
speed = np.sqrt(2*Tmax / (rho*CD0*S))
#speed = (2*Power/(rho*S*CD0))**(1./3.)
#print(speed)
distance = 9000.
flight_time = distance / speed

plt.plot(chord,speed,label='Max Speed (ft/s)')
plt.plot(chord,weight,label='Weight (lbs)')
plt.plot(chord,flight_time,label='Flight Time (sec)')
plt.plot(chord,AR,label='Aspect Ratio')
plt.legend()
plt.grid()
plt.xlabel('Chord (ft)')

M2 = weight / flight_time
plt.figure()
plt.plot(chord,M2)
plt.grid()
plt.xlabel('Chord (ft)')
plt.ylabel('Weight (lbs) / Flight Time (sec)')
plt.show()