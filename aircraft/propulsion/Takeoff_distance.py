import numpy as np
import matplotlib.pyplot as plt

def takeoff(ThrustN):
    Vstall_ms = 9.7186496 #m/s
    W_lbs = 15.0 #weight of aircraft
    m_kg = W_lbs/2.2 #convert weight in pounds to kg 
    a_ms2 = ThrustN/m_kg
    t = Vstall_ms/a_ms2
    xf_m = 0.5*a_ms2*t**2
    xf_ft = xf_ft = xf_m*3.28
    return xf_ft

Tlbs = np.linspace(1,67,1000) #0 to 15 lbs of thrust in Newtons
TN = (Tlbs/2.2)*9.81 #convert pounds to Newtons
xf_ft = takeoff(TN)

Tlbstar = 14.4
TNstar = Tlbstar/2.2*9.81 #pounds of thrust of DBF plane
xf_ftstar = takeoff(TNstar)

print("Thrust lbs = ",Tlbstar," Takeoff Distance = ",xf_ftstar)

plt.figure()
plt.grid()
plt.title(' Net Thrust vs Distance Eng Units')
plt.xlabel('Distance (ft)')
plt.ylabel('Thrust (lbf)')
plt.plot(xf_ft,Tlbs)
plt.plot(xf_ftstar,Tlbstar,'b*')
plt.show()