import numpy as np
import matplotlib.pyplot as plt

###This will be an example NACA 2412 aerodynamic analysis
##of a fixed wing aircraft

##Density and Viscosity of air
rho = 1.225 #kg/m^3
mu = 1.81*10**-5 #SI

##Conceptual Estimates of our aircraft
c = 0.25 #m
b = 1.0
S = b*c
AR = b**2/S
Vest = 10.0 #m/s - best guess right now

#compute the reynolds
Re = rho * Vest * c / mu
print('Estimated Reynolds Number = ',np.round(Re))

##Read the txt file data from airfoil tools into memory
data = np.loadtxt('naca_2412.txt')
alpha_deg = data[:,0]
Cl = data[:,1]
Cd = data[:,2]
LD = Cl/Cd

Cl0 = 0.285 #this is from the graph
Cla = (0.8 - 0.285) / (5*np.pi/180.0)
print('Cl0 , Cla = ',Cl0,Cla)
alpha_rad = alpha_deg*np.pi/180.0
Cl_linear = Cl0 + Cla*alpha_rad

##Generate 3 plots
plt.figure()
plt.plot(alpha_deg,Cl,'b*')
plt.plot(alpha_deg,Cl_linear,'r-')
plt.grid()
plt.xlabel('AoA (deg)')
plt.ylabel('Cl')

plt.figure()
plt.plot(alpha_deg,Cd,'b*')
plt.grid()
plt.xlabel('AoA (deg)')
plt.ylabel('Cd')

plt.figure()
plt.plot(alpha_deg,LD,'b*')
plt.grid()
plt.xlabel('AoA (deg)')
plt.ylabel('L/D')

##Compute Vstall and VL/D
Clmax = np.max(Cl)
CLmax = Clmax / (1 + Clmax/(np.pi*AR))
print('Clmax , CLmax = ',Clmax,CLmax)
W = 1.6*9.81 #weight in Newtons for a 1.6 kg aircraft
Vstall = np.sqrt(2*W/(rho*S*CLmax))
print('Vstall = ',Vstall)

alpha_rad_LDmax = 6.05*np.pi/180.0
ClmaxLD = Cl0 + Cla*alpha_rad_LDmax
CLmaxLD = ClmaxLD / (1 + ClmaxLD/(np.pi*AR))
print('ClmaxLD , CLmaxLD = ',ClmaxLD,CLmaxLD)
VLD = np.sqrt(2*W/(rho*S*CLmaxLD))
print('V max L/D = ',VLD)

plt.show()