import numpy as np
import matplotlib.pyplot as plt

#Tg = 3197
#Tkg = Tg/1000.0
#TN = Tkg*9.81
#print(TN)
P = 772  ##POWER FROM WEBSITE
P *= 1.50 ##Add 50% for induced power
rho = 1.225  #DENSITY AT SEA-LEVEL
R = (12./12.)/(3.28*2)  ##RADIUS OF PROP IN METERS
A = np.pi*R**2  ##AREA OF ROTOR
kv = 600.0  ##KV of MOTOR
V = 6*3.7  ##VOLTAGE OF BATTER
omega_rpm = kv*V
omega_rads = omega_rpm*np.pi/30.0
print('Omega RPM = ',omega_rpm)
#CT = 2*TN / (rho*A*(omega_rads*R)**2)
#Cp = CT**(3.0/2.0)/np.sqrt(2)
Cp = P / (rho*A*(omega_rads*R)**3)
print('Power Coefficient = ',Cp)
CT = (Cp*np.sqrt(2))**(2.0/3.0)
print('Thrust Coefficient=',CT)
TN = 0.5*rho*A*(omega_rads*R)**2*CT
P = rho*A*(omega_rads*R)**3*Cp
Tkg = TN/9.81
Tg = Tkg*1000
print('Thrust (g) = ',Tg)
