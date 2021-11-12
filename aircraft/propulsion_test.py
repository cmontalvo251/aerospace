import numpy as np
import matplotlib.pyplot as plt

##Throttle Position
throttle = 0.95
pressure= 101325 ##Pa
print('Pressure = ',pressure)
T = 20  ##Celsius
print('Temperature = ',T)

#Tg = 3197
#Tkg = Tg/1000.0
#TN = Tkg*9.81
#print(TN)
P = 300.0  ##POWER FROM WEBSITE
P *= 1.50 ##Add 50% for induced power
##Compute Density
Rmo = 287
## p = rho * R * T
rho = pressure / (Rmo*(T+273))
print('rho = ',rho)
#rho = 1.225  #DENSITY AT SEA-LEVEL
R = (4./12.)/(3.28*2)  ##RADIUS OF PROP IN METERS
A = np.pi*R**2  ##AREA OF ROTOR
kv = 600.0  ##KV of MOTOR
V = 3.0*3.7*throttle ##VOLTAGE OF BATTER + 90% of throttle
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
W = Tg*4
print('Weight (g) = ',W)

###Now recompute with a lower throttle position and diff pressure and temp
throttle = 1.0
pressure_new = pressure*0.9
Tnew = T*1.1


V = 3.0*3.7*throttle ##VOLTAGE OF BATTER + 90% of throttle
omega_rpm = kv*V
omega_rads = omega_rpm*np.pi/30.0
print('Omega RPM NEW = ',omega_rpm)
rhonew = pressure_new / (Rmo*(Tnew+273))
print('rho new = ',rhonew)
TN = 0.5*rhonew*A*(omega_rads*R)**2*CT
Tkg = TN/9.81
Tg = Tkg*1000
print('Thrust (g) NEW = ',Tg)
Wnew = Tg*4
print('Weight (g) (Full Throttle, diff P/T) = ',Wnew)
print('Weight (g) (90% Throttle) = ',W)