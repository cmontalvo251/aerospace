"""
Modification of Dr.Carlos Montalvo's propulsion code

Daniel Hodgins 17 Feb 2022
Notes:
    code was modified to account for pitch and remove velocity cancelation
    find thrust for going into or against wind
"""

import numpy as np

print('Computing Propulsion of a RC Rotor')

rho = 1.225 ##Density of air
P = 600.0 * 1.0 ##Extra 1.5 I think is for induced power; not used in modification
R = (12.5/(12*3.281*2)) ##Radius in meters for a 12.5 inch prop
pitch=5/39.37
A = np.pi*R**2
KV = 800  ##KV rpm per volt from the data sheet 
V = 11.1  #voltage in the battery
omega_rpm = KV*V ##angular velocity in rpm
omega_rads = omega_rpm * 2*np.pi / 60.
print('Omega (RPM) = ',omega_rpm)

#this section has the cancellation issue, omega_rads gets cancelled
CP = P / (rho * A * omega_rads**3 * R**3)
print('CP = ',CP)
CT = (CP*np.sqrt(2))**(2.0/3.0)
print('CT = ',CT)
TN = 0.5*rho*A*(omega_rads**2)*(R**2)*CT
Tkg = TN/9.81
Tlbs = 2.205*Tkg
Tg = Tkg*1000
print('Thrust (g) = ',Tg)
print('Thrust (lbs) =',Tlbs)


omega_rps=omega_rpm/60
acV=pitch*omega_rps#m/s aircraft speed
windV=0#windvelocity
delp=0.5*rho*(acV**2-windV**2)#kg*(m/s)**2/m**3=N/(m**2)
TN=A*delp#N
Tkg=TN/9.81
Tlbs=Tkg*2.205
Tg=Tkg*1000
print('Thrust (g) = ',Tg)
print('Thrust (lbs) =',Tlbs)