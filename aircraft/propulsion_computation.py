import numpy as np

print('Computing Propulsion of a RC Rotor')

# USing this motor - https://www.horizonhobby.com/product/avian-5055-650kv-outrunner-brushless-motor/SPMXAM4745.html
# Avian 650KV

rho = 1.225 ##Density of air
P = 4400.0 ##Extra 1.5 I think is for induced power
Diameter_Prop_Inches = 20.0
KV = 200.0  ##KV rpm per volt from the data sheet 
V = 44.4  #voltage in the battery
R = (Diameter_Prop_Inches/(12*3.28*2)) ##Radius in meters for a 12.5 inch prop
A = np.pi*R**2
omega_rpm = KV*V ##angular velocity in rpm
omega_rads = omega_rpm * 2*np.pi / 60.
print('Omega (RPM) = ',omega_rpm)
CP = P / (rho * A * omega_rads**3 * R**3)
print('CP = ',CP)
CT = (CP*np.sqrt(2))**(2.0/3.0)
print('CT = ',CT)
TN = 0.5*rho*A*(omega_rads**2)*(R**2)*CT
Tkg = TN/9.81
Tlbs = 2.2*Tkg
Tg = Tkg*1000
print('Thrust (g) = ',Tg)
print('Thrust (lbs) =',Tlbs)
