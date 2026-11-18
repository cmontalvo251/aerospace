import numpy as np

print('Computing Propulsion of a RC Rotor')

# Using this motor - https://store.tmotor.com/product/tmotor-u7-v2-motor-u-power.html?srsltid=AYCX9vJvzrj8Ble8QHrlJ0JDBNnm_6f7SHYYm095y9qJWkzlASpuIaU7
# U7 KV 280

#########INPUTS################
Diameter_Prop_Inches = 22.0
KV = 280.0 #KV rating
V = 24.0  #voltage in the battery
I = 35.0 #Maximum continuous current in motor
###############################

###COMPUTATIONS#########
rho = 1.225 ##Density of air
R = (Diameter_Prop_Inches/(12*3.28*2)) ##Radius in meters for a 12.5 inch prop
print('Radius of Prop (m) = ',R)
A = np.pi*R**2
print('Disk Area (m^2) = ',A)
P = V*I
print('Power (W) = ',P)
omega_rpm = KV*V ##angular velocity in rpm
omega_rads = omega_rpm * 2*np.pi / 60.
print('Omega (rads) = ',omega_rads)
print('Omega (RPM) = ',omega_rpm)
CP = P / (rho * A * omega_rads**3 * R**3)
print('CP = ',CP)
CT = (CP*np.sqrt(2))**(2.0/3.0)
print('CT = ',CT)
TN = rho*A*(omega_rads**2)*(R**2)*CT
Tkg = TN/9.81
Tlbs = 2.2*Tkg
Tg = Tkg*1000
print('Thrust (g) = ',Tg)
print('Thrust (lbs) =',Tlbs)
print('Remember the above is theoretical and a vast overestimate. My recommendation is to add 30% of losses if planning for actual flight tests')
print('Thrust with 30% losses (g) = ',Tg*0.7)
print('Thrust with 30% losses (lbs) =',Tlbs*0.7)
