import numpy as np
import matplotlib.pyplot as plt

###This code is going to compute how much power is used for a specific maneuever.
final_angle = 45.0*np.pi/180.0 #degrees
inertia = np.array([0.13,0.15,0.18]) #kg-m^2 - Get this number from that txt file I sent you
max_torque = 0.025 #N-m
power = 6.0 #W
angular_acceleration = max_torque/inertia  #rad/s
print('Alfa = ',angular_acceleration)

#factor can be anything from 0 to 0.5
for f in np.linspace(0,0.5,3):
    print('Factor = ',f)
    ###First maneuver
    #manuever_angle*factor =  0.5*angular_acceleration*t1^2
    t1 = np.sqrt(2*final_angle*f/angular_acceleration)
    print('T1 = ',t1)
    
    ###Coast Period
    angular_velocity = t1*angular_acceleration
    print('Angular Velocity = ',angular_velocity)
    first_angle = final_angle*f
    print('First Angle = ',first_angle)
    coast_angle = final_angle - first_angle*2.0
    print('Coast Angle = ',coast_angle)
    tcoast = coast_angle / angular_velocity
    print('Coast Time = ',tcoast)
    
    ##Put it all together
    t2 = t1
    total_time = t1 + tcoast + t2
    print('Total Time = ',total_time)
    total_time_deriv = 2*t1 + final_angle*(1-2*f)/(angular_acceleration*t1)
    print('Total Time (Deriv) = ',total_time_deriv)
    total_accel_time = t1 + t2
    print('Total Accel Time = ',total_accel_time)
    energy_storage = power*total_accel_time/3600.0 ##watt hours
    print('Energy Storage = ',energy_storage)
    energy_storage_derivation = 2*power*np.sqrt(2*final_angle*f/angular_acceleration)/3600.0
    print('Energy Storage (Deriv) = ',energy_storage_derivation)