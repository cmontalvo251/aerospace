#!/usr/bin/python

##Import toolboxes
import numpy as np
import sys


##Grab filename
if len(sys.argv) > 1:
    filename = sys.argv[1]
    print "Filename = ",filename
else:
    print "Filename not found"
    sys.exit();

##Open file name a parse through file
file_in = open(filename+'.csv','rb')

file_out = open(filename+'_processed.csv','wb')

for line in file_in:
    #Check for comment
    out_row = ''
    if line[0] != '#':
        #Split row
        row = line.split(' ')

        # % Time (s) Altitude (ft) Vertical velocity (ft/s) Vertical acceleration (ft/s²) Total velocity (ft/s)
        # Total acceleration (ft/s²) Position East of launch (ft) Position North of launch (ft) Lateral distance (ft)
        # Lateral velocity (ft/s) Lateral acceleration (ft/s²) Latitude (°) Longitude (°) Angle of attack (°) Roll rate (r/s)
        # Pitch rate (r/s) Yaw rate (r/s) Mass (oz) Longitudinal moment of inertia (lb·ft²) Rotational moment of inertia (lb·ft²)
        # % Thrust (N) Pitch moment coefficient (?) Yaw moment coefficient (?)  
        
        Time = row[0]
        z = row[1] #ft
        vertical_velocity = row[2] #ft/s
        Vertical_acceleration = row[3]
        Total_velocity = row[4] #ft/s
        Total_acceleration = row[5]
        Position_East_of_Launch = row[6] #ft
        Position_North_of_Launch = row[7] #ft
        lateral_distance = row[8] #ft
        lateral_direction = row[9] #deg
        lateral_velocity = row[10] #ft/s
        Lateral_acceleration = row[11]
        Latitude = row[12]
        Longitude = row[13]
        Gravitational_acceleration = row[14]
        alpha = row[15] #deg
        p = row[16] #rad/s
        q = row[17] #rad/s
        r = row[18] #rad/s
        Mass = row[19] #oz
        Propellant_mass = row[20]
        Longitudinal_moment_of_inertia = row[21] #lb-ft^2
        Rotational_moment_of_inertia = row[22] #lb-ft^2
        CP_location = row[23]
        CG_location = row[24]
        Stability_margin_calibers = row[25]
        Mach_number = row[26]
        Reynolds_number = row[27]
        Thrust = row[28] #Newtons
        Drag_force = row[29]
        Drag_coefficient = row[30]
        Axial_drag_coefficient = row[31]
        Friction_drag_coefficient = row[32]
        Pressure_drag_coefficient = row[33]
        Base_drag_coefficient = row[34]
        Normal_force_coefficient = row[35]
        Pitch_moment_coefficient = row[36]
        Yaw_moment_coefficient = row[37]
        Side_force_coefficient = row[38]
        Roll_moment_coefficient = row[39]
        Roll_forcing_coefficient = row[40]
        Roll_damping_coefficient = row[41]
        Pitch_damping_coefficient = row[42]
        Reference_length = row[43]
        Reference_area = row[44]
        theta = row[45] #deg
        psi = row[46] #deg
        Wind_velocity = row[47]
        Air_temperature = row[48]
        Air_pressure = row[49]
        Speed_of_sound = row[50]
        Simulation_time_step = row[51]
        Computation_time = row[52]
        c = ','
        out_row += str(Time) + c + Position_East_of_Launch + c + Position_North_of_Launch + c + lateral_distance + c + Latitude + c + Longitude + c + z + c + vertical_velocity + c + lateral_velocity + c + Total_velocity + c + theta + c + psi + c + p + c + q + c + r + c + lateral_direction + c + alpha + c + Mass + c + Longitudinal_moment_of_inertia + c + Rotational_moment_of_inertia + c + Thrust + '\n'
        #Time - 1 
        #EAST(x) - 2 
        #NORTH(y) - 3
        #lateral_distance(x) - 4
        #Latitude - 5
        #longitude - 6
        #z - 7
        #vertical_velocity - 8
        #lateral_velocity - 9
        #Total_velocity - 10
        #theta - 11
        #psi - 12
        #-p - 13
        #-q - 14
        #-r - 15
        #-lateral_direction - 16
        #-alpha - 17
        #-Mass - 18
        #-Iyy - 19
        #-Ixx - 20 
        #-Thrust - 21
        file_out.write(out_row)


file_in.close()
file_out.close()
        
        
    
    
