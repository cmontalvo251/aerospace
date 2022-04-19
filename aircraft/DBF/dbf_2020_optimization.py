import numpy as np
import matplotlib.pyplot as plt

###Mission 1 - Just fly
#M1 = 1.0 #Everyone gets 1 point for completing the mission

##Mission 2 - Fly with as many containers as possible
#M2 = 1 + (Ncontainers/time) / max2

##Mission 3 is to deploy during flight
#M3 = 2 + (Nlaps*sensor_length*sensor_weight)/ max3

###Ground Mission - Load payload with lots of payloads
#GM = Min_time / time

###In my opinion Mission 1 can't be optimized
###Ground mission is optimized when the number of payloads is minimized but I don't know how much time is
###saved when we increase the payload. Maybe add 2 seconds per payload?

##Anyway I think the best we can do is optimize Mission2 and Mission 3 or rather
#(Ncontainers/time) + (Nlaps*sensor_length) ##I forgot about sensor_weight
#Now of course time will change as we add more payload because we have to fly at a higher
#Angle of attack which will increase drag but let's just look at Ncontainers vs sensor_length

##So what I want to optimize is Ncontainers + sensor_length
#cost = Ncontainers + sensor_length

#Ok so for starters let's assume a volume
#volume = 10 #m^3 #We will non-dimensionalize this later

##Then Let's assume the size of the payload bay is volume*(1x1 m^2)
##However if you look at that volume equation really the cross sectional area has more to do with area of the sensor
##versus the length so let's just look at length of the payload bay then
#So let's say that the length of the payload is fixed
length_payload = 10 #meters

##Ok so now we can look at varying the sensor_length from 0 to length_payload
#But we can't divide by zero. So let's set the lower limit to 1 cm
sensor_length = np.linspace(0.01,length_payload,100)

##Once we know sensor_length we can compute number of containers that can fit in the length of the payload
Ncontainers = np.floor(length_payload/sensor_length) ##Has to be an integer number

##Ok so now we can compute cost
cost = sensor_length + Ncontainers

##Non dimensionalize first
cost_nd = cost/length_payload
sensor_nd = sensor_length/length_payload

##And then plot sensor_length on the x-axis and cost on the y-axis
#plt.plot(sensor_nd,cost_nd,'b-')
#^This plot looks kind of odd though so let's make the optimization routine a bit more sophisticated.

#Let's assume the vehicle has certain parameters
Vcruise = 20 #m/s
length_of_track = 2500.0/3.28 ##feet to meters
time_for_1_lap = length_of_track/Vcruise #seconds
print('1 Lap Time = ',time_for_1_lap)

###Mission 1 - Just fly
#M1 = 1.0 #Everyone gets 1 point for completing the mission
##Mission 2 - Fly with as many containers as possible
#M2 = 1 + (Ncontainers/time) / max2
##Mission 3 is to deploy during flight
#M3 = 2 + (Nlaps*sensor_length)/ max3
###Ground Mission - Load payload with lots of payloads
#GM = Min_time / time

#Let's try and optimize Ncontainers/time + Nlaps*sensor_length

##First Mission2
laps_Mission_2 = 3.0 ##Need to do 5 laps
time_M2 = time_for_1_lap*laps_Mission_2

#Now Mission3
max_time = 10.0*60 #you have 10 minutes
Nlaps = np.floor(max_time/time_for_1_lap)

#Ok now let's compute the new cost
cost2 = (Ncontainers/time_M2) + (Nlaps*sensor_length)

cost2_nd = cost2/length_payload

plt.figure()
plt.plot(sensor_nd,cost2_nd)
plt.grid()
plt.xlabel('Sensor Length (% of Payload Bay Length)')
plt.ylabel('Mission 2/3 Points (Non-Dimensionalized by Payload Bay Length)')

##For a final one let's make some assumptions about best times
M2_best = 10./100. ##10 containers over 100 seconds
M3_best = 20.*length_payload  #20 laps and X meter payload

cost3 = (Ncontainers/time_M2)/M2_best + (Nlaps*sensor_length)/M3_best

cost3_nd = cost3/length_payload

plt.figure()
plt.plot(sensor_nd,cost3_nd)
plt.grid()
plt.xlabel('Sensor Length (% of Payload Bay Length)')
plt.ylabel('Mission 2/3 Points (Non-Dimensionalized by Payload Bay Length & Other Teams)')
plt.show()