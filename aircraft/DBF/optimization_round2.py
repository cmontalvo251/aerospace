import numpy as np
import matplotlib.pyplot as plt

bvec = np.linspace(1,5,100)
M2vec = []
M3vec = []
MT = []
Nsensorvec = []

for b in bvec:
    #b = 5.0 #ft
    AR = 7.0
    c = b/AR
    rho = 0.00238
    S = b*c
    wf = 0.4*b
    SF = wf**2
    CD0_AF = 0.05
    CD0 = CD0_AF + SF/S*0.5
    CDA2 = 1.3
    AoA = 10.0*np.pi/180.0
    CL0_AF = 0.4
    CLA_AF = 5.6
    CLA = CLA_AF/(1+CLA_AF/(np.pi*AR))
    CL = CL0_AF + CLA*AoA
    K = np.pi/AR
    CD = (CD0+K*CL**2)
    ##Assum a weight
    Weight_max = b*2.0
    for x in range(0,100):
        Tmax = (CD/CL)*Weight_max
        #T = D = 0.5*rho*V^2*S*CD
        Vmax = np.sqrt(2.0*Tmax/(S*CD*rho))
        #L=W=0.5*rho*V^2*S*CL
        ##Assume a weight
        Weight_max = 0.5*rho*Vmax**2*S*CL
        
    rangeTrack = 2500.0
    time_1_lap = rangeTrack/Vmax
    print('Vmax = ',Vmax)
    print('Weight = ',Weight_max)
    print('Time for 1 lap = ',time_1_lap)

    ##Mission 3
    mission_time = 10*60.0 #seconds
    Nlaps = (mission_time/time_1_lap)
    print('NLaps = ',Nlaps)

    ###Sensor stuff
    length_of_fuselage = 0.8*b
    length_of_sensor = 0.8*length_of_fuselage
    
    print('Length of Sensor = ',length_of_sensor)
    
    diameter_of_sensor = 1.0/12.0
    
    print('Diameter of Sensor = ',diameter_of_sensor)

    Nsensors = ((wf/diameter_of_sensor)**2)

    print('Number of Sensors = ',Nsensors)

    #sensor_weight_total = 0.15*Weight_max

    #sensor_weight = sensor_weight_total/Nsensors
    
    sensor_weight = 200.0/1000.0*2.2

    print('Sensor Weight = ',sensor_weight)

    M3 = Nlaps*length_of_sensor*sensor_weight

    print('M3 = ',M3)

    M2 = Nsensors/(3*time_1_lap)

    print('M2 = ',M2)
    
    M2vec.append(M2)
    M3vec.append(M3)
    
    Nsensorvec.append(Nsensors)
    
    
M2vec = np.asarray(M2vec)
M3vec = np.asarray(M3vec)

plt.figure()
plt.plot(bvec,M2vec,label='M2')
plt.legend()

plt.figure()
plt.plot(bvec,M3vec,label='M3')
plt.legend()

plt.figure()
plt.plot(bvec,M2vec,label='M2')
plt.plot(bvec,M3vec,label='M3')
plt.legend()

plt.figure()
plt.plot(bvec,M2vec+M3vec)

plt.figure()
plt.plot(bvec,Nsensorvec)

#plt.figure()
#plt.plot(bvec,time_1_lap_vec)

plt.show()




