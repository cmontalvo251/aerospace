import numpy as np
import matplotlib.pyplot as plt

Nsyringes = np.arange(0,10,1)
We = 6.0 #lbf
Wsyringes = 1.0/8.0 #lbf
W = We + Nsyringes*Wsyringes
b = 7.5
c = 1.0
S = b*c
AR = b**2/S
Cla = 4.8
alfa_LD = 10*np.pi/180.0
CLa = Cla/(1+Cla/(np.pi*AR))
CL = CLa*alfa_LD
LD = 4.0
rho = 0.00238 ##slugs/ft^3
velocity = np.sqrt(2*W/(rho*S*CL))
distance = 2500*2+2*100 ##feet
time = distance / velocity



CD = CL/LD
T = 0.5*rho*(velocity**2)*S*CD
plt.figure()
plt.plot(Nsyringes,T)

Tmax = 11.0
m = W/32.2
accel = Tmax/m
time_to_accel = velocity / accel
delta_time = 6*time_to_accel
distance_to_accel = 0.5*accel*time_to_accel**2

plt.figure()
plt.plot(Nsyringes,accel)
plt.figure()
plt.plot(Nsyringes,delta_time-delta_time[0])
plt.ylabel('Delta Time')

plt.figure()
M2bar = Nsyringes / (time[0]+delta_time)
plt.plot(Nsyringes,M2bar)
plt.ylabel('M2bar')

plt.figure()
plt.plot(Nsyringes,distance_to_accel)
plt.show()