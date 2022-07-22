import numpy as np
import math
import matplotlib.pyplot as plt
import scipy.integrate as I

#CONSTANTS
m_r = 399.19 #mass of rocket, kg
r_p = 3396200 #radius of planet, m
m_p = 6.41712044416609e23 #mass of planet, kg
w_p = 7.0882359e-5 #spin rate of plannet
G = 6.67408e-11 #gravitational constant, m^3/(kg*s^2)
mu = G*(m_p+m_r) #gravitational parameter
LA_vec = np.linspace(0.1,np.pi/2,100) #launch angle, rad
V1 = np.zeros(len(LA_vec))
V2 = np.zeros(len(LA_vec))
P1 = np.zeros(len(LA_vec))
E1 = np.zeros(len(LA_vec))

#OPTIMIZE V1
i = 0
for launch_ang in LA_vec:
	r_a = 0.0
	r_d = 353.0 #km
	v = 1000.0 #starting guess
	while r_a < r_d:
		#INITIAL CONDITIONS
		v_vec = [v*np.cos(launch_ang), v*np.sin(launch_ang)+r_p*w_p]
		v_tot = np.sqrt(v_vec[0]**2+v_vec[1]**2)
		r_vec = [r_p, 0.0]

		#ANGULAR MOMENTUM
		h_vec = np.cross(r_vec,v_vec)
		#h = np.sqrt(h_vec[0]**2+h_vec[1]**2)
		h = h_vec

		#ECCENTRICITY
		e_vec = []
		e_vec.append(1.0/mu*((v_tot**2.0-mu/r_p)*r_vec[0]-np.dot(r_vec,v_vec)*v_vec[0]))
		e_vec.append(1.0/mu*((v_tot**2.0-mu/r_p)*r_vec[1]-np.dot(r_vec,v_vec)*v_vec[1]))
		e = np.sqrt(e_vec[0]**2+e_vec[1]**2)

		#PARAMETER
		p = h**2.0/mu

		#APOAPSIS
		r_a = (p/(1.0-e)-r_p)/1000.0

		#STEP UP VELOCITY
		v = v + 0.1
	#STORE EVERYTHING
	V1[i] = v - 0.1
	E1[i] = e
	P1[i] = p
	i = i + 1
	print([i,len(LA_vec)])

#OPTIMIZE V2
for i in range(len(LA_vec)):
	V2[i] = np.sqrt(mu/(r_d*1000+r_p))-(np.sqrt(mu/P1[i])*(1-E1[i]))
	print([i,len(LA_vec)])

#TOTAL V1 AND V2
VTOT = np.zeros(len(LA_vec))
for i in range(len(LA_vec)):
	VTOT[i] = abs(V1[i]) + abs(V2[i])

#PRINT MINIMUM V1 AND V2 BASED ON TOTAL
print(V1[np.argmin(VTOT)])
print(V2[np.argmin(VTOT)])
print(LA_vec[np.argmin(VTOT)]*180/np.pi)

###DELTA V's FROM OTIS###
V1_OTIS = 1990.555
V2_OTIS = 1661.2
LA_OTIS = 46.43

#PLOT EVERYTHING
plt.close('all')

plt.figure()
plt.plot(LA_vec*180/np.pi,V1,label='V1',color = 'red')
plt.plot(LA_vec*180/np.pi,V2,label='V2',color = 'blue')
plt.plot(LA_vec*180/np.pi, abs(V1)+abs(V2), label = 'Total',color = 'green')
plt.plot(LA_OTIS,V1_OTIS,'*',color = 'red')
plt.plot(LA_OTIS,V2_OTIS,'*',color = 'blue')
plt.xlabel('Launch Angle (deg)')
plt.ylabel(r'$\Delta$V (m/s)')
plt.title(r'$r_a = %.1f$ km and $e = %.1f$' % (r_d,0.0))
plt.xticks(np.arange(0, 95, 5.0))
plt.legend()
plt.grid()

plt.show()