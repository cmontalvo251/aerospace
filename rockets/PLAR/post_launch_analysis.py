import matplotlib.pyplot as plt
import numpy as np

def truncate(outvec,timevec,ts,te):
    outvec = outvec[timevec > ts]
    timevec_ = timevec[timevec > ts]
    outvec = outvec[timevec_ < te]
    return outvec

def filter(vec,tau):
    outvec = 0*vec
    outvec[0] = vec[0]
    for x in range(1,len(vec)-1):
        outvec[x] = (1-tau)*outvec[x-1] + tau*vec[x]
    return outvec

#1.) Open the file and read in all accel/gyro data streams
files = ['launch_data.txt','launch_data.txt']
Ax_avg = []
Ay_avg = []
Az_avg = []
gx_avg = []
gy_avg = []
gz_avg = []
for file in files:
    data = np.loadtxt(file)
    time = data[:,0]
    Ax = data[:,1]
    Ay = data[:,2]
    Az = data[:,3]
    A = np.sqrt(Ax**2 + Ay**2 + Az**2)
    gx = data[:,4]
    gy = data[:,5]
    gz = data[:,6]
    plt.figure()
    plt.grid()
    plt.xlabel('Time (sec)')
    plt.ylabel('Acceleration (m/s^2)')
    #plt.plot(time,A)
    #1.5) Truncate data to be within launch window
    time_events = time[A>35]
    A_events = A[A>35]
    #plt.plot(time_events,A_events,'b*')
    tlaunch = time_events[0]
    print('Time of Launch = ',tlaunch)
    tland = time_events[-1]
    print('Time of Landing = ',tland)
    Ax_t = truncate(Ax,time,tlaunch-1,tland+1)
    Ay_t = truncate(Ay,time,tlaunch-1,tland+1)
    Az_t = truncate(Az,time,tlaunch-1,tland+1)
    gx_t = truncate(gx,time,tlaunch-1,tland+1)
    gy_t = truncate(gy,time,tlaunch-1,tland+1)
    gz_t = truncate(gz,time,tlaunch-1,tland+1)
    time_t = truncate(time,time,tlaunch-1,tland+1)
    plt.plot(time_t,Ax_t,label='Ax_t')
    plt.plot(time_t,Ay_t,label='Ay_t')
    plt.plot(time_t,Az_t,label='Az_t')
    #2.) Low Pass Filter the data - remove all high frequency content (new_data = old_data*(1-s) + s*filtered_data
    tau = 0.5 #1.0 is no filter and 0.0 is overfiltering
    Ax_tf = filter(Ax_t,tau)
    Ay_tf = filter(Ay_t,tau)
    Az_tf = filter(Az_t,tau)
    gx_tf = filter(gx_t,tau)
    gy_tf = filter(gy_t,tau)
    gz_tf = filter(gz_t,tau)
    plt.plot(time_t,Ax_tf,label='Ax_tf')
    plt.plot(time_t,Ay_tf,label='Ay_tf')
    plt.plot(time_t,Az_tf,label='Az_tf')
    plt.legend()
    #3.) Translate the accel data streams to numerous center of masses of the projectile using 2D dynamics
    #Skip this one and assume the accelerometer is at the center of mass or we integrate the acceleration at the
    #Accelerometer point rather than Cg
    #3.5) Average all data streams together. (This will remove the bias)
    if len(Ax_avg) == 0:
        Ax_avg = Ax_tf
        Ay_avg = Ay_tf
        Az_avg = Az_tf
        gx_avg = gx_tf
        gy_avg = gy_tf
        gz_avg = gz_tf
        numvars = 1
    else:
        Ax_avg += Ax_tf
        Ay_avg += Ay_tf
        Az_avg += Az_tf
        gx_avg += gx_tf
        gy_avg += gy_tf
        gz_avg += gz_tf
        numvars+=1
        
##Compute the average
Ax_avg /= numvars
Ay_avg /= numvars
Az_avg /= numvars
gx_avg /= numvars
gy_avg /= numvars
gz_avg /= numvars

plt.figure()
plt.grid()
plt.xlabel('Time (sec)')
plt.ylabel('Acceleration Average (m/s^2)')
plt.plot(time_t,Ax_avg,label='Ax')
plt.plot(time_t,Ay_avg,label='Ay')
plt.plot(time_t,Az_avg,label='Az')

plt.figure()
plt.grid()
plt.xlabel('Time (sec)')
plt.ylabel('Gyro Average (rad/s)')
plt.plot(time_t,gx_avg,label='gx')
plt.plot(time_t,gy_avg,label='gy')
plt.plot(time_t,gz_avg,label='gz')

#4.) Need to run some sort of Kalman Filter (RK4 integrator + Linear Least Squared Estimator) to obtain Euler angles/quaternions for the entire duration of the flight
#Init quaternion
q0 = 0*Ax_avg
q1 = 0*Ax_avg
q2 = 0*Ax_avg
q3 = 0*Ax_avg
q0[0] = 1.0
GEARTH = 9.81
for i in range(0,len(Ax_avg)-1):
    ##Extract Ax,Ay,Az
    ax = Ax_avg[i]
    ay = Ay_avg[i]
    az = Az_avg[i]
    ##Convert to Gs
    ax /= GEARTH
    ay /= GEARTH
    az /= GEARTH

    #// Normalise accelerometer measurement
    recipNorm = np.sqrt(ax * ax + ay * ay + az * az)
    ax /= recipNorm
    ay /= recipNorm
    az /= recipNorm

    #// Estimated direction of gravity and vector perpendicular to magnetic flux
    halfvx = q1[i] * q3[i] - q0[i] * q2[i]
    halfvy = q0[i] * q1[i] + q2[i] * q3[i]
    halfvz = q0[i] * q0[i] - 0.5 + q3[i] * q3[i]

    #// Error is sum of cross product between estimated and measured direction of gravity
    halfex = (ay * halfvz - az * halfvy);
    halfey = (az * halfvx - ax * halfvz);
    halfez = (ax * halfvy - ay * halfvx);

    #// Apply proportional feedback
    gxi = gx_avg[i]
    gyi = gy_avg[i]
    gzi = gz_avg[i]
    twoKp = 2.0
    gxi += twoKp * halfex
    gyi += twoKp * halfey
    gzi += twoKp * halfez

    #// Integrate rate of change of quaternion
    elapsedTime = time_t[i+1] - time_t[i]
    gxi *= (0.5 * elapsedTime);		#// pre-multiply common factors
    gyi *= (0.5 * elapsedTime);
    gzi *= (0.5 * elapsedTime);
    qa = q0[i]
    qb = q1[i]
    qc = q2[i]
    q0[i+1] = q0[i] + (-qb * gxi - qc * gyi - q3[i] * gzi)
    q1[i+1] = q1[i] + (qa * gxi + qc * gzi - q3[i] * gyi)
    q2[i+1] = q2[i] + (qa * gyi - qb * gzi + q3[i] * gxi)
    q3[i+1] = q3[i] + (qa * gzi + qb * gyi - qc * gxi);

    #// Normalise quaternion
    Norm = np.sqrt(q0[i+1] * q0[i+1] + q1[i+1] * q1[i+1] + q2[i+1] * q2[i+1] + q3[i+1] * q3[i+1])
    q0[i+1] /= Norm
    q1[i+1] /= Norm
    q2[i+1] /= Norm
    q3[i+1] /= Norm

roll = np.arctan2(2.*(q0*q1 + q2*q3),1-2*(q1**2 + q2**2))
pitch = np.arcsin(2.*(q0*q2-q3*q1)) 
yaw = np.arctan2(2.*(q0*q3 + q1*q2),1-2.*(q2**2 + q3**2))

plt.figure()
plt.grid()
plt.xlabel('Time (sec)')
plt.ylabel('Euler Angles (rad)')
plt.plot(time_t,roll,label='Roll')
plt.plot(time_t,pitch,label='Pitch')
plt.plot(time_t,yaw,label='Yaw')
plt.legend()

#5.) Take the body frame acceleration vector and transform it to the inertial frame
Ax_I = 0*Ax_avg
Ay_I = 0*Ay_avg
Az_I = 0*Az_avg
for i in range(0,len(Ax_avg)):
    ax = Ax_avg[i]
    ay = Ay_avg[i]
    az = Az_avg[i]
    ab = np.array([ax,ay,az])
    phi = roll[i]
    theta = pitch[i]
    psi = yaw[i]
    ctheta = np.cos(theta)
    cpsi = np.cos(psi)
    sphi = np.sin(phi)
    stheta = np.sin(theta)
    cphi = np.cos(phi)
    spsi = np.sin(psi)
    ttheta = np.tan(theta)
    TIB = np.asarray([[ctheta*cpsi,sphi*stheta*cpsi-cphi*spsi,cphi*stheta*cpsi+sphi*spsi],[ctheta*spsi,sphi*stheta*spsi+cphi*cpsi,cphi*stheta*spsi-sphi*cpsi],[-stheta,sphi*ctheta,cphi*ctheta]])
    aI = np.matmul(TIB,ab)
    Ax_I[i] = aI[0]
    Ay_I[i] = aI[1]
    ##5.5) Substract Gravity inertia
    Az_I[i] = aI[2] - GEARTH
    
plt.figure()
plt.grid()
plt.xlabel('Time (sec)')
plt.ylabel('Acceleration Inertial (m/s^2)')
plt.plot(time_t,Ax_I,label='Ax')
plt.plot(time_t,Ay_I,label='Ay')
plt.plot(time_t,Az_I,label='Az')

#6.) double integrate to get position (somehow you're gonna have to remove bias)
vx = 0*Ax_I
vy = 0*Ay_I
vz = 0*Az_I

for i in range(0,len(Ax_I)-1):
    elapsedTime = time_t[i+1]-time_t[i]
    vx[i+1] = vx[i] + Ax_I[i]*elapsedTime
    vy[i+1] = vy[i] + Ay_I[i]*elapsedTime
    vz[i+1] = vz[i] + Az_I[i]*elapsedTime
    
plt.figure()
plt.grid()
plt.xlabel('Time (sec)')
plt.ylabel('Velocity Inertial (m/s)')
plt.plot(time_t,vx,label='Vx')
plt.plot(time_t,vy,label='Vy')
plt.plot(time_t,vz,label='Vz')

px = 0*Ax_I
py = 0*Ay_I
pz = 0*Az_I

for i in range(0,len(Ax_I)-1):
    elapsedTime = time_t[i+1]-time_t[i]
    px[i+1] = px[i] + vx[i]*elapsedTime
    py[i+1] = py[i] + vy[i]*elapsedTime
    pz[i+1] = pz[i] + vz[i]*elapsedTime
    
plt.figure()
plt.grid()
plt.xlabel('Time (sec)')
plt.ylabel('Position Inertial (m)')
plt.plot(time_t,px,label='x')
plt.plot(time_t,py,label='y')
plt.plot(time_t,pz,label='z')

#7.) Profit.

plt.show()