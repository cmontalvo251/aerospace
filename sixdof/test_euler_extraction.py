import numpy as np
import random as rand
import sixdof as dof

#%%%%Let's make some random roll pitch and yaw values
phi = rand.uniform(-np.pi,np.pi)
theta = rand.uniform(-np.pi,np.pi)
psi = rand.uniform(-np.pi,np.pi)
phi = 30.0*np.pi/180.0
theta = 60.0*np.pi/180.0
psi = -20*np.pi/180.0
ptp = np.asarray([phi,theta,psi])
print('ptp = ',ptp)
ptp_deg = ptp*180.0/np.pi
print('ptp (deg) = ',ptp_deg)

# %%%Create the Rotation Matrix
TIB = dof.R123(phi,theta,psi);
TBI = np.transpose(TIB);

ptp_ = dof.extract_Euler(TBI);
print('ptp extracted = ',ptp_)

# %%%Now let's convert the Euler Angles to Quaternions
q0123_e2q = dof.euler2quat(ptp_)
print('q0123 converted = ',q0123_e2q)

# %%%Let's test the norm component
# n = norm(q0123_e2q);

# %%%Now let's make the transformation matrix using quaternions
TIB_quat = dof.RQUAT(q0123_e2q);
TBI_quat = np.transpose(TIB_quat)

# %%%Let's extract the quaternion
# %%I guess theoretically I could just extract the Euler Angles and
# %then convert to quaternions?
# [phi_q,theta_q,psi_q] = extract_Euler(TBI_quat);
# q0123_extract_Euler = euler2quat([phi_q,theta_q,psi_q])
# %%But that just seems roundabout.

# %%This method below will give you one of the possible quaternion solutions
# %%The reason is because Euler angles are non-unique 
# %%For example take psi = pi and psi = -pi
# %They are both technically the same direction 
# %There is a discontinuity in the Euler angles causing the issue
# %here in 4Dimensional space. 
# %%What you are doing is mapping 4D space into 3D space and
# %unfortunately you lose information
q0123_e = dof.extract_quaternion(TBI_quat)

print('quat extracted = ',q0123_e)

ptp_e = dof.quat2euler(q0123_e)
print('ptp converted extracted = ',ptp_e)
ptp_e_deg = ptp_e*180.0/np.pi
print('ptp converted extracted (deg) = ',ptp_e_deg)

