from math import cos,sin
from numpy import pi

phi = 0.0 #rad
theta = 0.0  #rad
psi = 45.0*pi/180.0 #rad

q0 = cos(phi/2)*cos(theta/2)*cos(psi/2) + sin(phi/2)*sin(theta/2)*sin(psi/2);
q1 = sin(phi/2)*cos(theta/2)*cos(psi/2) - cos(phi/2)*sin(theta/2)*sin(psi/2);
q2= cos(phi/2)*sin(theta/2)*cos(psi/2) + sin(phi/2)*cos(theta/2)*sin(psi/2);
q3 = cos(phi/2)*cos(theta/2)*sin(psi/2) - sin(phi/2)*sin(theta/2)*cos(psi/2);

print(str(q0) + ' !Quaternion 0')
print(str(q1) + ' !Quaternion 1')
print(str(q2) + ' !Quaternion 2')
print(str(q3) + ' !Quaternion 3')
