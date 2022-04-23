###MODULES NEEDED######
import numpy as np
import matplotlib.pyplot as plt
import cmath as cm

#################INPUTS TO SOFTWARE########################
####EVERYTHING MUST BE IN ENGLISH UNITS###

###AIRCRAFT WEIGHT######
weight = 636600 #lbf

###MAIN WING GEOMETRIC PARAMETERS###
main_wing_S = 5500 ##ft^2
main_wing_chord = 27.3 #ft
main_wing_wingspan = 195.7 #ft

#######CONSTANTS#######
g0 = 32.2 #ft/s^2
density = 0.00238 #slugs/ft^3

###COEFFIICENTS
CL0 = 1.11
CLA = 5.7
CDA = 0.66
CMA = -1.26
CD0 = 0.102
CLADOT = 6.7
theta0 = 0.0
CMADOT = -3.2
CLQ = 5.4
CMQ = -20.8

CYBETA = -0.96
CYP = 0.0
CYR = 0.0
CLBETA = -0.221
CLP = -0.45
CLR = 0.101
CNBETA = 0.15
CNP = -0.121
CNR = -0.3

###Mass
mass = weight / g0 ##slugs

##########Inertias####
Ix = 18.2e6 #slugs-ft^2
Iy = 33.1e6 #slugs-ft^2
Iz = 49.7e6 #slugs-ft^2
Ixz = 0.97e6 #slugs-ft^2

###Trim ###
Vtrim = 279.1 #ft/s

###Dynamic Pressure###
Q = 0.5*density*Vtrim**2

####COMPUTE STABILITY DERIVATIVES
XU = -(Q*main_wing_S*(2*CD0))/(mass*Vtrim)
CXA = -CDA+CL0
XW = Q*main_wing_S*CXA/(mass*Vtrim)
XWDOT = 0
XQ = 0
CZ0 = -CL0
ZU = (2*CZ0)*Q*main_wing_S/(mass*Vtrim)
mu = Q*main_wing_S/(mass*Vtrim)
CZA = -CLA
ZW = mu*CZA
ZWDOT = 0.
muu = Q*main_wing_S*main_wing_chord/(2*mass*Vtrim)
CZQ = -CLQ
ZQ = muu*CZQ
muuy = Q*main_wing_chord*main_wing_S/(Iy*Vtrim)
CMU = 0.
MU = muuy*CMU
MW = muuy*CMA
MWDOT = 0.
muuy2 = Q*(main_wing_chord**2)*main_wing_S/(2*Iy*Vtrim)
MQ = muuy2*CMQ
YV = mu*CYBETA
mub = Q*main_wing_S*main_wing_wingspan/(2*mass*Vtrim)
YP = mub*CYP
YR = mub*CYR
muix = Q*main_wing_S*main_wing_wingspan/(Ix*Vtrim)
LV = muix * CLBETA
muix2 = Q*main_wing_S*main_wing_wingspan**2/(2*Ix*Vtrim)
LP = muix2*CLP
LR = muix2*CLR
muiz = Q*main_wing_S*main_wing_wingspan/(Iz*Vtrim)
NV = muiz*CNBETA
muiz2 = Q*main_wing_S*main_wing_wingspan**2/(2.*Iz*Vtrim)
NP = muiz2*CNP
NR = muiz2*CNR

###Create Longitudinal A-Matrix
first_row = [XU,XW,0,-g0*np.cos(theta0)]
second_row = [ZU,ZW,Vtrim,-g0*np.sin(theta0)]
third_row = [MU+MWDOT*ZU,MW+MWDOT*ZW,MQ+Vtrim*MWDOT,-MWDOT*g0*np.sin(theta0)]
fourth_row = [0.,0.,1.,0.]
ALON = [first_row,second_row,third_row,fourth_row]

##Approximation to Phugoid
Aph = [[XU,-g0],[-ZU/Vtrim,0]]
##Approximation to Short Period
Asp = [[ZW/(1-ZWDOT),(Vtrim+ZQ)/(1-ZWDOT)],[MW+(MWDOT*ZW)/(1-ZWDOT),MQ+MWDOT*(Vtrim+ZQ)/(1-ZWDOT)]]

###Create Lateral A-Matrix
first_row = [YV,YP,g0*np.cos(theta0),YR-Vtrim]
second_row = [LV,LP,0.,LR]
third_row = [0.,1.,0.,0.]
fourth_row = [NV,NP,0.,NR]
ALAT = [first_row,second_row,third_row,fourth_row]

##Approximation to Dutch Roll
a = 1.0
b = -((LP*NR+Vtrim*NV-LR*NP)/(LP+NR)+Vtrim*(LV*NP-LP*NV)/(LP+NR)**2)
c = Vtrim*(LP*NV-LV*NP)/(LP+NR)
dutch_roll1 = -b/(2*a) + 0.5*cm.sqrt(b**2-4*a*c)
dutch_roll2 = -b/(2*a) - 0.5*cm.sqrt(b**2-4*a*c)

##Approximation to Roll Mode
ix = Ixz/Ix
iz = Ixz/Iz
roll_mode = (LP + ix*NP)/(1-ix*iz)

##Approximation to Spiral
spiral_mode = (NR-LR*NV/LV)

###COMPUTE EIGENVALUES
LONeigenvalues,LONeigenvectors = np.linalg.eig(ALON)
LATeigenvalues,LATeigenvectors = np.linalg.eig(ALAT)
PHeigenvalues,PHeigenvectors = np.linalg.eig(Aph)
SPeigenvalues,PHeigenvectors = np.linalg.eig(Asp)

###PLOT AND PRINT EVERYTHING 
plt.figure()
plt.grid()
plt.xlabel('Real')
plt.ylabel('Imaginary')

print('Longitudinal Modes:')
l = 0
for i in range(0,len(LONeigenvalues)):
    print(LONeigenvalues[i])
    l = not l
    if l:
        plt.plot(np.real(LONeigenvalues[i]),np.imag(LONeigenvalues[i]),'bx',label='Longitudinal Mode '+str(i))
    else:
        plt.plot(np.real(LONeigenvalues[i]),np.imag(LONeigenvalues[i]),'bx')
print('Short Period Approx:')
for i in range(0,len(SPeigenvalues)):
    print(SPeigenvalues[i])
    l = not l
    if l:
        plt.plot(np.real(SPeigenvalues[i]),np.imag(SPeigenvalues[i]),'rx',label='Short Period Approximation '+str(i))
    else:
        plt.plot(np.real(SPeigenvalues[i]),np.imag(SPeigenvalues[i]),'rx')
print('Phugoid Approx:')
for i in range(0,len(PHeigenvalues)):
    print(PHeigenvalues[i])
    l = not l
    if l:
        plt.plot(np.real(PHeigenvalues[i]),np.imag(PHeigenvalues[i]),'gx',label='Phugoid Approximation '+str(i))
    else:
        plt.plot(np.real(PHeigenvalues[i]),np.imag(PHeigenvalues[i]),'gx')
print('Lateral Modes:')
for i in range(0,len(LATeigenvalues)):
    print(LATeigenvalues[i])
    l = not l
    if l:
        plt.plot(np.real(LATeigenvalues[i]),np.imag(LATeigenvalues[i]),'bo',label='Lateral Mode '+str(i))
    else:
        plt.plot(np.real(LATeigenvalues[i]),np.imag(LATeigenvalues[i]),'bo')
print('Dutch Roll Approx:')
print(dutch_roll1)
plt.plot(np.real(dutch_roll1),np.imag(dutch_roll1),'ro',label='Dutch Roll Approximation')
print(dutch_roll2)
plt.plot(np.real(dutch_roll2),np.imag(dutch_roll2),'ro')
print('Roll Mode Approx:')
print(roll_mode)
plt.plot(np.real(roll_mode),np.imag(roll_mode),'go',label='Roll Approximation')
print('Spiral Mode Approx:')
print(spiral_mode)
plt.plot(np.real(spiral_mode),np.imag(spiral_mode),'mo',label='Sprial Approximation')
plt.legend()
plt.show()

