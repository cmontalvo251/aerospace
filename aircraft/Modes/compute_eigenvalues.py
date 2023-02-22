#################INPUTS TO SOFTWARE########################
####EVERYTHING MUST BE IN ENGLISH UNITS###

###AIRCRAFT WEIGHT######
weight = 3.08 #lbf

###AIRCRAFT GEOMETRY####
NoseTail = 50.0/12 #ft #Distance from nose to tail
Depth = 7.25/12 #ft    #Height from bottom of fuselage to top of wing
xcg = (2.0)/12.0 #ft center of mass from leading edge of aircraft

###MAIN WING GEOMETRIC PARAMETERS###
main_wing_wingspan = 50.0/12.0 ##feet
main_wing_chord = 10.0/12.0 ##feet
main_wing_S = 500./(12.0**2) ##ft^2
gamma = 0 #radians of dihedral
sweep = 0 #radians of wing sweep

###HORIZONTAL TAIL GEOMETRIC PARAMETERS###
horizontal_wing_wingspan = 26.0/12.0 ##ft
horizontal_wing_S = (175.0)/(12.0**2) ##ft^2
horizontal_wing_chord = horizontal_wing_S / horizontal_wing_wingspan
horizontal_xac = 25.0/12.0 #ft  #distance from main wing aerodynamic chord to tail aerodynamic chord

###VERTICAL TAIL GEOMETRIC PARAMETERS###
vertical_wing_wingspan = 13.0/12.0 ##ft
vertical_wing_S = (87.5)/(12.0**2) ##ft^2
vertical_wing_chord = vertical_wing_S / vertical_wing_wingspan

###All taken from airfoiltools.com

###MAIN WING AERODYNAMICS PARAMETERS####
###Clark Y Airfoil
Cl0w = 0.34995 ##unitless
Cd0w = 0.015342 ##unitless
Cdaw = 0.713431 ## /rad^2
Claw = 4.8872 ## /rad
Cm_acw = -0.065 ## unitless (pitch moment of the airfoil at AoA=0.0)
alfa_maxLD = 8.0 ##degrees

###HORIZONTAL WING AERODYNAMICS PARAMETERS####
###NACA 0008 
Cl0h = 0.0 ##unitless
Cd0h = 0.006 ##unitless
Cdah = 1.09145 ## /rad^2
Clah = 6.28 ## /rad
Cm_ach = 0.0 ## unitless

###MAIN WING AERODYNAMICS PARAMETERS####
###NACA 0008 
Cl0v = 0.0 ##unitless
Cd0v = 0.006 ##unitless
Cdav = 1.09145 ## /rad^2
Clav = 6.28 ## /rad
Cm_acv = 0.0 ## unitless

#############################################################

###MODULES NEEDED######
import numpy as np
import matplotlib.pyplot as plt
import cmath as cm

#######CONSTANTS#######
g0 = 32.2 #ft/s^2
density = 0.00238 #slugs/ft^3

###########################################################

##Compute Aspect Ratios
def aspect_ratio(wingspan,S):
    return wingspan**2 / S
main_wing_aspect_ratio = aspect_ratio(main_wing_wingspan,main_wing_S)
horizontal_wing_aspect_ratio = aspect_ratio(horizontal_wing_wingspan,horizontal_wing_S)
vertical_wing_aspect_ratio = aspect_ratio(vertical_wing_wingspan,vertical_wing_S)
print('Main Wing, Horizontal Stab, Vertical Stab Aspect Ratios = ',main_wing_aspect_ratio,horizontal_wing_aspect_ratio,vertical_wing_aspect_ratio)

##TAIL VOLUME COEFFICIENTS
def volume_coeff(Si,xaci,S,c):
    return Si*xaci/(S*c)
VH = volume_coeff(horizontal_wing_S,horizontal_xac,main_wing_S,main_wing_chord)
VV = volume_coeff(vertical_wing_S,horizontal_xac,main_wing_S,main_wing_chord)
zvprime = -vertical_wing_wingspan/2.0

###Aerodynamic Center
xac = main_wing_chord / 4.0

##Stability Margin
sm = xcg - xac
smbar = sm / main_wing_chord

##Convert to wing body 3D parameters
def airfoil2wing(Cla,AR):
    return Cla/(1+Cla/(np.pi*AR))
CLAw = airfoil2wing(Claw,main_wing_aspect_ratio)
CLAh = airfoil2wing(Clah,horizontal_wing_aspect_ratio)
CLAv = airfoil2wing(Clav,vertical_wing_aspect_ratio)
CMACw = Cm_acw
CMACh = Cm_ach
CMACv = Cm_acv
CL0w = Cl0w
CD0w = Cd0w
def wing_body(coeffw,coeffh,coeffv,Sh,Sv,S):
    return coeffw + coeffh * (Sh / S) + coeffv * (Sv / S)
CD0wb = wing_body(Cd0w,Cd0h,Cd0v,horizontal_wing_S,vertical_wing_S,main_wing_S)
print('CD0wb = ',CD0wb)
CDAwb = wing_body(Cdaw,Cdah,0.0,horizontal_wing_S,vertical_wing_S,main_wing_S)
CL0wb = wing_body(Cl0w,Cl0h,0.0,horizontal_wing_S,vertical_wing_S,main_wing_S)
print('CL0wb = ',CL0wb)
CLAwb = wing_body(CLAw,CLAh,0.0,horizontal_wing_S,vertical_wing_S,main_wing_S)
CM0wb = CMACw + smbar*CL0w + CMACh * (horizontal_wing_S / main_wing_S)
CMAwb = CLAw*smbar - VH*CLAh

###LONGITUDINAL COEFFICIENTS####
CLQ = 2*VH*CLAh
CMQ = -horizontal_xac/horizontal_wing_chord*CLQ

###LATERAL COEFFICIENTS###
CYBETA = -1.0*vertical_wing_S/main_wing_S*CLAv
CNBETA = VV*CLAv
CLBETA = -2.0/(3.0*np.pi)*CLAw*np.sin(gamma) - 2.0/(3.0*np.pi)*CL0w*np.sin(2*sweep) - zvprime/main_wing_wingspan*CYBETA
CYR = 2.0*VV*CLAv
CLR = CL0w/4.0 - 2.0*zvprime/main_wing_wingspan*VV*CLAv
CNR = -CD0w/4.0 - 2.0*horizontal_xac/main_wing_wingspan*VV*CLAv
CYP = -8.0/(3.0*np.pi)*vertical_wing_wingspan*vertical_wing_S/(main_wing_wingspan*main_wing_S)*CLAv
CLP = -CLAw/6.0
CNP = -horizontal_xac/main_wing_wingspan*CYP - CL0w/8.0

###Mass
mass = weight / g0 ##slugs

##########Approximate Inertias###########
Ix = (mass/12.)*main_wing_wingspan**2
Iy = (mass/12.)*((NoseTail)**2+(Depth)**2)
Iz = (mass/12.)*((NoseTail)**2+main_wing_wingspan**2)

###FIND TRIM VELOCITY###
###Assume Max L/D is trim##
##First compute Lift coefficient at max LD
CL = CL0w + CLAw * alfa_maxLD * np.pi/180.0
##Then Assume Lift = Weight
Vtrim = np.sqrt(2*weight/(density*main_wing_S*CL))
print('Vtrim (ft/s) = ',Vtrim)

###Dynamic Pressure###
Q = 0.5*density*Vtrim**2

####COMPUTE STABILITY DERIVATIVES
XU = -(Q*main_wing_S*(2*CD0wb))/(mass*Vtrim)
CXA = -CDAwb #This is wrong 
XW = Q*main_wing_S*CXA/(mass*Vtrim)
XWDOT = 0
XQ = 0
CZ0 = -CL0wb
ZU = (2*CZ0)*Q*main_wing_S/(mass*Vtrim)
mu = Q*main_wing_S/(mass*Vtrim)
CZA = -CLAwb
ZW = mu*CZA
ZWDOT = 0.
muu = Q*main_wing_S*main_wing_chord/(2*mass*Vtrim)
CZQ = -CLQ
ZQ = muu*CZQ
muuy = Q*main_wing_chord*main_wing_S/(Iy*Vtrim)
CMU = 0.
MU = muuy*CMU
MW = muuy*CMAwb
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
theta0 = alfa_maxLD*np.pi/180.0
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
roll_mode = (LP + Ix*NP)/(1-Ix*Iz)

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
#print('Spiral Mode Approx:')
#print(spiral_mode)
#plt.plot(np.real(spiral_mode),np.imag(spiral_mode),'mo',label='Sprial Approximation')
plt.legend()
plt.show()

