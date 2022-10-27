import numpy as np
import matplotlib.pyplot as plt

##CONSTANTS
G = 6.6743e-11

##Mass Properties of Center Body
mDidymos = 5.2e11
RDidymos = 390.0
muDidymos = G*mDidymos

###########PRE-IMPACT#############

###Orbital Elements of Dimorphis Pre-Impact
aPRE = 1190.
ePRE = 0.0

##Orbital Period
TPRE = 2*np.pi/np.sqrt(muDidymos)*aPRE**(3./2.)
TPRE /= 60.0 #converts to minutes
TPRE /= 60.0 #converts to hours
print('Orbital Period Pre Impact = ',TPRE,' hrs')

##Velocity of Orbit
# h = a * v
# h = sqrt(p*mu)
# p = a*(1-e^2)
pPRE = aPRE*(1-ePRE**2)
hPRE = np.sqrt(pPRE*muDidymos)
vPRE = hPRE/aPRE
print('Velocity Pre Impact = ',vPRE,' m/s')

############POST IMPACT############
TPOST = (11.0 + 23./60.)*3600.0
#TPOST = 2*np.pi/np.sqrt(muDidymos)*aPOST**(3./2.)
aPOST = (TPOST*np.sqrt(muDidymos)/(2*np.pi))**(2./3.)
print('Semi Major Axis Post = ',aPOST)
# a = (ra + rp)/2
raPOST = aPRE
rpPOST = aPOST*2.0 - raPOST
print('Apoaps = ',raPOST,' Periaps = ',rpPOST)
ePOST = (raPOST - rpPOST)/(2*aPOST)
print('Eccentricity Post Impact = ',ePOST)
pPOST = aPOST*(1-ePOST**2)
hPOST = np.sqrt(pPOST*muDidymos)
vaPOST = hPOST / raPOST
vpPOST = hPOST / rpPOST
print('Velocity Apo = ',vaPOST,' Velocity Peri = ',vpPOST)

###########PLOT EVERYTHING#########

###PREIMPACT
vu = np.linspace(0,2*np.pi,1000)
xPRE = aPRE*np.cos(vu)
yPRE = aPRE*np.sin(vu)

##POSTIMPACT
rPOST = pPOST / ( 1 + ePOST*np.cos(vu))
xPOST = rPOST*np.cos(vu)
yPOST = rPOST*np.sin(vu)

plt.plot(xPRE,yPRE,label='Pre-Impact')
plt.plot(xPOST,yPOST,label='Post-Impact')
plt.legend()
plt.grid()

plt.show()