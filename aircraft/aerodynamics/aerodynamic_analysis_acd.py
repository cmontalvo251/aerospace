import numpy as np
import matplotlib.pyplot as plt

###This will be an example NACA 2412 aerodynamic analysis
##of a fixed wing aircraft

##Density and Viscosity of air (these are constants)
rho = 1.225 #kg/m^3
mu = 1.81*10**-5 #SI

##Conceptual Estimates of our aircraft (You must put these in depending on the size of aircraft you are designing)
c = 0.25 #m
b = 1.0
S = b*c
AR = b**2/S
Vest = 10.0 #m/s - best guess right now
W = 1.6*9.81 #weight in Newtons for a 1.6 kg aircraft - This will need to be changed based on your weight estimate

#compute the reynolds number
Re = rho * Vest * c / mu
print('Estimated Reynolds Number = ',np.round(Re))

##Read the txt file data from airfoil tools into memory (go to airfoiltools.com, find an airfoil and click "Details"
#for the Reynolds number you find from above)
data = np.loadtxt('naca_2412.txt')
#The text file you copy from the website must be in the same folder as this python script
#The first column is AoA in degrees
alpha_deg = data[:,0]
#The second and third columns are the sectional lift and drag coefficients.
Cl = data[:,1]
Cd = data[:,2]
#L/D is just Cl/Cd
LD = Cl/Cd

##Using the Cl vs AoA graph you can compute Cl0 and Cla
Cl0 = 0.285 #this is from the graph
Cla = (0.8 - 0.285) / (5*np.pi/180.0) #This is just a slope formula - m = Rise / Run
print('Cl0 , Cla = ',Cl0,Cla)

##This is your linear approximation to lift
alpha_rad = alpha_deg*np.pi/180.0
Cl_linear = Cl0 + Cla*alpha_rad

##Generate 3 plots
plt.figure()
plt.plot(alpha_deg,Cl,'b*')
plt.plot(alpha_deg,Cl_linear,'r-')
plt.grid()
plt.xlabel('AoA (deg)')
plt.ylabel('Cl')

plt.figure()
plt.plot(alpha_deg,Cd,'b*')
plt.grid()
plt.xlabel('AoA (deg)')
plt.ylabel('Cd')

plt.figure()
plt.plot(alpha_deg,LD,'b*')
plt.grid()
plt.xlabel('AoA (deg)')
plt.ylabel('L/D')

##Compute Vstall and VL/D
Clmax = np.max(Cl) #This will grab the maximum lift coefficient from the data obtained from airfoiltools
CLmax = Clmax / (1 + Clmax/(np.pi*AR)) #This is the aspect ratio equation that takes wingtip vortices into account
print('Clmax , CLmax = ',Clmax,CLmax)
Vstall = np.sqrt(2*W/(rho*S*CLmax))
print('Vstall = ',Vstall)

##Using the L/D vs AoA graph you can find the Angle of attack of maximum L/D
alpha_rad_LDmax = 6.05*np.pi/180.0
#Then use that AoA to compute the lift coefficient at that AoA
ClmaxLD = Cl0 + Cla*alpha_rad_LDmax
#and then convert that to 3D lift using the aspect ratio equation
CLmaxLD = ClmaxLD / (1 + ClmaxLD/(np.pi*AR))
print('ClmaxLD , CLmaxLD = ',ClmaxLD,CLmaxLD)
VLD = np.sqrt(2*W/(rho*S*CLmaxLD))
print('V max L/D = ',VLD)

plt.show()