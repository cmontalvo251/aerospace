import numpy as np
import matplotlib.pyplot as plt
import sys
import scipy.integrate as I

#Global variable
mu = 0.0
REarth = 0.0

def EOMs(state,t):
    global mu,REarth
    x = state[0]
    z = state[1]
    velx = state[2]
    velz = state[3]
    
    ###Kinematics
    xdot = velx
    zdot = velz
    
    ##Dynamics
    ##Gravitational Acceleration
    rSat = np.sqrt((x**2) +(z**2))
    if rSat < REarth:
        xdbldot = 0.0
        zdbldot = 0.0
    else:
        xdbldot = -mu*x/(rSat**3)
        zdbldot = -mu*z/(rSat**3)
        
    statedot = np.asarray([xdot,zdot,xdbldot,zdbldot])
    return statedot

class Earth_Orbit():
    def __init__(self,ra,rp):
        self.Numerical = False #Default
        if (ra < rp):
            print('Apogee is less than perigee...exiting')
            sys.exit()
        #Get Constants
        self.constants()
        #Ingest apogee and perigee
        self.ra = ra*1000.0 + self.REarth #Note that ra and rp are assumed to be AGL
        self.rp = rp*1000.0 + self.REarth #and in km
        #First Compute Orbital Elements that matter to compute the analytic solution
        self.Orbital_Elements()
        #Then compute the Analytic Orbit
        self.Analytic_Orbit()
        #Let's make a circle for Earth for plotting
        nuE = np.linspace(0,2*np.pi,100)
        self.xE = self.REarth*np.cos(nuE)
        self.yE = self.REarth*np.sin(nuE)
    
    def constants(self):
        global REarth,mu
        self.G = 6.6742*10**-11; #%%Gravitational constant
        self.MEarth = 5.9736e24 #earth mass
        self.REarth = 6357000.0 #earth radius
        REarth = self.REarth #for global variable
        ##Get rid of self to set global variable for numerical simulation
        mu = self.G*self.MEarth #gravitational parameters
        self.mu = mu #but also set one for self.
        
    def Orbital_Elements(self):
        print('Apogee = ',self.ra)
        print('Perigee = ',self.rp)
        #Semi-Major Axis
        self.a = (self.ra+self.rp)/2.0
        #Eccentricity
        self.ecc = (self.ra - self.rp)/(self.ra+self.rp)
        print('Eccentricity = ',self.ecc)
        #Parameter
        self.p = self.a*(1-self.ecc**2)
        #Period of Orbit for numerical simulation
        self.T = self.a**(3./2.)*2*np.pi/(np.sqrt(self.mu))
        print('Period of Orbit (hrs) = ',self.T/3600.0)
        #Angular Momentum
        self.h = np.sqrt(self.p*self.mu)
    
    def Analytic_Orbit(self):
        ##True Anamoly
        self.nu = np.linspace(0,2*np.pi,1000)
        #Radius from Earth
        self.r = self.p/(1.0+self.ecc*np.cos(self.nu))
        #Components in Cartesian Coordinates
        self.xp = -self.r*np.cos(self.nu)
        self.yq = self.r*np.sin(self.nu)
        #Compute Radius
        self.r = np.sqrt(self.xp**2+self.yq**2)
        #And Altitude of Earth
        self.alt = self.r - self.REarth
    
    def Numerical_Orbit(self,N):
        #Set Numerical to True for plotting
        self.Numerical = True
        ##Time vector and therfore timestep is set by the period of the orbit
        ##and the number of data points given by the user
        self.t = np.linspace(0,self.T,N)
        ##Start Simulating assuming we are at perigee
        stateinitial = np.asarray([self.rp,0.0,0.0,self.h/self.rp])
        ##Run the integrator
        stateout = I.odeint(EOMs,stateinitial,self.t)
        #Extract States
        self.xpN = -stateout[:,0]
        self.yqN = stateout[:,1]
        self.xdotN = stateout[:,2]
        self.ydotN = stateout[:,3]
        #Compute Total Velocity
        self.Velocity = np.sqrt(self.xdotN**2+self.ydotN**2)
        #Compute Radius
        self.rN = np.sqrt(self.xpN**2+self.yqN**2)
        #And Altitude of Earth
        self.altN = self.rN - self.REarth
    
    def make_plots(self,pp=None):
        ##First let's plot the 2D Orbit
        plt.figure()
        plt.plot(self.xE/1000.0,self.yE/1000.0,'b-',label='Earth')
        plt.plot(self.xp/1000.0,self.yq/1000.0,'r-',label='Analytic Orbit (Kepler)')
        if self.Numerical == True:
            plt.plot(self.xpN/1000.0,self.yqN/1000.0,'g-',label='Numerical Orbit (Newton)')
        plt.grid()
        #plt.title('2D Orbit')
        plt.xlabel('X axis (km)')
        plt.ylabel('Y axis (km)')
        plt.legend()
        if pp != None:
            pp.savefig()
        
        if self.Numerical == True:
            plt.figure()
            plt.plot(self.t,self.altN/1000.0,'b-')
            plt.grid()
            #plt.title('Altitude vs Time')
            plt.xlabel('Time (sec)')
            plt.ylabel('Altitude of Earth (km)')
            if pp != None:
                pp.savefig()

##You can import this module or just run this script and it will plot some defaults
if __name__ == '__main__':
    ##Set the perigee (lowest point) in km
    rp = 400.0
    ##and Apogee (highest point) in km
    ra = 160000.0
    ##The default here is to just make the analytic orbit which runs very quickly.
    orbit = Earth_Orbit(ra,rp)
    ##If you want the numerical simulation which essentially solves Keplers Time of
    ##Flight problem you can run self.numerical(ndatapts)
    #Increasing the number of data points increases accuracy but slows down computation
    orbit.Numerical_Orbit(1000)
    
    ##Then you can run make plots. If you don't run the numerical orbit you'll just
    #get a plot of the analytic orbit
    orbit.make_plots()
    
    #Show everything otherwise you won't see the plots
    plt.show()