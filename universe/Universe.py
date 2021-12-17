import numpy as np
import copy as C
import matplotlib.pyplot as plt
import plotting as P
import mio as fileIO
import sys

##This code will eventually simulate the entire solar system but for now I'd 
##like to go through all the chapters in Fundamentals of Astrodynamics and
##plot stuff so that I understand all the things


############## CHAPTER 1 #################

## Keplers 3 Laws
# 1 - The orbit of each planet is an ellipse
# 2 - The line joining the planet to the sun sweeps out equal areas in equal times
# 3 - The square of the period of a planet is proportional to the cube of its mean distance

#Newton's 3 Laws
#1 - Every body continues with equal velocity unless acted on by other forces
#2 - The rate of change of momentum is proportional to the force applied
#3 - Every action has an equal and opposite reaction

class JPL():
    def __init__(self,julian_day):
        #run the coordinate transformation
        self.planetsInit()
        print('Planets Initialized')
        satellites = self.computePlanetLocations(julian_day)
        print('============Satellites Created==============')
        ##Now that we've looped through the satellites it's time to create a solar system Class
        self.MilkyWay = SolarSystem(satellites,'The Solar System')

    def AnimateOrbits(self,pp,julian_day,day_skip,num_skips,pause_time):
        framenumber = 0
        print('Animating Orbits')
        ##Let's just try and plot today
        plt.close("all")

        plti = P.plottool(12,'X (AU)','Y (AU)','Skip = '+str(0))
        XTRACES = []
        YTRACES = []
        for j in range(0,num_skips):
            print('j=',j)
            #Recompute orbits? But how?
            #Ok moved some things around here we go
            self.MilkyWay.satellites = self.computePlanetLocations(julian_day+j*day_skip)
            #self.MilkyWay.Orbit()

            plt.cla()
            #Plot system in a top down view -- REally it'd be nice if we could plot the orbital plane somehow
            plt.title('Skip = '+str(j))         
            for i in range(0,self.MilkyWay.numsatellites):
                offsetx = self.MilkyWay.satellites[3].x0/self.MilkyWay.AU*0
                offsety = self.MilkyWay.satellites[3].y0/self.MilkyWay.AU*0
                x = self.MilkyWay.satellites[i].x0/self.MilkyWay.AU-offsetx
                y = self.MilkyWay.satellites[i].y0/self.MilkyWay.AU-offsety
                if j == 0:
                    print('J == 0')
                    tracex = []
                    tracey = []
                    XTRACES.append(tracex)
                    YTRACES.append(tracey)
                else:
                    XTRACES[i].append(x)
                    YTRACES[i].append(y)
                    plti.plot(XTRACES[i],YTRACES[i],color=self.MilkyWay.satellites[i].color,label=self.MilkyWay.satellites[i].name)
                #plti.plot(self.MilkyWay.satellites[i].x/self.MilkyWay.AU-offsetx,self.MilkyWay.satellites[i].y/self.MilkyWay.AU-offsety,label=self.MilkyWay.satellites[i].name,color=self.MilkyWay.satellites[i].color)
                plti.plot(x,y,marker='o',color=self.MilkyWay.satellites[i].color)
            plt.legend(loc='upper right')
            plt.grid()
            #plt.axis('equal')
            #if self.numsatellites < 7:
            plt.axis('square')
            plt.xlim([-4,4])
            plt.ylim([-4,4])
            if j == 0:
                plt.pause(1.0)
            plt.pause(pause_time)
            strnumber = str(framenumber)
            strnumber = '0'*(4-len(strnumber)) + strnumber
            filename = 'Frames/'+strnumber+'.png'
            framenumber+=1
            print(filename)
            plt.savefig(filename)

    def planetsInit(self):
        self.G = 6.67408e-11 #m3 kg-1 s-2
        #First we need to open the correction parameters
        self.correction_parameters = fileIO.dlmread('Outer_Planets_Corrections.txt',' ')
        #print correction_parameters
        file = open('Solar_System_Orbital_Elements.txt')
        ctr = 0
        self.names = ['Sun','Mercury','Venus','Earth','Mars','Jupiter','Saturn','Uranus','Neptune','Pluto']
        self.colorwheel = ['tan','orange','blue','red','orange','yellow','green','blue','grey']
        self.AU = 149597870700.0
        self.a0 = []
        self.e0 = []
        self.i0 = []
        self.L0 = []
        self.wbar0 = []
        self.OMEGA0 = []
        self.adot = []
        self.edot = []
        self.idot = []
        self.Ldot = []
        self.wbardot = []
        self.OMEGAdot = []
        for line in file:
            if len(line) > 0:
                numbers = line.split(',')
                #This means that the current row is the '0' numbers
                if ctr == 0:
                    #print(self.names[planet_number])
                    a0 = np.float(numbers[0])
                    self.a0.append(a0)
                    e0 = np.float(numbers[1])
                    self.e0.append(e0)
                    i0 = np.float(numbers[2])
                    self.i0.append(i0)
                    L0 = np.float(numbers[3])
                    self.L0.append(L0)
                    wbar0 = np.float(numbers[4])
                    self.wbar0.append(wbar0)
                    OMEGA0 = np.float(numbers[5])
                    self.OMEGA0.append(OMEGA0)
                if ctr == 1:
                    adot = np.float(numbers[0])
                    self.adot.append(adot)
                    edot = np.float(numbers[1])
                    self.edot.append(edot)
                    idot = np.float(numbers[2])
                    self.idot.append(idot)
                    Ldot = np.float(numbers[3])
                    self.Ldot.append(Ldot)
                    wbardot = np.float(numbers[4])
                    self.wbardot.append(wbardot)
                    OMEGAdot = np.float(numbers[5])                    
                    self.OMEGAdot.append(OMEGAdot)
                ctr+=1
                if ctr == 2:
                    ctr = 0

    def computePlanetLocations(self,julian_day):
        #This will import the text file Solar_System_Orbital_Elements and save the orbital elements
        #Sun is always assumed at the center of the universe
        self.Sun = Satellite(1.989e30,432169*5280./3.28,np.asarray([0,0,0]),np.asarray([0,0,0]),'Sun','yellow',0)
        satellites = [self.Sun]
        planet_number = 1
        for i in range(0,len(self.names)-2):
            name = self.names[planet_number]
            #print(name)
            this_planet = Satellite(0,432169*5280./3.28,np.asarray([0,0,0]),np.asarray([0,0,0]),name,self.colorwheel[planet_number-1],0)
            #First compute the orbital Elements for this particular Julian Day
            T = (julian_day - 2451545.)/36525.0
            this_planet.a = (self.a0[i] + T*self.adot[i])*self.AU
            this_planet.initial_enorm = self.e0[i] + T*self.edot[i] ##This is in radians
            this_planet.initial_i = self.i0[i] + T*self.idot[i]
            this_planet.initial_L = self.L0[i] + T*self.Ldot[i]
            this_planet.wbar = self.wbar0[i] + T*self.wbardot[i]
            this_planet.initial_OMEGA = self.OMEGA0[i] + T*self.OMEGAdot[i]
            #Compute the rest of the orbital elements using aprx_pos_planets.pdf
            if planet_number > 4:
                #print('Applying correction factors')
                b = self.correction_parameters[planet_number-5,0]
                c = self.correction_parameters[planet_number-5,1]
                s = self.correction_parameters[planet_number-5,2]
                f = self.correction_parameters[planet_number-5,3]
            else:
                b = 0
                c = 0
                s = 0
                f = 0
            #print(b,c,s,f)
            self.ComputeCoordinates(this_planet,T,b,c,s,f)
            satellites.append(this_planet)
            planet_number += 1
        return satellites

    def ComputeCoordinates(self,planet,T,b,c,s,f):
        #print('Computing Planet Coordinates')
        #print('i = ',planet.initial_i)
        #print('e = ',planet.initial_enorm)
        #Compute estar
        planet.initial_enormstar = planet.initial_enorm*180./np.pi
        #print('estar = ',planet.initial_enormstar)
        #Argument of the perihelion
        planet.initial_w = planet.wbar - planet.initial_OMEGA
        #Mean Anomaly
        planet.initial_M = planet.initial_L - planet.wbar + b*T**2 + c*np.cos(f*T) + s*np.sin(f*T)
        #print('M = ',planet.initial_M)
        #Need to modulus M
        while planet.initial_M > 180:
            planet.initial_M -= 360
        while planet.initial_M < -180:
            planet.initial_M += 360
        #print('M(modulus) = ',planet.initial_M)
        #Solve for Eccentric Anomaly
        #M = E - estar*sin(E)
        M = planet.initial_M
        estar = planet.initial_enormstar
        e = planet.initial_enorm
        E = M + estar*np.sin(M*np.pi/180.0)
        dM = 1
        while abs(dM) > 1e-6:
            dM = M - (E - estar*np.sin(E*np.pi/180.0))
            #print('dM = ',dM)
            dE = dM/(1.0-e*np.cos(E*np.pi/180.0))
            E += dE
        #print('E = ',E)#,'E-estar*sin(E)',E-estar*np.sin(E*np.pi/180.0)
        
        #Compute coordinate of planet in ecliptic plane of the planet
        xprime = planet.a*(np.cos(E*np.pi/180.0)-planet.initial_enorm)
        yprime = planet.a*np.sqrt(1-planet.initial_enorm**2)*np.sin(E*np.pi/180.0)
        zprime = 0.0

        #Compute the semi latus rectum
        planet.initial_semilatus = planet.a*(1-planet.initial_enorm**2)
        #print('p = ',planet.initial_semilatus)

        #Convert certain parameters to radians
        planet.initial_w *= np.pi/180.0
        planet.initial_OMEGA *= np.pi/180.0
        planet.initial_i *= np.pi/180.0

        #Compute coordinate of planet in the J2000 frame or the ecliptic plane of the sun
        w = planet.initial_w
        OMEGA = planet.initial_OMEGA
        I = planet.initial_i
        planet.x0 = (np.cos(w)*np.cos(OMEGA) - np.sin(w)*np.sin(OMEGA)*np.cos(I))*xprime + (-np.sin(w)*np.cos(OMEGA)-np.cos(w)*np.sin(OMEGA)*np.cos(I))*yprime
        planet.y0 = (np.cos(w)*np.sin(OMEGA)+np.sin(w)*np.cos(OMEGA)*np.cos(I))*xprime + (-np.sin(w)*np.sin(OMEGA)+np.cos(w)*np.cos(OMEGA)*np.cos(I))*yprime
        planet.z0 = (np.sin(w)*np.sin(I))*xprime + (np.cos(w)*np.sin(I))*yprime

        ##Convert to numpy array
        planet.initial_pos = np.asarray([planet.x0,planet.y0,planet.z0])

        #print('Current Position of Planet = ',planet.x0,planet.y0,planet.z0)

        #Compute Period of Planet
        T = planet.a**(3./2.)*2*np.pi/(np.sqrt(self.G*self.Sun.M))
        #print('Period (days) = ',T/86400.)

class UniverseParameters():
    def __init__(self):
        print('Creating Standard Celestial Bodies')
        self.Earth = Satellite(5.972e24,6357000,np.asarray([0.,0.,0.]),np.asarray([0.,0.,0.]),'Earth','blue',(92.96e6)*5280/3.28)
        self.Sun = Satellite(1.989e30,432169*5280./3.28,np.asarray([0,0,0]),np.asarray([0,0,0]),'Sun','yellow',0)
        self.Moon = Satellite(7.34767309e22,1079*5280./3.28,np.asarray([0,0,0]),np.asarray([0,0,0]),'Moon','gray',238900*5280/3.28)
        self.G = 6.67408e-11 #m3 kg-1 s-2
        print('Universe Set')

##Let's create a class so we can easily get gravitational accelerations and other properties
class Satellite():
    def __init__(self,M,r,pos,vel,name,color,toCenter):
        self.M = M
        self.r = r ##assume here that r is a scalar
        self.G = 6.67408e-11 #m3 kg-1 s-2
        self.mu = self.G*self.M
        self.name = name
        self.color = color
        self.toCenter = toCenter

        self.initial_pos = pos #Never changes. Always a constant
        self.x0 = pos[0]
        self.y0 = pos[1]
        self.z0 = pos[2]
        self.initial_vel = vel

        self.tracex = []
        self.tracey = []
        
        self.nominal_pos = pos #Used for rk4 step
        self.nominal_vel = vel

        self.current_pos = pos #Used for rk4 step
        self.current_vel = vel
        self.current_accel = np.asarray([0.,0.,0.]) ##Can't compute this until we know what other celestial bodies are in the universe

        self.k1pos = []
        self.k1vel = []
        self.k2pos = []
        self.k2vel = []
        self.k3pos = []
        self.k3vel = []
        self.k4pos = []
        self.k4vel = []


        ##Computed Quantities
        self.xyz = []
        self.xyzdot = []
        self.h = []
        self.hnorm = []
        self.vnorm = []
        self.rnorm = []
        self.phi = []

        #print('Created Satellite = ',name)

    def Gravity(self,r):
        ##Assume that r is a vector in IJK frame
        rnorm = np.linalg.norm(r)
        if abs(rnorm) < 1e-2:
            accel = np.asarray([0,0,0])
        else:
            accel = -(self.G*self.M/rnorm**3)*r ##This is just the acceleration
        #in order to get force you need to multiply my mass of the satellite
        return accel
    
    def CircularVelocity(self,altitude):
        #Compute the velocity required to achieve circular orbit at an altitude
        vel_circ = np.sqrt(self.mu/(self.r+altitude))
        print('Computing Orbit around ',self.name)
        print('Orbit Velocity (m/s) = ',vel_circ,' at an altitude of (m)',altitude)
        #Compute Orbital Period for kicks
        self.CircularPeriod(altitude)
        return vel_circ
    
    def CircularPeriod(self,altitude):
        #The square of the period of a planet is proportional to the cube of its mean distance
        period = 2*np.pi*np.sqrt(((self.r+altitude)**3)/self.mu)
        print('Orbital Period (sec) = ',period,'at an altitude of (m)',altitude)

    def plotgravity(self,pp):
        r = np.linspace(self.r/2.0,3*self.r,100)
        a = []
        for ri in r:
            rvec = np.asarray([ri,0,0])
            accel = self.Gravity(rvec)
            anorm = np.linalg.norm(accel)
            a.append(anorm)
        plti = P.plottool(12,'Distance from Center of Planet (radii)','Gravitational Acceleration (m/s^2)',self.name)
        plti.plot(r/self.r,a)
        pp.savefig()
        
class SolarSystem():
    def __init__(self,satellites,name):
        self.AU = 149597870700.0
        self.satellites = satellites
        self.numsatellites = len(satellites)
        self.name = name
        print('Created Solar System = ',self.name)
        print('Number of Satellites = ',self.numsatellites)
        print('Satellite Names:')
        for i in range(0,self.numsatellites):
            print(self.satellites[i].name)

        #Assume the first satellite is the central planet (Sun/Earth etc and then compute orbital elements of all others)
        #I wonder if we should compute have mu be a parameter of the satellite itself. But maybe this is ok?
        self.ComputeOrbitalElements()

    def Derivatives(self):
        ##Need to get acceleration from all other satellites except itself
        #print('Derivatives')
        for i in range(0,self.numsatellites):
            self.satellites[i].current_accel = np.asarray([0.,0.,0.])
            #print(self.satellites[i].name)
            #print('pos = ',self.satellites[i].current_pos)
            #print('vel = ',self.satellites[i].current_vel)
            for j in range(0,self.numsatellites):
                dr = self.satellites[i].current_pos - self.satellites[j].current_pos
                #print('dr=',dr)
                self.satellites[i].current_accel += self.satellites[j].Gravity(dr)
            #print('accel = ',self.satellites[i].current_accel)
            if i == 0 and self.fixed == True:
                self.satellites[i].current_accel = np.asarray([0.,0.,0.])
    def Simulate(self,tfinal,timestep,tnext,fixed):
        print('Simulating System = ',self.name)
        self.fixed = fixed
        #This will kick off an RK4 routine to the position of all the planets and satellites in the simulation
        t = 0.
        tthresh = 0.
        self.time = []
        while t <= tfinal:
            ##Save Current Positions and Velocities
            for i in range(0,self.numsatellites):
                self.satellites[i].nominal_pos = self.satellites[i].current_pos
                self.satellites[i].nominal_vel = self.satellites[i].current_vel
                #print(self.satellites[i].name,self.satellites[i].nominal_pos)
                #Print to home for user
                if t >= tthresh:
                    self.satellites[i].xyz.append(self.satellites[i].nominal_pos)
                    rtnorm = np.linalg.norm(self.satellites[i].nominal_pos)
                    self.satellites[i].rnorm.append(rtnorm)
                    
                    self.satellites[i].xyzdot.append(self.satellites[i].nominal_vel)
                    vtnorm = np.linalg.norm(self.satellites[i].nominal_vel)
                    self.satellites[i].vnorm.append(vtnorm)
                    
                    ht = np.cross(self.satellites[i].nominal_pos,self.satellites[i].nominal_vel)
                    htnorm = np.linalg.norm(ht)
                    self.satellites[i].h.append(ht)
                    self.satellites[i].hnorm.append(htnorm)

                    self.satellites[i].phi.append(np.arccos(htnorm/(rtnorm*vtnorm)))

            if t>=tthresh:
                #print('Simulation Time = ',round(t))
                self.time.append(t)
                tthresh += tnext

            ##Call Derivatives 4 times
            self.Derivatives()

            ##Save First Step
            for i in range(0,self.numsatellites):
                self.satellites[i].k1pos = self.satellites[i].current_vel
                self.satellites[i].k1vel = self.satellites[i].current_accel

            #Step State
            for i in range(0,self.numsatellites):
                self.satellites[i].current_pos = self.satellites[i].nominal_pos + (timestep/2.0)*self.satellites[i].k1pos
                self.satellites[i].current_vel = self.satellites[i].nominal_vel + (timestep/2.0)*self.satellites[i].k1vel            
                
            #Call Derivatives Again
            self.Derivatives()
            
            #Save Second Step
            for i in range(0,self.numsatellites):
                self.satellites[i].k2pos = self.satellites[i].current_vel
                self.satellites[i].k2vel = self.satellites[i].current_accel
                
            #Step State
            for i in range(0,self.numsatellites):
                self.satellites[i].current_pos = self.satellites[i].nominal_pos + (timestep/2.0)*self.satellites[i].k2pos
                self.satellites[i].current_vel = self.satellites[i].nominal_vel + (timestep/2.0)*self.satellites[i].k2vel                        
                          

            #Call Derivatives a third time
            self.Derivatives()
            
            #Save Third Step
            for i in range(0,self.numsatellites):
                self.satellites[i].k3pos = self.satellites[i].current_vel
                self.satellites[i].k3vel = self.satellites[i].current_accel
                
            #Step State
            for i in range(0,self.numsatellites):
                self.satellites[i].current_pos = self.satellites[i].nominal_pos + (timestep)*self.satellites[i].k3pos
                self.satellites[i].current_vel = self.satellites[i].nominal_vel + (timestep)*self.satellites[i].k3vel                                                            
                    
            #Call Derivatives a third time
            self.Derivatives()
            
            #Save Second Step
            for i in range(0,self.numsatellites):
                self.satellites[i].k4pos = self.satellites[i].current_vel
                self.satellites[i].k4vel = self.satellites[i].current_accel
                
            #Step Final State
            for i in range(0,self.numsatellites):
                #print(self.satellites[i].k1pos,self.satellites[i].k2pos)
                rk4pos = (1./6.)*(self.satellites[i].k1pos + 2.0*self.satellites[i].k2pos + 2.0*self.satellites[i].k3pos + self.satellites[i].k4pos)
                rk4vel = (1./6.)*(self.satellites[i].k1vel + 2.0*self.satellites[i].k2vel + 2.0*self.satellites[i].k3vel + self.satellites[i].k4vel)
                self.satellites[i].current_pos = self.satellites[i].nominal_pos + rk4pos*timestep
                self.satellites[i].current_vel = self.satellites[i].nominal_vel + rk4vel*timestep           
            
            ##Step Time
            t+=timestep

        #Convert Everything to numpy arrays?
        for i in range(0,self.numsatellites):
            self.satellites[i].xyz = np.asarray(self.satellites[i].xyz)
            self.satellites[i].rnorm = np.asarray(self.satellites[i].rnorm)
            self.satellites[i].xyzdot = np.asarray(self.satellites[i].xyzdot)
            self.satellites[i].vnorm = np.asarray(self.satellites[i].vnorm)
            self.satellites[i].h = np.asarray(self.satellites[i].h)
            self.satellites[i].hnorm = np.asarray(self.satellites[i].hnorm)
            self.satellites[i].phi = np.asarray(self.satellites[i].phi)
        self.time = np.asarray(self.time)

        #Return the Satellites so the User Can manipulate them?
        #return self.satellites
        print('Simulation Complete')
        
    def PlotSystem(self,pp,zoomed):
        # draw sphere
        u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
        xsph = np.cos(u)*np.sin(v)/self.AU
        ysph = np.sin(u)*np.sin(v)/self.AU
        zsph = np.cos(v)/self.AU
        fig = plt.figure(self.name)
        ax = fig.add_subplot(111,projection='3d')
        amax = 0
        for i in range(0,self.numsatellites):
            #Plot Entire Trajectory
            x = self.satellites[i].xyz[:,0]/self.AU
            y = self.satellites[i].xyz[:,1]/self.AU
            z = self.satellites[i].xyz[:,2]/self.AU
            if np.max(x) > amax:
                amax = np.max(x)
            if np.max(y) > amax:
                amax = np.max(y)
            if np.max(z) > amax:
                amax = np.max(z)
            ax.plot(x,y,z, color = self.satellites[i].color, linestyle = 'solid',label=self.satellites[i].name)
            #Plot a sphere as the beginning of the orbit
            #ax.plot_wireframe(x[0]+self.satellites[i].r*xsph,y[0]+self.satellites[i].r*ysph,z[0]+self.satellites[i].r*zsph,color=self.satellites[i].color)            
            #Or just plot a marker instead
            #ax.scatter(x[0],y[0],z[0],c=self.satellites[i].color,marker='o',facecolor=self.satellites[i].color,edgecolor=self.satellites[i].color)
            this_color = self.satellites[i].color
            ax.scatter(x[0],y[0],z[0],c=this_color,marker='o',facecolor=this_color,edgecolor=this_color)
            #ax.scatter(0,0,0,c='r',marker='o',facecolor='r',edgecolor='r')
        plt.title(self.name)
        ax.set_xlabel('X (1000*km)')
        ax.set_ylabel('Y (1000*km)')
        ax.set_zlabel('Z (1000*km)')
        ax.set_zlim([-amax,amax])
        plt.legend()
        if zoomed >= 0:
            xmin = np.min(self.satellites[1].xyz[:,0])/self.AU
            xmax = np.max(self.satellites[1].xyz[:,0])/self.AU
            ymin = np.min(self.satellites[1].xyz[:,1])/self.AU
            ymax = np.max(self.satellites[1].xyz[:,1])/self.AU
            zmin = np.min(self.satellites[1].xyz[:,2])/self.AU
            zmax = np.max(self.satellites[1].xyz[:,2])/self.AU
            ax.set_xlim([xmin,xmax])
            ax.set_ylim([ymin,ymax])
            ax.set_zlim([zmin,zmax])
        pp.savefig()

        #Plot system in a top down view -- REally it'd be nice if we could plot the orbital plane somehow
        plti = P.plottool(12,'X (1000*km)','Y (1000*km)','Orbital Plane View')
        for i in range(0,self.numsatellites):
            plti.plot(self.satellites[i].xyz[:,0]/self.AU,self.satellites[i].xyz[:,1]/self.AU,label=self.satellites[i].name,color=self.satellites[i].color,marker='o')
        plt.legend()
        #plt.axis('equal')
        plt.xlim([-20,20])
        plt.ylim([-20,20])
        pp.savefig()

    def PlotOrbit(self,pp,zoomed):
        # draw sphere
        u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
        xsph = np.cos(u)*np.sin(v)/self.AU
        ysph = np.sin(u)*np.sin(v)/self.AU
        zsph = np.cos(v)/self.AU
        fig = plt.figure(self.name)
        ax = fig.add_subplot(111,projection='3d')
        amax = 0
        for i in range(0,self.numsatellites):
            print(self.satellites[i].name)
            #Plot Entire Trajectory
            x = self.satellites[i].x/self.AU
            #print(x)
            y = self.satellites[i].y/self.AU
            z = self.satellites[i].z/self.AU
            if np.max(x) > amax:
                amax = np.max(x)
            if np.max(y) > amax:
                amax = np.max(y)
            if np.max(z) > amax:
                amax = np.max(z)
            ax.plot(x,y,z, color = self.satellites[i].color, linestyle = 'solid',label=self.satellites[i].name)
            #Plot a sphere as the beginning of the orbit
            #ax.plot_wireframe(x[0]+self.satellites[i].r*xsph,y[0]+self.satellites[i].r*ysph,z[0]+self.satellites[i].r*zsph,color=self.satellites[i].color)            
            #Or just plot a marker instead
            #ax.scatter(x[0],y[0],z[0],c=self.satellites[i].color,marker='o',facecolor=self.satellites[i].color,edgecolor=self.satellites[i].color)
            this_color = self.satellites[i].color
            ax.scatter(self.satellites[i].x0/self.AU,self.satellites[i].y0/self.AU,self.satellites[i].z0/self.AU,c=this_color,marker='o',facecolor=this_color,edgecolor=this_color)
            
        plt.title(self.name)
        ax.set_xlabel('X (AU)')
        ax.set_ylabel('Y (AU)')
        ax.set_zlabel('Z (AU)')
        ax.set_zlim([-amax,amax])
        plt.legend()
        if zoomed >= 0:
            xmin = np.min(self.satellites[i].x/self.AU)
            xmax = np.max(self.satellites[i].x/self.AU)
            ymin = np.min(self.satellites[i].y/self.AU)
            ymax = np.max(self.satellites[i].y/self.AU)
            zmin = np.min(self.satellites[i].z/self.AU)
            zmax = np.max(self.satellites[i].z/self.AU)           
            ax.set_xlim([xmin,xmax])
            ax.set_ylim([ymin,ymax])
            ax.set_zlim([zmin,zmax])
        pp.savefig()

        #Plot system in a top down view -- REally it'd be nice if we could plot the orbital plane somehow
        plti = P.plottool(12,'X (AU)','Y (AU)','Orbital Plane View')
        for i in range(0,self.numsatellites):
            plti.plot(self.satellites[i].x/self.AU,self.satellites[i].y/self.AU,label=self.satellites[i].name,color=self.satellites[i].color)
            plti.plot(self.satellites[i].x0/self.AU,self.satellites[i].y0/self.AU,marker='o',color=self.satellites[i].color)
        plt.legend()
        #plt.axis('equal')
        #if self.numsatellites < 7:
        plt.axis('square')
            #plt.xlim([-2.4,2.4])
            #plt.ylim([-2,2])
        pp.savefig()

    def PlotMayavi(self):
        from mayavi import mlab
        import webcolors as WC
        #create a sphere
        u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
        xsph = np.cos(u)*np.sin(v)
        ysph = np.sin(u)*np.sin(v)
        zsph = np.cos(v)

        #fig = plt.figure(self.name)
        #ax = fig.add_subplot(111,projection='3d')

        f = self.AU
        
        amax = 0
        print('-----------MAYAVI--------------')
        for i in range(0,self.numsatellites):
        #for i in range(0,2):
            print(self.satellites[i].name)
            #Plot Entire Trajectory
            x = self.satellites[i].x/f
            #print(x)
            y = self.satellites[i].y/f
            z = self.satellites[i].z/f
            if np.max(x) > amax:
                amax = np.max(x)
            if np.max(y) > amax:
                amax = np.max(y)
            if np.max(z) > amax:
                amax = np.max(z)

            #Plot orbit in Mayavi
            print(self.satellites[i].color)
            colRGB = WC.name_to_rgb(self.satellites[i].color)
            print(colRGB)
            col_tuple = (colRGB.red/255,colRGB.green/255,colRGB.blue/255)
            print(col_tuple)
            s=mlab.plot3d(x,y,z,color=col_tuple)
            #ax.plot(x,y,z, color = self.satellites[i].color, linestyle = 'solid',label=self.satellites[i].name)
            
            print(self.satellites[i].x0,self.satellites[i].y0,self.satellites[i].z0)

            ##Plot the current location of the planet
            r = self.satellites[i].r/self.AU*100
            xx = (r*xsph)+self.satellites[i].x0/f
            yy = (r*ysph)+self.satellites[i].y0/f
            zz = (r*zsph)+self.satellites[i].z0/f
            s=mlab.mesh(xx,yy,zz,color=col_tuple)
            print(np.max(np.max(xx)))
            #ax.plot_wireframe(xx,yy,zz,color=self.satellites[i].color)            

            dist = np.sqrt(self.satellites[i].x0**2 + self.satellites[i].y0**2 + self.satellites[i].z0**2)/f
            print(dist)

        #plt.show()
        mlab.view(azimuth=-0, elevation=45, distance=2*dist)
        mlab.orientation_axes()
        mlab.show()

    def Orbit(self):
        print('Computing Orbit based on Orbital Elements')
        #Alright so here we're going to make a vector for the true anomaly
        nu = np.linspace(0,2*np.pi,1000)
        for i in range(0,self.numsatellites):
            print(self.satellites[i].name)
            if i == 0:
                self.satellites[i].xecl = 0*nu
                self.satellites[i].yecl = 0*nu
                self.satellites[i].zecl = 0*nu
                self.satellites[i].x = 0*nu
                self.satellites[i].y = 0*nu
                self.satellites[i].z = 0*nu
            else:
                self.satellites[i].rorbit = self.satellites[i].initial_semilatus/(1.0 + self.satellites[i].initial_enorm*np.cos(nu))
                #self.satellites[i].xecl = self.satellites[i].rorbit*np.cos(nu)
                #self.satellites[i].yecl = self.satellites[i].rorbit*np.sin(nu)  ##Perhaps this is wrong?
                self.satellites[i].xecl = self.satellites[i].a*(np.cos(nu)-self.satellites[i].initial_enorm)
                self.satellites[i].yecl = self.satellites[i].a*np.sqrt(1-self.satellites[i].initial_enorm**2)*np.sin(nu)
                self.satellites[i].zecl = 0.0*nu
                OMEGA = self.satellites[i].initial_OMEGA
                I = self.satellites[i].initial_i
                w = self.satellites[i].initial_w
                ##So I think this is wrong here. Need to fix this.
                ##So go back and look at the computation of the orbit based on the orbital elements using aprx*.pdf
                #because I think something is definitely wrong here.
                #xprime = self.satellites[i].initial_semilatus
                xprime = self.satellites[i].xecl
                yprime = self.satellites[i].yecl                
                self.satellites[i].x = (np.cos(w)*np.cos(OMEGA) - np.sin(w)*np.sin(OMEGA)*np.cos(I))*xprime + (-np.sin(w)*np.cos(OMEGA)-np.cos(w)*np.sin(OMEGA)*np.cos(I))*yprime
                self.satellites[i].y = (np.cos(w)*np.sin(OMEGA) + np.sin(w)*np.cos(OMEGA)*np.cos(I))*xprime + (-np.sin(w)*np.sin(OMEGA)+np.cos(w)*np.cos(OMEGA)*np.cos(I))*yprime
                self.satellites[i].z = (np.sin(w)*np.sin(I))*xprime + (np.cos(w)*np.sin(I))*yprime
                #self.satellites[i].x = (np.cos(nu)*np.cos(OMEGA) - np.sin(nu)*np.sin(OMEGA)*np.cos(I))*xprime
                #self.satellites[i].y = (np.cos(nu)*np.sin(OMEGA) + np.sin(nu)*np.cos(OMEGA)*np.cos(I))*xprime
                #self.satellites[i].z = (np.sin(nu)*np.sin(I))*xprime

    def PlotMechanicalEnergy(self,pp):
        ##This assumes you have a two body system. The 0 body is the central body and 1 is the orbiting satellite
        plti = P.plottool(12,'Time (sec)','Mechanical Energy (J)','Orbiting Satellite Energy')
        orbiting_satellite = self.satellites[1]
        center_satellite = self.satellites[0]
        ###Equation from Chapter 1
        mech_energy = orbiting_satellite.vnorm**2/2.0 - center_satellite.mu/orbiting_satellite.rnorm
        plti.plot(self.time,mech_energy)
        #plti.get_yaxis().get_major_formatter().set_useOffset(False)
        plt.gcf().subplots_adjust(left=0.25)
        pp.savefig()

    def plotPositionVelocity(self,pp):
        plti = P.plottool(12,'Time (sec)','Position (m)','Position of Satellite')
        for i in range(0,self.numsatellites):
            plti.plot(self.time,self.satellites[i].xyz)
            plti.plot(self.time,self.satellites[i].rnorm,color='black')
        pp.savefig()
        plti = P.plottool(12,'Time (sec)','Velocity (m/s)','Velocity of Satellite')
        for i in range(0,self.numsatellites):
            plti.plot(self.time,self.satellites[i].xyzdot)
            plti.plot(self.time,self.satellites[i].vnorm,color='black')
        pp.savefig()


    def plotAngularMomentum(self,pp):
        plti = P.plottool(12,'Time (sec)','Flight-Path Angle (deg)','Flight Path Angle')
        orbiting_satellite = self.satellites[1]
        plti.plot(self.time,orbiting_satellite.phi*2.0*np.pi/180.)
        pp.savefig()

        plti2 = P.plottool(12,'Time (sec)','Angular Momentum (m^2/s)','Angular Momentum')
        plti2.plot(self.time,orbiting_satellite.hnorm)
        plt.gcf().subplots_adjust(left=0.25)
        pp.savefig()

    def ComputeOrbitalElements(self):
        #This will compute a satellites orbital elements based on position self.initial_pos and velocity self.initial_vel
        #In order to compute this you really need to have an body on the inside.
        central_satellite = self.satellites[0] #It's assume your solar system is set up this way                
        for i in range(1,self.numsatellites):
            if np.linalg.norm(self.satellites[i].initial_pos) > 1e-2 and np.linalg.norm(self.satellites[i].initial_vel) > 1e-2:               
                print('Computing Orbital Elements....')
                print(self.satellites[i].name)
                #R is already known from the setup
                print('r = ',self.satellites[i].initial_pos)
                #Thus the x0,y0,z0 coordinates are already known
                self.satellites[i].x0 = self.satellites[i].initial_pos[0]
                self.satellites[i].y0 = self.satellites[i].initial_pos[1]
                self.satellites[i].z0 = self.satellites[i].initial_pos[2]
                ##Compute rnorm
                self.satellites[i].initial_rnorm = np.linalg.norm(self.satellites[i].initial_pos)
                print('rnorm = ',self.satellites[i].initial_rnorm)
                #v is also already known
                print('v = ',self.satellites[i].initial_vel)
                #Compute norm of v
                self.satellites[i].initial_vnorm = np.linalg.norm(self.satellites[i].initial_vel)
                print('vnorm = ',self.satellites[i].initial_vnorm)
                #Compute angular momentum
                self.satellites[i].initial_h = np.cross(self.satellites[i].initial_pos,self.satellites[i].initial_vel)
                print('h = ',self.satellites[i].initial_h)
                #Compute norm of angular momentum
                self.satellites[i].initial_hnorm = np.linalg.norm(self.satellites[i].initial_h)
                print('hnorm = ',self.satellites[i].initial_hnorm)
                #Compute Eccentricity vector
                self.satellites[i].initial_e = (1.0/(central_satellite.mu))*((self.satellites[i].initial_vnorm**2-central_satellite.mu/self.satellites[i].initial_rnorm)*self.satellites[i].initial_pos-np.dot(self.satellites[i].initial_pos,self.satellites[i].initial_vel)*self.satellites[i].initial_vel)
                print('e = ',self.satellites[i].initial_e)
                #Compute Eccentricity itself
                self.satellites[i].initial_enorm = np.linalg.norm(self.satellites[i].initial_e)
                print('enorm = ',self.satellites[i].initial_enorm)
                #Compute the semi-latus rectum
                self.satellites[i].initial_semilatus = self.satellites[i].initial_hnorm**2/central_satellite.mu
                print('p = ',self.satellites[i].initial_semilatus)
                #Compute the inclination
                self.satellites[i].initial_i = np.arccos(self.satellites[i].initial_h[2]/self.satellites[i].initial_hnorm)
                print('i = ',self.satellites[i].initial_i)
                #Compute the line of nodes
                self.satellites[i].initial_n = np.asarray([-self.satellites[i].initial_h[1],self.satellites[i].initial_h[0],0])
                print('n = ',self.satellites[i].initial_n)
                #Compute norm of lines of nodes
                self.satellites[i].initial_nnorm = np.linalg.norm(self.satellites[i].initial_n)
                print('nnorm = ',self.satellites[i].initial_nnorm)
                #Compute the longitude of the ascending node OMEGA
                #Compute the argument of the periapsis
                if abs(self.satellites[i].initial_nnorm) < 1e-2:
                    self.satellites[i].initial_OMEGA = 0.0
                    self.satellites[i].initial_w = 0.0
                else:
                    self.satellites[i].initial_OMEGA = np.arccos(self.satellites[i].initial_n[0]/self.satellites[i].initial_nnorm)
                    self.satellites[i].initial_w = np.arccos(np.dot(self.satellites[i].initial_n,self.satellites[i].initial_e)/(self.satellites[i].initial_nnorm*self.satellites[i].initial_enorm))
                print('OMEGA = ',self.satellites[i].initial_OMEGA)
                print('w = ',self.satellites[i].initial_w)
                #Compute the true anomaly
                if abs(self.satellites[i].initial_enorm) < 1e-2:
                    self.satellites[i].initial_v0 = 0.0
                else:
                    self.satellites[i].initial_v0 = np.arccos(np.dot(self.satellites[i].initial_e,self.satellites[i].initial_pos)/(self.satellites[i].initial_enorm*self.satellites[i].initial_rnorm))
                print('v0 = ',self.satellites[i].initial_v0)
