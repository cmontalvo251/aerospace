########################################################
#   NAME OF SOFTWARE
#   WHO WROTE IT
#   LAST EDIT DATE
#   LIST OF INPUTS
#   LIST OF OUTPUTS
#   DIRECTIONS ON HOW TO RUN
#   DESCRIPTION OF SOFTWARE
#########################################################

###Modules
import numpy as np
import Orbit as O
import os
import matplotlib.pyplot as plt
import mag_field as igrf
from matplotlib.backends.backend_pdf import PdfPages

###Functions
def RiemannSum(x,t):
    out = 0.0
    for i in range(1,len(t)):
        out += x[i]*(t[i]-t[i-1])
    return out

def Solar(height,width,Surface_Area,Ixx,time,Tor):
    rad_pressure = 4.5e-6 # Pa
    Surface_Area = height*width #m^2
    solar_rad_force = Surface_Area*rad_pressure #N
    solar_rad_torque = solar_rad_force*np.max([height/2.0,width/2.0]) #Nm
    Tor = solar_rad_torque + 0*time
        
    #print(solar_rad_torque)
    plt.figure()
    plt.plot(time,Tor)
    plt.grid()
    plt.xlabel('Time (sec)')
    plt.ylabel('Solar Radiation Torque (Nm)')
    pdfhandle.savefig()
    
    ang_accel_solar = solar_rad_torque/Ixx #rad/s^2
    gamma = ang_accel_solar*time
    return Tor,gamma

def Grav(x,Mass,Grav_Coeff,Mass_Earth,Rad_Earth,height,width):
    Grav_Accel_Bot = Grav_Coeff*Mass_Earth/(Rad_Earth + x - height/2.0*np.sqrt(2.0)/2.0)**2
    Grav_Accel_Top = Grav_Coeff*Mass_Earth/(Rad_Earth + x + height/2.0*np.sqrt(2.0)/2.0)**2
    #print(Grav_Accel_Bot) 
    Grav_Force_Bot = Grav_Accel_Bot*Mass
    Grav_Force_Top = Grav_Accel_Top*Mass
    #print(Grav_Force_Bot)
    
    plt.figure()
    plt.plot(x/1000.0,Grav_Force_Bot)
    plt.grid()
    plt.xlabel('Altitude (km)')
    plt.ylabel('Gravitational Force (N)')
    pdfhandle.savefig()
    
    Del_Grav = Grav_Force_Top - Grav_Force_Bot
    Grav_Torque = abs(Del_Grav*np.max([height/2.0,width/2.0]))*np.sqrt(2.0)/2.0
    return Grav_Torque

def Aerodynamics(x,Rad_Earth,Mass_Earth,Grav_Coeff,CD,Surface_Area,height,width,time,velocity):
    beta = 0.1354 #inverse km
    rhos = 1.225 #kg/m^3
    rho = rhos*np.exp(-beta*(x/1000.0))
    
    plt.figure()
    plt.plot(time,rho)
    plt.grid()
    plt.xlabel('Time (sec)')
    plt.ylabel('Density (kg/m^3)')
    pdfhandle.savefig()
    
    Disturbance = 0.5*rho*velocity**2*Surface_Area*CD
    #print(Disturbance)
    
    plt.figure()
    plt.plot(time,Disturbance)
    plt.grid()
    plt.xlabel('Time (sec)')
    plt.ylabel('Aerodynamic Disturbances (N)')
    pdfhandle.savefig()
    
    Aero = Disturbance*np.max([height/2.0,width/2.0])/2.0
    return Aero

def Dipole(Solar_Torque,btotal_dipole,btotal_arr):
    Dipole_Moment = Solar_Torque[0]/(btotal_dipole*1e-09)
    #print(Solar_Torque[0],btotal_dipole,Dipole_Moment)
    MRD = Dipole_Moment*(btotal_arr*1e-09)
    return MRD

###Class
class CubeSat():
    def __init__(self,FS,Ixx,Iyy,Izz,wmax,length,width,height,Mission_Duration,CD,rp,ra,Mass_Sat,Mag_Moment,Mag_Moment_np,Mag_Mass_np,Chosen_Momentum,pdfhandle):
        self.FS = FS
        self.Ixx = Ixx
        self.Iyy = Iyy
        self.Izz = Izz
        self.wmax = wmax
        self.length = length
        self.width = width
        self.height = height
        self.Surface_Area = self.height*self.width
        self.Mission_Duration = Mission_Duration
        self.CD = CD
        self.rp = rp
        self.ra = ra
        self.Mass_Sat = Mass_Sat
        self.Mag_Moment = Mag_Moment
        self.Mag_Moment_np = Mag_Moment_np
        self.Mag_Mass_np = Mag_Mass_np
        self.Chosen_Momentum = Chosen_Momentum
        self.Orbit_Analysis(pdfhandle)
        self.Magnetic_Field_Model() #Using IGRF
        self.Disturbance_Torques()
        self.Reaction_Wheel_Analysis()
        self.Magnetorquer_Analysis()
        
    def Orbit_Analysis(self,pdfhandle):
        #Make sure this is working properly
        orbit = O.Earth_Orbit(self.ra,self.rp)
        orbit.Numerical_Orbit(1000)
        #Get time as a vector to use later on
        self.time = orbit.t #sec
        #print(self.time)
        self.altitude = orbit.altN #m
        self.Grav_Coeff = orbit.G
        self.Mass_Earth = orbit.MEarth #kg
        self.Rad_Earth = orbit.REarth #m
        self.nu = orbit.nu
        self.velocity = orbit.Velocity
        #Print plots of orbit
        orbit.make_plots(pdfhandle)
        
    def Magnetic_Field_Model(self):
        lat_arr = 0*self.altitude #latitude is constant
        lon_arr = self.nu*180.0/np.pi ##longitude is the same as True Anamoly for planar orbits
        date = 2021
        self.bx_arr,self.by_arr,self.bz_arr,self.btotal_arr = igrf.bfield_arr(lat_arr,lon_arr,self.altitude/1000.0,date)
        
        result = igrf.Mag_Dipole()
        self.bx_dipole = result[3]
        self.by_dipole = result[4]
        self.bz_dipole = result[5]
        self.btotal_dipole = result[6]
        print('B-field = ',self.bx_dipole,self.by_dipole,self.bz_dipole,self.btotal_dipole)
        
        #print(btotal_arr)
        plt.figure()
        plt.plot(self.altitude/1000.0,self.btotal_arr)
        plt.grid()
        plt.xlabel('Altitude (km)')
        plt.ylabel('Magnetic Field Strength (nT)')
        pdfhandle.savefig()
        
        #plt.figure()
        #plt.plot(self.nu,self.btotal_arr)
        #plt.xlabel('True Anamoly (rad)')
        #plt.ylabel('B-Field (nT)')
        #plt.grid()
        #pdfhandle.savefig()
        
        plt.figure()
        plt.plot(self.time,self.btotal_arr)
        plt.grid()
        plt.xlabel('Time (sec)')
        plt.ylabel('Magnetic Field Strength (nT)')
        pdfhandle.savefig()
        
    def Disturbance_Torques(self):
        ##Compute the disturbance torques based on orbit        
        ##solar radiation pressure
        self.Tor = 0.0
        self.Solar_Torque,self.Ang_Vel_Solar = Solar(self.height,self.width,self.Surface_Area,self.Ixx,self.time,self.Tor)
        Solar_Momentum = RiemannSum(self.Ang_Vel_Solar,self.time)
        #print(Ang_Vel_Solar)
        #plt.figure()
        #plt.plot(self.time,self.Ang_Vel_Solar)
        #plt.grid()
        #plt.xlabel('Time (sec)')
        #plt.ylabel('Solar Angular Velocity (rad/s)')
        #pdfhandle.savefig()
        
        plt.figure()
        plt.plot(self.time,self.Solar_Torque)
        plt.grid()
        plt.xlabel('Time (sec)')
        plt.ylabel('Solar Radiation Pressure Torque (Nm)')
        pdfhandle.savefig()
        
        ##aerodynamic torques
        self.Aero_Torque = Aerodynamics(self.altitude,self.Rad_Earth,self.Mass_Earth,self.Grav_Coeff,self.CD,self.Surface_Area,self.height,self.width,self.time,self.velocity)
        #print(self.Aero_Torque)
        plt.figure()
        plt.plot(self.time,self.Aero_Torque)
        plt.grid()
        plt.xlabel('Time (sec)')
        plt.ylabel('Aerodynamic Torque (Nm)')
        pdfhandle.savefig()
        
        ##gravity gradient torque
        self.Gravity_Torque = Grav(self.altitude,self.Mass_Sat,self.Grav_Coeff,self.Mass_Earth,self.Rad_Earth,self.height,self.width)
        #print(self.Gravity_Torque)
        plt.figure()
        plt.plot(self.time,self.Gravity_Torque)
        plt.grid()
        plt.xlabel('Time (sec)')
        plt.ylabel('Gravitational Gradient Torque (Nm)')
        pdfhandle.savefig()
        
        ##magnetic resonance dipole
        self.Mag_Res_Dipole = Dipole(self.Solar_Torque,self.btotal_dipole,self.btotal_arr)
        #print(self.Mag_Res_Dipole)
        plt.figure()
        plt.plot(self.time,self.Mag_Res_Dipole)
        plt.grid()
        plt.xlabel('Time (sec)')
        plt.ylabel('Magnetic Resonance Dipole Torque (Nm)')
        pdfhandle.savefig()
        
        #Total Disturbances
        self.Tot_Dis_Tor = self.Solar_Torque + self.Aero_Torque + self.Gravity_Torque + self.Mag_Res_Dipole
        self.Total_Momentum = RiemannSum(self.Tot_Dis_Tor,self.time)
        print('Total Momentum Change from Disturbances (Nms) = ',self.Total_Momentum)
        
        plt.figure()
        plt.plot(self.time,self.Solar_Torque,'b-',label='Solar Radiation Pressure')
        plt.plot(self.time,self.Aero_Torque,'r-',label='Aerodynamics')
        plt.plot(self.time,self.Gravity_Torque,'g-',label='Gravity Gradient')
        plt.plot(self.time,self.Mag_Res_Dipole,'y-',label='Magnetic Resonance Dipole')
        plt.plot(self.time,self.Tot_Dis_Tor,'k-',label='Total Disturbance Torque')
        plt.xlabel('Time (sec)')
        plt.ylabel('Disturbance Torques (Nm)')
        plt.grid()
        plt.legend()
        pdfhandle.savefig()
        
        plt.figure()
        plt.plot(self.altitude/1000.0,self.Solar_Torque,'b-',label='Solar Radiation Pressure')
        plt.plot(self.altitude/1000.0,self.Aero_Torque,'r-',label='Aerodynamics')
        plt.plot(self.altitude/1000.0,self.Gravity_Torque,'g-',label='Gravity Gradient')
        plt.plot(self.altitude/1000.0,self.Mag_Res_Dipole,'y-',label='Magnetic Resonance Dipole')
        plt.plot(self.altitude/1000.0,self.Tot_Dis_Tor,'k-',label='Total Disturbance Torque')
        plt.xlabel('Altitude (km)')
        plt.ylabel('Disturbance Torques (Nm)')
        plt.grid()
        plt.legend()
        pdfhandle.savefig()
        
        #plt.figure()
        #plt.plot(self.time,self.Solar_Torque,'b-',label='Solar Radiation Pressure')
        #plt.plot(self.time,self.Aero_Torque,'r-',label='Aerodynamics')
        #plt.plot(self.time,self.Gravity_Torque,'g-',label='Gravity Gradient')
        #plt.plot(self.time,self.Mag_Res_Dipole,'y-',label='Magnetic Resonance Dipole')
        #plt.plot(self.time,self.Tot_Dis_Tor,'k-',label='Total Disturbance Torque')
        #plt.yscale('log')
        #plt.xlabel('Time (sec)')
        #plt.ylabel('Disturbance Torques (Nm)')
        #plt.grid()
        #plt.legend()
        #pdfhandle.savefig()
        
        #plt.figure()
        #plt.plot(self.altitude/1000.0,self.Solar_Torque,'b-',label='Solar Radiation Pressure')
        #plt.plot(self.altitude/1000.0,self.Aero_Torque,'r-',label='Aerodynamics')
        #plt.plot(self.altitude/1000.0,self.Gravity_Torque,'g-',label='Gravity Gradient')
        #plt.plot(self.altitude/1000.0,self.Mag_Res_Dipole,'y-',label='Magnetic Resonance Dipole')
        #plt.plot(self.altitude/1000.0,self.Tot_Dis_Tor,'k-',label='Total Disturbance Torque')
        #plt.yscale('log')
        #plt.xlabel('Altitude (km)')
        #plt.ylabel('Disturbance Torques (Nm)')
        #plt.grid()
        #plt.legend()
        #pdfhandle.savefig()
        
        ##disturbance torques per axis
        #Disturbance_Torques = 1.0
        #self.Disturbance_Torques_per_axis = self.Total_Momentum/np.sqrt(3)
        #self.d_torques = np.asarray([self.Disturbance_Torques_per_axis,self.Disturbance_Torques_per_axis,self.Disturbance_Torques_per_axis])
        
    def Reaction_Wheel_Analysis(self):
        ##Compute reaction wheels that will provide the necessary momentum storage
        wx = self.wmax*np.pi/180.0
        wy = self.wmax*np.pi/180.0
        wz = self.wmax*np.pi/180.0      
        I = np.asarray([self.Ixx,self.Iyy,self.Izz])
        w = np.asarray([wx,wy,wz])
        H = I*w
        self.Hreq = H*self.FS
        print('Reaction Wheel Momentum Requirement (Nms) = ',np.max(self.Hreq))
        ##Based on disturbance torques
        ##Can compute number of times during mission to desaturate?
        Total_Orbits = (self.Mission_Duration * 30 * 24 * 60 * 60) / self.time[-1]
        print('Total Number of Orbits per Mission =  ', Total_Orbits)
        # Compute number of orbits before we need to desaturate
        number_of_orbits_before_desat_required = np.max(self.Hreq) / self.Total_Momentum
        print('Number of Orbits before Desat Required (Req) = ', number_of_orbits_before_desat_required)
        ##Compute the number of times we need to desaturate
        self.Desat_Mans = Total_Orbits / number_of_orbits_before_desat_required
        print('Number of Desaturation Maneuvers (Req) = ', self.Desat_Mans)

        number_of_orbits_before_desat_required_chosen = self.Chosen_Momentum / self.Total_Momentum
        print('Number of Orbits before Desat Required (Chosen) = ',number_of_orbits_before_desat_required_chosen)
        ##Compute the number of times we need to desaturate
        self.Desat_Mans_chosen = Total_Orbits / number_of_orbits_before_desat_required_chosen
        print('Number of Desaturation Maneuvers (Chosen) = ',self.Desat_Mans_chosen)

    def Magnetorquer_Analysis(self):
        Mag_Tor = self.btotal_arr*(1e-09)*self.Mag_Moment
        #print(Mag_Tor)
        
        plt.figure()
        plt.plot(self.altitude/1000.0,Mag_Tor)
        plt.grid()
        plt.xlabel('Altitude (km)')
        plt.ylabel('Torque from Magnetorquer (Nm)')
        pdfhandle.savefig()
        
        plt.figure()
        plt.plot(self.time,Mag_Tor)
        plt.grid()
        plt.xlabel('Time (sec)')
        plt.ylabel('Torque from Magnetorquer (Nm)')
        pdfhandle.savefig()
        
        plt.figure()
        plt.plot(self.time,Mag_Tor,'b-',label='Magnetorquer')
        plt.plot(self.time,self.Tot_Dis_Tor,'k-',label='Total Disturbance')
        plt.xlabel('Time (sec)')
        plt.ylabel('Torque (Nm)')
        plt.grid()
        plt.legend()
        pdfhandle.savefig()
        
        ###Now let's compute the delta torque
        self.Del_Tor = Mag_Tor - self.Tot_Dis_Tor
        #print(self.Del_Tor)
        Needed_Alt = []
        for x in range(len(self.Del_Tor)):
            if(self.Del_Tor[x] <= 0):
                #print(self.altitude[x])
                Needed_Alt.append(self.altitude[x])
        print('Altitude MagTs Turn Off (m from center of Earth):',Needed_Alt[0])
        Needed_Alt_AGL = Needed_Alt[0] - self.Rad_Earth
        print('Altitude MagTs Turn Off (m-AGL):',Needed_Alt_AGL)
        print('Altitude MagTs Turn Off (km-AGL):',Needed_Alt_AGL/1000.0)

        MagT_Alt = []
        for x in range(len(self.Tot_Dis_Tor)):
            #Check for Magnetorquer Effectiveness
            Delta_Tor = Mag_Tor[x] - self.Tot_Dis_Tor[x]
            ##If our magnetorquers are inneffective
            if(Delta_Tor <= 0):  ##Either disturbance torques go negative
            #if (self.altitude[x]/1000.0 > 20000): ##or I just turn off at aspecific altitude
                ##We need to turn the magnetorquers off
                Mag_Tor[x] = 0.0
                ##And then recompute delta_torque
                Delta_Tor = Mag_Tor[x] - self.Tot_Dis_Tor[x]
            ###Then here we append the DeltaTorque without magnetorquers on
            MagT_Alt.append(Delta_Tor)
        #print(MagT_Alt)
        Momentum_Needed = RiemannSum(MagT_Alt,self.time)
        print('Total Momentum Absorbed Per Orbit with Magnetorquers turning off (N-m-s) = ',Momentum_Needed)
        #self.Chosen_Momentum = np.max(self.Hreq)
        num_orbits_desat_magTs_off = self.Chosen_Momentum/Momentum_Needed
        print('Momentum Needed from RW Datasheet (Nms) = ',self.Chosen_Momentum)
        print('Number of Orbits Required to Desaturate RWs while turning off = ',num_orbits_desat_magTs_off)
        
        plt.figure()
        plt.plot(self.time,MagT_Alt)
        plt.xlabel('Time (sec)')
        plt.ylabel('Delta Torque (Magnetorquer - Disturbance Torques) (Nm)')
        plt.title('Magnetorquers Off at Alt (km) = '+str(Needed_Alt_AGL/1000.0))
        plt.grid()
        pdfhandle.savefig()

        plt.figure()
        plt.plot(self.time,self.Del_Tor)
        plt.xlabel('Time (sec)')
        plt.ylabel('Delta Torque (Magnetorquer - Disturbance Torques) (Nm)')
        plt.title('Magnetorquers Always ON')
        plt.grid()
        pdfhandle.savefig()

        ##Compute total momentum dump capability of mag Torquers
        Momentum_Diff = RiemannSum(self.Del_Tor,self.time)
        print('Total Momentum Absorbed Per Orbit with Magnetorquers (N-m-s) = ',Momentum_Diff)

        #Determine number of orbits required to desaturate rws
        num_orbits_desat_magTs = self.Chosen_Momentum/Momentum_Diff
        print('Number of Orbits Required to Desaturate RWs = ',num_orbits_desat_magTs)
        
        Norbits = []
        Tot_Mom = []
        for x in range(0,10):
            Mag_Tor_np = self.btotal_arr*(1e-09)*self.Mag_Moment_np[x]
            self.Del_Tor_np = Mag_Tor_np - self.Tot_Dis_Tor
            #print(self.Del_Tor_np[214],self.altitude[214],self.Del_Tor_np[215],self.altitude[215])
            Momentum_Diff_np = RiemannSum(self.Del_Tor_np,self.time)
            #print('Total Momeuntum Absorbed Per Orbit with Magnetorquers (N-m-s) = ',Momentum_Diff_np)
            Tot_Mom.append(Momentum_Diff_np)
            num_orbits_desat_magTs_np = self.Chosen_Momentum/Momentum_Diff_np
            #print('Number of Orbits Required to Desaturate RWs = ',num_orbits_desat_magTs_np)
            Norbits.append(num_orbits_desat_magTs_np)
            
        Norbits_np = np.asarray(Norbits)
        #print(Norbits_np)
        Tot_Mom_np = np.asarray(Tot_Mom)
        threshold = 0.0
        Norbits_clip = Norbits_np[Norbits_np>threshold]
        Tot_Mom_clip = Tot_Mom_np[Norbits_np>threshold]
        Mag_Mass_clip = self.Mag_Mass_np[Norbits_np>threshold]
        print('Number of Orbits to Desat RWs = ',Norbits_clip)
        print('Total Momentum Absorbed Per Orbit with Magnetorquers (N-m-s) = ',Tot_Mom_clip)
        plt.figure()
        plt.plot(Mag_Mass_clip,Norbits_clip,'.')
        plt.xlabel('Mass of Magnetorquers (g)')
        plt.ylabel('Number of Orbits to Desaturate RWs')
        plt.grid()
        pdfhandle.savefig()

##Inputs
example_inputs = np.loadtxt('ABEX_GNC_Data_File.txt')

FS = example_inputs[0] #factor of safety

Ixx = example_inputs[1] #kg-m^2
Iyy = example_inputs[2]
Izz = example_inputs[3]

wmax = example_inputs[4] #deg/s - Not sure where this comes from

length = example_inputs[5]/100.0 #meters
width = example_inputs[6]/100.0
height = example_inputs[7]/100.0

Mission_Duration = example_inputs[8] #months

CD = example_inputs[9] #Drag

rp = example_inputs[10] #perigee
#print('Perigee (km) = ',rp)
ra = example_inputs[11] #apogee

Mass_Sat = example_inputs[12]

Chosen_Momentum = example_inputs[13]

Mag_Moment = example_inputs[14]
#print('Magnetic Moment (Amp-m^2) = ',Mag_Moment)

mag = []
for x in range(15,25):
    rows = example_inputs[x]
    mag.append(rows)
    
Mag_Moment_np = np.asarray(mag)
#print('Mag Moment = ',Mag_Moment_np)

mag_mass = []
for x in range(25,35):
    row = example_inputs[x]
    mag_mass.append(row)
    
Mag_Mass_np = np.asarray(mag_mass)
#print('Mag Mass = ',Mag_Mass_np)

##Run the function above
pdfhandle = PdfPages('GNC_Analysis_Tool.pdf')
GNC = CubeSat(FS,Ixx,Iyy,Izz,wmax,length,width,height,Mission_Duration,CD,rp,ra,Mass_Sat,Mag_Moment,Mag_Moment_np,Mag_Mass_np,Chosen_Momentum,pdfhandle)
pdfhandle.close()
print('Plotting Routine Complete for Python')
#os.system('evince GNC_Analysis_Tool.pdf &')
#plt.show()