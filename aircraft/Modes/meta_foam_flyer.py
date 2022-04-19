import numpy as np
import sys
import matplotlib.pyplot as plt
import scipy.optimize as RF
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.patches as patches

#Pre-amble
#plt.close("all") #close all figures. 

##So I'm going to make a class called wing. because well we have two and I want to make this code absurd
class WING():
    def __init__(self,name,AR,LE):
        print 'Created ',name,' wing'
        self.name = name
        self.AR = AR
        self.LE = LE

    def compute_wingspan(self,S):
        self.S = S
        self.wingspan = np.sqrt(self.AR*self.S)

    def compute_chord(self):
        self.chord = self.S/self.wingspan
        #If we compute the chord we may as well compute the cg location and aero center
        self.compute_centers()

    def compute_centers(self):
        self.xac = self.LE + self.chord/4.0
        self.xcg = self.LE + self.chord/2.0

    def display_wing_geometry(self):
        print '------------',self.name,' wing properties-------------'
        print 'Wingspan (ft) = ',self.wingspan
        print 'Chord (ft) = ',self.chord
        print 'Area (ft^2) = ',self.S
        print 'AR = ',self.AR
        print 'Aerodynamic Center = ',self.xac

    def choose_airfoil(self,airfoil):
        if airfoil == '0008':
            self.Cl0 = 0.0
            self.Cd0 = 0.006
            self.Cda = 1.09145
            self.Cla = 2.0*np.pi
            self.CMAC = 0.0
        else:
            print 'Airfoil Not found. You need to code it in yourself in the choose_airfoil routine'
            sys.exit()

        #Now that we have the airfoil parameters we can get the wing parameters
        self.CL0 = self.AReqn(self.Cl0)
        self.CD0 = self.AReqn(self.Cd0)
        self.CDA = self.AReqn(self.Cda)
        self.CLA = self.AReqn(self.Cla)
        
        #Print Airfoil parameters
        print '-------------',self.name,' airfoil coefficients'
        print 'Cl0 = ',self.Cl0
        print 'Cd0 = ',self.Cd0
        print 'Cda = ',self.Cda
        print 'Cla = ',self.Cla
        #Print Wing parameters
        print '-------------',self.name,' wing coefficients'
        print 'CL0 = ',self.CL0
        print 'CD0 = ',self.CD0
        print 'CDA = ',self.CDA
        print 'CLA = ',self.CLA

    def AReqn(self,coeff):
        return coeff/(1. + coeff/(np.pi*self.AR))

    def Cl(self,alpha):
        return self.Cl0 + self.Cla*alpha

    def CL(self,alpha):
        return self.CL0 + self.CLA*alpha

    def Cd(self,alpha):
        return self.Cd0 + self.Cda*alpha**2

    def CD(self,alpha):
        return self.CD0 + self.CDA*alpha**2

##I also want a class called Airplane
class AIRPLANE():
    def __init__(self,AR,main_LE,ARt,tail_LE,de_ratio,dr_ratio,da_ratio,aileron_span_ratio):
        self.main_wing = WING('main',AR,main_LE)
        self.tail_wing = WING('tail',ARt,tail_LE)
        ARv = ARt/2.0
        self.vertical_wing = WING('vertical',ARv,tail_LE)
        self.de_ratio = de_ratio/100.
        self.dr_ratio = dr_ratio/100.
        self.da_ratio = da_ratio/100.
        self.ail2 = aileron_span_ratio[0]/100.
        self.ail1 = aileron_span_ratio[1]/100.
        self.pdfhandle = PdfPages('plots.pdf') #This code will make a pdf of figures called plots.
        #^Pretty slick.
        self.density = 0.00238
        self.sos = 1116.
        self.mu = 3.77e-7

    def no_thrust(self,alpha):
        y = self.CD0 + self.CDA*alpha**2 - np.tan(alpha)*self.CL0 - np.tan(alpha)*self.CLA*alpha
        return y

    def Reynolds(self,V):
        self.Re = self.density*V*self.main_wing.chord/self.mu
        print '--------------------------------'
        print 'Reynolds Number = ',self.Re
        print '--------------------------------'

    def Mach_No(self,V):
        self.M = V/self.sos

    def draw_plane(self):
        #Create a Figure
        fig1 = plt.figure()
        ax1 = fig1.add_subplot(111)
        #Draw the main wing
        ax1.add_patch(patches.Rectangle((self.main_wing.LE,-self.main_wing.wingspan/2.0),self.main_wing.chord,self.main_wing.wingspan,facecolor="blue"))
        #Tail wing
        ax1.add_patch(patches.Rectangle((self.tail_wing.LE,-self.tail_wing.wingspan/2.0),self.tail_wing.chord,self.tail_wing.wingspan,facecolor="red"))
        #Draw the center of mass
        plt.plot(self.xcg,0,'ko',markersize=20,label='xCG')
        #Draw the ac of main wing and tail and overall
        plt.plot(self.main_wing.xac,0,'ro',markersize=10,label='xAC_wing')
        plt.plot(self.tail_wing.xac,0,'bo',markersize=10,label='xAC_tail')
        plt.plot(self.xac,0,'go',markersize=20,label='xAC')
        #Set the axis limits
        ax1.set_xlim([0,self.tail_wing.LE+self.tail_wing.chord])
        ax1.set_ylim([-self.main_wing.wingspan/2.0,self.main_wing.wingspan/2.0])
        #Make aspect ratio equal
        plt.axes().set_aspect('equal', 'datalim')
        #Legend
        plt.legend()
        self.pdfhandle.savefig()

    def choose_airfoil(self,main_wing_airfoil,tail_wing_airfoil,vertical_wing_airfoil):
        self.main_wing.choose_airfoil(main_wing_airfoil)
        self.tail_wing.choose_airfoil(tail_wing_airfoil)
        self.vertical_wing.choose_airfoil(vertical_wing_airfoil)
        ### COMPUTE ALL LONGITUDINAL COEFFICIENTS
        #Now combine all the coefficients to get wing and tail combo
        #Note that this assumes the tail is a flat plate
        self.CM0 = self.main_wing.CMAC + self.smbar*self.main_wing.CL0
        self.CL0 = self.main_wing.CL0
        self.CD0 = self.main_wing.CD0 + (self.tail_wing.S/self.main_wing.S)*self.tail_wing.CD0
        self.CDA = self.main_wing.CDA + (self.tail_wing.S/self.main_wing.S)*self.tail_wing.CDA
        self.CLA = self.main_wing.CLA + (self.tail_wing.S/self.main_wing.S)*self.tail_wing.CLA
        self.CMA = self.CLA*self.smbar - self.VH*self.tail_wing.CLA
        self.CLQ = 2*self.VH*self.tail_wing.CLA
        self.CMQ = -self.lt/self.tail_wing.chord*self.CLQ
        self.tail_wing.CLDE = self.de_ratio*self.tail_wing.CLA
        self.CLDE = self.tail_wing.CLDE*(self.tail_wing.S/self.main_wing.S)
        self.CMDE = -(self.smbar + self.lt/self.main_wing.chord)*self.CLDE
        #Now print everything
        print '----------------Longitudinal Coefficients---------------'
        print 'CM0 = ',self.CM0
        print 'CL0 = ',self.CL0
        print 'CD0 = ',self.CD0
        print 'CDA = ',self.CDA
        print 'CLA = ',self.CLA
        print 'CMA = ',self.CMA
        print 'CLQ = ',self.CLQ
        print 'CMQ = ',self.CMQ
        print 'CLDE = ',self.CLDE
        print 'CMDE = ',self.CMDE
        self.Lift_Drag_Curves()

        ###COMPUTE ALL LATERAL COEFFICIENTS
        self.CYBETA = -1.0*self.vertical_wing.S/self.main_wing.S*self.vertical_wing.CLA
        self.CNBETA = self.VV*self.vertical_wing.CLA
        gamma = 0 #unless you have dihedral
        sweep = 0 #unless you have
        zvprime = -self.vertical_wing.wingspan/2.0
        self.CLBETA = -2.0/(3.0*np.pi)*self.main_wing.CLA*np.sin(gamma) - 2.0/(3.0*np.pi)*self.CL0*np.sin(2*sweep) - zvprime/self.main_wing.wingspan*self.CYBETA
        self.CYR = 2.0*self.VV*self.vertical_wing.CLA
        self.CLR = self.CL0/4.0 - 2.0*zvprime/self.main_wing.wingspan*self.VV*self.vertical_wing.CLA
        self.CNR = -self.CD0/4.0 - 2.0*self.lt/self.main_wing.wingspan*self.VV*self.vertical_wing.CLA
        self.CYP = -8.0/(3.0*np.pi)*self.vertical_wing.wingspan*self.vertical_wing.S/(self.main_wing.wingspan*self.main_wing.S)*self.vertical_wing.CLA
        self.CLP = -self.main_wing.CLA/6.0
        self.CNP = -self.lt/self.main_wing.wingspan*self.CYP - self.CL0/8.0
        self.CYDA = 0.0 ###typically zero for most a/c
        self.vertical_wing.CLDR = self.dr_ratio*self.vertical_wing.CLA
        self.CYDR = self.vertical_wing.CLDR*(self.vertical_wing.S/self.main_wing.S)
        self.CNDR = -self.VV*self.vertical_wing.CLDR
        self.y2 = self.main_wing.wingspan*self.ail2
        self.y1 = self.main_wing.wingspan*self.ail1
        self.CLDA = (2*self.main_wing.CLA*self.da_ratio*self.main_wing.chord)/(self.main_wing.S*self.main_wing.wingspan)*((self.y2**2)/2.-self.y1**2/2.)
        self.CLDR = -self.vertical_wing.S*zvprime/(self.main_wing.S*self.main_wing.wingspan)*self.vertical_wing.CLDR
        self.CNDA = 0.0 ##Assumed to be zero for this aircraft

        print 'CYBETA = ',self.CYBETA
        print 'CNBETA = ',self.CNBETA
        print 'CLBETA = ',self.CLBETA
        print 'CYR = ',self.CYR
        print 'CLR = ',self.CLR
        print 'CNR = ',self.CNR
        print 'CYP = ',self.CYP
        print 'CLP = ',self.CLP
        print 'CNP = ',self.CNP
        print 'CYDA = ',self.CYDA
        print 'CYDR = ',self.CYDR
        print 'CNDR = ',self.CNDR
        print 'CLDA = ',self.CLDA
        print 'CLDR = ',self.CLDR
        print 'CNDA = ',self.CNDA

    def compute_Eigenvalues(self,Vtrim):
        self.Vtrim = Vtrim[0]
        self.Q = 0.5*self.density*self.Vtrim**2
        self.XU = -(self.Q*self.main_wing.S*(2*self.CD0))/(self.mass*self.Vtrim)
        self.CXA = -self.CDA #This is wrong
        self.XW = self.Q*self.main_wing.S*self.CXA/(self.mass*self.Vtrim)
        self.XWDOT = 0
        self.XQ = 0
        self.CZ0 = -self.CL0
        self.ZU = (2*self.CZ0)*self.Q*self.main_wing.S/(self.mass*self.Vtrim)
        self.mu = self.Q*self.main_wing.S/(self.mass*self.Vtrim)
        self.CZA = -self.CLA
        self.ZW = self.mu*self.CZA
        self.ZWDOT = 0.
        self.muu = self.Q*self.main_wing.S*self.main_wing.chord/(2*self.mass*self.Vtrim)
        self.CZQ = -self.CLQ
        self.ZQ = self.muu*self.CZQ
        self.muuy = self.Q*self.main_wing.chord*self.main_wing.S/(self.Iy*self.Vtrim)
        self.CMU = 0.
        self.MU = self.muuy*self.CMU
        self.MW = self.muuy*self.CMA
        self.MWDOT = 0.
        self.muuy2 = self.Q*(self.main_wing.chord**2)*self.main_wing.S/(2*self.Iy*self.Vtrim)
        self.MQ = self.muuy2*self.CMQ
        self.YV = self.mu*self.CYBETA
        self.mub = self.Q*self.main_wing.S*self.main_wing.wingspan/(2*self.mass*self.Vtrim)
        self.YP = self.mub*self.CYP
        self.YR = self.mub*self.CYR
        self.muix = self.Q*self.main_wing.S*self.main_wing.wingspan/(self.Ix*self.Vtrim)
        self.LV = self.muix * self.CLBETA
        self.muix2 = self.Q*self.main_wing.S*self.main_wing.wingspan**2/(2*self.Ix*self.Vtrim)
        self.LP = self.muix2*self.CLP
        self.LR = self.muix2*self.CLR
        self.muiz = self.Q*self.main_wing.S*self.main_wing.wingspan/(self.Iz*self.Vtrim)
        self.NV = self.muiz*self.CNBETA
        self.muiz2 = self.Q*self.main_wing.S*self.main_wing.wingspan**2/(2.*self.Iz*self.Vtrim)
        self.NP = self.muiz2*self.CNP
        self.NR = self.muiz2*self.CNR
        print '------------------------'
        print 'Trim Velocity (ft/s) = ',Vtrim
        print 'Dynamic Pressure = ',self.Q
        print 'XU = ',self.XU
        print 'XWDOT = ',self.XWDOT
        print 'XQ = ',self.XQ
        print 'ZU = ',self.ZU
        print 'ZW = ',self.ZW
        print 'ZWDOT = ',self.ZWDOT
        print 'ZQ = ',self.ZQ
        print 'MU = ',self.MU
        print 'MW = ',self.MW
        print 'MQ = ',self.MQ
        print '------------------------'
        print 'LONGITUDINAL MATRIX'
        theta0 = 0.
        g0 = 32.2
        self.ALON = [[self.XU,self.XW,0,-g0*np.cos(theta0)],[self.ZU,self.ZW,self.Vtrim,-g0*np.sin(theta0)],[self.MU+self.MWDOT*self.ZU,self.MW+self.MWDOT*self.ZW,self.MQ+self.Vtrim*self.MWDOT,-self.MWDOT*g0*np.sin(theta0)],[0.,0.,1.,0.]]
        eigenvalues,eigenvectors = np.linalg.eig(self.ALON)
        self.lon_modes = eigenvalues
        print self.ALON
        for i in range(0,len(self.lon_modes)):
            print self.lon_modes[i]
        #print np.shape(self.ALON)
        print '-------------------------'
        #Atest = [[-0.0188,0.0415,0,-32.174],[-0.1987,-0.5193,264.312,0.],[0.,-0.0018,-0.4898,0.],[0.,0.,1.,0.]]
        #e,v = np.linalg.eig(Atest)
        #print e
        print 'LATERAL MATRIX'
        self.ALAT = [[self.YV,self.YP,g0*np.cos(theta0),self.YR-self.Vtrim],[self.LV,self.LP,0.,self.LR],[0.,1.,0.,0.],[self.NV,self.NP,0.,self.NR]]
        eigenvalues,eigenvectors = np.linalg.eig(self.ALAT)
        self.lat_modes = eigenvalues
        print self.ALAT
        for i in range(0,len(self.lat_modes)):
            print self.lat_modes[i]
        print '--------------------------'
    def Lift_Drag_Curves(self):
        print 'Generating Plots....'
        #LIFT vs. ALPHA
        plt.figure()
        alpha = np.linspace(-10.,10.,100)*np.pi/180.
        #Airfoil Lift
        Cl_wing_airfoil = self.main_wing.Cl(alpha)
        Cl_tail_airfoil = self.tail_wing.Cl(alpha)
        plt.plot(alpha*180./np.pi,Cl_wing_airfoil,'b--',label=self.main_wing.name+' airfoil')
        plt.plot(alpha*180./np.pi,Cl_tail_airfoil,'r--',label=self.tail_wing.name+' airfoil')
        plt.grid()
        plt.xlabel('Angle of Attack (deg)')
        plt.ylabel('Lift Coefficient (nd)')
        #Wing lift
        CL_wing = self.main_wing.CL(alpha)
        CL_tail = self.tail_wing.CL(alpha)
        plt.plot(alpha*180./np.pi,CL_wing,'b-',label=self.main_wing.name+' wing')
        plt.plot(alpha*180./np.pi,CL_tail,'r-',label=self.tail_wing.name+' wing')
        #Total Aircraft Lift
        CL_total = self.CL(alpha,0)
        CL_plus30 = self.CL(alpha,30*np.pi/180.)
        CL_minus30 = self.CL(alpha,-30*np.pi/180.)
        plt.plot(alpha*180./np.pi,CL_total,'k-',label='Total Aircraft')
        plt.plot(alpha*180./np.pi,CL_plus30,'k--',label='+30 deg Elev')
        plt.plot(alpha*180./np.pi,CL_minus30,'k-.',label='-30 deg Elev')
        plt.legend(loc='best')
        self.pdfhandle.savefig()
        print 'CL vs Alpha'

        #DRAG vs. ALPHA
        plt.figure()
        #Airfoil Drag
        Cd_wing_airfoil = self.main_wing.Cd(alpha)
        Cd_tail_airfoil = self.tail_wing.Cd(alpha)
        plt.plot(alpha*180./np.pi,Cd_wing_airfoil,'b--',label=self.main_wing.name+' airfoil')
        plt.plot(alpha*180./np.pi,Cd_tail_airfoil,'r--',label=self.tail_wing.name+' airfoil')
        plt.grid()
        plt.xlabel('Angle of Attack (deg)')
        plt.ylabel('Drag Coefficient (nd)')
        #Wing lift
        CD_wing = self.main_wing.CD(alpha)
        CD_tail = self.tail_wing.CD(alpha)
        plt.plot(alpha*180./np.pi,CD_wing,'b-',label=self.main_wing.name+' wing')
        plt.plot(alpha*180./np.pi,CD_tail,'r-',label=self.tail_wing.name+' wing')
        #Total Aircraft Lift
        CD_total = self.CD(alpha)
        plt.plot(alpha*180./np.pi,CD_total,'k-',label='Total Aircraft')
        plt.legend(loc='best')
        self.pdfhandle.savefig()
        print 'CD vs Alpha'

        #CL vs. DRAG
        plt.figure()
        plt.plot(Cd_wing_airfoil,Cl_wing_airfoil,'b--',label=self.main_wing.name+' airfoil')
        plt.plot(Cd_tail_airfoil,Cl_tail_airfoil,'r--',label=self.tail_wing.name+' airfoil')
        plt.grid()
        plt.ylabel('Lift Coefficient (nd)')
        plt.xlabel('Drag Coefficient (nd)')
        #Wing lift
        plt.plot(CD_wing,CL_wing,'b-',label=self.main_wing.name+' wing')
        plt.plot(CD_tail,CL_tail,'r-',label=self.tail_wing.name+' wing')
        #Total Aircraft Lift
        plt.plot(CD_total,CL_total,'k-',label='Total Aircraft')
        plt.legend(loc='best')
        self.pdfhandle.savefig()
        print 'CL vs CD'

        alpha_deg = alpha*180./np.pi

        #L/D
        plt.figure()
        plt.plot(alpha_deg,Cl_wing_airfoil/Cd_wing_airfoil,'b--',label=self.main_wing.name+' airfoil')
        plt.plot(alpha_deg,Cl_tail_airfoil/Cd_tail_airfoil,'r--',label=self.tail_wing.name+' airfoil')
        plt.grid()
        plt.ylabel('Lift Coefficient (nd)')
        plt.xlabel('Angle of Attack (deg)')
        #Wing lift
        plt.plot(alpha*180./np.pi,CL_wing/CD_wing,'b-',label=self.main_wing.name+' wing')
        plt.plot(alpha*180./np.pi,CL_tail/CD_tail,'r-',label=self.tail_wing.name+' wing')
        #Total Aircraft Lift
        plt.plot(alpha*180./np.pi,CL_total/CD_total,'k-',label='Total Aircraft')
        plt.legend(loc='best')
        self.pdfhandle.savefig()
        print 'L/D vs Alpha'

        #CMA vs. AOA
        plt.figure()
        plt.grid()
        plt.xlabel('Angle of Attack (deg)')
        plt.ylabel('Pitch Moment Coefficient (nd)')
        #Total Pitch Moment
        CM_total = self.CM(alpha,0)
        CM_plus30 = self.CM(alpha,30*np.pi/180.)
        CM_minus30 = self.CM(alpha,-30*np.pi/180.)
        plt.plot(alpha*180./np.pi,CM_total,'k-',label='Total Aircraft')
        plt.plot(alpha*180./np.pi,CM_plus30,'k--',label='+30 deg')
        plt.plot(alpha*180./np.pi,CM_minus30,'k-.',label='-30 deg')
        plt.legend(loc='best')
        self.pdfhandle.savefig()
        print 'CMA vs Alpha'

    def CL(self,alpha,de):
        return self.CL0 + self.CLA*alpha + self.CLDE*de

    def CD(self,alpha):
        return self.CD0 + self.CDA*alpha**2

    def CM(self,alpha,de):
        return self.CM0 + self.CMA*alpha + self.CMDE*de
        
    def compute_wingspan(self,S,St):
        self.main_wing.compute_wingspan(S)
        self.tail_wing.compute_wingspan(St)
        self.vertical_wing.compute_wingspan(St/2.0)

    def compute_VH(self):
        self.VH = self.tail_wing.S*self.lt/(self.main_wing.S*self.main_wing.chord)

    def compute_VV(self):
        self.VV = self.vertical_wing.S*self.lt/(self.main_wing.S*self.main_wing.wingspan)

    def compute_chord(self):
        self.main_wing.compute_chord()
        self.tail_wing.compute_chord()
        self.vertical_wing.compute_chord()
        #Once we compute the chord we also need to compute lt the distance between aero centers
        self.lt = self.tail_wing.xac - self.main_wing.xac
        #We may as well compute the aero center as well
        self.compute_aero_center()
        #Probably need the tail volume coefficient
        self.compute_VH()
        #Need the vertical tail volume coefficient too
        self.compute_VV()
        #And the elevator chord length
        self.elevator_chord = self.de_ratio*self.tail_wing.chord

    def compute_aero_center(self):
        self.xac = (self.main_wing.S*self.main_wing.xac + self.tail_wing.S*self.tail_wing.xac)/(self.main_wing.S + self.tail_wing.S)

    def display_geometry(self):
        self.main_wing.display_wing_geometry()
        self.tail_wing.display_wing_geometry()
        self.vertical_wing.display_wing_geometry()
        print '-------------Wing and Tail Combo----------'
        print 'Total Weight (lbf) = ',self.weight
        print 'Distance from Main Wing to Tail Aero Centers (ft) = ',self.lt
        print 'Inertias = ',self.Ix,self.Iy,self.Iz
        print 'Aerodynamic Center (ft) = ',self.xac
        print 'Center of Mass (ft) = ',self.xcg
        print 'Stability Margin (nd) = ',self.smbar
        print 'Tail Volume Coefficient (nd) = ',self.VH
        print 'Elevator Chord (ft) = ',self.elevator_chord

    def compute_wing_body_geometry(self,rho_foam):
        ##All weights in lbf and position of items in feet to make sure we have correct ENGLISH units

        ##Servos
        m_servo = 0.0308
        ##Also need position of every item
        x_servo1 = 0.6666666667 ##Aileron (feet)
        x_servo2 = 1.6666666667 ##Elevator (feet)
        #This aircraft has no rudder
        #ESC
        m_esc = 0.0572
        x_esc = 0.25
        #Motor
        m_motor = 0.2068
        x_motor = -0.0416666667
        #Receiver
        m_rec = 0.0264
        x_rec = 0.875
        #This airplane has an Arduino on it
        m_ard = 0.176
        x_ard = -0.9583333333
        #And thus two batteries
        m_2cell = 0.055
        m_3cell = 0.121
        x_2cell = -1.2916666667
        x_3cell = -0.375

        #In order to get this we need to determine the size of the aircraft. So. the plan is to
        #make this a function and compute the weight of the aircraft depending on area of the wing
        #Alright since S and St are inputs to this function we can just compute the centers ourselves
        m_main_wing = rho_foam*self.main_wing.S
        m_tail_wing = rho_foam*self.tail_wing.S
        print 'Weight of Main Wing = ',m_main_wing
        print 'Weight of Tail Wing = ',m_tail_wing

        #Ok to make this code modular we will put everything into arrays here
        weights = np.array([m_servo,m_servo,m_esc,m_motor,m_rec,m_ard,m_2cell,m_3cell,m_main_wing,m_tail_wing])
        positions = np.array([x_servo1,x_servo2,x_esc,x_motor,x_rec,x_ard,x_2cell,x_3cell,self.main_wing.xcg,self.tail_wing.xcg])

        #Then our total weight is just
        self.weight = np.sum(weights)
        self.mass = self.weight/32.2
        #Approximate Inertias
        self.Ix = (self.mass/12.)*self.main_wing.wingspan**2
        self.Iy = (self.mass/12.)*self.lt**2
        self.Iz = (self.mass/12.)*(self.lt**2 + self.main_wing.wingspan**2)

        #and center of mass is
        self.xcg = np.sum(weights*positions)/self.weight

        #so we can also compute the stability margin
        self.SM = self.xcg - self.xac
        self.smbar = self.SM/self.main_wing.chord


##################THIS IS THE ACTUAL START OF THE ROUTINE########################
##################EVERYTHING ABOVE IS JUST SUBROUTINES AND CLASSES##############
        
##One of the things we need to design around is our wing loading
##For this aircraft we want
WS = 3.0 #lbf/ft^2 How do I know this? Because I've built this airplane before. You'll need to change this
#for your airplane.

##We also need our Tail Area Ratio
tail_area_ratio = 33.48837209 #How do I know this? Cuz it looks good.
#Also need to choose how big of an elevator we want
elevator_ratio = 50.
#How big of a rudder?
rudder_ratio = 50.
#Aileron? Need chord ratio and span ratio
aileron_chord_ratio = 50.
aileron_span_ratio = [100.,20.]

#We also need to choose an AR of the main wing and tail
AR = 4.0 #Why? Cuz this is a fighter and we kind of built this airplane without any calculations
ARt = 2.0 #Why? Same reason as above
#WE also need the positions of the main wing and tail
main_LE = 5.75/12.0
tail_LE = 27./12.0
#So I have some fancy classes above so I'm going to initialize those here
meta = AIRPLANE(AR,main_LE,ARt,tail_LE,elevator_ratio,rudder_ratio,aileron_chord_ratio,aileron_span_ratio)

###The First step in designing an airplane is getting an estimate
##of your weight. Some of these things we know other we don't
#In order to make this code really slick we're going to have S and St be inputs
#to our routine here as well as the position of the tail
##Alright so now we need to iterate since we only have an estimate of weight

##So the hard part. Weight of the fuselage. This depends on the area of the wing and tail among other things
density_of_foam = 0.1412905923 #lbf/ft^3

meta.weight = 1.0  #First we need a guess of our weight. I'm going to assume the airplane is 1 lb
N = 100
for x in range(0,N): #So we will iterate N times
    #Since we have the weight of our airplane we can compute the Area of both wings put together
    print 'Iteration = ',x
    print 'Current Weight = ',meta.weight,' (lbf)'
    S = meta.weight*WS
    print 'Area of Wing (ft^2) = ',S
    #Based on that area we need to compute the area of the tail
    St = tail_area_ratio/100.*S
    #Now that we know that we can compute wingspan of main wing and tail
    meta.compute_wingspan(S,St)
    #At this point we can get the chord assuming rectangular wings. You'll need
    #Different code if you want swept or something different
    meta.compute_chord()
    #Once we have the size of the main wing, tail, wingspan and chord we can compute the mass and center of mass
    meta.compute_wing_body_geometry(density_of_foam)
    #^So this routine computes the weight of the aircraft
    #Which means when we iterate we should hopefully converge to the actual value. Not sure but hey let's see
    #If forget what this optimization method is called but it's pretty stout

#Once the routine Kicks out we can kick out alot of stuff.
meta.display_geometry()

#With the geometry displayed we can make a figure with rectangles and draw all pertinent things on the axes.
meta.draw_plane()

#Let's get a ball park Reynolds number
meta.Reynolds(49.)

#And then it's time to design our airfoils
#The first argument is the airfoil of the main wing and the second is the airfoil of the second
#Note that I've hard coded 0008 from XFLR5 but you'll have to get your own stuff
meta.choose_airfoil('0008','0008','0008') ##so this is main wing, tail, and vertical stabilizer
#^^This routine also makes a bunch of plots such as CLvsalpha,cdvsalpha,etc.....

#Alright so now let's make some plots that are a function of velocity.
V = np.linspace(20,50,1000)

##With V set we can get CL
CL = 2*meta.weight/(meta.density*meta.main_wing.S*V**2)

#With CL we can get alpha and de but that requires some matrix algebra
# CL - CL0 = [CLA CLDE]*[alpha ; de]
# 0        = [CMA CMDE]*[alpha ; de]
# and a loop
alpha = []
de = []
mat = np.asarray([[meta.CLA,meta.CLDE],[meta.CMA,meta.CMDE]])
invmat = np.linalg.inv(mat)
for CLi in CL:
    sol = np.asarray([CLi-meta.CL0,0])
    x = np.matmul(invmat,sol)
    alpha.append(x[0])
    de.append(x[1])
alpha = np.asarray(alpha)
de = np.asarray(de)

#Now we can get CD and CM
CD = meta.CD(alpha)
CM = meta.CM0 + meta.CMA*alpha + meta.CMDE*de #<<---plotting this is pointless. This is always zero

#Total Drag and of course Thrust
D = 0.5*meta.density*V**2*meta.main_wing.S*CD

#Find Angle of Attack for No control
l = np.where(np.abs(de) < 1e-2)
print '-----------NO CONTROL------------------'
if len(l[0]) == 0:
    print 'No control Flight Mode Impossible'
else:
    print '--------No Control Flight Mode-------'
    print 'Angle of Attack = ',alpha[l[0]]
    print 'CL (nd) = ',CL[l[0]]
    print 'V (ft/s) = ',V[l[0]]
    print 'Elevator (rad) = ',de[l[0]]
    print 'CD (nd) = ',CD[l[0]]
    print 'L/D (nd) = ',CL[l[0]]/CD[l[0]]
    print 'D (lbf) = ',D[l[0]]

#Find Angle of Attack for MAX L/D or Min Drag
l_md = np.where(D==np.min(D))
print '-------------MINIMUM DRAG / MAX L/D -------------'
alpha_md = alpha[l_md[0]]
print 'Angle of Attack (rad) = ',alpha_md
CL_md = CL[l_md[0]]
print 'CL (nd) = ',CL_md
V_md = V[l_md[0]]
print 'V (ft/s) = ',V_md
de_md = de[l_md[0]]
print 'Elevator (rad) = ',de_md
CD_md = CD[l_md[0]]
print 'CD (nd) = ',CD_md
print 'L/D (nd) = ',CL_md/CD_md
D_md = D[l_md[0]]
print 'D (lbf) = ',D_md

meta.compute_Eigenvalues(V_md)

#Find No thrust Angle of Attack
alpha_no_thrust = np.linspace(-10,10,100)*np.pi/180.0
y = meta.no_thrust(alpha_no_thrust)
sol = RF.bisect(meta.no_thrust,0.0,0.15)

print '-------------No Thrust Solution----------------'
alpha_nt = sol
CD_nt = meta.CD(sol)
de_nt = (-meta.CM0 - meta.CMA*alpha_nt)/meta.CMDE
CL_nt = meta.CL(alpha_nt,de_nt)
V_nt = np.sqrt(2*meta.weight/(meta.density*meta.main_wing.S*(CL_nt*np.cos(alpha_nt) + CD_nt*np.sin(alpha_nt))))
print 'Alpha (rad) = ',alpha_nt
print 'CL (nd) = ',CL_nt
print 'Velocity (ft/s) = ',V_nt
print 'Elevator (rad) = ',de_nt
print 'CD (nd) = ',CD_nt
print 'L/D (nd) = ',CL_nt/CD_nt
print '-------------------------------------------------'

#Now we can plot stuff
print 'Generating more plots.....'

plt.figure()
plt.plot(alpha_no_thrust*180.0/np.pi,y)
plt.plot(sol*180.0/np.pi,0,'rx')
plt.xlabel('Alpha (deg)')
plt.ylabel('Root Finding Function')
plt.grid()
meta.pdfhandle.savefig()
print 'No Thrust'

plt.figure()
plt.plot(V,CL)
plt.plot(V_md,CL_md,'rx')
plt.grid()
plt.xlabel('Velocity (ft/s)')
plt.ylabel('Lift Coefficient (nd)')
#plt.ylim([0,2.0]) #Max CL is 2.0 so no need to plot that
meta.pdfhandle.savefig()
print 'CL vs V'

plt.figure()
plt.plot(V,alpha*180./np.pi)
plt.plot(V_md,alpha_md*180./np.pi,'rx')
plt.grid()
plt.xlabel('Velocity (ft/s)')
plt.ylabel('Alpha (deg)')
#plt.ylim([0,2.0]) #Max CL is 2.0 so no need to plot that
meta.pdfhandle.savefig()
print 'Alpha vs V'

plt.figure()
plt.plot(V,de*180./np.pi)
plt.plot(V_md,de_md*180./np.pi,'rx')
plt.grid()
plt.xlabel('Velocity (ft/s)')
plt.ylabel('Evelator Angle (deg)')
#plt.ylim([0,2.0]) #Max CL is 2.0 so no need to plot that
meta.pdfhandle.savefig()
print 'DE vs V'

plt.figure()
plt.plot(V,CD)
plt.plot(V_md,CD_md,'rx')
plt.grid()
plt.xlabel('Velocity (ft/s)')
plt.ylabel('Drag Coefficient (nd)')
#plt.ylim([0,2.0]) #Max CL is 2.0 so no need to plot that
meta.pdfhandle.savefig()
print 'CD vs V'

plt.figure()
plt.plot(V,CL/CD)
plt.plot(V_md,CL_md/CD_md,'rx')
plt.grid()
plt.xlabel('Velocity (ft/s)')
plt.ylabel('L/D (nd)')
#plt.ylim([0,2.0]) #Max CL is 2.0 so no need to plot that
meta.pdfhandle.savefig()
print 'L/D vs V'

plt.figure()
plt.plot(V,D)
plt.plot(V_md,D_md,'rx')
plt.grid()
plt.xlabel('Velocity (ft/s)')
plt.ylabel('D (lbf)')
#plt.ylim([0,2.0]) #Max CL is 2.0 so no need to plot that
meta.pdfhandle.savefig()
print 'D vs V'

#Close the pdf
meta.pdfhandle.close()

