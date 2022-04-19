#!/usr/bin/python

import numpy as np
import sys
import matplotlib.pyplot as plt
import scipy.optimize as RF
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.patches as patches
import os

#Pre-amble
#plt.close("all") #close all figures. 

##So I'm going to make a class called wing. because well we have two and I want to make this code absurd
class WING():
    def __init__(self,name,AR,LE):
        print('Created ',name,' wing')
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
        print('------------',self.name,' wing properties-------------')
        print('Wingspan (ft) = ',self.wingspan)
        print('Chord (ft) = ',self.chord)
        print('Area (ft^2) = ',self.S)
        print('AR = ',self.AR)
        print('Aerodynamic Center = ',self.xac)

    def choose_airfoil(self,airfoil):
        if airfoil == '0008':
            self.Cl0 = 0.0
            self.Cd0 = 0.006
            self.Cda = 1.09145
            self.Cla = 2.0*np.pi
            self.CMAC = 0.0
        elif airfoil == 'ClarkY':
            self.Cl0 = 0.34995
            self.Cd0 = 0.015342
            self.Cda = 0.713431
            self.Cla = 4.8872
            self.CMAC = -0.065
        else:
            print('Airfoil Not found. You need to code it in yourself in the choose_airfoil routine')
            sys.exit()

        #Now that we have the airfoil parameters we can get the wing parameters
        self.CL0 = self.AReqn(self.Cl0)
        self.CD0 = self.AReqn(self.Cd0)
        self.CDA = self.AReqn(self.Cda)
        self.CLA = self.AReqn(self.Cla)
        
        #Print Airfoil parameters
        print('-------------',self.name,' airfoil coefficients')
        print('Cl0 = ',self.Cl0)
        print('Cd0 = ',self.Cd0)
        print('Cda = ',self.Cda)
        print('Cla = ',self.Cla)
        #Print Wing parameters
        print('-------------',self.name,' wing coefficients')
        print('CL0 = ',self.CL0)
        print('CD0 = ',self.CD0)
        print('CDA = ',self.CDA)
        print('CLA = ',self.CLA)

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
    def __init__(self,AR,main_LE,ARt,tail_LE,de_ratio,ARv):
        self.main_wing = WING('main',AR,main_LE)
        self.tail_wing = WING('tail',ARt,tail_LE)
        self.vertical_wing = WING('vertical',ARv,tail_LE)
        self.de_ratio = de_ratio/100.
        self.pdfhandle = PdfPages('plots.pdf') #This code will make a pdf of figures called plots.
        self.density = 0.00238 #Density of air
        self.sos = 1116. #Speed of sound
        self.mu = 3.77e-7 #Free viscocity of air

    def no_thrust(self,alpha):
        y = self.CD0 + self.CDA*alpha**2 - np.tan(alpha)*self.CL0 - np.tan(alpha)*self.CLA*alpha
        return y

    def Reynolds(self,V):
        self.Re = self.density*V*self.main_wing.chord/self.mu
        print('--------------------------------')
        print('Reynolds Number = ',self.Re)
        print('--------------------------------')

    def Mach_No(self,V):
        self.M = V/self.sos

    def draw_plane(self):
        #Create a Figure
        fig1 = plt.figure()
        plt.xlabel('Feet from nose of the airplane to tail wing trailing edge')
        plt.ylabel('Feet from centerline of the fuselage to wing tip')
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
        print('----------------Longitudinal Coefficients---------------')
        print('CM0 = ',self.CM0)
        print('CL0 = ',self.CL0)
        print('CD0 = ',self.CD0)
        print('CDA = ',self.CDA)
        print('CLA = ',self.CLA)
        print('CMA = ',self.CMA)
        print('CLQ = ',self.CLQ)
        print('CMQ = ',self.CMQ)
        print('CLDE = ',self.CLDE)
        print('CMDE = ',self.CMDE)
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
        
        #This adds the last few coefficients
        self.Cyv = (Cr/Cv)*self.vertical_wing.CLA
        self.Cydelr = (Sv/Sm)*self.Cyv
        self.CLdela = ((self.CLA*Ca)/((Bm**2)*self.main_wing.chord)*((Y2**2)-(Y1**2)))
        self.CLdelr = -(zvprime/Bm)*self.Cydelr
        self.Cndela = 0.216 #This number was a guess and may need adjusted
        self.Cndelr = -(self.smbar+self.lt/self.main_wing.chord)*self.Cydelr

        print('CYBETA = ',self.CYBETA)
        print('CNBETA = ',self.CNBETA)
        print('CLBETA = ',self.CLBETA)
        print('CYR = ',self.CYR)
        print('CLR = ',self.CLR)
        print('CNR = ',self.CNR)
        print('CYP = ',self.CYP)
        print('CLP = ',self.CLP)
        print('CNP = ',self.CNP)
        print('Cydelr = ',self.Cydelr)
        print('CLdela = ',self.CLdela)
        print('CLdelr = ',self.CLdelr)
        print('Cndela = ',self.Cndela)
        print('Cndelr = ',self.Cndelr)

    def Lift_Drag_Curves(self):
        print('Generating Plots....')
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
        print('CL vs Alpha')

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
        print('CD vs Alpha')

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
        print('CL vs CD')

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
        print('L/D vs Alpha')

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
        print('CMA vs Alpha')

    def CL(self,alpha,de):
        return self.CL0 + self.CLA*alpha + self.CLDE*de

    def CD(self,alpha):
        return self.CD0 + self.CDA*alpha**2

    def CM(self,alpha,de):
        return self.CM0 + self.CMA*alpha + self.CMDE*de
        
    def compute_wingspan(self,S,St,Sv):
        self.main_wing.compute_wingspan(S)
        self.tail_wing.compute_wingspan(St)
        self.vertical_wing.compute_wingspan(Sv)

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
        print('Wing Aero Distance = ',self.lt)
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
        print('-------------Wing and Tail Combo----------')
        print('Total Weight (lbf) = ',self.weight)
        print('Distance from Main Wing to Tail Aero Centers (ft) = ',self.lt)
        print('Inertias = ',self.Ix,self.Iy,self.Iz)
        print('Aerodynamic Center (ft) = ',self.xac)
        print('Center of Mass (ft) = ',self.xcg)
        print('Stability Margin (nd) = ',self.smbar)
        print('Tail Volume Coefficient (nd) = ',self.VH)
        print('Elevator Chord (ft) = ',self.elevator_chord)

    def compute_wing_body_geometry(self):
        #Then our total weight is just
        self.weight = 2.72 #lbf
        self.mass = self.weight/32.2 #slugs
        #Approximate Inertias
        self.Ix = (self.mass/12.)*self.main_wing.wingspan**2
        NoseTail = 42.5/12 #ft
        Depth = 7.25/12 #ft
        self.Iy = (self.mass/12.)*((NoseTail)**2+(Depth)**2)
        #self.Iz = (self.mass/12.)*(self.lt**2 + self.main_wing.wingspan**2)
        self.Iz = (self.mass/12.)*((NoseTail)**2+self.main_wing.wingspan**2)

        #and center of mass is
        self.xcg = 13.0/12.0 #ft from leading edge

        #so we can also compute the stability margin
        self.SM = self.xcg - self.xac
        self.smbar = self.SM/self.main_wing.chord

#Measurements taken from META 1
Weight = 2.72 #lbf       
Bm = 59.0/12.0 #Main wingspan (ft)
Bt = 23.0/12.0 #Tail wingspan (ft)
Cm = 9.0/12.0 #Main wing chord (ft)
Ct = 7.5/12.0 #Tail wing chord (ft)
Ce = 2.5/12.0 #Elevator chord (ft)
Sm = Bm*Cm #Area of main wing (ft**2)
St = Bt*Ct #Area of tail wing (ft**2)
WS = Weight/Sm #weight of plane divide by area of main wing (lbf/ft**2)
tail_area_ratio = (St/Sm)*100
elevator_ratio = (Ce/Ct)*100
AR = ((Bm)**2)/(Sm) #Ratio of main wing's wingspan divided by main wing's area
ARt = ((Bt)**2)/(St) #Ratio of tail wing's wingspan divided by tail wing's area
main_LE = 10.0/12.0 #Main wing leading edge (ft)
tail_LE = 32.0/12.0 #Tail wing leading edge (ft)

#Vertical Tail Parameters
Bv = 8./12.
Cv = 7./12.
Sv = Bv*Cv
ARv = Bv**2/Sv
Y1 = 26.75/12. #Start of Aleron
Y2 = 14.25/12. #End of Aleron
Ca = 2.25/12. #Chord of Aleron

#Measurements needed for last few coefficients
Cr = 2./12. #Chord of the rudder (ft)
Br = 7.38/12. #Wingspan of the rudder (ft)

#So I have some fancy classes above so I'm going to initialize those here
meta = AIRPLANE(AR,main_LE,ARt,tail_LE,elevator_ratio,ARv)

#Now that we know that we can compute wingspan of main wing and tail
meta.compute_wingspan(Sm,St,Sv)

#At this point we can get the chord assuming rectangular wings. You'll need
#Different code if you want swept or something different
meta.compute_chord()

#Once we have the size of the main wing, tail, wingspan and chord we can compute the mass and center of mass
meta.compute_wing_body_geometry()

#Once the routine Kicks out we can kick out alot of stuff.
meta.display_geometry()

#With the geometry displayed we can make a figure with rectangles and draw all pertinent things on the axes.
meta.draw_plane()

#Let's get a ball park Reynolds number
meta.Reynolds(42.64) 

#And then it's time to design our airfoils
#The first argument is the airfoil of the main wing and the second is the airfoil of the second
#Note that I've hard coded 0008 from XFLR5 but you'll have to get your own stuff
meta.choose_airfoil('ClarkY','0008','0008') #so this is main wing, tail, and vertical stabilizer
#^^This routine also makes a bunch of plots such as CLvsalpha,cdvsalpha,etc.....

#Alright so now let's make some plots that are a function of velocity.
#Close the pdf
meta.pdfhandle.close()

os.system('evince plots.pdf &')
