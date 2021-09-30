#!/usr/bin/python

import numpy as np 
import matplotlib.pyplot as plt
import sys
import os
from matplotlib.backends.backend_pdf import PdfPages
#from matplotlib.ticker import FormatStrFormatter

##Let's make a class to utilize object oriented programming
class Vehicle():
    def __init__(self):
        #Set up Simdata
        T0 = 0.0
        TFINAL = 10
        TIMESTEP = 0.01
        self.simdata = np.asarray([T0,TFINAL,TIMESTEP])
        ##Set up the massdata
        self.mass = 5.6
        Ixx = 1.0
        Iyy = 2.0
        Izz = 3.0
        self.I = np.asarray([[Ixx,0,0],[0,Iyy,0],[0,0,Izz]])
        self.Iinv = np.linalg.inv(self.I)
        ##Set up initial conditions
        x0 = 0.0
        y0 = 0.0
        z0 = -100
        phi0 = 0.0
        theta0 = 0.0
        psi0 = 0.0
        u0 = 20.0
        v0 = 0.0
        w0 = 0.0
        p0 = 0.0
        q0 = 0.0
        r0 = 0.0
        self.stateinitial = np.asarray([x0,y0,z0,phi0,theta0,psi0,u0,v0,w0,p0,q0,r0])
        
    def ForceMoment(self,t,state):
        #Force model is just 9.81 m/s^2 for now.
        Fgrav = self.mass*np.asarray([0,0,9.81])
        Faero = np.asarray([0,0,-100])
        F = Fgrav + Faero
        #Moment is zero for now
        M = np.asarray([0,0,0])
        return F,M
        
    def Derivatives(self,t,state):
        #Need to compute statedot
        #x = state[0]
        #y = state[1]
        #z = state[2]
        phi = state[3]
        theta = state[4]
        psi = state[5]
        u = state[6]
        v = state[7]
        w = state[8]
        p = state[9]
        q = state[10]
        r = state[11]

        #Set up vectors
        ptp = np.asarray([phi,theta,psi])
        uvw = np.asarray([u,v,w])
        pqr = np.asarray([p,q,r])
        
        #Kinematics
        ctheta = np.cos(theta)
        cpsi = np.cos(psi)
        sphi = np.sin(phi)
        stheta = np.sin(theta)
        cphi = np.cos(phi)
        spsi = np.sin(psi)
        ttheta = np.tan(theta)
        TIB = np.asarray([[ctheta*cpsi,sphi*stheta*cpsi-cphi*spsi,cphi*stheta*cpsi+sphi*spsi],[ctheta*spsi,sphi*stheta*spsi+cphi*cpsi,cphi*stheta*spsi-sphi*cpsi],[-stheta,sphi*ctheta,cphi*ctheta]])
        xyzdot = np.matmul(TIB,uvw)
        H = np.asarray([[1.0,sphi*ttheta,cphi*ttheta],[0.0,cphi,-sphi],[0.0,sphi/ctheta,cphi/ctheta]])
        ptpdot = np.matmul(H,pqr)
        
        #Force and Moment Model 
        F,M = self.ForceMoment(t,state) ##THESE ARE IN THE INERTIAL FRAME!!!
        
        #Dynamics
        uvwdot = np.matmul(np.transpose(TIB),F)/self.mass - np.cross(pqr,uvw)
        pqrdot = np.matmul(self.Iinv,M-np.cross(pqr,np.matmul(self.I,pqr)))

        dxdt = np.concatenate([xyzdot,ptpdot,uvwdot,pqrdot])

        return dxdt

    def Integrate(self):
        #In order to run the RK4 engine we need the simdata block
        #to create a time vector from the simdata block
        tinitial = self.simdata[0]
        tfinal = self.simdata[1]
        timestep = self.simdata[2]
        #Then we can run the RK4 engine
        t = tinitial

        #Although scripting languages like python can plot natively. I will generate an outputfile so if anyone wants
        #to plot in a different language they can.
        outfile = open('Python.OUT','w')

        #Run the integrator
        print('Begin RK4 Integrator')
        state = self.stateinitial
        while t <= tfinal:
            #Output Contents to File
            list = state.tolist()
            s = ", ".join(map(str,list))
            outfile.write(str(t)+','+s+'\n')
            #RK4 Call
            k1 = self.Derivatives(t,state)
            k2 = self.Derivatives(t+timestep/2.0,state+k1*timestep/2.0)
            k3 = self.Derivatives(t+timestep/2.0,state+k2*timestep/2.0)
            k4 = self.Derivatives(t+timestep,state+k3*timestep)
            phi = (1.0/6.0)*(k1 + 2*k2 + 2*k3 + k4)
            #Step State
            state += phi*timestep
            t+=timestep
            print('Time =',t)
        outfile.close()
        print('RK4 Integration Complete')

def euler2quat(ptp):
    #%%%Input is a 3x1 vector and output is a 4x1 vector

    phi = ptp[0];
    theta = ptp[1]
    psi = ptp[2]

    q0 = np.cos(phi/2)*np.cos(theta/2)*np.cos(psi/2) + np.sin(phi/2)*np.sin(theta/2)*np.sin(psi/2);
    q1 = np.sin(phi/2)*np.cos(theta/2)*np.cos(psi/2) - np.cos(phi/2)*np.sin(theta/2)*np.sin(psi/2);
    q2 = np.cos(phi/2)*np.sin(theta/2)*np.cos(psi/2) + np.sin(phi/2)*np.cos(theta/2)*np.sin(psi/2);
    q3 = np.cos(phi/2)*np.cos(theta/2)*np.sin(psi/2) - np.sin(phi/2)*np.sin(theta/2)*np.cos(psi/2);

    return np.asarray([q0,q1,q2,q3]);
        
def quat2euler(q0123):
    q0 = q0123[0]
    q1 = q0123[1]
    q2 = q0123[2]
    q3 = q0123[3]

    phi = (np.arctan2(2*(q0*q1 + q2*q3),1-2*(q1**2 + q2**2)))
    theta = np.arcsin(2*(q0*q2-q3*q1))
    psi = np.arctan2(2*(q0*q3 + q1*q2),1-2*(q2**2 + q3**2))

    return np.asarray([phi,theta,psi])

def extract_Euler(T):
    #%%%Assuming R is a 3x3 matrix extract phi,theta,psi Euler angles
    #%%%assuming a 3-2-1 transformation sequence
    #%%%Let R be defined such that v(body) = T v(inertial)
    #Using Mark Costello's notation this would be TBI

    theta = -np.arcsin(T[0][2]);
    sphi  = T[1][2]/np.cos(theta);
    cphi = T[2][2]/np.cos(theta);
    phi = np.arctan2(sphi,cphi);
    spsi = T[0][1]/np.cos(theta);
    cpsi = T[0][0]/np.cos(theta);
    psi = np.arctan2(spsi,cpsi);

    ptp = np.zeros(3)

    ptp[0] = phi
    ptp[1] = theta
    ptp[2] = psi

    return ptp

def extract_quaternion(TBI):
    #%%%Assuming T is a 3x3 matrix, extract the quaternion vector (q0,q1,q2,q3)
    #%%%assuming a 3-2-1 transformation sequence
    #%%%Let T be defined such that v(body) = T v(inertial)
    #%%In Mark Costello's notation, T would TBI
    alfa = TBI[0,1] + TBI[1,0]
    bata = TBI[2,0] + TBI[0,2]
    gama = TBI[2,1] + TBI[1,2]
    q1squared = alfa*bata/(4.0*gama);

    #%%%Solutions split into two different solutions here
    q1a = np.sqrt(q1squared);
    q1b = -np.sqrt(q1squared);
    q2a = alfa/(4.0*q1a);
    q2b = alfa/(4.0*q1b);
    q3a = bata/(4.0*q1a);
    q3b = bata/(4.0*q1b);

    #%%%These are the same so just pick one
    #%q0asquared = TBI(1,1) + q2a^2 + q3a^2 - q1a^2
    #%q0bsquared = TBI(1,1) + q2b^2 + q3b^2 - q1b^2
    q0squared = TBI[0,0] + q2a**2 + q3a**2 - q1a**2

    #%%%Solution However still splits into 4 possible solutions
    #You can get around this though by enforcing q0a to be positive
    #See this article here %%http://planning.cs.uiuc.edu/node151.html
    q0a = np.sqrt(q0squared);
    #q0b = -sqrt(q0squared);

    #%%%Here are my 2 possible solutions
    q0123aa = np.asarray([q0a,q1a,q2a,q3a])
    q0123ab = np.asarray([q0a,q1b,q2b,q3b])
    
    #%According to this website there are multiple solutions 
    #%that yield the same euler angles. 
    #%http://planning.cs.uiuc.edu/node151.html
    #%So what we want to do is compute the euler angles from the matrix
    ptp = extract_Euler(TBI)
    #ptp = [phi,theta,psi];
    #%Then check and see which have the same euler angles
    #because you restricted q0 to be positive above though.
    #This will only result in one unique solution from the 2.
    #Thus you can break out of this loop as soon as you find it.
    quats = [q0123aa,q0123ab]
    for q0123j in quats:
        ptpj = quat2euler(q0123j)
        val = abs(sum(ptpj-ptp))
        if val < 1e-10:
            return q0123j

def RQUAT(q0123):
    #%compute R such that v(inertial) = R v(body)
    #%Using Mark Costello's notation this would be TIB

    q0 = q0123[0]
    q1 = q0123[1]
    q2 = q0123[2]
    q3 = q0123[3]

    R = np.asarray([[q0**2+q1**2-q2**2-q3**2,2*(q1*q2-q0*q3),2*(q0*q2+q1*q3)],[2*(q1*q2+q0*q3),(q0**2-q1**2+q2**2-q3**2),2*(q2*q3-q0*q1)],[2*(q1*q3-q0*q2),2*(q0*q1+q2*q3),q0**2-q1**2-q2**2+q3**2]])

    return R

def R123(phi,theta,psi):
    #%compute R such that v(inertial) = R v(body)
    #%Compute sines and cosines
    ctheta = np.cos(theta);
    stheta = np.sin(theta);
    sphi = np.sin(phi);
    cphi = np.cos(phi);
    spsi = np.sin(psi);
    cpsi = np.cos(psi);
    #%Kinematics
    R = np.array([[ctheta*cpsi,sphi*stheta*cpsi-cphi*spsi,cphi*stheta*cpsi+sphi*spsi],[ctheta*spsi,sphi*stheta*spsi+cphi*cpsi,cphi*stheta*spsi-sphi*cpsi],[-stheta,sphi*ctheta,cphi*ctheta]]);
    return R

if __name__ == '__main__':

    #Create a Vehicle class which also sets all the standard parameters
    mav = Vehicle()

    #Here we go ahead and integrate the equations of motion
    mav.Integrate()

    print('Python Module Complete')

    #Once the integration is complete. It's time to read the entire file and plot all states
    #Although it's possible to simply create vectors like in MATLAB I elected to do it this way
    #to emulate how it would be done in C++ or Fortran.
    file = open('Python.OUT','r')
    data = []
    for line in file:
        row = line.split(',')
        row_np = []
        for r in row:
            row_np.append(np.float(r))
        data.append(row_np)

    data_np = np.asarray(data)
    time = data_np[:,0]
    size = np.shape(data_np)
    NOSTATES = size[1]

    #sys.exit()

    ##Save Figures natively to PDF
    print('Generating Plots')

    pdfhandle = PdfPages('python_plots.pdf')
    ylabels = ['Time (sec)','X (m)','Y (m)','Z (m)','PHI (rad)','THETA (rad)','PSI (rad)','U (m/s)','V (m/s)','W (m/s)','P (rad/s)','Q (rad/s)','R (rad/s)']
    for idx in range(0,NOSTATES):
        plt.figure()
        #print(time)
        #print(data_np[:,idx])
        plt.plot(time,data_np[:,idx])
        plt.grid()
        plt.xlabel('Time (sec)')
        plt.ylabel(ylabels[idx])
        pdfhandle.savefig()
        print(ylabels[idx])
    
    pdfhandle.close()
    print('Plotting Routine Complete for Python')
    os.system('evince python_plots.pdf &')
    #plt.show()

# Copyright - Carlos Montalvo 2016
# You may freely distribute this file but please keep my name in here
# as the original owner
