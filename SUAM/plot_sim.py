#!/usr/bin/python3

import numpy as np 
import matplotlib.pyplot as pyplot
import sys
sys.path.append('/home/gypsiewanderer13/Dropbox/pdf')
from pdf import *
import mymath as M 
import os
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import cm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import mpl_toolkits.mplot3d as a3
import matplotlib.patches as patch
import copy as CPY

def plottool3(title,x,y,z,xlabel,ylabel,zlabel):
    fig = plt.figure('3-D')
    ax = fig.add_subplot(111,projection='3d')
    #ax.plot_wireframe(x,y,z, color = 'blue', linestyle = 'solid')
    ax.plot(x,y,z)
    plt.title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_zlabel(zlabel)
    return ax #just in case they want access to the ax variable


def get_state_data(inputfilename):
    print('Loading File...')
    statedata = open(inputfilename)
    time = []
    x = []
    y = []
    z = []
    phi = []
    theta = []
    psi = []
    u = []
    v = []
    w = []
    p = []
    q = []
    r = []
    T1 = []
    T2 = []
    T3 = []
    T4 = []
    T5 = []
    T6 = []
    T7 = []
    T8 = []
    T9 = []
    gx = []
    gy = []
    gz = []
    lenfile = 0
    datastate = []

    for line in statedata:
        lenfile+=1
        if len(line) > 2:
            row = line.split(' ')
        else:
            print ('whoops')

        if len(row) > 1:
            #print row
            row_np = []
            for x in row:
                try:
                    val = np.float(x)
                except:
                    pass
                row_np.append(val)
            datastate.append(np.asarray(row_np))

    print('datastate made')
    return np.asarray(datastate)

def create_state_plots(datastate,pp):
        #create sense hat plots
        time = datastate[:,0]
        x = datastate[:,1]
        y = datastate[:,2]
        z = datastate[:,3]
        phi = datastate[:,4] * (180/np.pi)
        theta = datastate[:,5] * (180/np.pi)
        psi = datastate[:,6] * (180/np.pi)
        u =  datastate[:,7]
        v = datastate[:,8]
        w = datastate[:,9]
        p = datastate[:,10]
        q = datastate[:,11]
        r = datastate[:,12]
        gammax = datastate[:,13]
        gammay = datastate[:,14]
        gammaz = datastate[:,15]
        T1 = datastate[:,16]
        T2 = datastate[:,17]
        T3 = datastate[:,18]
        T4 = datastate[:,19]
        T5 = datastate[:,20]
        T6 = datastate[:,21]
        T7 = datastate[:,22]
        T8 = datastate[:,23]
        T9 = datastate[:,24]
        Thrusts = datastate[:,16:24]
        mu1 = datastate[:,25]
        mu2 = datastate[:,26]
        mu3 = datastate[:,27]
        mu4 = datastate[:,28]
        mu5 = datastate[:,29]
        mu6 = datastate[:,30]
        mu7 = datastate[:,31]
        mu8 = datastate[:,32]
        mu9 = datastate[:,33]
        Torques = datastate[:,34:42]

        #plot x
        plt.figure()
        plt.plot(time,x,label = 'x')
        plt.xlabel('Time (s)')
        plt.ylabel('x (ft)')
        plt.grid()
        plt.legend()
        plt.title('x v. t')
        pp.savefig()
        

        #plot y

        plt.figure()
        plt.plot(time,y,label = 'y')
        plt.xlabel('Time (s)')
        plt.ylabel('y (ft)')
        plt.grid()
        plt.legend()
        plt.title('y vs. t')
        pp.savefig()
        #plot z
        plt.figure()
        plt.plot(time,z,label = 'z')
        plt.xlabel('Time (s)')
        plt.ylabel('z (ft)')
        plt.grid()
        plt.legend()
        plt.title('z vs. t')
        pp.savefig()
        #plot phi
        plt.figure()
        plt.plot(time,phi,label = '$\\phi$')
        plt.xlabel('Time (s)')
        plt.ylabel('$\\phi$ (Degrees)') #this way of adding greek letters works!!
        plt.grid()
        plt.legend()
        plt.title('$\\phi$ v. t')
        pp.savefig()

                #plot phi
        plt.figure()
        plt.plot(time,theta,label = '$\\theta$')
        plt.xlabel('Time (s)')
        plt.ylabel('$\\theta$ (Degrees)')
        plt.grid()
        plt.legend()
        plt.title('$\\theta$ v. t')
        pp.savefig()

                #plot phi
        plt.figure()
        plt.plot(time,psi,label = '$\\psi$')
        plt.xlabel('Time (s)')
        plt.ylabel('$\\psi$ (Degrees)')
        plt.grid()
        plt.legend()
        plt.title('$\\psi$ v. t')
        pp.savefig()

                #plot phi
        plt.figure()
        plt.plot(time,u,label = 'u')
        plt.xlabel('Time (s)')
        plt.ylabel('u (ft/s)')
        plt.grid()
        plt.legend()
        plt.title('u v. t')
        pp.savefig()

                #plot phi
        plt.figure()
        plt.plot(time,v,label = 'v')
        plt.xlabel('Time (s)')
        plt.ylabel('v (ft/s)')
        plt.grid()
        plt.legend()
        plt.title('v v. t')
        pp.savefig()

                #plot phi
        plt.figure()
        plt.plot(time,w,label = 'w')
        plt.xlabel('Time (s)')
        plt.ylabel('w (ft/s)')
        plt.grid()
        plt.legend()
        plt.title('w v. t')
        pp.savefig()

                #plot phi
        plt.figure()
        plt.plot(time,p,label = 'p')
        plt.xlabel('Time (s)')
        plt.ylabel('p (rad/s)')
        plt.grid()
        plt.legend()
        plt.title('p v. t')
        pp.savefig()

                #plot phi
        plt.figure()
        plt.plot(time,q,label = 'q')
        plt.xlabel('Time (s)')
        plt.ylabel('q (rad/s)')
        plt.grid()
        plt.legend()
        plt.title('q v. t')
        pp.savefig()

                #plot phi
        plt.figure()
        plt.plot(time,r,label = 'r')
        plt.xlabel('Time (s)')
        plt.ylabel('r (rad/s)')
        plt.grid()
        plt.legend()
        plt.title('r v. t')
        pp.savefig()

                #plot phi
        plt.figure()
        plt.plot(time,T1,label = 'Thrust 1')
        plt.xlabel('Time (s)')
        plt.ylabel('T1(lbf)')
        plt.grid()
        plt.legend()
        plt.title('Thrust 1 v. t')
        pp.savefig()

                #plot phi
        plt.figure()
        plt.plot(time,T2,label = 'Thrust 2')
        plt.xlabel('Time (s)')
        plt.ylabel('T2 (lbf)')
        plt.grid()
        plt.legend()
        plt.title('Thrust 2 v. t')
        pp.savefig()

        #plot phi
        plt.figure()
        plt.plot(time,T3,label = 'Thrust 3')
        plt.xlabel('Time (s)')
        plt.ylabel('T3 (lbf)')
        plt.grid()
        plt.legend()
        plt.title('Thrust 3 v. t')
        pp.savefig()
                #plot phi
        plt.figure()
        plt.plot(time,T4,label = 'Thrust 4')
        plt.xlabel('Time (s)')
        plt.ylabel('T4 (lbf)')
        plt.grid()
        plt.legend()
        plt.title('Thrust 4 v. t')
        pp.savefig()    

        plt.figure()
        plt.plot(time,T5,label = 'Thrust 5')
        plt.xlabel('Time (s)')
        plt.ylabel('T5  (lbf)')
        plt.grid()
        plt.legend()
        plt.title('Thrust 5 v. t')
        pp.savefig()

        plt.figure()
        plt.plot(time,T6,label = 'Thrust 6')
        plt.xlabel('Time (s)')
        plt.ylabel('T6 (lbf)')
        plt.grid()
        plt.legend()
        plt.title('Thrust 6 v. t')
        pp.savefig()

        plt.figure()
        plt.plot(time,T7,label = 'Thrust 7')
        plt.xlabel('Time (s)')
        plt.ylabel('T7  (lbf)')
        plt.grid()
        plt.legend()
        plt.title('Thrust 7 v. t')
        pp.savefig()

        plt.figure()
        plt.plot(time,T8,label = 'Thrust 8')
        plt.xlabel('Time (s)')
        plt.ylabel('T8  (lbf)')
        plt.grid()
        plt.legend()
        plt.title('Thrust 8 v. t')
        pp.savefig()

        plt.figure()
        plt.plot(time,T9,label = 'Thrust 9')
        plt.xlabel('Time (s)')
        plt.ylabel('T9 (lbf)')
        plt.grid()
        plt.legend()
        plt.title('Thrust 9 v. t')
        pp.savefig()

        ##Plot mu1-8
        plt.figure()
        plt.plot(time,mu1,label='1')
        plt.plot(time,mu2,label='2')
        plt.plot(time,mu3,label='3')
        plt.plot(time,mu4,label='4')
        plt.plot(time,mu5,label='5')
        plt.plot(time,mu6,label='6')
        plt.plot(time,mu7,label='7')
        plt.plot(time,mu8,label='8')
        plt.grid()
        plt.legend()
        plt.xlabel('Time (s)')
        plt.ylabel('PWM Signal (us)')
        plt.title('Control Signal')
        pp.savefig()

        plt.figure()
        plt.plot(time,mu9)
        plt.grid()
        plt.xlabel('Time (s)')
        plt.ylabel('PWM Signal (us)')
        plt.title('Pusher Signal')
        pp.savefig()

        ##Plot delta mu1-8
        mulist = [mu1,mu2,mu3,mu4,mu5,mu6,mu7,mu8]
        muavg = 0*mu1
        nummotorsrunning = 0.0
        for mu in mulist:
            mucopy = CPY.deepcopy(mu)
            mucopy[mucopy<1100] = 0
            if sum(mucopy) != 0:
                muavg += mu
                nummotorsrunning+=1
        muavg /= nummotorsrunning
        plt.figure()
        munumber = 0
        for mu in mulist:
            munumber+=1
            mucopy = CPY.deepcopy(mu)
            mucopy[mucopy<1100] = 0
            if sum(mucopy) != 0:
                plt.plot(time,mu-muavg,label=str(munumber))
        #plt.plot(time,mu2-muavg,label='2')
        #plt.plot(time,mu3-muavg,label='3')
        #plt.plot(time,mu4-muavg,label='4')
        #plt.plot(time,mu5-muavg,label='5')
        #plt.plot(time,mu6-muavg,label='6')
        #plt.plot(time,mu7-muavg,label='7')
        #plt.plot(time,mu8-muavg,label='8')
        plt.grid()
        plt.legend()
        plt.xlabel('Time (s)')
        plt.ylabel('Delta PWM Signal (us)')
        plt.title('Control Signal')
        pp.savefig()

        plt.figure()
        munumber = 0
        r,c = np.shape(Thrusts)
        for idx in range(0,c):
            munumber+=1
            plt.plot(time,Thrusts[:,idx],label=str(munumber))
        plt.grid()
        plt.legend()
        plt.xlabel('Time (s)')
        plt.ylabel('Thrust (N)')
        plt.title('All thrust')
        pp.savefig()

        plt.figure()
        munumber = 0
        r,c = np.shape(Torques)
        for idx in range(0,c):
            munumber+=1
            plt.plot(time,Torques[:,idx],label=str(munumber))
        plt.grid()
        plt.legend()
        plt.xlabel('Time (s)')
        plt.ylabel('Torques (N-m)')
        plt.title('All Torques')
        pp.savefig()

        z3 = np.asarray([-z,time])
        ##3D plots
        #plt.figure()
        plottool3('Flight Path',x,y,-z,'x','y','z')
        pp.savefig()
 
if __name__ == '__main__':

        print('Processing State output File')

        inputfilename = 'Output_Files/C++.OUT'

        SHOWPLOTS = 0
        pp = PDF(SHOWPLOTS,plt)

        data_all = []
        #inputfilename = inputfilenames[0]
        data = get_state_data(inputfilename)
        #data_all.append(data)
        
        #datastate = data_all
        #print len(datastate[0])

        create_state_plots(data,pp)

        pp.close()
