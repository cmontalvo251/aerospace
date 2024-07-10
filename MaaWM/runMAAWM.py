#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import csv
import sys
from gps import *
from pdf import *
from plotting import * #This is in BlackBox.git
import pitot as P

class AERIAL():
    def __init__(self):
        self.t = []
        self.x = []
        self.y = []
        self.z = []
        self.phi = []
        self.theta = []
        self.psi = []
        self.u = []
        self.v = []
        self.w = []

#Figure out how many quads and airplanes there are
ics = open('Input_Files/Aircraft.ICS')
num_aircraft = int(ics.readline().split()[0])
num_airplanes = int(ics.readline().split()[0])
num_quads = num_aircraft - num_airplanes

pp = PDF(0,plt)

aircraft = []
for a in range(0,num_aircraft):
    aircraft.append(AERIAL())

#Read in state.out file
print 'Reading State.OUT File'
fid = open('Output_Files/State.OUT')
s = 0
nostates = 12
for line in fid:
    states = line.split()
    s = 0
    for a in range(0,num_aircraft):
        aircraft[a].t.append(float(states[0]))
        aircraft[a].x.append(float(states[s+1]))
        aircraft[a].y.append(float(states[s+2]))
        aircraft[a].z.append(float(states[s+3]))
        aircraft[a].phi.append(float(states[s+4]))
        aircraft[a].theta.append(float(states[s+5]))
        aircraft[a].psi.append(float(states[s+6]))
        aircraft[a].u.append(float(states[s+7]))
        aircraft[a].v.append(float(states[s+8]))
        aircraft[a].w.append(float(states[s+9]))
        s += nostates

#Convert all to numpy arrays
for a in range(0,num_aircraft):
    print 'Aircraft = ',a
    aircraft[a].t = np.asarray(aircraft[a].t)
    aircraft[a].x = np.asarray(aircraft[a].x)
    aircraft[a].y = np.asarray(aircraft[a].y)
    aircraft[a].z = np.asarray(aircraft[a].z)
    aircraft[a].phi = np.asarray(aircraft[a].phi)
    aircraft[a].theta = np.asarray(aircraft[a].theta)
    aircraft[a].psi = np.asarray(aircraft[a].psi)
    
    #And plot everything
    print 'X'
    plti = plottool(12,'Time (sec)','Range (ft)','X Coordinate')
    plti.plot(aircraft[a].t,aircraft[a].x)
    plt.gcf().subplots_adjust(left=0.15)
    pp.savefig()
    
    print 'Y'
    plti = plottool(12,'Time (sec)','CrossRange (ft)','Y Coordinate')
    plti.plot(aircraft[a].t,aircraft[a].y)
    plt.gcf().subplots_adjust(left=0.15)
    pp.savefig()
    
    print 'XvsY'
    plti = plottool(12,'Range (ft)','CrossRange (ft)','X vs Y Coordinate')
    plti.plot(aircraft[a].x,aircraft[a].y)
    plt.gcf().subplots_adjust(left=0.15)
    pp.savefig()
    
    print 'Z'
    plti = plottool(12,'Time (sec)','Altitude (ft)','Z Coordinate')
    plti.plot(aircraft[a].t,aircraft[a].z)
    #print aircraft[a].z
    plt.gcf().subplots_adjust(left=0.15)
    pp.savefig()
    
    print 'RP'
    plti = plottool(12,'Time (sec)','Roll Pitch (rad)','RP')
    plti.plot(aircraft[a].t,aircraft[a].phi,label='Roll')
    plti.plot(aircraft[a].t,aircraft[a].theta,label='Theta')
    plti.legend()
    plt.gcf().subplots_adjust(left=0.15)
    pp.savefig()

    print 'Yaw'
    plti = plottool(12,'Time (sec)','Yaw (rad)','Yaw')
    plti.plot(aircraft[a].t,aircraft[a].psi)
    plt.gcf().subplots_adjust(left=0.15)
    pp.savefig()

    print 'U'
    plti = plottool(12,'Time (sec)','Forward Speed (ft/s)','U')
    plti.plot(aircraft[a].t,aircraft[a].u)
    plt.gcf().subplots_adjust(left=0.15)
    pp.savefig()

    print 'VvsW'
    plti = plottool(12,'Time (sec)','Side and Vertical Velocity (ft/s)','V and W')
    plti.plot(aircraft[a].t,aircraft[a].v,label='v')
    plti.plot(aircraft[a].t,aircraft[a].w,label='w')
    plti.legend()
    plt.gcf().subplots_adjust(left=0.15)
    pp.savefig()

print 'TopDown Overall'
plti = plottool(12,'Range (ft)','CrossRange (ft)','X vs Y Coordinate')
for a in range(0,num_aircraft):
    plti.plot(aircraft[a].x,aircraft[a].y)
plt.gcf().subplots_adjust(left=0.15)
pp.savefig()

#Open the Winds File
file = open('Output_Files/Winds.OUT')
t = []
x = []
y = []
z = []
vx = []
vy = []
vz = []
print 'Reading Winds.OUT'
for line in file:
    row = line.split()
    x.append(float(row[0]))
    y.append(float(row[1]))
    z.append(float(row[2]))
    t.append(float(row[3]))
    vx.append(float(row[4]))
    vy.append(float(row[5]))
    vz.append(float(row[6]))

print 'Vx'
scatter(x,y,vx,'Range (ft)','CrossRange (ft)','Vx (ft/s)','Measured Winds')
pp.savefig()

print 'Vy'
scatter(x,y,vy,'Range (ft)','CrossRange (ft)','Vy (ft/s)','Measured Winds')
pp.savefig()

print 'Vz'
scatter(x,y,vz,'Range (ft)','CrossRange (ft)','Vz (ft/s)','Measured Winds')
pp.savefig()


pp.close()

    


