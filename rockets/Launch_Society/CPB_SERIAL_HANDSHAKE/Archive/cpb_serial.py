import serial
import sys
import time
import numpy as np

class CPX:
##############################################################################
#
def __init__(self,portname,index):
print('Attempting to open serial port')
self.index = index
self.filename = 'CPX'+str(index)+'.txt'
self.file = open(self.filename,'w')
try:
self.serialPort = serial.Serial(portname,115200)
print('if no errors port opened')
print('portname: ',portname)
except IOError:
print('Error opening serial port')
#print('portname: ',portname)
sys.exit()

81

self.norm = 0

##############################################################################
####
def convert(self,input):
input = input.strip('b')
input = input.strip("'")
input = input.strip('\\n')
input = input.strip('\\r')
return np.float(input)

##############################################################################
#####
def read(self):
inBuffer = self.serialPort.readline()
instr = str(inBuffer)
buffer_list = instr.split(' ')
#print(self.index,":",buffer_list)
self.t = self.convert(buffer_list[0])
self.x = self.convert(buffer_list[1])
self.y = self.convert(buffer_list[2])
self.z = self.convert(buffer_list[3])
self.ax = self.convert(buffer_list[4])
self.ay = self.convert(buffer_list[5])
self.az = self.convert(buffer_list[6])
self.gx = self.convert(buffer_list[7])

82

self.gy = self.convert(buffer_list[8])
self.gz = self.convert(buffer_list[9])

##########################################################################

self.norm1 = np.sqrt(self.x**2 + self.y**2 + self.z**2)
self.norm2 = np.sqrt(self.ax**2 + self.ay**2 + self.az**2)
self.norm = (self.norm1 + self.norm2)/2

##########################################################################
#print('CPX = ',self.index,":",self.x, self.y, self.z, self.ax, self.ay, self.az,

self.gx,self.gy,self.gz,'advrage = ',self.norm)

####EDIT THIS PART TO OUTPUT AXYZ and GXYZ ########### DONE
outstring = str(self.t) + " " + str(self.x) + " " + str(self.y) + " " + str(self.z) + " " +
str(self.ax) + " " + str(self.ay) + " " + str(self.az) + " " + str(self.gx) + " " + str(self.gy) + " " +
str(self.gz) + "\n"

self.file.write(outstring)
self.file.flush()
def extractdata(self):
###EDIT THIS TO EXTRACT ALL DATA STREAMS #### DONE
data = np.loadtxt(self.filename,delimiter=' ')
self.tvec = data[:,0]
self.xvec = data[:,1]
self.yvec = data[:,2]
self.zvec = data[:,3]

83

self.axvec = data[:,4]
self.ayvec = data[:,5]
self.azvec = data[:,6]
self.gxvec = data[:,7]
self.gyvec = data[:,8]
self.gzvec = data[:,9]
self.tvec -= self.tvec[0]

##############################################################################
####
def close(self):
self.serialPort.close()
self.file.close()

##############################################################################
####
##############################################################################
####
time.sleep(30)
startTime = time.time()
CPX_list = []
counter = 0
for x in range(25):
try:
addr = "/dev/ttyACM" + str(x)

84

cpx = CPX(addr,counter)
print(counter)
CPX_list.append(cpx)
counter += 1
except:
pass
CPX1 = CPX_list[0]
CPX2 = CPX_list[1]
CPX3 = CPX_list[2]
CPX4 = CPX_list[3]
#CPX5 = CPX_list[4]
#CPX6 = CPX_list[5]
#CPX7 = CPX_list[6]

##########################################
#CPX1.norm = 0
#spike = CPX1.norm
counter_new = 0
while True:
##########################################
try:
CPX1.read()
spike = CPX1.norm

85

except:
#pass
spike = CPX2.norm
try:
CPX2.read()
except:
#pass
spike = CPX3.norm
try:
CPX3.read()
except:
#pass
spike = CPX4.norm
try:
CPX4.read()
except:
pass
# spike = CPX5.norm
# try:
# CPX5.read()
# except:
#pass
# spike = CPX6.norm

86

# try:
# CPX6.read()
# except:
# pass
#spike = CPX7.norm
#try:
#CPX7.read()
#except:
#break
elapsedTime = time.time() - startTime
if(spike > 15):
break
if(elapsedTime > 10):
startTime = time.time()
print('I am working')
filename = 'CPX1/' + 'CPX' + str(counter_new) + '.txt'
CPX1.file.close()
CPX1.file = open(filename,'w')
filename = 'CPX2/' + 'CPX' + str(counter_new) + '.txt'
CPX2.file.close()
CPX2.file = open(filename,'w')
filename = 'CPX3/' + 'CPX' + str(counter_new) + '.txt'
CPX3.file.close()

87

CPX3.file = open(filename,'w')
filename = 'CPX4/' + 'CPX' + str(counter_new) + '.txt'
CPX4.file.close()
CPX4.file = open(filename,'w')
counter_new += 1

#print(spike,CPX1.norm,CPX2.norm,CPX3.norm,CPX4.norm,CPX5.norm,CPX6.norm)
time.sleep(0.01)
##############################################################################
####
print('spike detected')
CPX1.close()
CPX2.close()
CPX3.close()
CPX4.close()
#CPX5.close()
#CPX6.close()
#CPX7.close()

print('Closed serial')
print('program end')

CPX1.extractdata()
CPX2.extractdata()

88

CPX3.extractdata()
CPX4.extractdata()
#CPX5.extractdata()
#CPX6.extractdata()
#CPX7.extractdata()