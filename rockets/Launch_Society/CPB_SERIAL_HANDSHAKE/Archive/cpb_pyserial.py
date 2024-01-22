#!/usr/bin/python3

import serial
import sys
import time
import numpy as np

print('Attempting to open Serial port')
try:
    serialPort = serial.Serial("/dev/ttyACM0",115200)
    print('If no errors port opened')
except IOError:
    print('Error opening port exiting')
    sys.exit()

def convert(input):
    input = input.strip('b')
    input = input.strip("'")
    input = input.strip('\\n')
    input = input.rstrip('\\r')
    return np.float(input)

###DATA COLLECTION BRANCH
norm = 0
while norm < 25:
    ##Read
    inBuffer = serialPort.readline()
    instr = str(inBuffer)
    buffer_list = instr.split(' ')
    x = convert(buffer_list[0])
    y = convert(buffer_list[1])
    z = convert(buffer_list[2])
    norm = np.sqrt(x**2 + y**2 + z**2)
    print(x,y,z,norm)
    time.sleep(0.1)

print('Large Accel Spike Detected')
serialPort.close()
print('Closed Serial Port')
print('Program End')



