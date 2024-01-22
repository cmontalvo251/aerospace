#!/usr/bin/python3

import serial
import sys
import time
import numpy as np

print('Attempting to open Serial port')
try:
    serialPort = serial.Serial("COM11",115200)
    print('If no errors port opened')
except IOError:
    print('Error opening port exiting')
    sys.exit()

###SEND DF COMMAND
df = 0.0
while True:
    df += 0.5
    if df > 1.0:
        df = 0.0
    print('Sending df command = ',df)
    #serialPort.write(str(df).encode('ascii')  )
    time.sleep(1.0)
    output = ''
    while True:
        a = serialPort.read()
        output += str(a).strip('b').strip("'")
        if a == b'\n':
            print(output)
            output = ''
    
##THIS IS TO READ
#while True:
#    a = serialPort.read()
#    print(a)
#    time.sleep(0.1)
    
serialPort.close()
print('Closed Serial Port')
print('Program End')



