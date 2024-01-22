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


###SEND DF COMMAND
df = 0.5
serialPort.write(df)

serialPort.close()
print('Closed Serial Port')
print('Program End')



