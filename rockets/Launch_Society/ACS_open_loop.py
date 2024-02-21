# SPDX-FileCopyrightText: 2018 Kattni Rembor for Adafruit Industries
#
# SPDX-License-Identifier: MIT

import time
import board
import pwmio
import digitalio
import adafruit_adxl37x
from analogio import AnalogOut
from analogio import AnalogIn
import adafruit_bmp3xx
#import neopixel
import math
import busio


##Setup the BMP390
i2c = board.I2C()
bmp = adafruit_bmp3xx.BMP3XX_I2C(i2c)
print("Pressure: {:6.1f}".format(bmp.pressure))
print("Temperature: {:5.2f}".format(bmp.temperature))

accel = adafruit_adxl37x.ADXL375(i2c)
print("%f %f %f m/s^2" % accel.acceleration)
#print(dir(adafruit_adxl37x))
#print(adafruit_adxl37x._STANDARD_GRAVITY)
#print(adafruit_adxl37x.DataRate)
#print(adafruit_adxl37x.Range)
#print(adafruit_adxl37x._ADXL347_MULTIPLIER)

xcal = 6.7
ycal = 3.5
zcal = -6.0

#for x in range(0,100):
norm = 0
threshold = 20
while norm < threshold:
    #p = bmp.pressure
    #t = bmp.temperature
    x,y,z = accel.acceleration
    x-=xcal
    y-=ycal
    z-=zcal
    norm = math.sqrt(x*x + y*y + z*z)
    #print(p,t,x,y,z)
    #print((x,y,z))
    print((norm,))
    time.sleep(0.1)

# Initialize PWM output for the servo (on pin A2):
servo1 = pwmio.PWMOut(board.A3, frequency=50)
servo2 = pwmio.PWMOut(board.A1, frequency=50)

# Create a function to simplify setting PWM duty cycle for the servo:
def servo_duty_cycle(x, frequency=50):
    period_ms = 1.0 / frequency * 1000.0
    duty_cycle = int(x / (period_ms / 65535.0))
    return duty_cycle

pulse1 = [2.57,1.15,1.5,1.85,2.2,2.57]
pulse2 = [1.83,0.4,0.7675,1.145,1.523,1.83]

print('Waiting for main motor burnout....')
time.sleep(10)

run = True
# Main loop will run forever moving between 1.0 and 2.0 mS long pulses:
while run==True:
    for x in range(0,6):
        servo1.duty_cycle = servo_duty_cycle(pulse1[x])
        servo2.duty_cycle = servo_duty_cycle(pulse2[x])
        print("S1 = ", pulse1[x],"S2 = ", pulse2[x],x)
        time.sleep(2)
        if x==5:
            run = False


