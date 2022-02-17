###Import 
import numpy as np
import matplotlib.pyplot as plt

###Inputs
##BARFOOT WING
cr = 7.5 ##root chord (use consistent units)
ct = 5.2 ##tip chord
sweep_angle = np.arctan2(8.5,18)*180.0/np.pi ##in degrees
wing_span = 2.0*18.0 ##wingspan

##Random Wing
'''
cr = 8.0
ct = 8.0
wing_span = 18.0
sweep_angle  = 0.0
'''

###Debug Prints
print('Root Chord = ',cr)
print('Tip Chord = ',ct)
print('Sweep Angle (deg) = ',sweep_angle)
print('Wingspan = ',wing_span)
##Compute the Area as well
b1 = cr
b2 = ct
h = wing_span/2.0
S = (b1+b2)*h
print('Area = ',S)

###Create a function to get the chord vs wingspan
def chord(y):
  m = (ct - cr)/(wing_span/2.0)
  return m*y + cr

###Create a function to get the quarter chord line
def quarter_chord_line(y):
  c = chord(y)
  ly = y*np.tan(sweep_angle*np.pi/180.0)
  return ly + c/4.0

###Plot Chord vs Wingspan
y = np.linspace(0,wing_span/2.0,1000)
plt.figure()
plt.plot(y,chord(y))
plt.xlabel('y')
plt.ylabel('Chord Length')
plt.grid()

##Plot Quater chord line as a function of wing_span
plt.figure()
plt.plot(y,quarter_chord_line(y))
plt.xlabel('y')
plt.ylabel('x_(c/4)')
plt.grid()

##Alright now integrate equation 2.29 of Caughey's notes using a Reimann Sum
xMAC = 0.0
dy = y[1]-y[0]
for yi in y:
  xMAC += quarter_chord_line(yi)*chord(yi)*dy
xMAC*=2.0/S
print('xMAC = ',xMAC)

###MEAN AERODYNAMIC CHORD
MAC = 0.0
for yi in y:
    MAC += chord(yi)
MAC/=len(y)
print('MAC = ',MAC)

###Plot the Wing for kicks to make sure it looks right
plt.figure()
plt.plot([0,wing_span/2.0],[0,np.tan(sweep_angle*np.pi/180.0)*wing_span/2.0]) ##Leading edge of right side
plt.plot([0,-wing_span/2.0],[0,np.tan(sweep_angle*np.pi/180.0)*wing_span/2.0]) ##Leading edge of left side
plt.plot([0,wing_span/2.0],[cr,ct+np.tan(sweep_angle*np.pi/180.0)*wing_span/2.0]) ##trailing edge of right side
plt.plot([0,-wing_span/2.0],[cr,ct+np.tan(sweep_angle*np.pi/180.0)*wing_span/2.0]) ##trailing edge of right side
plt.plot([wing_span/2.0,wing_span/2.0],[ct+np.tan(sweep_angle*np.pi/180.0)*wing_span/2.0,np.tan(sweep_angle*np.pi/180.0)*wing_span/2.0]) #end cap right side
plt.plot([-wing_span/2.0,-wing_span/2.0],[ct+np.tan(sweep_angle*np.pi/180.0)*wing_span/2.0,np.tan(sweep_angle*np.pi/180.0)*wing_span/2.0]) #end cap right side
plt.plot([0,wing_span/2.0],[quarter_chord_line(0),quarter_chord_line(wing_span/2.0)]) #quarter chord line of right side
plt.plot([0,-wing_span/2.0],[quarter_chord_line(0),quarter_chord_line(wing_span/2.0)]) #quarter chord line of left side
plt.plot([0,0],[xMAC,xMAC],'b*')
plt.grid()

plt.show()
