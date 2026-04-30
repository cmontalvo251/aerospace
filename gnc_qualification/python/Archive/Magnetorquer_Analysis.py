#####Import Extra Modules#####
import numpy as np
import matplotlib.pyplot as plt
import xlrd

#####Function to get values from excel sheet#####
def cell_contents(sheet,row_x):
    results = []
    for col_x in range(5,sheet.ncols-6):
        cell = sheet.cell(row_x,col_x)
        #print(cell)
        results.append(np.float(cell.value))
    return results

#####Inputs to Model#####
###Orbit
Apogee = 35000.0 #km
Perigee = 185.0 #km (AGL)

##Drag Characteristics 
CD = 1.0

##Size and mass of CubeSat
height = 34.0/100. 
width =  22.6/100. #meters 6U
mass = 13.0 #kg

#####Compute Disturbance Torques#####

###Equation to compute solar radiation pressure 
rad_pressure = 4.5e-6 # Pa
##use this to compute torque on body

Surface_Area = height*width
solar_rad_force = Surface_Area*rad_pressure
solar_rad_torque = solar_rad_force*np.max([height/2.0,width/2.0])
print('Solar Rad Torque = ',solar_rad_torque)

#Earth Parameters
REarth = 6.371e6 #meters
MEarth = 5.972e24 #kg
G = 6.6742*10**-11; #%%Gravitational constant
mu = G*MEarth

###Gravity Gradient
## We can use G*m/r^2 using AEGIS Backup

##At perigee
accel_top = G*MEarth/(REarth + Perigee*1000.0 + height)**2
accel_bottom = G*MEarth/(REarth + Perigee*1000.0)**2
force_top = accel_top*mass
force_bottom = accel_bottom*mass
delf = force_top-force_bottom
grav_torque_perigee = abs(delf*np.max([height/2.0,width/2.0]))

accel_top = G*MEarth/(REarth + Apogee*1000.0 + height)**2
accel_bottom = G*MEarth/(REarth + Apogee*1000.0)**2
force_top = accel_top*mass
force_bottom = accel_bottom*mass
delf = force_top-force_bottom
grav_torque_apogee = abs(delf*np.max([height/2.0,width/2.0]))

grav_torque = abs(0.5*(grav_torque_perigee + grav_torque_apogee))

print('Gravity Gradient Torque = ',grav_torque)

##Aerodynamic Torque
#ra = Apogee+Rearth
rp = Perigee*1000.0+REarth
#nu = np.linspace(0,2*np.pi,1000)
#a = (ra+rp)/2.
#e = (ra - rp)/(ra+rp)
#p = a*(1-e**2)

## Velocity at perigee - Use Equation 2.5-4 to get velocity
Vp = np.sqrt(MEarth*G/rp)

##Need to compute density as a function of altitude 
## h is the norm of position - Rearth (AGL)
beta = 0.1354 #inverse km
rhos = 1.225 #kg/m^3
rho = rhos*np.exp(-beta*Perigee)

D = 0.5*rho*Vp**2*Surface_Area*CD
aero_torque = D*np.max([height/2.0,width/2.0])/2.0

print('Aerodynamic Torque = ',aero_torque)

#Max Disturbance Torque (Nm)
Disturbance_Torques = aero_torque+solar_rad_torque+grav_torque
Disturbance_Torques_per_axis = Disturbance_Torques/np.sqrt(3)

print('Total Disturbance Torque = ',Disturbance_Torques)
print('Total Disturbance Torque Per Axis = ',Disturbance_Torques_per_axis)

#####Plot Empirical Fit of MagT (This was done using IGRF Model)#####
d = np.linspace(Perigee,Apogee,1000)*1000
r = REarth + d
C = 8.3203e24
B_field_nT = C/r**3
plt.figure()
plt.plot(d,B_field_nT)
plt.ylabel('nT')
plt.xlabel('km')
plt.grid()

#####Get data from sheets#####
loc = ('GNC_Main.xlsx')
wb = xlrd.open_workbook(loc)
sheet = wb.sheet_by_index(4)

#####Get Names of all MagTs#####
names = []
for y in range(3,13):
    for x in range(1,sheet.ncols-14):
        cell = sheet.cell(y,x)
        #print(cell)
        names.append(cell.value)
print(names)

#####Convert Spreadsheet into a numpy array#####
data = []
for x in range(3,13):
    rows = cell_contents(sheet,x)
    data.append(rows)
#print(data)
data_np = np.asarray(data)
#print(data_np[:,0])

#####Extract Relevant Data#####
mag_moment = data_np[:,3]
#print(mag_moment)

plt.figure()
ctr = 0
for moment in mag_moment:
  r = ((moment*C*1e-9)/Disturbance_Torques_per_axis)**(1.0/3.0)
  limit = r - REarth
  print('MagT: ',names[ctr],' Effectiveness Limit (km): ',limit/1000.0)
  Torque = moment*B_field_nT*1e-9
  plt.plot(d/1000.0,Torque,label=names[ctr])
  ctr+=1
plt.xlabel('Distance AGL (km)')
plt.ylabel('Torque (N-m)')
plt.grid()
plt.legend()
plt.plot(d/1000.0,0*d+Disturbance_Torques_per_axis)

plt.show()