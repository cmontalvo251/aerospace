#####Import Extra Modules#####
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import xlrd

#####Function to get values from excel sheet#####
def cell_contents(sheet,row_x):
    results = []
    for col_x in range(6,sheet.ncols-8):
        cell = sheet.cell(row_x,col_x)
        if cell.value == '':
            val = 0.0
        else:
            val = np.float(cell.value)
        results.append(val)
    return results

#####Get data from sheets#####
loc = ('GNC_Main.xlsx')
wb = xlrd.open_workbook(loc)
sheet = wb.sheet_by_index(6)

#####Get Names of all Thrusters#####
names = []
for y in range(3,42):
    for x in range(1,sheet.ncols-20):
        cell = sheet.cell(y,x)
        #print(cell)
        names.append(cell.value)
#print(names)

data = []
for x in range(3,42):
    rows = cell_contents(sheet,x)
    data.append(rows)
#print(data)
data_np = np.asarray(data)
#print(data_np[:,0])

#####Extract Data#####
###Thrust from RCS
Thrust = data_np[:,3]
#print(Thrust.shape)
#print(Thrust)
###Isp
Isp = data_np[:,4]
###Power
Power = data_np[:,7]
###Mass of RCS System
Mass = data_np[:,6]
###Volume of RCS System
Volume = data_np[:,0]
###Max Propellant Mass
Prop_Mass = data_np[:,5]
#print(Prop_Mass)

############################INPUTS#####################################

#Length of Mission (Months)
Mission_Duration = 12.0
#Inertia of Satellite (kg-m^2)
Ixx = 0.2159   #kg-m^2
Iyy = 0.1679
Izz = 0.0713
mass = 13.0 #kg
#I = np.asarray([Ixx,Iyy,Izz])
#Reaction Wheel Storage (Nms) - Calculated with RWS analysis tool
Hreq_X = 0.07536332 
Hreq_Y = 0.05860816 
Hreq_Z = 0.024888
#Hreq = np.asarray([Hreq_X,Hreq_Y,Hreq_Z])
#Factor of Safety - Do we need this since our Reaction Wheel storage has a factor of safety?
FS = 1.0 #This assumes we have just enough propellant

###Orbit
Apogee = 35000.0 #km
Perigee = 185.0 #km (AGL)

##Drag Characteristics 
CD = 1.0

##Size of CubeSat
height = 34.0/100. 
width =  22.6/100. #meters 6U

#######################################################################

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

#####Calculations for Each Axis#####
Thrust_Transpose = np.transpose(Thrust)
#print(Thrust_Transpose.shape)
#print(Thrust_Transpose)
Moment_Arm = np.cbrt(Volume) ##This is a pretty terrible assumption but it's all we can do right now.
#print(Volume_Comp)
Torque_RCS = Thrust_Transpose*Moment_Arm
#print(Torque_RCS)

###X-axis
Angular_Accel_X = Disturbance_Torques_per_axis/Ixx #rad/s
Max_Angular_Velocity_X = Hreq_X/Ixx #rad/s
Time_to_Sat_sec_X = Max_Angular_Velocity_X/Angular_Accel_X #sec
Time_to_Sat_Months_X = ((Time_to_Sat_sec_X/3600.0)/24.0)/30.0 #months
Number_of_Desat_Man_X = Mission_Duration/Time_to_Sat_Months_X
Number_of_Desat_Man_Estimation_X = np.ceil(Number_of_Desat_Man_X)

Accel_RCS_X = Torque_RCS/Ixx
Time_to_Desat_X = Max_Angular_Velocity_X/Accel_RCS_X
Exit_Vel_X = Isp*9.81
Mdot_X = Thrust/Exit_Vel_X
Mass_per_Desat_X = Mdot_X*Time_to_Desat_X
Mass_per_Mission_X = Mass_per_Desat_X*Number_of_Desat_Man_Estimation_X
Mass_per_Mission_X_grams = Mass_per_Mission_X*1000.0
#print(Mass_per_Mission_X_grams)

###Y-axis
Angular_Accel_Y = Disturbance_Torques_per_axis/Iyy #rad/s
Max_Angular_Velocity_Y = Hreq_Y/Iyy #rad/s
Time_to_Sat_sec_Y = Max_Angular_Velocity_Y/Angular_Accel_Y #sec
Time_to_Sat_Months_Y = ((Time_to_Sat_sec_Y/3600.0)/24.0)/30.0 #months
Number_of_Desat_Man_Y = Mission_Duration/Time_to_Sat_Months_Y
Number_of_Desat_Man_Estimation_Y = np.ceil(Number_of_Desat_Man_Y)

Accel_RCS_Y = Torque_RCS/Iyy
Time_to_Desat_Y = Max_Angular_Velocity_Y/Accel_RCS_Y
Exit_Vel_Y = Isp*9.81
Mdot_Y = Thrust/Exit_Vel_Y
Mass_per_Desat_Y = Mdot_Y*Time_to_Desat_Y
Mass_per_Mission_Y = Mass_per_Desat_Y*Number_of_Desat_Man_Estimation_Y
Mass_per_Mission_Y_grams = Mass_per_Mission_Y*1000.0
#print(Mass_per_Mission_Y_grams)

###Z-axis
Angular_Accel_Z = Disturbance_Torques_per_axis/Izz #rad/s
Max_Angular_Velocity_Z = Hreq_Z/Izz #rad/s
Time_to_Sat_sec_Z = Max_Angular_Velocity_Z/Angular_Accel_Z #sec
Time_to_Sat_Months_Z = ((Time_to_Sat_sec_Z/3600.0)/24.0)/30.0 #months
Number_of_Desat_Man_Z = Mission_Duration/Time_to_Sat_Months_Z
Number_of_Desat_Man_Estimation_Z = np.ceil(Number_of_Desat_Man_Z)

Accel_RCS_Z = Torque_RCS/Izz
Time_to_Desat_Z = Max_Angular_Velocity_Z/Accel_RCS_Z
Exit_Vel_Z = Isp*9.81
Mdot_Z = Thrust/Exit_Vel_Z
Mass_per_Desat_Z = Mdot_Z*Time_to_Desat_Z
Mass_per_Mission_Z = Mass_per_Desat_Z*Number_of_Desat_Man_Estimation_Z
Mass_per_Mission_Z_grams = Mass_per_Mission_Z*1000.0
#print(Mass_per_Mission_Z_grams)

###Total Prop Mass for Mission
Tot_Prop_Mass_for_Mission = Mass_per_Mission_X + Mass_per_Mission_Y + Mass_per_Mission_Z
Tot_Prop_Mass_FS = FS*Tot_Prop_Mass_for_Mission
#print(Tot_Prop_Mass_FS)

#Power for each axis
Power_X = Power*2
Power_Y = Power*2
Power_Z = Power*2
Tot_Power = (Power_X + Power_Y + Power_Z)
#print(Tot_Power)

#######################OUTPUTS##########################################

##Which Thrusters are viable given the mass constraint
##If mass required of propellant for all axes combined is greather than max propellant of system, we need to remove it from the output
Table = {'Name': [],'Power X-axis (W)': [],'Power Y-axis (W)': [],'Power Z-axis (W)': [],
         'Total Power (W)': [],'Propellant Mass Used X-axis (kg)': [],
         'Propellant Mass Used Y-axis (kg)': [],
         'Propellant Mass Used Z-axis (kg)': [],
         'Total Propellant Mass Used (kg)': [],
         'Total Propellant Mass (kg)': [],
         'Mass of RCS (kg)': [],
         'Volume of RCS (m^3)': []}

for x in range(0,len(Tot_Prop_Mass_FS)):
  if Tot_Prop_Mass_FS[x] <= Prop_Mass[x]:
    Table['Name'].append(names[x])
    Table['Power X-axis (W)'].append(Power[x])
    Table['Power Y-axis (W)'].append(Power[x])
    Table['Power Z-axis (W)'].append(Power[x])
    Table['Total Power (W)'].append(Tot_Power[x])
    Table['Propellant Mass Used X-axis (kg)'].append(Mass_per_Mission_X[x])
    Table['Propellant Mass Used Y-axis (kg)'].append(Mass_per_Mission_Y[x])
    Table['Propellant Mass Used Z-axis (kg)'].append(Mass_per_Mission_Z[x])
    Table['Total Propellant Mass Used (kg)'].append(Tot_Prop_Mass_FS[x])
    Table['Total Propellant Mass (kg)'].append(Prop_Mass[x])
    Table['Mass of RCS (kg)'].append(Mass[x])
    Table['Volume of RCS (m^3)'].append(Volume[x])

Final_Table = pd.DataFrame(Table)
pd.set_option('display.max_columns', None)
print(Final_Table)

########################################################################