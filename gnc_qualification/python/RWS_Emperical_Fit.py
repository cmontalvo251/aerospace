#####Import Extra Models#####
import numpy as np
import matplotlib.pyplot as plt
import xlrd

#####Function to get values from excel sheet#####
def cell_contents(sheet,row_x):
    results = []
    for col_x in range(5,sheet.ncols-5):
        cell = sheet.cell(row_x,col_x)
        #print(cell)
        results.append(np.float(cell.value))
    return results

#####Get data from sheets#####
loc = ('GNC_Main.xlsx')
wb = xlrd.open_workbook(loc)
sheet = wb.sheet_by_index(5)

#print(cell_contents(sheet,2))
data = []
for x in range(3,30):
    rows = cell_contents(sheet,x)
    data.append(rows)
#print(data)
data_np = np.asarray(data)
print(data_np[:,0]) #This should match Column F of GNC_Main

storage = data_np[:,3]
volume = data_np[:,0]
mass = data_np[:,2]
power = data_np[:,5]
#print(power)

#####Clip data to desired range#####
threshold_high = 1.0
threshold_low = 1.0E-09

storage_clip = storage[storage<threshold_high]
#print(storage_clip)
volume_clip = volume[storage<threshold_high]
mass_clip = mass[storage<threshold_high]
power_clip = power[storage<threshold_high]
#print(power_clip)
power_fixed = power_clip[threshold_low<power_clip]
#print(power_fixed)
storage_fixed = storage_clip[threshold_low<power_clip]
#print(storage_fixed)

#####Fit Volume Data#####
coeff_volume = np.polyfit(storage_clip,volume_clip,1)
storage_fit = np.linspace(np.min(storage_clip),np.max(storage_clip),100)
volume_fit = np.polyval(coeff_volume,storage_fit)
print('Volume = '+str(coeff_volume[1])+'(Momenum_Storage) + '+str(coeff_volume[0]))

#####Fit Mass Data#####
coeff_mass = np.polyfit(storage_clip,mass_clip,1)
storage_fit = np.linspace(np.min(storage_clip),np.max(storage_clip),100)
mass_fit = np.polyval(coeff_mass,storage_fit)
print('Mass = '+str(coeff_mass[1])+'(Momenum_Storage) + '+str(coeff_mass[0]))

#####Fit Power Data#####
coeff_power = np.polyfit(storage_fixed,power_fixed,1)
storage_fit = np.linspace(np.min(storage_fixed),np.max(storage_fixed),100)
power_fit = np.polyval(coeff_power,storage_fit)
print('Power = '+str(coeff_power[1])+'(Momenum_Storage) + '+str(coeff_power[0]))

plt.figure()
plt.plot(storage_clip,volume_clip,'b*',label='Data')
plt.plot(storage_fit,volume_fit,'r-',label='Empirical Fit')
plt.grid()
plt.xlabel('Momentum Storage (Nms)')
plt.ylabel('Volume (m^3)')
plt.legend()

plt.figure()
plt.plot(storage_clip,mass_clip,'b*',label='Data')
plt.plot(storage_fit,mass_fit,'r-',label='Empirical Fit')
plt.grid()
plt.xlabel('Momentum Storage (Nms)')
plt.ylabel('Mass (kg)')
plt.legend()

plt.figure()
plt.plot(storage_fixed,power_fixed,'b*',label='Data')
plt.plot(storage_fit,power_fit,'r-',label='Empirical Fit')
plt.grid()
plt.xlabel('Momentum Storage (Nms)')
plt.ylabel('Power (W)')
plt.legend()

#####Create the Input Output Tool#####
###########################INPUTS###########################
###6U
Ixx = 0.2159   #kg-m^2
Iyy = 0.1679
Izz = 0.0713
wmax = 10.0 #deg/s
wx = wmax*np.pi/180.0  #rad/s
wy = wmax*np.pi/180.0
wz = wmax*np.pi/180.0
safety_factor = 2.0
#################################################################

#####Make arrays for simple computation#####
I = np.asarray([Ixx,Iyy,Izz])
w = np.asarray([wx,wy,wz])

#####First computation is Angular Momentum#####
H = I*w

#####Then Compute the Required Momentum Storage#####
Hreq = safety_factor*H
print('Momentum Required (Nms) = ',Hreq)

#####Volume Required for Reaction Wheels Assuming the Empirical Fit#####
volume_required = np.polyval(coeff_volume,Hreq)
print('Volume of Reaction Wheel Assuming Empirical Model (m^3) = ',volume_required)

#####Mass Required for Reaction Wheels Assuming the Empirical Fit#####
mass_required = np.polyval(coeff_mass,Hreq)
print('Mass of Reaction Wheel Assuming Empirical Model (kg) = ',mass_required)

#####Power Required for Reaction Wheels Assuming the Empirical Fit#####
power_required = np.polyval(coeff_power,Hreq)
print('Power of Reaction Wheel Assuming Empirical Model (W) = ',power_required)

plt.show()