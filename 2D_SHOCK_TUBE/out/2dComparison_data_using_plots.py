'''Results visualization of shock tube case 1 '''

import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import numpy as np
import math 
import csv
# To run this plot file we need the results from SHOCK_TUBE_case1_L_0.5_0.2 case.
# Inside the out folder the exact Riemann solver solutions from the https://cococubed.com/code_pages/exact_riemann.shtml and the .out files can be found.
# The simualted data extracted from the line y = ypoints/2

rho2  = []
vel2x = []
vel2y = []
temp2 = []

# Read exact results from output1273 to the following list
simulated_out_name = "ST_2_"
exact_results_file = "output310"
exactpoints      = []
densityexact     = []
temperatureexact = []
velexact         = []
presexact        = []
nr               = []

points     = 600            # Number of points along x-axis
total_time = 310            # defined simulation time
ypoints    = 4               # Number of points along Y-axis        
ITERATION  = str(total_time) # defined simulation time
sphr       = 1.4           # specif heat ratio when calculating the Mach number

# Start reading exact results files
with open(exact_results_file,"r") as file:
    lines = file.readlines()[2:]
    for i in lines:
        p = i.split()
        exactpoints     .append(float(p[0]))
        densityexact    .append(float(p[2]))
        temperatureexact.append(float(p[6]))
        presexact       .append(float(p[3]))
        velexact        .append(float(p[4]))

# Start reading simulated results files
with open( simulated_out_name + ITERATION + ".out","r") as file:
    lines = file.readlines()[1:]
    for i in lines:
        p = i.split()
        # print(p)
        rho2 .append(float(p[0]))
        vel2x.append(float(p[1]))
        vel2y.append(float(p[2]))
        temp2.append(float(p[3]))
        nr.append(float(p[4]))

# Creating numpy arrays and reshaping  for ease of data handling
rho2  = np.array(rho2)
vel2x = np.array(vel2x)
vel2y = np.array(vel2y)
T2    = np.array(temp2)
nr    = np.array(nr)

vel2x = np.reshape(vel2x,(ypoints,points))
vel2y = np.reshape(vel2y,(ypoints,points))
rho2  = np.reshape(rho2,(ypoints,points))
T2    = np.reshape(T2,(ypoints,points))
nr    = np.reshape(nr,(ypoints,points))

# Points along x-axis
x = np.linspace(0, points, points, endpoint=False) 
# fig1, [[ax1,ax2],[ax3,ax4]] = plt.subplots(2,2)
fig1, ([ax1,ax3]) = plt.subplots(1,2)
fig2, ([ax4,ax2]) = plt.subplots(1,2)

# Compare temperature
x = x/points
ax4.plot(x,T2[ypoints//2,:],label='2D CLBM',color='k',linestyle='-')
ax4.plot(x,temperatureexact,label='Exact', color='red',linestyle='--')
# ax4.grid()
ax4.grid(which='major', color='#DDDDDD', linewidth=0.5)
ax4.grid(which='minor', color='#EEEEEE', linewidth=0.25)
ax4.minorticks_on()
ax4.set_ylabel(r'T',fontsize = 18)
ax4.set_xlabel(r'x/L',fontsize = 18)
ax4.tick_params(axis='x', labelsize=15)
ax4.tick_params(axis='y', labelsize=15)
ax4.legend(bbox_to_anchor =(0.75, 1.15), ncol = 2)

# Compare density
ax2.plot(x,rho2[ypoints//2,:],label='2D CLBM',color='k',linestyle='-')
# ax2.grid()
ax2.plot(x,densityexact, color='red',linestyle='--')
ax2.grid(which='major', color='#DDDDDD', linewidth=0.5)
ax2.grid(which='minor', color='#EEEEEE', linewidth=0.25)
ax2.minorticks_on()
ax2.set_ylabel(r'$\rho$',fontsize = 18)
ax2.set_xlabel(r'x/L',fontsize = 18)
ax2.tick_params(axis='x', labelsize=15)
ax2.tick_params(axis='y', labelsize=15)

# Compare velocity
ax3.plot(x,vel2x[ypoints//2,:],label='2D CLBM',color='k',linestyle='-')
# ax3.grid()
ax3.plot(x,velexact, color='red',linestyle='--')
ax3.grid(which='major', color='#DDDDDD', linewidth=0.5)
ax3.grid(which='minor', color='#EEEEEE', linewidth=0.25)
ax3.minorticks_on()
ax3.set_ylabel(r'$u_x$',fontsize = 18)
ax3.set_xlabel(r'x/L',fontsize = 18)
ax3.tick_params(axis='x', labelsize=15)
ax3.tick_params(axis='y', labelsize=15)


# Compare pressure
pressure2 = rho2[ypoints//2,:]*T2[ypoints//2,:]                             # Calculating the Pressure using the state equation mentioned in the thesis text
ax1.plot(x,pressure2,label='2D CLBM',color='k',linestyle='-')
# ax1.grid()
ax1.plot(x,presexact, color='red',linestyle='--',label="Exact")
ax1.grid(which='major', color='#DDDDDD', linewidth=0.5)
ax1.grid(which='minor', color='#EEEEEE', linewidth=0.25)
ax1.minorticks_on()
ax1.set_ylabel('p',fontsize = 18)
ax1.tick_params(axis='x', labelsize=15)
ax1.tick_params(axis='y', labelsize=15)
ax1.set_xlabel(r'x/L',fontsize = 18)
ax1.legend(bbox_to_anchor =(0.75, 1.15), ncol = 2)

# Calculating the exact local Mach number
ma2              = vel2x[ypoints//2,:]/np.sqrt(sphr*1*T2[ypoints//2,:])     # Calculating the Mach number using the Ma = velocity / sqrt(sphr * temperature)
max_value = np.amax(ma2)
print("Maximum value using numpy.amax():", max_value)
velexact         = np.array(velexact)
temperatureexact = np.array(temperatureexact)
exactmach        = velexact/np.sqrt(sphr*1*temperatureexact)

# Plotting the local Mach number in a seperate figure
fig = plt.figure("Ma")
# fig, ([ax5]) = plt.subplots(1,1)
plt.plot(x,ma2,label="2D CLBM",color='k',linestyle='-')
plt.plot(x,exactmach,color='red',linestyle='--',label="Exact")
plt.grid(which='major', color='#DDDDDD', linewidth=0.5)
plt.grid(which='minor', color='#EEEEEE', linewidth=0.25)
plt.minorticks_on()
# plt.grid()
plt.ylabel('Ma',fontsize = 18)
plt.xlabel(r'x/L',fontsize = 18)
plt.tick_params(axis='x', labelsize=15)
plt.tick_params(axis='y', labelsize=15)
plt.legend(bbox_to_anchor =(0.75, 1.15), ncol = 2)
# plt.legend()
# Plotting the newton raphson in a seperate figure
fig = plt.figure("nr")
# fig, ([ax5]) = plt.subplots(1,1)
plt.scatter(x,nr[ypoints//2,:],marker='x')
plt.grid(which='major', color='#DDDDDD', linewidth=0.5)
plt.grid(which='minor', color='#EEEEEE', linewidth=0.25)
plt.minorticks_on()
# plt.grid()
plt.ylabel('Newton-Raphson iterations',fontsize = 18)
plt.xlabel(r'x/L',fontsize = 18)
plt.tick_params(axis='x', labelsize=15)
plt.tick_params(axis='y', labelsize=15)
# plt.legend(bbox_to_anchor =(0.75, 1.15), ncol = 2)
# plt.legend()
plt.show()
