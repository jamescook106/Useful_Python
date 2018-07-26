# master.py
# James Widdicombe
# Last Updated 25/07/2018
# A master analysis script that include L2M, L2H and mass measurments
# for GRChombo HDF5 files

# Load the modules
import yt
import os
import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import fsolve
from yt import derived_field
import time
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# Enable Parallelism
yt.enable_parallelism()

# Plot parameters
line = 2.5 # Line width
alp = 0.6 # alpha

# CHANGE ME
# Loading dataset (Load from one folder up)
ds = yt.load('../../outMatterSF_*.3d.hdf5')

# Define H2
def _H2(field, data):
        return data["Ham"]*data["Ham"]

# Define M2
def _M2(field, data):
        return data["Mom1"]*data["Mom1"] + data["Mom2"]*data["Mom2"] + data["Mom3"]*data["Mom3"]

# Define Modified Rho
def _RhoJam(field, data):
        return data["rho"]/(data["chi"] * np.sqrt(data["chi"]))

# mkdir the plot directories
if not os.path.exists('rho_size'):
    os.mkdir('rho_size')

# Arrays for output data
time_data = []
L2H_data = []
L2M_data = []
RhoJam_Average_data = []
Mass_Total_data=[]
star_size_initial_x_data=[]
SizeyData = []
max_rho_pos_x_data = []
max_rho_pos_y_data = []
star_mass_center_x=[]
star_mass_center_y=[]
star_mass_center_average=[]

# Other factors
total_box_size = float(ds[0].domain_right_edge[0])
center = total_box_size/2.
Quality = 'cubic'

# Plot Start Time
plot_start_time = time.time()
program_runtime_data = []
program_frames=[]
counter=0

for i in ds:

    #Timings
    i_start = time.time()
    counter+=1
    program_frames.append(counter)

    # Current Simulation Time
    time_data.append(i.current_time)

    # Add the M2 and L2 Fields
    i.add_field("H2",_H2, units = "")
    i.add_field("M2",_M2, units = "")

    # Add the RhoJam Field to the data
    i.add_field("RhoJam",_RhoJam, units = "")

    # Exclude the edges of the simulation
    #ad = i.r[8:504,8:504,8:504]

    # All Data
    ad = i.all_data()

    # L2H
    meanH2 = ad.mean("H2", weight="cell_volume")
    L2H = np.sqrt(meanH2)
    L2H_data.append(L2H)

    # L2M
    meanM2 = ad.mean("M2", weight="cell_volume")
    L2M = np.sqrt(meanM2)
    L2M_data.append(L2M)

    # Average RhoJam
    Rho_Jam_Average = ad.mean("RhoJam", weight="cell_volume")
    RhoJam_Average_data.append(Rho_Jam_Average)

    # Total RhoJam
    Mass_Total = Rho_Jam_Average*(total_box_size**3)
    Mass_Total_data.append(Mass_Total)

    program_runtime_data.append(time.time()-i_start)
    np.savetxt('program_runtime_data.out',program_runtime_data)
    np.savetxt('L2M.out',L2M_data)
    np.savetxt('L2H.out',L2H_data)
    np.savetxt('RhoJam_Average.out',RhoJam_Average_data)
    np.savetxt('Mass_Total.out',Mass_Total_data)

    # Look at initial Stars
    # Line out from 0 to center (x), at the center of (y and z)
    c = i.ray((0,center,center),(center,center,center))

    # Name x coordinate and rho from lineout
    x = c["x"]
    rho = c["rho"]

    # Find maximum rho value within lineout
    rho_max = float(rho.max())

    # Find the coordinate of where the maximum value is
    x_max = np.where(rho == rho_max)[0][0]
    x_max_val = x[x_max]
    max_rho_pos_x_data.append(x_max_val)

    # Create a function that goes to negative when we are below 95% of rho max
    rhofun = interp1d(x, rho-0.05*rho_max, kind=Quality, fill_value="extrapolate")

    # Solve the function for where it goes to zero using a best guess past and
    # before the the max value of x
    x_1 = fsolve(rhofun,float(x[x_max])-1)[0]
    x_2 = fsolve(rhofun,float(x[x_max])+1)[0]

    # Size of the 95% rho in center
    size_x = (x_2-x_1)/2.
    star_size_initial_x_data.append(size_x)

    np.savetxt('star_size_inimax_rho_postial_x.out',star_size_initial_x_data)

    cy = i.ray((x_max_val,0,center),(x_max_val,total_box_size,center))

    y = cy["y"]
    rhoy = cy["rho"]

    rhoy_max = float(rhoy.max())

    rhoyfun = interp1d(y, rhoy-0.05*rhoy_max, kind=Quality)
    y_1 = fsolve(rhoyfun,center-1)[0]
    y_2 = fsolve(rhoyfun,center+1)[0]

    size_y = (y_2-y_1)/2.
    SizeyData.append(size_y)

    plt.figure(figsize=(20,20))
    ax1 = plt.subplot(211)
    plt.scatter(x,rho, alpha = alp,label = "t ="+str(i.current_time))
    plt.xlim([x_1-10,x_2+10])
    #plt.ylim([-0.0005,0.008])
    plt.legend(loc = "upper right")
    plt.xlabel(r'$t~[1/m]$')
    plt.ylabel(r'$\rho~[M_{pl}^2 m^2]$')

    ax2 = plt.subplot(212)
    plt.plot(time_data,star_size_initial_x_data, linewidth = line,color = "red", alpha = alp,label = "size = 5% maxRho x dir ")
    plt.plot(time_data,SizeyData, linewidth = line,linestyle = "--",color = "blue", alpha = alp,label = "size = 5% maxRho y dir ")
#    plt.xlim([0,400])
    plt.legend(loc = "lower right")
    plt.xlabel(r'$t~[1/m]$')
    plt.ylabel(r'$x~[1/m]$')
    plt.savefig(("rho_size/rho%05d.png" % counter),bbox_inches = 'tight')
    plt.close()

    if x_max_val>center-5 and x_max_val<center+5:
        if size_x<1.0:
            print size_x
            star_mass_center_x.append(0.0)
            star_mass_center_y.append(0.0)
            star_mass_center_average.append(0.0)
        elif size_y<1.0:
            print size_y
            star_mass_center_x.append(0.0)
            star_mass_center_y.append(0.0)
            star_mass_center_average.append(0.0)
        else:
            sp = i.sphere([center, center, center], size_x)
            RhoJamTotalx = float(sp.mean("RhoJam", weight="cell_volume"))*(4./3.)*np.pi*np.power((size_x),3)
            star_mass_center_x.append(RhoJamTotalx)
            sp2 = i.sphere([center, center, center], size_y)
            RhoJamTotaly = float(sp2.mean("RhoJam", weight="cell_volume"))*(4./3.)*np.pi*np.power((size_y),3)
            star_mass_center_y.append(RhoJamTotaly)
            star_mass_center_average.append(np.mean([star_mass_center_y,star_mass_center_x]))
    else:
        star_mass_center_x.append(0.0)
        star_mass_center_y.append(0.0)
        star_mass_center_average.append(0.0)

# Plots
# Update rcParams so that we get all axis labelling
rcParams.update({'figure.autolayout': True})
rcParams['axes.formatter.limits'] = [-5,5]

# Timing Plot
plt.figure(1)
plt.title('Program Run Time')
plt.plot(program_frames,program_runtime_data)
plt.xlabel('Data Files')
plt.ylabel('Time $[s]$')
plt.savefig('pictures_program_time.png')
plt.close()

# L2H
plt.figure(2)
plt.plot(time_data,L2H_data)
plt.title('$\\mathcal{H}$ vs time')
plt.ylabel('$\\mathcal{H}$')
plt.xlabel('Time $[1/m]$')
plt.savefig('H.png')
plt.close()

# L2M
plt.figure(3)
plt.plot(time_data,L2M_data)
plt.title('$\\mathcal{M}$ vs time')
plt.ylabel('$\\mathcal{M}$')
plt.xlabel('Time $[1/m]$')
plt.savefig('M.png')
plt.close()

# RhoJam_Average
plt.figure(4)
plt.plot(time_data,RhoJam_Average_data)
plt.title('$\\tilde{\\rho}$ vs time')
plt.ylabel('$\\tilde{\\rho} \\, [M_{pl}^2m^2]$')
plt.xlabel('Time $[1/m]$')
plt.savefig('RhoJam_Average.png')
plt.close()

# Total Mass
plt.figure(5)
plt.plot(time_data,Mass_Total_data)
plt.plot(time_data,star_mass_center_y)
plt.plot(time_data,star_mass_center_x)
plt.plot(time_data,star_mass_center_average)
plt.title('$M$ vs time')
plt.ylabel('$M \\, [M_{pl}^2/m]$')
plt.xlabel('Time $[1/m]$')
plt.savefig('Mass_Total.png')
plt.close()

# Max rho pos
plt.figure(6)
plt.plot(time_data,max_rho_pos_x_data)
plt.title('$x,y$ vs time')
plt.xlabel('Time $[1/m]$')
plt.savefig('max_rho_pos.png')
plt.close()
