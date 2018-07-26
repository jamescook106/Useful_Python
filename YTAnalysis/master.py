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

# Arrays for output data
time_data = []
L2H_data = []
L2M_data = []

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

    # Add the M2 and L2 Fields
    i.add_field("H2",_H2, units = "")
    i.add_field("M2",_M2, units = "")

    # Exclude the edges of the simulation
    ad = i.r[8:504,8:504,8:504]

    # L2H
    meanH2 = ad.mean("H2", weight="cell_volume")
    L2H = np.sqrt(meanH2)
    L2H_data.append(L2H)

    # L2M
    meanM2 = ad.mean("M2", weight="cell_volume")
    L2M = np.sqrt(meanM2)
    L2M_data.append(L2M)

    program_runtime_data.append(time.time()-i_start)
    np.savetxt('program_runtime_data.out',program_runtime_data)
    np.savetxt('L2M.out',L2M_data)
    np.savetxt('L2H.out',L2H_data)

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
plt.close

# L2H
plt.figure(2)
plt.plot(time_data,L2H_data)
plt.title('$\\mathcal{H}$ vs time')
plt.ylabel('$\\mathcal{H}$')
plt.xlabel('Time $[s]$')
plt.savefig('H.png')

# L2M
plt.figure(3)
plt.plot(time_data,L2M_data)
plt.title('$\\mathcal{M}$ vs time')
plt.ylabel('$\\mathcal{M}$')
plt.xlabel('Time $[s]$')
plt.savefig('M.png')
