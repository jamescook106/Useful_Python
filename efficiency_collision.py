# efficiency_collision.py
# James Widdicombe
# Last updated 02/08/2018
# A script to calculate the efficiency of star
# collision, M_f/M_i

# Load the modules
import os
import time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams

# Matplotlib Settings
rcParams.update({'figure.autolayout': True})
rcParams['axes.formatter.limits'] = [-5,5]

# Input Directory for data files
input_directory = 'YTAnalysis'

# Output Directory
output_directory = 'plots_efficiency'

# Box Size
box_size=512

# Create Output Directory
if not os.path.exists(output_directory):
    os.mkdir(output_directory)

# Data
time_data = np.genfromtxt(input_directory + '/' + 'time_data.out')
total_mass_data = np.genfromtxt(input_directory + '/' + 'Mass_Total.out')
star_mass_data = np.genfromtxt(input_directory + '/' + 'star_mass_ellipse.out')

# System Mass Plot
plt.figure(1)
plt.plot(time_data,total_mass_data, label='Total Mass')
plt.plot(time_data,star_mass_data, label='Star Mass')
plt.ylabel('$M \\, [M_{pl}^2/m]$')
plt.xlabel('Time $[1/m]$')
plt.legend()
plt.savefig(output_directory+'/Mass_Total.png')
plt.close()

# Gather Mf
if time_data[0]==0:
    initial_mass=total_mass_data[0]
else:
    print 'Warning: Time 0 in data file is not 0'
    initial_mass=total_mass_data[0]

# Collision useful information
collision_time = 0.
collision_index = 0
collision_recorder = 0.
post_collision_star_mass=[]
collided_time = []
collided_mass = []

# Detect first non-zero central star mass
for i in range(0,len(time_data)):
    if star_mass_data[i]!=0:
        collided_time.append(time_data[i])
        collided_mass.append(star_mass_data[i])
        if collision_recorder==0:
            collision_time=time_data[i]
            collision_index=i
            collision_recorder=1
        if time_data[i]-collision_time>100 and time_data[i]-collision_time<box_size:
            post_collision_star_mass.append(star_mass_data[i])

# Calculate Efficiency
star_mass_mean = np.mean(post_collision_star_mass)
star_mass_std = np.std(post_collision_star_mass)
efficiency = star_mass_mean/initial_mass
efficiency_error = star_mass_std/initial_mass

# Print useful information
print collision_index
print collision_time
print star_mass_mean
print star_mass_std
print efficiency
print efficiency_error

# Plotting tools
constant_mass = []
constant_mass_plus = []
constant_mass_minus = []

for i in range(collision_index,len(time_data)):
    constant_mass.append(star_mass_mean)
    constant_mass_plus.append(star_mass_mean+star_mass_std)
    constant_mass_minus.append(star_mass_mean-star_mass_std)

# Sanity Check on Mass
plt.figure(2)
plt.plot(collided_time,collided_mass)
plt.plot(collided_time,constant_mass,color='orange')
plt.plot(collided_time,constant_mass_plus,linestyle='--',alpha=0.5,color='orange')
plt.plot(collided_time,constant_mass_minus,linestyle='--',alpha=0.5,color='orange')
plt.ylabel('$M \\, [M_{pl}^2/m]$')
plt.xlabel('Time $[1/m]$')
plt.legend()
plt.savefig(output_directory+'/Star_Mass.png')
plt.close()
