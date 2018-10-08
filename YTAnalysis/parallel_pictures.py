# parallel_pictures.py
# James Widdicombe
# Last Updated 17/08/2018
# Plotting script for GRChombo HDF5 files

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

# Timings
start_time = time.time()

# Enable Parallelism
yt.enable_parallelism()

# mkdir the plot directories
if yt.is_root():
	if not os.path.exists('rho'):
		os.mkdir('rho')
	#if not os.path.exists('K'):
		#os.mkdir('K')
	#if not os.path.exists('Pi'):
		#os.mkdir('Pi')

# Plot parameters
line = 2.5 # Line width
alp = 0.6 # alpha

# CHANGE ME
# Loading dataset (Load from one folder up)
data_location = '../*.3d.hdf5' # Data file location

#Loading dataset
ts = yt.load(data_location)

# Plot Start Time
plot_start_time = time.time()
program_runtime_data = []
program_frames=[]
counter=0

for i in ts.piter():

    # Plot rho
    slc = yt.SlicePlot(i,'z','rho',width = (90.0, 'cm'))
    #slc.annotate_grids()
    slc.set_log('rho',False)
    slc.set_buff_size(1024)
    #slc.set_zlim('rho', 0.0007, 0.0000)
    slc.set_cmap(field="rho", cmap='dusk')
    slc.set_xlabel(r'x $\left[\frac{1}{m}\right]$')
    slc.set_ylabel(r'y $\left[\frac{1}{m}\right]$')
    slc.set_colorbar_label("rho", r'$\rho \ \left[\frac{M_{pl}^2}{m}\right]$')
    slc.annotate_text((0.13, 0.92),( 'time = '+str(float(i.current_time))+' 1/m' ) , coord_system='figure',text_args={'color':'white'})
    slc.set_window_size(10)
    slc.save('rho/')

    # Plot Pi
    #slc = yt.SlicePlot(i,'z','Pi',width = (90.0, 'cm'))
    #slc.annotate_grids()
    #slc.set_log('Pi',False)
    #slc.set_buff_size(1024)
    #slc.set_zlim('Pi', 0.0007, 0.0000)
    #slc.set_cmap(field="Pi", cmap='dusk')
    #slc.set_xlabel(r'x $\left[\frac{1}{m}\right]$')
    #slc.set_ylabel(r'y $\left[\frac{1}{m}\right]$')
    #slc.set_colorbar_label("Pi",r'$\Pi$')
    #slc.annotate_text((0.13, 0.92),( 'time = '+str(float(i.current_time))+' 1/m' ) , coord_system='figure',text_args={'color':'white'})
    #slc.set_window_size(10)
    #slc.save('Pi/')

    # Plot K
    #slc = yt.SlicePlot(i,'z','K',width = (90.0, 'cm'))
    #slc.annotate_grids()
    #slc.set_log('K',False)
    #slc.set_buff_size(1024)
    #slc.set_zlim('K', 0.0007, 0.0000)
    #slc.set_cmap(field="K", cmap='dusk')
    #slc.set_xlabel(r'x $\left[\frac{1}{m}\right]$')
    #slc.set_ylabel(r'y $\left[\frac{1}{m}\right]$')
    #slc.set_colorbar_label("K",r'$K$')
    #slc.annotate_text((0.13, 0.92),( 'time = '+str(float(i.current_time))+' 1/m' ) , coord_system='figure',text_args={'color':'white'})
    #slc.set_window_size(10)
    #slc.save('K/')
