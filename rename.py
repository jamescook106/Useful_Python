# rename.py
# James Widdicombe
# A script that copies a series on png and reassigns them into a standard order
# that can be used with ffmpeg. Then creates a movie from those files
# Last Updated June 2018

# Import Modules
import os, sys, shutil
import glob

# Function to count png files
def PngCounter(str_input_directory):
    return len(glob.glob1(str_input_directory,"*.png"))

# Function with custom file output
def OutputName(int_file_number,str_output_path):
    if int_file_number<10:
        return str_output_path+'/00000'+str(int_file_number)+'.png'
    elif int_file_number<100 and int_file_number>9:
        return str_output_path+'/0000'+str(int_file_number)+'.png'
    elif int_file_number<1000 and int_file_number>99:
        return str_output_path+'/000'+str(int_file_number)+'.png'
    elif int_file_number<10000 and int_file_number>999:
        return str_output_path+'/00'+str(int_file_number)+'.png'
    elif int_file_number<100000 and int_file_number>9999:
        return str_output_path+'/0'+str(int_file_number)+'.png'

# Input and Output Directories
Variable = 'rho'
Input_Directory = "YTAnalysis/"+Variable
Output_Path = 'Output'

# Make an Output_Path
if not os.path.exists(Output_Path):
    os.mkdir(Output_Path)

Total_Files = PngCounter(Input_Directory)

while_counter = 0
name_counter = 0

# Rename the png
while while_counter<Total_Files:
    if name_counter<10:
        j = Input_Directory + '/outMatterSF_00000' + str(name_counter) + '.3d.hdf5_Slice_z_'+Variable+'.png'
        if os.path.isfile(j):
            print j
            shutil.copyfile(j,OutputName(while_counter,Output_Path))
            while_counter+=1
    elif name_counter>9 and name_counter<100:
        j = Input_Directory + '/outMatterSF_0000' + str(name_counter) + '.3d.hdf5_Slice_z_'+Variable+'.png'
        if os.path.isfile(j):
            print j
            shutil.copyfile(j,OutputName(while_counter,Output_Path))
            while_counter+=1
    elif name_counter>99 and name_counter<1000:
        j = Input_Directory + '/outMatterSF_000' + str(name_counter) + '.3d.hdf5_Slice_z_'+Variable+'.png'
        if os.path.isfile(j):
            print j
            shutil.copyfile(j,OutputName(while_counter,Output_Path))
            while_counter+=1
    elif name_counter>999 and name_counter<10000:
        j = Input_Directory + '/outMatterSF_00' + str(name_counter) + '.3d.hdf5_Slice_z_'+Variable+'.png'
        if os.path.isfile(j):
            print j
            shutil.copyfile(j,OutputName(while_counter,Output_Path))
            while_counter+=1
    elif name_counter>9999 and name_counter<100000:
        j = Input_Directory + '/outMatterSF_0' + str(name_counter) + '.3d.hdf5_Slice_z_'+Variable+'.png'
        if os.path.isfile(j):
            print j
            shutil.copyfile(j,OutputName(while_counter,Output_Path))
            while_counter+=1
    name_counter+=1

# Create the movie
os.chdir(Output_Path)
os.system('ffmpeg -framerate 12 -i %06d.png -c:v libx264 -r 30 -pix_fmt yuv420p movie.mp4')
os.chdir('..')
shutil.copyfile(Output_Path+'/movie.mp4',os.getcwd()+'/'+Variable+'.mp4')
shutil.rmtree(Output_Path)
