# data_compare.py
# Last updated 15/08/2018
# Last Formatted Dec 2019
# James Widdicombe
# Simple script to compare data from a serial and parallel
# run of the same program

# Load the modules
import numpy as np

# Input Directory for data files
serial_directory = "serial"
parallel_directory = "parallel"

# Time Data
# Serial
serial_time_data = np.genfromtxt(serial_directory + "/" + "time.out")

# Parallel
parallel_time_data = np.genfromtxt(parallel_directory + "/" + "time.out")

# Input data have different file lengths
time_difference_index = min(len(serial_time_data), len(parallel_time_data))

test_1 = 0
for i in range(time_difference_index):
    if parallel_time_data[i] - serial_time_data[i] != 0:
        print(i)
        test_1 = 1
if test_1 == 0:
    print("Test 1 Passed - Time files align")
else:
    print("Test 1 FAILED - Time files do not align")
# All input files
input_name_array = [
    "Weyl4_l2_m0_data.out",
    "Weyl4_l2_m1_data.out",
    "Weyl4_l2_m2_data.out",
    "Weyl4_l2_m1n_data.out",
    "Weyl4_l2_m2n_data.out",
    "Weyl4_l3_m0_data.out",
    "Weyl4_l3_m1_data.out",
    "Weyl4_l3_m2_data.out",
    "Weyl4_l3_m3_data.out",
    "Weyl4_l3_m1n_data.out",
    "Weyl4_l3_m2n_data.out",
    "Weyl4_l3_m3n_data.out",
    "Weyl4_l4_m0_data.out",
    "Weyl4_l4_m1_data.out",
    "Weyl4_l4_m2_data.out",
    "Weyl4_l4_m3_data.out",
    "Weyl4_l4_m4_data.out",
    "Weyl4_l4_m1n_data.out",
    "Weyl4_l4_m2n_data.out",
    "Weyl4_l4_m3n_data.out",
    "Weyl4_l4_m4n_data.out",
]

serial_master = []
parallel_master = []

# Load all spherical harmonics
for i in input_name_array:
    serial_master.append(
        np.genfromtxt(
            serial_directory + "/" + i,
            dtype=np.complex128,
            converters={0: lambda s: complex(s.decode().replace("+-", "-"))},
        )
    )
    parallel_master.append(
        np.genfromtxt(
            parallel_directory + "/" + i,
            dtype=np.complex128,
            converters={0: lambda s: complex(s.decode().replace("+-", "-"))},
        )
    )
# Check all output numbers are the same
test_2 = 0
for i in range(len(serial_master)):
    for j in range(min(len(serial_master[i]), len(parallel_master[i]))):
        if serial_master[i][j] - parallel_master[i][j] != 0:
            print(serial_master[i][j] - parallel_master[i][j])
            test_2 = 1
if test_2 == 0:
    print("Test 2 Passed - All output files are the same")
else:
    print("Test 2 FAILED - All output files are NOT the same")
