# Example file for converting outputs of numerical simulation for Gyoto
#
# Copyright 2026 Nicolas Aimar
#
# This file is part of Gyoto.
#
# Gyoto is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Gyoto is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Gyoto.  If not, see <http://www.gnu.org/licenses/>.

### Explanation:
#
# This script converts the outputs of a numerical simulation (GRMHD or PIC)
# into FITS files readable by Gyoto: one FITS file per simulation snapshot.
#
# HOW IT WORKS
#   1. SETUP: set the input/output paths and the physical options
#      (fluid or kinetic, velocity, magnetic field).
#   2. READING: for each snapshot, read the grid, the time and the physical
#      quantities from your simulation file into numpy arrays.
#   3. WRITING: the arrays are automatically written to a FITS file.
#
# WHAT YOU MUST DO
#   - Edit the SETUP section.
#   - Fill the parts marked ## TO BE ADAPTED ##, ## TO BE READ FROM FILE ##
#     and ## TO BE ADDED ACCORDING TO YOUR FILES ##.
#   - Provide at least:
#       fluid   : density and temperature, shape (nx1, nx2, nx3)
#       kinetic : j_I, shape (nx1, nx2, nx3, nfreq, npitch)
#   - Optionally provide: vel (3, nx1, nx2, nx3), Bfield (4, nx1, nx2, nx3)
#     and, for kinetic runs, the other radiative coefficients (same shape as j_I).
#     Set the matching flags (velocity, magnetic_field) accordingly.
#   - For a 2D simulation, set the unused dimension length to 1.
#
# WHAT YOU MUST NOT CHANGE
#   - Variable names (nx1, nx2, nx3, nfreq, npitch, ntot, time_array,
#     x1_array, ..., density, temperature, j_I, vel, Bfield).
#   - The shape and axis order of the arrays.
#   - Everything below the line "DO NOT CHANGE ANYTHING BELOW"
#     (conversion and FITS writing).


# ========================================================================
# SETUP
# 
import numpy as np
import gyoto
import math
import os

# 1- Directory and files managment

inputdir        = "examples_Sim_Files/blk/"  # directory of input file(s) from the simulation
fname_prefix    = "GRMHD-example_"           # prefix of the name of file(s) (without the number and extension)
fname_suffix    = ".blk"                     # extension of the input file(s)
nfile           = 11                         # number of file(s) to read and convert

outputdir       = "fits/"                    # output directoy where the fits file(s) will be generated
fitsname_prefix = "data"

if not os.path.exists(outputdir):
    os.makedirs(outputdir)

# 2- Physical parameters
fluid           = True                       # boolean that specifiy if the simulations are fluid (GRMHD; True) or kinetic (PIC; False) 

magnetic_field  = 0                          # int that specify if the magnetic field is in the simulation outputs and given to Gyoto (0 = no; 1 = yes)
velocity        = 1                          # boolean that specify if the velocity is in the simulation outputs and given to Gyoto (0 = no; 1 = yes)


# ========================================================================
# READING and CONVERTING
#
# 1- The first step is to fill the time array of the simulation output files.
# This is because in each FITS file created we must include the time of this file
# as well as the overall time array to allow GYOTO to read the correct files.
#
# 2- Read the content of each file
#
# 3- Create the associated FITS file readable by GYOTO
#
# NOTE: This section MUST be filled according to your simulation outputs.
# The array that must be filled are : density and temperature for fluid or at least j_I for kinetics
# To which you can add vel and Bfield (3D-velocity and 4-magnetic field respectively) if provided (flags are True)

# 1- Fill time_array
print("Filling time_array...")
time_array = np.zeros(nfile)
for ii in range(nfile):
    curnum=str(ii)
    with open(inputdir+fname_prefix+curnum.zfill(4)+fname_suffix,"r") as datafile:

        time_array[ii] = float(ii) # example for the script to not crash. ## TO BE ADAPTED ##
print("Done.")


# Looping on all the simulation output files
for ii in range(nfile):
    print("File index ", ii, " over ", nfile-1)
    
    # 2- Reading simulation outputs files
    curnum=str(ii)
    datafile=open(inputdir+fname_prefix+curnum.zfill(4)+fname_suffix,"r") 
    
    
    # Reading grid size from file. DO NOT CHANGE NAMES
    # If simulation is 2D put the correct dimension length to 1, it will be automatically treated.
    nx1 = 1 # example for the script to not crash. ## TO BE READ FROM FILE ##
    nx2 = 1 # example for the script to not crash. ## TO BE READ FROM FILE ##
    nx3 = 1 # example for the script to not crash. ## TO BE READ FROM FILE ##
    if not fluid:
        nfreq = 1 # example for the script to not crash. ## TO BE READ FROM FILE ##
        npitch = 1 # example for the script to not crash. ## TO BE READ FROM FILE ##
    
    # defining the grid arrays
    x1_array = np.ndarray(shape=(nx1), dtype=float)
    x2_array = np.ndarray(shape=(nx2), dtype=float)
    x3_array = np.ndarray(shape=(nx3), dtype=float)
    if not fluid:
        frequency_array   = np.ndarray(shape=(nfreq), dtype=float)
        pitch_angle_array = np.ndarray(shape=(npitch), dtype=float)
    
    # Fill the grid arrays from file
    ## TO BE ADDED ACCORDING TO YOUR FILES ##
    #
    
    
    
    
    #
    if (fluid):
        ntot = nx1*nx2*nx3 # DO NOT CHANGE THIS
        print("nx1, nx2, nx3, ntot= ", nx1, nx2, nx3, ntot)
    else:
        ntot = nx1*nx2*nx3*nfreq*npitch # DO NOT CHANGE THIS
        print("nx1, nx2, nx3, nfreq, npitch, ntot= ", nx1, nx2, nx3, nfreq, npitch, ntot)
    
    
    # Defining the arrays that will be filled with simulation data
    # DO NOT CHANGE THE STRUCTURE OF THE ARRAYS
    if (fluid):
        density     = np.ndarray(shape=(nx1, nx2, nx3), dtype=float)
        temperature = np.ndarray(shape=(nx1, nx2, nx3), dtype=float)
    else:
        j_I         = np.ndarray(shape=(nx1, nx2, nx3, nfreq, npitch), dtype=float)
        
        # Uncomment the following array if provided in simulation output files
        j_Q         = np.ndarray(shape=(nx1, nx2, nx3, nfreq, npitch), dtype=float)
        #j_U         = np.ndarray(shape=(nx1, nx2, nx3, nfreq, npitch), dtype=float)
        #j_V         = np.ndarray(shape=(nx1, nx2, nx3, nfreq, npitch), dtype=float)
        
        #a_I         = np.ndarray(shape=(nx1, nx2, nx3, nfreq, npitch), dtype=float)
        #a_Q         = np.ndarray(shape=(nx1, nx2, nx3, nfreq, npitch), dtype=float)
        #a_U         = np.ndarray(shape=(nx1, nx2, nx3, nfreq, npitch), dtype=float)
        #a_V         = np.ndarray(shape=(nx1, nx2, nx3, nfreq, npitch), dtype=float)
        
        #r_Q         = np.ndarray(shape=(nx1, nx2, nx3, nfreq, npitch), dtype=float)
        #r_U         = np.ndarray(shape=(nx1, nx2, nx3, nfreq, npitch), dtype=float)
        #r_V         = np.ndarray(shape=(nx1, nx2, nx3, nfreq, npitch), dtype=float)
        
    if (velocity):
        vel=np.ndarray(shape=(3, nx1, nx2, nx3), dtype=float)
    
    if (magnetic_field):
        Bfield=np.ndarray(shape=(4, nx1, nx2, nx3), dtype=float)
    
    
    # Filling the arrays from simulation data file
    print("Reading data...") 
    ## TO BE ADDED ACCORDING TO YOUR FILES ##
    #
    
    
    
    
    #
            

    datafile.close()
    print("Done.")

    # At this point you should have in any case the time and the spatial 1D-arrays
    # For fluid simulations, you must have at least the DENSITY and TEMPERATURE 2D-arrays
    # For kinetic simulations, you must have at least the TOTAL_EMISSION (j_I) 4D-array (2 spatial dimensions + frequency + pitch angle)
    # The other radiative transfer coefficients are optional
    # Optional arrays are : VELOCITY and MAGNETIC_FIELD

    # 3- Create the associated FITS file readable by GYOTO
    #
    # DO NOT CHANGE ANYTHING BELOW
    #

    # Rearranging stuff for gyoto
    timegyoto = gyoto.core.array_double_fromnumpy1(time_array)
    x1gyoto = gyoto.core.array_double_fromnumpy1(x1_array)
    x2gyoto = gyoto.core.array_double_fromnumpy1(x2_array)
    x3gyoto = gyoto.core.array_double_fromnumpy1(x3_array)
    if not fluid:
        freqgyoto = gyoto.core.array_double_fromnumpy1(frequency_array)
        pitchAnglegyoto = gyoto.core.array_double_fromnumpy1(pitch_angle_array)
    
    
    if (velocity):
        v11D = gyoto.core.array_double_fromnumpy3(vel[0])
        v21D = gyoto.core.array_double_fromnumpy3(vel[1])
        v31D = gyoto.core.array_double_fromnumpy3(vel[2])
    
    if (magnetic_field):
        B01D = gyoto.core.array_double_fromnumpy3(Bfield[0])
        B11D = gyoto.core.array_double_fromnumpy3(Bfield[1])
        B21D = gyoto.core.array_double_fromnumpy3(Bfield[2])
        B31D = gyoto.core.array_double_fromnumpy3(Bfield[3])
    
    fits = gyoto.FitsRW()
    filename = outputdir+fitsname_prefix+str(ii).zfill(4)+'.fits'
    
    fptr = fits.fitsCreate(filename)
    
    fits.fitsWriteKey(fptr, "NB_X0", nfile)
    fits.fitsWriteKey(fptr, "NB_X1", nx1)
    fits.fitsWriteKey(fptr, "NB_X2", nx2)
    fits.fitsWriteKey(fptr, "NB_X3", nx3)
    fits.fitsWriteKey(fptr, "TIME", time_array[ii])
    fits.fitsWriteKey(fptr, "BINFILE", magnetic_field)
    if not fluid:
        fits.fitsWriteKey(fptr, "NB_FREQ", nfreq)
        fits.fitsWriteKey(fptr, "NB_PITCH", npitch)        
    
    fits.fitsWriteHDUData(fptr, "X0", timegyoto, nfile)
    fits.fitsWriteHDUData(fptr, "X1", x1gyoto, nx1)
    fits.fitsWriteHDUData(fptr, "X2", x2gyoto, nx2)
    fits.fitsWriteHDUData(fptr, "X3", x3gyoto, nx3)
    if not fluid:
        fits.fitsWriteHDUData(fptr, "FREQ", freqgyoto, nfreq)
        fits.fitsWriteHDUData(fptr, "PITCH", pitchAnglegyoto, npitch)
    
    if (velocity):
        fits.fitsWriteHDUData(fptr, "VELOCITY1", v11D, ntot)
        fits.fitsWriteHDUData(fptr, "VELOCITY2", v21D, ntot)
        fits.fitsWriteHDUData(fptr, "VELOCITY3", v31D, ntot)
    if (magnetic_field):
        fits.fitsWriteHDUData(fptr, "B0", B01D, ntot)
        fits.fitsWriteHDUData(fptr, "B1", B11D, ntot)
        fits.fitsWriteHDUData(fptr, "B2", B21D, ntot)
        fits.fitsWriteHDUData(fptr, "B3", B31D, ntot)
    
    if (fluid):
        density1D = gyoto.core.array_double_fromnumpy3(density)
        temperature1D = gyoto.core.array_double_fromnumpy3(temperature)
        fits.fitsWriteHDUData(fptr, "NUMBERDENSITY", density1D, ntot)
        fits.fitsWriteHDUData(fptr, "TEMPERATURE", temperature1D, ntot)
    else:
        J_I1D = gyoto.core.array_double_fromnumpy1(np.ascontiguousarray(j_I, dtype=np.float64).ravel(order='C'))
        fits.fitsWriteHDUData(fptr, "J_I", J_I1D, ntot)
        
        try:
            J_Q1D = gyoto.core.array_double_fromnumpy1(np.ascontiguousarray(j_Q, dtype=np.float64).ravel(order='C'))
            fits.fitsWriteHDUData(fptr, "J_Q", J_Q1D, ntot)
        except NameError:
            print("Coefficient j_Q not defined, skipping...")
        try:
            J_U1D = gyoto.core.array_double_fromnumpy1(np.ascontiguousarray(j_U, dtype=np.float64).ravel(order='C'))
            fits.fitsWriteHDUData(fptr, "J_U", J_U1D, ntot)
        except NameError:
            print("Coefficient j_U not defined, skipping...")
        try:
            J_V1D = gyoto.core.array_double_fromnumpy1(np.ascontiguousarray(j_V, dtype=np.float64).ravel(order='C'))
            fits.fitsWriteHDUData(fptr, "J_V", J_V1D, ntot)
        except NameError:
            print("Coefficient j_V not defined, skipping...")
        try:
            A_I1D = gyoto.core.array_double_fromnumpy1(np.ascontiguousarray(a_I, dtype=np.float64).ravel(order='C'))
            fits.fitsWriteHDUData(fptr, "ALPHA_I", A_I1D, ntot)
        except NameError:
            print("Coefficient a_I not defined, skipping...")
        try:
            A_Q1D = gyoto.core.array_double_fromnumpy1(np.ascontiguousarray(a_Q, dtype=np.float64).ravel(order='C'))
            fits.fitsWriteHDUData(fptr, "ALPHA_Q", A_Q1D, ntot)
        except NameError:
            print("Coefficient a_Q not defined, skipping...")
        try:
            A_U1D = gyoto.core.array_double_fromnumpy1(np.ascontiguousarray(a_U, dtype=np.float64).ravel(order='C'))
            fits.fitsWriteHDUData(fptr, "ALPHA_U", A_U1D, ntot)
        except NameError:
            print("Coefficient a_U not defined, skipping...")
        try:
            A_V1D = gyoto.core.array_double_fromnumpy1(np.ascontiguousarray(a_V, dtype=np.float64).ravel(order='C'))
            fits.fitsWriteHDUData(fptr, "ALPHA_V", A_V1D, ntot)
        except NameError:
            print("Coefficient a_V not defined, skipping...")
        try:
            R_Q1D = gyoto.core.array_double_fromnumpy1(np.ascontiguousarray(r_Q, dtype=np.float64).ravel(order='C'))
            fits.fitsWriteHDUData(fptr, "R_Q", R_Q1D, ntot)
        except NameError:
            print("Coefficient r_Q not defined, skipping...")
        try:
            R_U1D = gyoto.core.array_double_fromnumpy1(np.ascontiguousarray(r_U, dtype=np.float64).ravel(order='C'))
            fits.fitsWriteHDUData(fptr, "R_U", R_U1D, ntot)
        except NameError:
            print("Coefficient r_U not defined, skipping...")
        try:
            R_V1D = gyoto.core.array_double_fromnumpy1(np.ascontiguousarray(r_V, dtype=np.float64).ravel(order='C'))
            fits.fitsWriteHDUData(fptr, "R_V", R_V1D, ntot)
        except NameError:
            print("Coefficient r_V not defined, skipping...")
        
    fits.fitsClose(fptr)
    print("")
