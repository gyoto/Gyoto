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
# 
# Fixed grid (adaptive mesh not supported yet)


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


# 2- Physical parameters
fluid           = True                       # boolean that specifiy if the simulations are fluid (GRMHD; True) or kinetic (PIC; False) 

magnetic_field  = 0                          # int that specify if the magnetic field is in the simulation outputs and given to Gyoto (0 = no; 1 = yes)
velocity        = 1                          # boolean that specify if the velocity is in the simulation outputs and given to Gyoto (0 = no; 1 = yes)

#

######

# ========================================================================
# READING
#
# This section MUST be adapted to your simulation outputs

### Parameters for FLUID simulations ###
#
# Here we define the temperature from density.
# This can be changed according to your simulations.

numberDensityfactor = 1e7 # scaling factor for the number density into SI units 
Tinner = 3e5 # Temperature at inner radius
gamma = 2./3.
# T(rho) = Tinner*(rho/rhoinner)^gamma*(r/rinner)
rhoinner = 1.8887412394 # rhoinner is the density at the max, at t=tmin, phi=phimin, r=rmin
rinner   = 9.0
#

### Parameters for KINETIC simulations ###
#
# TO BE ADDED

### Reading grid details
# In this example, the grid is equatorial (no vertical axis)  
datafile=open(inputdir+fname_prefix+'0000'+fname_suffix,"r")
line1=datafile.readline() # first line comments 
line2=datafile.readline() # grid details
columns=line2.split()

ntot=int(np.asarray(columns)[0])
nx1=int(np.asarray(columns)[1])
nx3=int(np.asarray(columns)[2])
nx2=1


line3=datafile.readline() # time info

alllines=datafile.readlines()

## These arrays MUST be defined with these names and filled with the grid
time_array = np.zeros(nfile)
x1_array   = np.zeros(nx1) # In this example x1 = radius
x2_array   = np.array([math.pi/2.]) # x2 = polar angle; only 1 value at equatorial plane because 2D case
x3_array   = np.zeros(nx3) # x3 = azimuthal angle

linesradius = alllines[0:nx1]
for i in range(nx1):
    x1_array[i] = float(linesradius[i].split()[0])

for j in range(nx3):
    x3_array[j] = float(alllines[j*nx1].split()[1])

print("nx1, nx3, ntot= ",nx1,nx3,ntot)
assert(nx1*nx2*nx3==ntot)

print("rhoinner= ",rhoinner)

datafile.seek(0) # back at first line of file

# Fill time_array # IMPORTANT TO BE HERE
for ii in range(nfile):
    curnum=str(ii)
    datafile=open(inputdir+fname_prefix+curnum.zfill(4)+fname_suffix,"r") 

    line1=datafile.readline() # first line comments
    line2=datafile.readline() # grid details
    line3=datafile.readline() # time
    columns=line3.split()
    time_array[ii]=float(np.asarray(columns)[0])

print("")

if (fluid):
    rho=np.ndarray(shape=(nfile, nx1, nx2, nx3), dtype=float) # this structure must not be changed
    temperature=np.ndarray(shape=(nfile, nx1, nx2, nx3), dtype=float) # this structure must not be changed
    vel=np.ndarray(shape=(nfile, 3, nx1, nx2, nx3), dtype=float) # this structure must not be changed
    if (magnetic_field):
        Bfield=np.ndarray(shape=(nt, 4, nx1, nx2, nx3), dtype=float) # this structure must not be changed
    
for ii in range(nfile):

    print("File index ",ii," over ",nfile-1)
    
    if (fluid):
        minrad=mintheta=minphi=minrho=minv1=minv2=minv3=minT=minB0=minB1=minB2=minB3=1e30
        maxrad=maxtheta=maxphi=maxrho=maxv1=maxv2=maxv3=maxT=maxB0=maxB1=maxB2=maxB3=-minrad
    else:
        print("TO BE ADDED")
    
    
    curnum=str(ii)
    datafile=open(inputdir+fname_prefix+curnum.zfill(4)+fname_suffix,"r") 

    line1=datafile.readline() # first line comments
    line2=datafile.readline() # grid details
    columns=line2.split()

    ntotc=int(np.asarray(columns)[0])
    nx1c=int(np.asarray(columns)[1])
    nx3c=int(np.asarray(columns)[2])
    
    assert(ntotc==ntot)
    assert(nx1c==nx1)
    assert(nx3c==nx3)

    line3=datafile.readline() # time
    columns=line3.split()
    
    alllines=datafile.readlines()
    datafile.close()
    for x1 in range(nx1):
        for x2 in range(nx2):
            for x3 in range(nx3):
                linecur = alllines[x2*nx1*nx3+x3*nx1+x1] # because there is no theta variation, it is the one which evolve the slowest, then it is x3 (phi) and finally x1 (radius) which evolves the fastest; TO BE ADAPTED
                columns = linecur.split()
                if (fluid):
                    rho[ii, x1, x2, x3] = rhocur = float(columns[2])*numberDensityfactor
                    if rhocur<minrho:
                        minrho=rhocur
                    if rhocur>maxrho:
                        maxrho=rhocur
                    
                    temperature[ii, x1, x2, x3] = Tcur = Tinner*((rhocur/rhoinner)**gamma)*(x1/rinner)
                    if Tcur<minT:
                        minT=Tcur
                    if Tcur>maxT:
                        maxT=Tcur
                else:
                    print("TO BE ADDED")
                
                
                if (velocity):
                    v1cur = float(columns[4]) # Adapt the column number
                    v2cur = float(columns[5]) # Adapt the column number
                    v3cur = 0. # here we have vr and vphi only
                    vel[ii, 0, x1, x2, x3] = v1cur
                    vel[ii, 1, x1, x2, x3] = v2cur
                    vel[ii, 2, x1, x2, x3] = v3cur
                    if v1cur<minv1:
                        minv1=v1cur
                    if v1cur>maxv1:
                        maxv1=v1cur
                    if v2cur<minv2:
                        minv2=v2cur
                    if v2cur>maxv2:
                        maxv2=v2cur
                    if v3cur<minv3:
                        minv3=v3cur
                    if v3cur>maxv3:
                        maxv3=v3cur
                
                if (magnetic_field):
                    B0cur = float(columns[7]) # Adapt the column number 
                    B1cur = float(columns[8]) # Adapt the column number
                    B2cur = float(columns[9]) # Adapt the column number
                    B3cur = float(columns[10]) # Adapt the column number
                    Bfield[ii, 0, x1, x2, x3] = B0cur
                    Bfield[ii, 1, x1, x2, x3] = B1cur
                    Bfield[ii, 2, x1, x2, x3] = B2cur
                    Bfield[ii, 3, x1, x2, x3] = B3cur
                    if B0cur<minB0:
                        minB0=B0cur
                    if B0cur>maxB0:
                        maxB0=B0cur
                    if B1cur<minB1:
                        minB1=B1cur
                    if B1cur>maxB1:
                        maxB1=B1cur
                    if B2cur<minB2:
                        minB2=B2cur
                    if B2cur>maxB2:
                        maxB2=B2cur
                    if B3cur<minB3:
                        minB3=B3cur
                    if B3cur>maxB3:
                        maxB3=B3cur
                
    
    print("Min, max x1= ",x1_array[0],x1_array[-1])
    print("Min, max x2= ",x2_array[0],x2_array[-1])
    print("Min, max x3= ",x3_array[0],x3_array[-1])
    if (velocity):
        print("Min, max v1= ",minv1,maxv1)
        print("Min, max v2= ",minv2,maxv2)
        print("Min, max v3= ",minv3,maxv3)
    if (magnetic_field):
        print("Min, max B0= ",minB0,maxB0)
        print("Min, max B1= ",minB1,maxB1)
        print("Min, max B2= ",minB2,maxB2)
        print("Min, max B3= ",minB3,maxB3)
    
    if (fluid):
        print("Min, max rho= ",minrho,maxrho)
        print("Min, max T= ",'{:e}'.format(minT),'{:e}'.format(maxT))
    else:
        print("TO BE ADDED")



# At this point you should have in any case the time and the two spatial 1D-arrays
# For fluid simulations, you must have at least the DENSITY and TEMPERATURE 2D-arrays
# For kinetic simulations, you must have at least the TOTAL_EMISSION (J_I) 4D-array (2 spatial dimensions + frequency + pitch angle)
# The other radiative transfer coefficients are optional
# Optional arrays are : VELOCITY and MAGNETIC_FIELD


# ========================================================================
# CONVERTING
#
# DO NOT CHANGE ANYTHING BELOW
#

if not os.path.exists(outputdir):
    os.makedirs(outputdir)
    
for ii in range(nfile):
    # Rearranging stuff for gyoto
    timegyoto = gyoto.core.array_double_fromnumpy1(time_array)
    x1gyoto = gyoto.core.array_double_fromnumpy1(x1_array)
    x2gyoto = gyoto.core.array_double_fromnumpy1(x2_array)
    x3gyoto = gyoto.core.array_double_fromnumpy1(x3_array)
    
    
    if (velocity):
        v11D = gyoto.core.array_double_fromnumpy3(vel[ii,0])
        v21D = gyoto.core.array_double_fromnumpy3(vel[ii,1])
        v31D = gyoto.core.array_double_fromnumpy3(vel[ii,2])
    
    if (magnetic_field):
        B01D = gyoto.core.array_double_fromnumpy3(Bfield[ii,0])
        B11D = gyoto.core.array_double_fromnumpy3(Bfield[ii,1])
        B21D = gyoto.core.array_double_fromnumpy3(Bfield[ii,2])
        B31D = gyoto.core.array_double_fromnumpy3(Bfield[ii,3])
    
    fits = gyoto.FitsRW()
    filename = outputdir+fitsname_prefix+str(ii).zfill(4)+'.fits' # MUST NOT BE CHANGED
    
    fptr = fits.fitsCreate(filename)
    
    fits.fitsWriteKey(fptr, "NB_X0", nfile)
    fits.fitsWriteKey(fptr, "NB_X1", nx1)
    fits.fitsWriteKey(fptr, "NB_X2", nx2)
    fits.fitsWriteKey(fptr, "NB_X3", nx3)
    fits.fitsWriteKey(fptr, "TIME", time_array[ii])
    fits.fitsWriteKey(fptr, "BINFILE", magnetic_field)
    
    fits.fitsWriteHDUData(fptr, "X0", timegyoto, nfile)
    fits.fitsWriteHDUData(fptr, "X1", x1gyoto, nx1)
    fits.fitsWriteHDUData(fptr, "X2", x2gyoto, nx2)
    fits.fitsWriteHDUData(fptr, "X3", x3gyoto, nx3)
    
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
        density1D = gyoto.core.array_double_fromnumpy3(rho[ii])
        temperature1D = gyoto.core.array_double_fromnumpy3(temperature[ii])
        fits.fitsWriteHDUData(fptr, "NUMBERDENSITY", density1D, ntot)
        fits.fitsWriteHDUData(fptr, "TEMPERATURE", temperature1D, ntot)
    else:
        print("TO BE ADDED")
        #J_I1D = gyoto.core.array_double_fromnumpy5(J_I[ii])
        #fits.fitsWriteHDUData(fptr, "J_I", density1D, ntot)
    fits.fitsClose(fptr)
    print("")
