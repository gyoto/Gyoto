#/bin/env python
# -*- coding: utf-8 -*-
# Example file for gyoto with MPI
#
# Copyright 2015 Thibaut Paumard
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
# The script examplifies how to use MPI parallelization with Gyoto
# commanded from Python. It performs the ray-tracing of one of the
# Gyoto examples in doc/examples/ using MPI parallelization. It does
# nothing with the result, uncomment the last few lines if you want to
# display them.
#
# The parallel computing model of Gyoto is MIMD (Multiple Instruction,
# Multiple Data) which allows not caring too much about the details in
# the high level program. We have a single 'manager' process that
# takes care of handling input/output and distributes the work load
# (ray-tracing of individual geodesics) to 'worker' processes,
# instances of gyoto-mpi-worker.VERSION where VERSION is the ABI
# version of Gyoto.
#
# There are two basic scenarios supported in Gyoto to start the
# manager and the workers, plus a cheat that works only under Open
# MPI:
#
# 1- mpirun starts only the rank 0, manager process, which is
#    responsible for starting the workers processes using
#    MPI_Comm_spawn();
#
# 2- mpirun starts one instance of the manager process and NWORKER
#    instances of the worker process;
#
# 3- only under Open MPI: mpirun/orterun starts NWORKERS+1 instances
#    of the manager program, and all these instances except rank 0
#    immediately transform themselves into worker processes using the
#    execlp() system call, so we are back in scenario 2 above.
#
### Synopsis:
#
# This script automatically detects how it was started and supports
# the three scenarios above:
#
# 1- mpirun -np 1 python example-mpi.py
# 2- mpirun -np 1 python example-mpi.py : -np NWORKERS gyoto-mpi-worker.VERSION
# 3- orterun -np <NWORKERS+1> python example-mpi.py
#
# where NWORKERS is the desired number of workers to use (e.g. 4) and
# VERSION is the gyoto ABI version (e.g. 5 or 5-unreleased). In the
# first form, the script will spawn 4 workers (instances of
# gyoto-mpi-workers.VERSION). In the second form, the script will use
# the NWORKERS instances launched by mpirun.
#
# The third version works only if the MPI environment is Open MPI
# (hence the choice to use the command orterun instead of mpirun in
# the example). In this case, the NWORKERS process with rank > 0 will
# respawn themselves as gyoto-mpi-worker.VERSION using os.execlp().
#
### Environment:
#
# The system must be able to find libgyoto, the plugins, and the
# command-line utilities gyoto and gyoto-mpi-workers.VERSION,
# therefore you may need to set PATH and LD_LIBRARY_PATH appropriately
# or the like. Make sure the version of gyoto that is found first in
# your path is the right version, i.e. that LD_LIBRARY_PATH and PATH
# are set consistently.
#
###

# 1- Support scenario 3 above
#
# Check whether we are running under an OpenMPI environment and have a
# rank > 0. In this case, respawn as gyoto, which will itself respawn
# as the right gyoto-mpi-worker.VERSION using exactly the same
# heuristics. Respawning as gyoto instead of gyoto-mpi-worker.VERSION
# saves us from hard-coding VERSION, but is slightly unsafe. Make sure
# the right (ABI-compatible) version of gyoto is first in your PATH!
import os
if os.getenv('OMPI_COMM_WORLD_RANK', '0') != '0':
    os.execlp("gyoto", "")

# 2- Let mpi4py initialize the MPI environment:
import mpi4py.MPI

# 3- Prepare Gyoto::Scenery to ray-trace
import numpy
import matplotlib as ml
import matplotlib.pyplot as plt
import gyoto
gyoto.loadPlugin("stdplug")
import gyoto_std

a=gyoto.Factory("../doc/examples/example-moving-star.xml")
sc=a.getScenery()
sc.nThreads(1)
sc.astrobj().opticallyThin(False)

# 4- Prepare storage for ray-traced quantities

# Prepare array for holding results
res=sc.screen().resolution()
intensity=numpy.zeros((res, res), dtype=float)
time=numpy.zeros((res, res), dtype=float)
distance=numpy.zeros((res, res), dtype=float)

# Store array pointers in AstrobjProperties
aop=gyoto.AstrobjProperties()
aop.intensity=gyoto.array_double.fromnumpy2(intensity)
aop.time=gyoto.array_double.fromnumpy2(time)
aop.distance=gyoto.array_double.fromnumpy2(distance)

# 5- Prepare Coord2dSet to select what Photons to launch
ii=gyoto.Range(1, res, 1)
jj=gyoto.Range(1, res, 1)
grid=gyoto.Grid(ii, jj)

# 6- Autodetect scenario 1 vs. scenarios 2/3
# Spawn processes and clone scenery into them:
world = mpi4py.MPI.COMM_WORLD
if (world.size > 1 and world.rank == 0):
    # This is scenario 2/3
    # We assume that the other ranks are instances of
    # gyoto-mpi-worker.VERSION and enrole them
    sc.mpiSpawn(-1)
else:
    # This is scenario 1
    # We spawn 4 instances of gyoto-mpi-worker.VERSION
    sc.mpiSpawn(4)
sc.mpiClone()

# 7- Ray-trace
sc.rayTrace(grid, aop)

# 8- Terminate MPI workers
sc.mpiTerminate()

# 9- Do something with the data
#plt.imshow(intensity)
#plt.show()
#plt.imshow(time)
#plt.show()
#plt.imshow(distance)
#plt.show()

print("All done, exiting")
