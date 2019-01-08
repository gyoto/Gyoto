#/bin/env python
# -*- coding: utf-8 -*-
# Example file for gyoto with MPI
#
# Copyright 2015-2018 Thibaut Paumard
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
# (ray-tracing of individual geodesics) to 'worker' processes. The
# workers are processes that run the static method
# Gyoto::Scenery::mpiWorker(). This can be done directly from Python
# or by involing instances of gyoto-mpi-worker.VERSION where VERSION
# is the ABI version of Gyoto.
#
# There are three basic scenarios supported in Gyoto to start the
# manager and the workers:
#
# 1- mpirun starts only the rank 0, manager process, which is
#    responsible for starting the workers processes using
#    MPI_Comm_spawn();
#
# 2- mpirun starts one instance of the manager process and NWORKER
#    instances of the worker process;
#
# 3- mpirun starts NP instances of a process; rank 0 becomes the
#    manager while rank > 0 invoke Scenery.mpiWorker().
#
### Synopsis:
#
# This script automatically detects how it was started and supports
# the three scenarios above:
#
# 1- mpirun -np 1 python example-mpi.py
# 2- mpirun -np 1 python example-mpi.py : -np NWORKERS gyoto-mpi-worker.VERSION
# 3- mpirun -np <NWORKERS+1> python example-mpi.py
#
# where NWORKERS is the desired number of workers to use (e.g. 4) and
# VERSION is the gyoto ABI version (e.g. 5 or 5-unreleased).
#
### Environment:
#
# The system must be able to find libgyoto, the plugins, and (in
# scenario 1) the command-line utility gyoto-mpi-workers.VERSION.
# Therefore you may need to set PATH and LD_LIBRARY_PATH appropriately
# or the like.
#
###

# 1- Let mpi4py initialize the MPI environment:
import mpi4py.MPI

# 2- Prepare Gyoto::Scenery to ray-trace
import numpy
import matplotlib as ml
import matplotlib.pyplot as plt
import gyoto.core
import gyoto.std

sc=gyoto.core.Factory("../doc/examples/example-moving-star.xml").scenery()
sc.nThreads(1)
sc.astrobj().opticallyThin(False)

# 3- Autodetect scenario
# Spawn processes and clone scenery into them:
world = mpi4py.MPI.COMM_WORLD
if (world.size == 1):
    # This is scenario 1
    # We spawn 4 instances of gyoto-mpi-worker.VERSION
    sc.mpiSpawn(4)
else:
    if (world.rank == 0):
        # This is scenario 2 or 3, this process is the manager
        print("Rank ", world.rank, " becoming manager")
        sc.mpiSpawn(-1)
    else:
        # This is scenario 3, this process is a worker
        print("Rank ", world.rank, " becoming worker")
        sc.mpiWorker()
        print("Rank ", world.rank, " terminating")
        exit()

# 4- Prepare storage for ray-traced quantities

# Prepare array for holding results
res=sc.screen().resolution()
intensity=numpy.zeros((res, res), dtype=float)
time=numpy.zeros((res, res), dtype=float)
distance=numpy.zeros((res, res), dtype=float)

# Store array pointers in AstrobjProperties
aop=gyoto.core.AstrobjProperties()
aop.intensity=gyoto.core.array_double.fromnumpy2(intensity)
aop.time=gyoto.core.array_double.fromnumpy2(time)
aop.distance=gyoto.core.array_double.fromnumpy2(distance)

# 5- Prepare Coord2dSet to select what Photons to launch
ii=gyoto.core.Range(1, res, 1)
jj=gyoto.core.Range(1, res, 1)
grid=gyoto.core.Grid(ii, jj)

sc.mpiClone()

# 6- Ray-trace
sc.rayTrace(grid, aop)

# 7- Terminate MPI workers
#
# If the workers were spawned, this is done automatically when the
# Scenery object sc is deleted (but it does not harm to do it by hand).
#
# If the processes were created by mpirun, this needs to be done
# exactly once. Don't try to use the processes after calling
# mpiTerminate().
#
# Note that the same workers, if not spawned, can be reused in
# anorther Scenery. They should really be terminated only after the
# last MPI ray-tracing has been perfomed.
sc.mpiTerminate()

# 8- Do something with the data
# plt.imshow(intensity)
# plt.show()
# plt.imshow(time)
# plt.show()
# plt.imshow(distance)
# plt.show()

print("All done, exiting")
