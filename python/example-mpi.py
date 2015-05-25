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


# Let mpi4py initilize the MPI environment:
import mpi4py.MPI

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

# Ray-trace scenery

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

# Prepare Coord2dSet to select what Photons to launch
ii=gyoto.Range(1, res, 1)
jj=gyoto.Range(1, res, 1)
grid=gyoto.Grid(ii, jj)

# Spawn processes and clone scenery into them:
sc.mpiSpawn(4)
sc.mpiClone()

# Ray-trace
sc.rayTrace(grid, aop)

# Terminate workers
sc.mpiTerminate()

#plt.imshow(intensity)
#plt.show()
#plt.imshow(time)
#plt.show()
#plt.imshow(distance)
#plt.show()

print("All done, exiting")
