#/bin/env python
# -*- coding: utf-8 -*-
# Example file for gyoto with MPI

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
aop.Intensity(intensity)
aop.EmissionTime(time)
aop.MinDistance(distance)

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
