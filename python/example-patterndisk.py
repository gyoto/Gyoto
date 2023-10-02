#!/usr/bin/env python3
#
# In this example script, we create a Gyoto Scnenery with a
# PatternDisk, save it, read it and check that the re-read scenery is
# identical to the saved one.
# This is the same example as yorick/check-patterndisk.i

import gyoto.core, gyoto.std
import numpy
import os
import matplotlib.pyplot as plt

### Create a metric
metric = gyoto.std.KerrBL()
metric.mass(4e6, "sunmass");

### Create PatternDisk
# Create opacity and intensity grids as numpy arrays.
# Get pointers in a format that Gyoto undestands.
# Warning: here we assume that size_t is the same as uint64.
gridshape=numpy.asarray( (1, 3, 11) , numpy.uint64)
pgridshape=gyoto.core.array_size_t.fromnumpy1(gridshape)

opacity=numpy.zeros(gridshape)
popacity=gyoto.core.array_double.fromnumpy3(opacity)
opacity[:, 0::2, 0::2]=100.
opacity[:, 1::2, 1::2]=100.

intensity=opacity*0.+1.;
pintensity=gyoto.core.array_double.fromnumpy3(intensity)

# Create PatternDisk, attach grids, set some parameters
pd=gyoto.std.PatternDisk()
pd.copyIntensity(pintensity, pgridshape)
pd.copyOpacity  (popacity, pgridshape)
pd.innerRadius(3)
pd.outerRadius(28)
pd.repeatPhi(8)
pd.metric(metric)
pd.rMax(50)

### Create screen
screen=gyoto.core.Screen()
screen.metric(metric)
screen.resolution(64)
screen.time(1000., "geometrical_time")
screen.distance(100., "geometrical")
screen.fieldOfView(30./100.)
screen.inclination(110., "degree")
screen.PALN(180., "degree")

### Create Scenery
sc=gyoto.core.Scenery()
sc.metric(metric)
sc.screen(screen)
sc.astrobj(pd)

### Save Scenery
pd.fitsWrite("!check-patterndisk.fits.gz")
gyoto.core.Factory(sc).write("check-patterndisk.xml")

### Read Scenery
sc2=gyoto.core.Factory("check-patterndisk.xml").scenery()

### Check
# Compare Sceneries
assert sc2.screen().dMax() == sc.screen().dMax(), "dmax was not conserved when writing and reading XML"
assert sc2.tMin() == sc.tMin(), "tmin was not conserved when writing and reading XML"

# Delete temporary files
os.unlink("check-patterndisk.xml")
os.unlink("check-patterndisk.fits.gz")

# Compare PatternDisks
# compare shape
pd2 = gyoto.std.PatternDisk(sc2.astrobj())
pgridshape2=gyoto.core.array_size_t(3)
pd2.getIntensityNaxes(pgridshape2)
for k in range (3):
    assert pgridshape2[k]==pgridshape[k], "shape of grid changed"
bufsize=gridshape.prod()
# compare intensity
buf=gyoto.core.array_double.frompointer(pd2.getIntensity())
for k in range(bufsize):
    assert buf[k] == pintensity[k], "Intensity changed"
# compare opacity
buf=gyoto.core.array_double.frompointer(pd2.opacity())
for k in range(bufsize):
    assert buf[k] == popacity[k], "Opacity changed"

### Ray-trace
ii=gyoto.core.Range(1, screen.resolution(), 1)
jj=gyoto.core.Range(1, screen.resolution(), 1)
grid=gyoto.core.Grid(ii, jj)
aop=gyoto.core.AstrobjProperties()
frame=numpy.zeros((screen.resolution(), screen.resolution()))
pframe=gyoto.core.array_double.fromnumpy2(frame)
aop.intensity=pframe
sc.rayTrace(grid, aop)
plt.imshow(frame, origin='lower')
plt.show()
