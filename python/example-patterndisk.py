#!/usr/bin/env python3
#
# In this example script, we create a Gyoto Scnenery with a
# PatternDisk, save it, read it and check that the re-read scenery is
# identical to the saved one.
# This is the same example as yorick/check-patterndisk.i

import os
import sys
import numpy
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import gyoto.core, gyoto.std

# Parse command line and optionally switch to PDF output
pdfname=None
dir_path = os.path.dirname(os.path.realpath(__file__))
examples_dir=dir_path+"/../doc/examples/"

for param in sys.argv:
    sparam=param.split("=")
    if os.path.basename(sparam[0])==os.path.basename(__file__):
        pass
    elif sparam[0]=="--pdf":
        if len(sparam)==2:
            pdfname=sparam[1]
        else:
            raise ValueError('--pdf argument expects a filename, e.g. --pdf=output.pdf')
    elif sparam[0]=="--examples-dir":
        if len(sparam)==2:
            examples_dir=sparam[1]
        else:
            raise ValueError('--examples_dir argument expects a directory, e.g. --examples-dir=../doc/examples')
    else:
        raise ValueError(f'unknown argument: {sparam[0]}')

pdf=None if pdfname is None else PdfPages(pdfname)
if len(examples_dir) > 0 and examples_dir[-1] != "/":
    examples_dir += "/"

### Create a metric
metric = gyoto.std.KerrBL()
metric.Mass = 4e6, "sunmass"

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
pd.repeatPhi(8)
pd.InnerRadius =  3
pd.OuterRadius = 28
pd.Metric      = metric
pd.RMax        = 50

### Create screen
screen=gyoto.core.Screen()
screen.Metric      = metric
screen.Resolution  =   64
screen.Time        = 1000., "geometrical_time"
screen.Distance    =  100., "geometrical"
screen.FieldOfView =   30./100.
screen.Inclination =  110., "degree"
screen.PALN        =  180., "degree"

### Create Scenery
sc=gyoto.core.Scenery()
sc.Metric  = metric
sc.Screen  = screen
sc.Astrobj = pd

### Save Scenery
pd.fitsWrite("!check-patterndisk.fits.gz")
gyoto.core.Factory(sc).write("check-patterndisk.xml")

### Read Scenery
sc2=gyoto.util.readScenery("check-patterndisk.xml")

### Check
# Compare Sceneries
assert sc2.Screen.DMax == sc.Screen.DMax, "dmax was not conserved when writing and reading XML"
assert sc2.MinimumTime == sc.MinimumTime, "tmin was not conserved when writing and reading XML"

# Delete temporary files
os.unlink("check-patterndisk.xml")
os.unlink("check-patterndisk.fits.gz")

# Compare PatternDisks
# compare shape
pd2 = gyoto.std.PatternDisk(sc2.Astrobj)
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
if pdf is None:
    plt.show()
else:
    pdf.savefig()
    plt.close()

if pdf:
    pdf.close()
