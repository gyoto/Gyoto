#/bin/env python
# -*- coding: utf-8 -*-
# Example file for gyoto
#
# Copyright 2014-2018 Thibaut Paumard
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

import sys
import os
import numpy
import matplotlib as ml
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import gyoto


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

# Simple stuff

scr=gyoto.core.Screen()
scr.Metric = gyoto.metric.KerrBL()
pos=scr.getObserverPos()

# Load Scenery

sc = gyoto.util.readScenery(examples_dir+"example-moving-star.xml")
sc.NThreads = 8
sc.Astrobj.OpticallyThin = False

scr=sc.Screen
dest=numpy.zeros(8, float)
scr.getRayTriad(1,1,dest)
dest=numpy.ndarray(3, float)
scr.coordToSky((0., 5., numpy.pi/2, 0), dest)
    
# Trace and plot NULL geodesic:

ph=gyoto.core.Photon()
ph.setInitialCondition(sc.Metric, sc.Astrobj, sc.Screen, 0., 0.)
ph.hit()
n=ph.get_nelements()

# We try to map Gyoto arrays to NumPy arrays wherever possible.

# Create NumPy arrays
t=numpy.ndarray(n)
r=numpy.ndarray(n)
theta=numpy.ndarray(n)
phi=numpy.ndarray(n)

# Call Gyoto method that takes these arrays as argument:
ph.get_t(t)
ph.getCoord(t, r, theta, phi)

plt.plot(t, r)
if pdf is None:
    plt.show()
else:
    pdf.savefig()
    plt.close()

# Trace and plot timelike geodesic
# We need to cast the object to a gyoto.std.Star:

wl=gyoto.astrobj.Star(sc.Astrobj)
wl.xFill(1000)

n=wl.get_nelements()

x=numpy.ndarray(n)
y=numpy.ndarray(n)
z=numpy.ndarray(n)

wl.get_xyz(x, y, z)

plt.plot(x, y)
if pdf is None:
    plt.show()
else:
    pdf.savefig()
    plt.close()


# Ray-trace scenery

# For that, we can use the short-hand:
sc.Quantities = 'Intensity EmissionTime MinDistance'
results=sc[:,:] # or: sc.rayTrace()

plt.imshow(results['Intensity'])
if pdf is None:
    plt.show()
else:
    pdf.savefig()
    plt.close()

plt.imshow(results['EmissionTime'])
if pdf is None:
    plt.show()
else:
    pdf.savefig()
    plt.close()

plt.imshow(results['MinDistance'])
if pdf is None:
    plt.show()
else:
    pdf.savefig()
    plt.close()


# Or we can do it manually to understand how the Gyoto API works:

res=sc.Screen.Resolution
intensity=numpy.zeros((res, res), dtype=float)
time=numpy.zeros((res, res), dtype=float)
distance=numpy.zeros((res, res), dtype=float)
aop=gyoto.core.AstrobjProperties()

# Here we will use the low-level AstrobjProperties facilities. This is
# one of a few Gyoto functionalities where NumPy arrays are not
# directly supported. We use lower-level C-like arrays through the
# gyoto.core.array_double and gyoto.core.array_unsigned_long classes. Beware
# that this type does not provide any safeguards, it is quite easy to
# get it to SEGFAULT. As we develop Gyoto, we try to remove the need
# for the gyoto.core.array_* classes in favor of NumPy arrays. Code that
# uses this... ``feature'' therefore may break in future releases.
#
# To (indirectly) use NumPy arrays with a functionality that requires
# gyoto.core.array_* arguments, create the arrays using numpy (see above:
# `intensity', `time' and `distance' arrays) , then cast them using
# the fromnumpyN static methods, where the digit N indicates the
# dimensionality of the NumPy array. The underlying storage belongs to
# the NumPy variable and will be deleted with it: don't use the
# array_double() variable (for anyting else that destroying it) past
# the destruction of the corresponding NumPy variable.

aop.intensity=gyoto.core.array_double.fromnumpy2(intensity)
aop.time=gyoto.core.array_double.fromnumpy2(time)
aop.distance=gyoto.core.array_double.fromnumpy2(distance)

ii=gyoto.core.Range(1, res, 1)
jj=gyoto.core.Range(1, res, 1)
grid=gyoto.core.Grid(ii, jj, "\rj = ")

sc.rayTrace(grid, aop)

plt.imshow(intensity)
if pdf is None:
    plt.show()
else:
    pdf.savefig()
    plt.close()

plt.imshow(time)
if pdf is None:
    plt.show()
else:
    pdf.savefig()
    plt.close()

plt.imshow(distance)
if pdf is None:
    plt.show()
else:
    pdf.savefig()
    plt.close()


# Another Scenery, with spectrum

sc=gyoto.util.readScenery(examples_dir+"example-polish-doughnut.xml")
sc.Screen.Resolution = 32
res=sc.Screen.Resolution
ns=sc.Screen.Spectrometer.NSamples
spectrum=numpy.zeros((ns, res, res), dtype=float)

ii=gyoto.core.Range(1, res, 1)
jj=gyoto.core.Range(1, res, 1)
grid=gyoto.core.Grid(ii, jj, "\rj = ")

aop=gyoto.core.AstrobjProperties()
aop.spectrum=gyoto.core.array_double.fromnumpy3(spectrum)
aop.offset=res*res

sc.rayTrace(grid, aop)

plt.imshow(spectrum[1,:,:])
if pdf is None:
    plt.show()
else:
    pdf.savefig()
    plt.close()


# Another Scenery, with impact coords, created from within Python

met=gyoto.metric.KerrBL()
met.Mass = 4e6, "sunmass"
ao=gyoto.astrobj.PageThorneDisk()
ao.Metric        = met
ao.OpticallyThin = False
ao.RMax          = 100
screen=gyoto.core.Screen()
screen.Distance    = 8, "kpc"
screen.Time        = 8, "kpc"
screen.Resolution  = 64
screen.Inclination = numpy.pi/4
screen.PALN        = numpy.pi
screen.Time        = 8, "kpc"
screen.FieldOfView = 100, "µas"
sc=gyoto.core.Scenery()
sc.Metric   = met
sc.Astrobj  = ao
sc.Screen   = screen
sc.Delta    = 1, "kpc"
sc.Adaptive = True
sc.NThreads = 8

res=sc.Screen.Resolution

ii=gyoto.core.Range(1, res, 1)
jj=gyoto.core.Range(1, res, 1)
grid=gyoto.core.Grid(ii, jj, "\rj = ")

ipct=numpy.zeros((res, res, 16), dtype=float)

aop=gyoto.core.AstrobjProperties()
aop.impactcoords=gyoto.core.array_double.fromnumpy3(ipct)
aop.offset=res*res

sc.rayTrace(grid, aop)

plt.imshow(ipct[:,:,0], interpolation="nearest", vmin=-100, vmax=0)
if pdf is None:
    plt.show()
else:
    pdf.savefig()
    plt.close()


# Trace one line of the above using alpha and delta

N=10

buf=numpy.linspace(screen.fieldOfView()*-0.5, screen.fieldOfView()*0.5, N)
a=gyoto.core.Angles(buf)
d=gyoto.core.RepeatAngle(screen.fieldOfView()*-0.5, N)
bucket=gyoto.core.Bucket(a, d)

ipct=numpy.zeros((N, 16), dtype=float)

aop=gyoto.core.AstrobjProperties()
aop.impactcoords=gyoto.core.array_double.fromnumpy2(ipct)
aop.offset=N

sc.rayTrace(bucket, aop)
plt.plot(buf, ipct[:,0])
if pdf is None:
    plt.show()
else:
    pdf.savefig()
    plt.close()


# Trace the diagonal of the above using i and j. The Range and Indices
# definitions below are equivalent.  Range is more efficient for a
# range, Indices can hold arbitrary indices.

ind=numpy.arange(1, res+1, dtype=numpy.uintp) # on 64bit arch...
ii=gyoto.core.Indices(ind)

# Or:
# ind=gyoto.core.array_size_t(res)
# for i in range(0, res):
#   ind[i]=i+1
# ii=gyoto.core.Indices(ind, res)

jj=gyoto.core.Range(1, res, 1)
bucket=gyoto.core.Bucket(ii, jj)

ipct=numpy.zeros((res, 16), dtype=float)

aop=gyoto.core.AstrobjProperties()
aop.impactcoords=gyoto.core.array_double.fromnumpy2(ipct)
aop.offset=res

sc.rayTrace(bucket, aop)

t=numpy.clip(ipct[:,0], a_min=-200, a_max=0)
plt.plot(t)
if pdf is None:
    plt.show()
else:
    pdf.savefig()
    plt.close()


# Any derived class can be instantiated from its name, as soon as the
# corresponding plug-in has been loaded into Gyoto. The standard
# plug-in is normally loaded automatically (and is always loaded when
# gyoto.std is imported), but this can also be forced with
# gyoto.core.requirePlugin():
gyoto.core.requirePlugin('stdplug')
tt=gyoto.core.Astrobj('Torus')
kerr=gyoto.core.Metric('KerrBL')

# If the name of the plug-in is suitable as a Python identifier
# (i.e. composed only of letters, digits, the underscore character,
# and does not start with a digit), it can be imported as a Python
# module and classes there-in can be used as expected, as long as they
# have unique names that are valid Python identifiers themselves:
import gyoto.stdplug
tt=gyoto.stdplug.Torus()
kerr=gyoto.stdplug.KerrBL()

# Most properties that can be set in an XML file can also be accessed
# from Python as object attributes:
kerr.Spin = 0.95
assert kerr.Spin == 0.95
# this include setting with unit:
kerr.Mass = 4e6, "sunmass"
# but to specify a unit for getting, we need to use the "get()"
# method:
assert kerr.get('Mass', 'sunmass') == 4e6

# Note that here we have loaded the stdplug plug-in as a user-provided
# plug-in, with access only to the generic API of metrics, astrobjs
# etc. By using the gyoto.std module instead, it becomes possible to
# access the methods specific to derived classes.
tr2=gyoto.std.Torus()
# and we can cast a generic pointer (from the gyoto extension) to a
# derived class:
tr=gyoto.std.Torus(tt)
assert tt.SmallRadius == tr.smallRadius()

# Another example: using a complex (i.e. compound) Astrobj. Note that
# Astrobj::Complex is renamed to ComplexAstrobj in the standard
# plug-in.
cplx=gyoto.std.ComplexAstrobj()
cplx.append(tr)
cplx.append(sc.astrobj())
sc.Astrobj = cplx
# the append() method is accessible when instantiating cplx from
# gyoto.std, not from gyoto.stdplug. Also, since stdplug contains both
# a metric and an astrobj that are called "Complex" and they have not
# been renamed, we need to be more specific about what we want:
cplx_generic = gyoto.stdplug.Astrobj('Complex')
try:
    cplx_generic.append(tr)
    worked = True
except AttributeError:
    worked=False
if worked: raise gyoto.core.Error("This should not have worked!")

print("All done, exiting")
if pdf is not None:
    pdf.close()
