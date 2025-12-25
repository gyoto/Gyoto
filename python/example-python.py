#/bin/env python
# -*- coding: utf-8 -*-
#
# Example file for gyoto: writing new Metric/Astrobj objects in Python
#
# Copyright 2025 Thibaut Paumard
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

# standard imports
import sys
import os
import numpy
import matplotlib as ml
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# Importing gyoto also imports gyoto.python but ignores failures
# silently. Here we want an early explicit message is gyoto.python
# cannot be imported, so let's import it explicitely.
import gyoto.python


# Parse command line and optionally switch to PDF output

pdfname=None
examples_dir="/../doc/examples/"
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

# Let's see how to load a class from an external module. The module
# has to be somewhere in sys.path. The examples that we will load are
# in doc/examples/gyoto_sample_{standard,thindisk,,metrics,spectra}.py
# relative to the top source directory. Read these files for a
# detailed overview of how to create a custom Gyoto class in Python,
# here we will cover the basics. First lest's add this directory to
# the path:
sys.path.insert(0, examples_dir)

## Instanciate a metric
# First create the C++ wrapper:
gg= gyoto.python.PythonMetric()

# Let this wrapper import our external module:
gg.Module = "gyoto_sample_metrics"

# Get the right class from this module:
gg.Class = "Minkowski"
# It is important to understand that "gg" is *not* an instance of the
# Python class gyoto_sample_metrics.Minkowski. It is an instance of
# gyoto.python.PythonMetric, which internally stores an instance of
# gyoto_sample_metrics.Minkowski. This instance is not directly
# accessible.

# This metric supports both spherical and Cartesian coordinates. We
# can use either Cartesian coordinates:
gg.Spherical = False
# or Spherical coordinates:
gg.Spherical = True

# Let's instanciate a FixedStar from gyoto_sample_standard.py:
# first create the C++ wrapper:
ao = gyoto.python.PythonStandard()

# Let this wrapper import our external module:
ao.Module = "gyoto_sample_standard"

# Get the right class from this module:
ao.Class = "FixedStar"
# Again, it is important to understand that "ao" is *not* an instance
# of the Python class gyoto_sample_standard.FixedStar. It is an
# instance of gyoto.python.PythonStandard, which internally stores an
# instance of gyoto_sample_standard.FixedStar. This instance is not
# directly accessible.

# Instead of relying on an external module, we can provide an inline
# module as text:
ao.Class = ""
ao.InlineModule = '''
import math
class FlaredDisk:
    opening=0.2
    rin=4
    rout=15
    def __call__(self, coord):
        r=math.sin(coord[2])*coord[1]
        h_r=abs(math.cos(coord[2]))
        return max(h_r-self.opening, self.rin-r, r-self.rout)
    def getVelocity(self, coord, vel):
        self.this.metric().circularVelocity(coord, vel)
'''
ao.Class = "FlaredDisk"

# And instead of providing a module (external or inline) and a class
# name, we can provide a class instance:
import math
class FlaredDisk:
    properties = {"Opening": "double",
                  "Rin": "double",
                  "Rout": "double"}
    opening=0.2
    rin=4
    rout=15
    def set (self, key, val):
        if key in self.properties:
            setattr(self, key.lower(), val)
    def get (self, key):
        if key in self.properties:
            return getattr(self, key.lower())
    def __call__(self, coord):
        r=math.sin(coord[2])*coord[1]
        h_r=abs(math.cos(coord[2]))
        return max(h_r-self.opening, self.rin-r, r-self.rout)
    def getVelocity(self, coord, vel):
        self.this.metric().circularVelocity(coord, vel)
instance = FlaredDisk()
ao.Instance = gyoto.core.gyotoid(instance)

print(ao)

# Attach the metric to the star. This will telle the star that we are
# working in Cartesian coordinates.
ao.Metric = gg

# Since this class declares properties, we can set them:
ao.Opening = 0.15 # set disk opening angle
ao.Rin     = 3.   # set disk inner radius

# It is a "Standard" object, defined by some function fcn such that
# the interior of the object is less than a critical value inside the
# object. This function is implemented in __call__ and the critical
# value lives in the C++ member variable accessible as
# ao.CriticalValue.
assert (ao((0, 4, 1.57, 0)) < ao.CriticalValue)
assert (ao((0, 4, 0   , 0)) > ao.CriticalValue)

## Lets try with base class
from gyoto_sample_standard import FlaredDisk

ao=FlaredDisk()

## A complete example, identical to
# example-python-flared-disk-derived-kerrbl.xml:

sc=gyoto.Scenery()

sc.Metric=gyoto.metric.KerrBL()
sc.Metric.Mass = 4e6, "sunmass"
sc.Metric.Spin = 0.52

sc.Screen = gyoto.core.Screen()
sc.Screen.Distance = 8, "kpc"
sc.Screen.Inclination = 90, "degree"
sc.Screen.PALN = 180, "degree"
sc.Screen.Time = 8, "kpc"
sc.Screen.FieldOfView = 200, "microas"
sc.Screen.Resolution = 32
sc.Screen.Spectrometer = gyoto.spectrometer.wave()
sc.Screen.Spectrometer.NSamples = 1
sc.Screen.Spectrometer.Band = 2.0e-6, 2.4e-6

sc.Astrobj = FlaredDisk()
sc.Astrobj.RMax = 50
sc.Astrobj.OpticallyThin = True

sc.Quantities = "Intensity"

# print class and module name:
print(ao.Module)
print(ao.Class)

# Ray-trace for the entire screen:
data=sc[:,:]

# Plot:
plt.imshow(data["Intensity"])
if pdf is None:
    plt.show()
else:
    pdf.savefig()
    plt.close()

## Try a thindisk

sc = gyoto.util.readScenery(examples_dir+"example-thin-disk.xml")

from gyoto_sample_thindisks import ThinDisk

sc.Screen.Time = 1000, "geometrical_time"
sc.Screen.PALN (0.)
sc.Astrobj = ThinDisk()
sc.Astrobj.InnerRadius = 10
sc.Astrobj.RMax = 1000
sc.Astrobj.Spectrum = gyoto.spectrum.BlackBody()
sc.Astrobj.Spectrum.Temperature = 6000
sc.Astrobj.OpticallyThin = True

# print class and module name:
print(sc.Astrobj.Module)
print(sc.Astrobj.Class)
print(sc)

# ray-trace
data=sc[:,:]

# Plot:
plt.imshow(numpy.sqrt(data["Intensity"]))
if pdf is None:
    plt.show()
else:
    pdf.savefig()
    plt.close()

# Cleanup and exit
print("All done, exiting")
if pdf is not None:
    pdf.close()
