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
examples_dir="../doc/examples/"
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

# This metric supports both spherical and Cartesian coordinates. Let's
# use Cartesian coordinates.
gg.Spherical = False

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

# Attach the metric to the star. This will telle teh star that we are
# working in Cartesian coordinates.
ao.Metric = gg

# Since this class implements set() and get(), we can use them:
ao.set("Position", (1., 2., 3.)) # set the star's position
ao.set("Radius", 0.5)            # set the star's position
# Unfortunately we cannot yet simply do ao.Position = (1.,2.,3.)

# It is a "Standard" object, defined by some function fcn such that
# the interior of the object is less than a critical value inside the
# object. This function is implemented in __call__ and the critical
# value lives in the C++ member variable accessible as
# ao.CriticalValue. To get access to the __call__ function, we first
# need to cast ao to a StandardAstrobj.
ao=gyoto.core.StandardAstrobj(ao)
assert (ao((0, 1.2, 2.1, 3.)) < ao.CriticalValue)
assert (ao((0, 1.5, 3, 2.)) > ao.CriticalValue)


## A complete example, identical to example-python-flared-disk.xml:

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

sc.Astrobj = gyoto.python.PythonStandard()
sc.Astrobj.RMax = 50
sc.Astrobj.CriticalValue = 0.
sc.Astrobj.SafetyValue = 0.3
sc.Astrobj.OpticallyThin = True
# This time the module is inline, not in an external file:
sc.Astrobj.InlineModule = '''
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
sc.Astrobj.Class = "FlaredDisk"

sc.Quantities = "Intensity"

# Ray-trace for the entire screen:
data=sc[:,:]

# Plot:
plt.imshow(data["Intensity"])
if pdf is None:
    plt.show()
else:
    pdf.savefig()
    plt.close()

# Cleanup and exit
print("All done, exiting")
if pdf is not None:
    pdf.close()
