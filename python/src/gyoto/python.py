'''Implement Gyoto Metric, Astrobj and Spectrum and objects in Python

Using this module, Gyoto can be extended directly in Python, without
the need to build a custom plugin. Plenty of examples are given in the
source of Gyoto under plugins/python/doc/examples:
https://github.com/gyoto/Gyoto/tree/master/plugins/python/doc/examples

Example:
========

import gyoto

metric = gyoto.metric.KerrBL()
metric.Mass = 4e6, "sunmass"
metric.Spin = 0.52

flareddisk = gyoto.python.PythonStandard()
flareddisk.Metric = metric
flareddisk.RMax = 50.
flareddisk.CriticalValue = 0.
flareddisk.SafetyValue = 0.3
flareddisk.OpticallyThin = True
flareddisk.InlineModule="""
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
"""
flareddisk.Class = "FlaredDisk"

spectrometer=gyoto.spectrometer.wave()
spectrometer.NSamples = 1
spectrometer.Band = 2.0e-6, 2.4e-6

screen = gyoto.core.Screen()
screen.Metric = metric
screen.Distance = 8., "kpc"
screen.Inclination = 90., "degree"
screen.PALN = 180., "degree"
screen.Time = 8., "kpc"
screen.FieldOfView = 200., "microas"
screen.Resolution = 32
screen.Spectrometer = spectrometer

scenery=gyoto.core.Scenery()
scenery.Metric = metric
scenery.Astrobj = flareddisk
scenery.Screen = screen
scenery.Quantities = "Intensity"

image = scenery[:,:]["Intensity"]

'''

from . import core, util, metric, astrobj, spectrum, spectrometer
import sys
import os.path


# First check whether the plug-in is already loaded, which is the case
# if the user has added a python plug-in to GYOTO_PLUGINS, or if this
# instance of Python is actually running inside gyoto.

pluglist = core.vector_string(())
sctr = core.getMetricSubcontractor("Python", pluglist, 1)
if sctr:
    plugin=pluglist[1]
else:
    plugin=None
del sctr

# If not already loaded, the name of the Gyoto plug-in that can be
# loaded be the same as the name of the Python executable. Let's try
# it, as well as python3 and python as fallbacks.
if (plugin is None):
    pluglist=(os.path.basename(os.path.realpath(sys.executable)),
              os.path.basename(sys.executable),
              "python3", "python")
    for plugin in pluglist:
        try:
            core.requirePlugin(plugin)
            break
        except core.Error:
            plugin = None

if plugin is None:
    raise core.Error("Could not load Python plugin, tried: "
                     +repr(pluglist)+", python3 and python")

del pluglist

__all__=[]

# Wrap the classes implemented in the plugin. Also add them to
# gyoto.metric, gyoto.astrobj and gyoto.spectrum.
for namespace, clsname, identifier in ((metric, "Python", "PythonMetric"),
                                       (astrobj, "Python::Standard",
                                        "PythonStandard"),
                                       (astrobj, "Python::ThinDisk",
                                        "PythonThinDisk"),
                                       (spectrum, "Python", "PythonSpectrum")):
    klass = util.make_class(namespace, clsname, plugin, identifier, __name__)
    setattr(sys.modules[__name__], identifier, klass)
    __all__.append(identifier)
    setattr(namespace, identifier, klass)
    namespace.__all__.append(identifier)
