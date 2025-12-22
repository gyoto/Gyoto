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
    plugin=pluglist[0]
else:
    plugin=None
del sctr

# If not already loaded, the name of the Gyoto plug-in that can be
# loaded should be the same as the name of the Python
# executable. Let's try it, as well as 'python3' and 'python' as
# fallbacks.
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
for namespace, clsname, identifier, base in (
      (metric, "Python", "PythonMetric", None),
      (astrobj, "Python::Standard", "PythonStandard", astrobj.StandardAstrobj),
      (astrobj, "Python::ThinDisk", "PythonThinDisk", astrobj.ThinDisk),
      (spectrum, "Python", "PythonSpectrum", None)):
    klass = util.make_class(namespace, clsname, plugin, identifier, __name__,
                            base)
    setattr(sys.modules[__name__], identifier, klass)
    __all__.append(identifier)
    setattr(namespace, identifier, klass)
    namespace.__all__.append(identifier)

class PythonBase():
    '''Base class for Gyoto object implemented in Python

    This class is meant as a base class for implementing Gyoto objects
    in the Python language using. It handles properties in a
    consistent manner.

    '''

    properties = dict()

    def __init__(self, base, *args):
        '''Initialize instance

        Set default values of properties including 'Instance'.

        '''
        self.__dict__["_PythonPluginClass"] = base
        self._PythonPluginClass.__init__(self, *args)
        super(self._PythonPluginClass, self).set("Instance", core.gyotoid(self))
        for key in self.properties:
            if hasattr(type(self), key):
                self.set(key, getattr(type(self), key))
            else:
                raise(core.Error(f"Please provide default value for '{key}'"))

    def _key2attr(self, key):
        '''Transform an property key to an attibute name

        return '_'+key.lower()

        This is not meant to be reimplemented in derived class.

        '''
        return '_'+key.lower()

    def set (self, key, *args):
        '''Set a Gyoto property

        If key is listed in the self.properties dictionary, set the
        Python attribute named '_'+key.lower() in self (direclty in
        self.__dict__, bypassing __setattr__). Else, forward the call
        to the underlying C++ instance.

        Derived classes may reimplement this method but should not
        depart too far from this logic.

        '''
        if key in self.properties:
            self.__dict__[self._key2attr(key)] = args[0]
        else:
            super(self._PythonPluginClass, self).set(*args)

    def get (self, key):
        '''Get a Gyoto property

        If key is listed in the self.properties dictionary, get the
        Python attribute named '_'+key.lower() in self. Else, forward
        the call to the underlying C++ instance.

        Derived classes may reimplement this method but should not
        depart too far from this logic.

        '''
        if key in self.properties:
            return self.__dict__[self._key2attr(key)]
        else:
            return super(self._PythonPluginClass, self).get(key)

class StandardBase(PythonBase, PythonStandard):
    '''Base class for Gyoto Standard Astrobjs implemented in Python

    This class is meant as a base class for implementing Gyoto
    astronomical objects in the Python language using
    gyoto.python.PythonStandard (from which it derives):

    In a nutshell:
    class MyAstrobj(StandardBase):
        properties={ 'Property1': 'type1',
                     'property2': 'type2', ... }

        Property1 = deafult1
        Property2 = deafult2

        def __call__(self, position):
            # position is inside the object if and only if
            # retval < self.CriticalValue
            return retval

    ao = MyAstrobj()

    sc=gyoto.core.Scenery()

    sc.Astrobj = ao
    ...

    Such classes can also be used from XML files:
    <Astrobj kind = "PythonStandard" plugin="fallback:python*"
      <Module>module_where_MyAstrobj_is_defined</Module>
      <Class>MyAstrobj</Class>
      any other PythonStandard Property, e.g. CriticalValue, RMax etc.
    </Astrobj>

    Derived classes must implement __call__(), should generally
    implement getVelocity(), and may implement emission(),
    transmission(), giveDelta() and integrateEmission(). __call__()
    should work in both spherical and Cartesian coordinates, the
    helper method coordKindIsSpherical() is provided to help determine
    this.

    StandardBase.set and StandardBase.get simply accept any key in the
    properties dict and set/get the corresponding attribute with name
    '_'+key.lower() (directly in self.__dict__, by-passing
    setattr). Note that these attributes are protected and cannot
    be set directly through setattr.

    '''

    properties = dict()

    def __init__(self, *args):
        '''Initialize instance

        1- initialize the underlying PythonStandard instance;
        2- set the Instance Property in the underlying PythonStandard object.

        '''
        PythonBase.__init__(self, PythonStandard, *args)

    def getVelocity(self, coord, vel):
        '''Get astronomical object velocity field at a given position

        The StandardBase implementation returns the circularVelocity
        as implemented in the underlying metric, which is reasonable
        for equatorial disks but not so much for any other object.

        Most derived classes should reimplement this.

        The StandardBase implementation can be used in two ways:

        instance.getVelocity(coord, vel)
        vel = instance.getVelocity(coord)

        Derived classes that reimplement it should support at least
        the first one.

        Parameters:

        coord: a numpy array of doubles of length at least 4
               containing the position at which to evaluate the
               velocity.

        vel: the four-velocity, in a numpy array which in the first
               form must be pre-allocated.

        '''
        self.this.Metric.circularVelocity(coord, vel)

    def coordKindIsSpherical(self):
        '''True if Metric.coordKind() is gyoto.core.GYOTO_COORDKIND_SPHERICAL
        '''
        return self.this.Metric.coordKind() == core.GYOTO_COORDKIND_SPHERICAL

class ThinDiskBase(PythonBase, PythonThinDisk):
    '''Base class for Gyoto ThinDisk Astrobjs implemented in Python

    This class is meant as a base class for implementing Gyoto
    astronomical objects in the Python language using
    gyoto.python.PythonThinDisk (from which it derives):

    In a nutshell:
    class MyAstrobj(ThinDiskBase):
        properties={ 'Property1': 'type1',
                     'property2': 'type2', ... }

        Property1 = deafult1
        Property2 = deafult2

        def __call__(self, position):
            # position is inside the object if and only if
            # retval < self.CriticalValue
            return retval

    ao = MyAstrobj()

    sc=gyoto.core.Scenery()

    sc.Astrobj = ao
    ...

    Such classes can also be used from XML files:
    <Astrobj kind = "PythonThinDisk" plugin="fallback:python*"
      <Module>module_where_MyAstrobj_is_defined</Module>
      <Class>MyAstrobj</Class>
      any other PythonThinDisk Property, e.g. InnerRadius, RMax etc.
    </Astrobj>

    Derived classes may implement __call__(), getVelocity(),
    emission(), transmission() and integrateEmission(). __call__()
    should work in both spherical and Cartesian coordinates, the
    helper method coordKindIsSpherical() is provided to help determine
    this.

    '''

    properties = dict()

    def __init__(self, *args):
        '''Initialize instance

        1- initialize the underlying PythonStandard instance;
        2- set the Instance Property in the underlying PythonStandard object.

        '''
        PythonBase.__init__(self, PythonThinDisk, *args)

    def coordKindIsSpherical(self):
        '''True if Metric.coordKind() is gyoto.core.GYOTO_COORDKIND_SPHERICAL
        '''
        return self.this.Metric.coordKind() == core.GYOTO_COORDKIND_SPHERICAL
