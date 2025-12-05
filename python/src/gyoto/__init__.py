"""The General relativitY Orbit Tracer of paris Observatory

Note that importing "gyoto" is deprecated and may cease to work in a
future release. Please update your code to import gyoto.core instead.

"""
import sys
import types
import importlib

# For backwards compatibility, expose gyoto.core as gyoto
from gyoto.core import *

# Provide a Pythonic wrapper around Scenery.rayTrace.
# The underlying C++-like interface remains accessible.
from gyoto import core, util
core.Scenery.rayTrace = util.rayTrace
core.Scenery.rayTrace.__doc__ += core._core.Scenery_rayTrace.__doc__
core.Scenery.__getitem__ = util.Scenery_getitem
core.Scenery.__getitem__.__doc__ += core.Scenery.rayTrace.__doc__

core.Worldline.getCartesian =  util._Worldline_getCartesian
core.Worldline.getCartesian.__doc__ = core._core.Worldline_getCartesian.__doc__
core.Worldline.getCoord =  util._Worldline_getCoord
core.Worldline.getCoord.__doc__ = core._core.Worldline_getCoord.__doc__

# import plugins as modules

def __getattr__(name):
    f'''Autoload submodules

    {__name__}.__getattr__ is called when Python encounters a call to
    an attribute that does not yet exist.

    It first tries to load an existing submodule. If it fails, it
    tries to load a Gyoto plugin by that name and expose it as a
    minimal submodule.

    '''
    print(f"in __getattr__")
    # __getattr__ shouldn't be called in that case, but still the
    # right answer:
    if name in sys.modules[__name__].__dict__:
        return sys.modules[__name__].__dict__[name]

    # Take care of standard attributes
    if name.startswith('__') and name.endsswith('__'):
        raise AttributeError(f"module '{__name__}' has no attribute '{name}'")

    # Now try to load a submodule, or build a submodule for a plugin
    fullname=f"{__name__}.{name}"
    try:
        # If a submodule by that name exists, load it
        mod = importlib.import_module(fullname)
        setattr(sys.modules[__name__], name, mod)
        return mod
    except ModuleNotFoundError:
        try:
            # Else try to load a plugin by that name 
            core.requirePlugin(name)
        except core.Error:
            # If still not successful, raise the correct error
            raise ModuleNotFoundError(f"No module named '{fullname}'")

    # We managed to load a plugin. Expose it as a module.
    def __getattr__(clsname):
       try:
           # check whether this is a Metric
           obj = core.Metric(clsname, (name,))
           return lambda : core.Metric(clsname, (name,))
       except core.Error:
           # check whether this is an Astrobj
           obj = core.Astrobj(clsname, (name,))
           return lambda : core.Astrobj(clsname, (name,))
       except core.Error:
           # check whether this is an Spectrum
           obj = core.Spectrum(clsname, (name,))
           return lambda : core.Spectrum(clsname, (name,))
       except core.Error:
           # check whether this is an Spectrometer
           obj = core.Spectrometer(clsname, (name,))
           return lambda : core.Spectrometer(clsname, (name,))
       except core.Error:
           raise AttributeError(f"module '{fullname}' has no attribute '{clsname}'")

    mod = types.ModuleType(fullname)
    mod.__getattr__=__getattr__
    sys.modules[fullname]=mod
    setattr(sys.modules[__name__], name, mod)
    
    return mod
