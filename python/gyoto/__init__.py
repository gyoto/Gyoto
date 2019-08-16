"""The General relativitY Orbit Tracer of paris Observatory

Note that importing "gyoto" is deprecated and may cease to work in a
future release. Please update your code to import gyoto.core instead.

"""
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
