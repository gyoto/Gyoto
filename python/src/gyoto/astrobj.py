'''Gyoto::Astrobj namespace

In order to emulate the C++ Gyoto::Astrobj namespace, this module will
load gyoto.std and gyoto.lorene (if available) and expose all Astrobjs
in here.

'''

import sys

import gyoto._namespaces as _namespaces
from gyoto.core import Astrobj as Generic
__all__ = _namespaces.make_namespace(Generic, globals())
del _namespaces
Complex=ComplexAstrobj

import gyoto.core, gyoto.util

def __getattr__(name):
    '''Allows instanciating any astrobj kind

    Calling
      gyoto.astrobj.Kind()
    is equivalent to
      gyoto.astrobj.Generic("Kind")
    
    '''
    # __getattr__ shouldn't be called in that case, but still the
    # right answer:
    if name in sys.modules[__name__].__dict__:
        return sys.modules[__name__].__dict__[name]

    # Take care of standard attributes
    if name.startswith('__') and name.endsswith('__'):
        raise AttributeError(f"module '{__name__}' has no attribute '{name}'")

    # Check that a class by that name is registered
    try:
        obj = Generic(name)
    except gyoto.core.Error:
        raise AttributeError(f"module '{__name__}' has no attribute '{name}'")

    # Make a constructor and cache it in the namespace
    constructor = gyoto.util.make_constructor(sys.modules[__name__], name)
    setattr(sys.modules[__name__], name, constructor)
    sys.modules[__name__].__all__.append(name)

    # Also return it
    return constructor
