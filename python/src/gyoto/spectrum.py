'''Gyoto::Spectrum namespace

In order to emulate the C++ Gyoto::Spectrum namespace, this module will
load gyoto.std and gyoto.lorene (if available) and expose all Spectra
in here.

'''
import sys

from . import _namespaces
from .core import Spectrum as Generic
__all__ = _namespaces.make_namespace(Generic, globals())
del _namespaces

from . import core, utils

def __getattr__(name):
    '''Allows instanciating any spectrum kind

    Calling
      gyoto.spectrum.Kind()
    is equivalent to
      gyoto.spectrum.Generic("Kind")
    
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
    except core.Error:
        raise AttributeError(f"module '{__name__}' has no attribute '{name}'")

    # Make a class and cache it in the namespace
    klass = utils.make_class(sys.modules[__name__], name, None, None, __name__)
    setattr(sys.modules[__name__], name, klass)
    sys.modules[__name__].__all__.append(name)

    # Also return it
    return klass
