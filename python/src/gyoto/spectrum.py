'''Gyoto::Spectrum namespace

In order to emulate the C++ Gyoto::Spectrum namespace, this module will
load gyoto.std and gyoto.lorene (if available) and expose all Spectra
in here.

'''

import gyoto._namespaces as _namespaces
from gyoto.core import Spectrum as Generic
__all__ = _namespaces.make_namespace(Generic, globals())
del _namespaces

def __getattr__(name):
    '''Allows instanciating any spectrum kind

    Calling
      gyoto.spectrum.Kind([pluglist])
    is equivalent to
      gyoto.spectrum.Generic("Kind" [, pluglist])
    
    '''
    return lambda *pluglist : Generic(name, *pluglist)
