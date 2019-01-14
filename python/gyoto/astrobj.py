'''Gyoto::Astrobj namespace

In order to emulate the C++ Gyoto::Astrobj namespace, this module will
load gyoto.std and gyoto.lorene (if available) and expose all Astrobjs
in here.

'''

import gyoto._namespaces as _namespaces
from gyoto.core import Astrobj as Generic
__all__ = _namespaces.make_namespace(Generic, globals())
del _namespaces
Complex=ComplexAstrobj
