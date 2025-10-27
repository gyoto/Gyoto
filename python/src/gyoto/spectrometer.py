'''Gyoto::Spectrometer namespace

In order to emulate the C++ Gyoto::Spectrometer namespace, this module will
load gyoto.std and gyoto.lorene (if available) and expose all Spectrometers
in here.

'''

import gyoto._namespaces as _namespaces
from gyoto.core import Spectrometer as Generic
__all__ = _namespaces.make_namespace(Generic, globals())
del _namespaces
Complex=ComplexSpectrometer
Uniform=UniformSpectrometer
