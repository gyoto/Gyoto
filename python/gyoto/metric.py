'''Gyoto::Metric namespace

In order to emulate the C++ Gyoto::Metric namespace, this module will
load gyoto.std and gyoto.lorene (if available) and expose all Metrics
in here.

'''

import gyoto._namespaces as _namespaces
from gyoto.core import Metric as Generic
__all__ = _namespaces.make_namespace(Generic, globals())
del _namespaces
