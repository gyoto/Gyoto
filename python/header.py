# -*- coding: utf-8 -*-
"""
The General relativitY Orbit Tracer of paris Observatory
"""
# This is necessary at least on Linux to let the libgyoto symbols be
# visible from the libgyoto-stdplug (and other plug-ins) symbols.
from sys import getdlopenflags as _gyoto_tmp_getdlopenflags, setdlopenflags as _gyoto_tmp_setdlopenflags
from ctypes import RTLD_GLOBAL as _gyoto_tmp_RTLD_GLOBAL
_gyoto_tmp_setdlopenflags(_gyoto_tmp_getdlopenflags() | _gyoto_tmp_RTLD_GLOBAL)
del _gyoto_tmp_getdlopenflags, _gyoto_tmp_setdlopenflags, _gyoto_tmp_RTLD_GLOBAL
