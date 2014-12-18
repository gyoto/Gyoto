"""
The General relativitY Orbit Tracer of paris Observatory
"""
# This is necessary at least on Linux to let the libgyoto symbols be
# visible from the libgyoto-stdplug (and other plug-ins) symbols.
import sys, ctypes
sys.setdlopenflags(sys.getdlopenflags() | ctypes.RTLD_GLOBAL)
