'''Backward compatibility alias for gyoto.utils'''
from warnings import warn
warn('gyoto.util was renamed gyoto.utils, update your code')
from .utils import *
