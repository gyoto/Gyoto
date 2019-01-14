"""The Gyoto standard plug-in

Note that the gyoto modules have been renamed:
  - gyoto -> gyoto.core
  - gyoto_std -> gyoto.std
  - gyoto_lorene -> gyoto.lorene

Please update your code accordingly.

"""

from gyoto.std import *

import warnings

warnings.warn ("The gyoto modules have been renamed, please see help(gyoto_std)", DeprecationWarning)
