"""The Gyoto Lorene plug-in

Note that the gyoto modules have been renamed:
  - gyoto -> gyoto.core
  - gyoto_std -> gyoto.std
  - gyoto_lorene -> gyoto.lorene

Please update your code accordingly.

"""

from gyoto.lorene import *

import warnings

warnings.warn ("The gyoto modules have been renamed, please see help(gyoto_lorene)", DeprecationWarning)
