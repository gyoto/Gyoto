import gi
gi.require_version("Gtk", "4.0")
from gi.repository import GLib

GLib.set_prgname("Gyoto")

from .widgets import *
from .utils import *
from .apps.gyotoy import gyotoy

from .widgets import __all__ as widgets__all__
from .utils import __all__ as utils__all__

__all__ = widgets__all__ + utils__all__ + ['gyotoy']

from .apps.gyoto_object_editor import GyotoObjectEditor

### The following should be achieved using Swig's extend mechanism
def edit(self, blocking=True):
    """Display a GUI to edit this object's Properties"""
    ObjectEditor.run(self, blocking=blocking)

import gyoto
gyoto.core.Object.edit=edit
