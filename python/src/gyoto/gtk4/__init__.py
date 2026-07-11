import gi
gi.require_version("Gtk", "4.0")
from gi.repository import GLib

GLib.set_prgname("Gyoto")

from .apps.object_editor import ObjectEditor

### The following should be achieved using Swig's extend mechanism
def edit(self, blocking=True):
    """Display a GUI to edit this object's Properties"""
    ObjectEditor.run(self, blocking=blocking)

import gyoto
gyoto.core.Object.edit=edit
