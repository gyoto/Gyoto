"""Gyoto GTK4: GTK4-based GUI Components for Gyoto

This module provides GTK4-based graphical user interface components
for the Gyoto library, enabling interactive visualization and editing
of general relativistic simulations.

Module Structure
---------------
- **widgets**: Custom GTK4 widgets for Gyoto-specific data types
  (ScientificSpin, VectorScientificSpin, Viewer3D, etc.)
- **apps**: High-level application windows
  (GyotoyApplicationWindow, GyotoObjectEditor)
- **utils**: Utility functions for GTK4 applications
  (show_error_dialog)

Main Entry Points
----------------
- gyotoy(): Launch the main Gyotoy application window for geodesic
  simulation

- gyoto.core.Object.edit(): Window for editing Gyoto object properties,
  installed as a gyoto.core.Object method

Note:
    Requires GTK4 and Matplotlib with GTK4 backend.

"""

import gi
gi.require_version("Gtk", "4.0")
from gi.repository import GLib

# Set application name for GTK
GLib.set_prgname("Gyoto")

# Import submodules
from . import widgets, utils

# Re-export main entry point for convenience
from .apps.gyotoy import gyotoy

# Public API
__all__ = ['widgets', 'apps', 'utils', 'gyotoy']

# Import for monkey-patching
from .apps.gyoto_object_editor import GyotoObjectEditor

# Add edit() method to gyoto.core.Object for convenient property editing
# Note: This should ideally be achieved using SWIG's extend mechanism
def edit(self, blocking=True):
    GyotoObjectEditor.run(self, blocking=blocking)
edit.__doc__ = GyotoObjectEditor.run.__doc__

# Monkey-patch the edit method onto gyoto.core.Object
import gyoto
gyoto.core.Object.edit = edit
