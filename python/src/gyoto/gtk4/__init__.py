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

# Public API
__all__ = ['widgets', 'apps', 'utils', 'gyotoy']

import gi
gi.require_version("Gtk", "4.0")
from gi.repository import GLib

# Set application name for GTK
GLib.set_prgname("Gyoto")

# Import submodules
from . import widgets, utils

# Re-export main entry point for convenience

# Interactive-session entry point:
def gyotoy(particle=None, blocking=True):
    """Gyotoy: GTK4 Application for Gyoto Geodesic Integration

    This application provides a graphical interface for simulating and
    visualizing geodesics (time-like or null) in spacetimes supported
    by the Gyoto library.

    Layout
    ------
    ┌────────────────────────────────────────────────────────────────────┐
    │ MyApp                                                     ☰        │
    ├────────────────────────────────────────────────────────────────────┤
    │ ┌──────────────────────────────┬────────────────────────────────┐  │
    │ │                              │  ○ Star                        │  │
    │ │                              │  ○ Photon                      │  │
    │ │                              │                                │  │
    │ │      Matplotlib canvas       │  ┌──────────────────────────┐  │  │
    │ │                              │  │ PropertyEditorBox        │  │  │
    │ │                              │  └──────────────────────────┘  │  │
    │ └──────────────────────────────┴────────────────────────────────┘  │
    ├────────────────────────────────────────────────────────────────────┤
    │ ████████████████████████████████████────────────────────────────── │
    │                                                                    │
    │ Status...                    N: [100  ]      ⏮        ▶        ⏹   │
    │                                                                    │
    └────────────────────────────────────────────────────────────────────┘

    Description
    -----------
    - **Left Panel**: 3D Matplotlib viewer displaying the particle
        trajectory.
    - **Right Panel**: Property editor for adjusting metric and
        particle parameters.
    - **Bottom Controls**: Play/pause/stop buttons, interpolation
        settings, and status display.
    - **Worker Process**: Heavy computations run in a background
        process to keep the UI responsive.

    Usage
    -----
    Run as a standalone application:
        python3 -m gyoto.gtk4.apps.gyotoy

    Or import and use programmatically:
        from gyoto.gtk4 import gyotoy
        window = gyotoy([particle])
    From ipython3, the application can be involed non-blocking:
        %gui gtk4
        from gyoto.gtk4 import gyotoy
        window = gyotoy([particle, ] blocking=False)
    An optional particle (gyoto.std.Star or gyoto.core.Photon) can be
    provided.

    """
    # lazy import to not get in the way of stand-alone execution
    from .apps.gyotoy import GyotoyApplicationWindow
    return GyotoyApplicationWindow.run(particle, blocking)

# Import for monkey-patching
from .apps.gyoto_object_editor import GyotoObjectEditor

# Add edit() method to gyoto.core.Object for convenient property editing
# Note: This should ideally be achieved using SWIG's extend mechanism
def edit(self, blocking=True):
    GyotoObjectEditor.run(self, blocking=blocking)
edit.__doc__ = GyotoObjectEditor.run.__doc__

# Monkey-patch the edit method onto gyoto.core.Object
from ..core import Object
Object.edit = edit
