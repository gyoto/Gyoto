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

from .utils import gui_launcher

# Set application name for GTK
GLib.set_prgname("Gyoto")

# Re-export main entry point for convenience

# Interactive-session entry point:
def gyotoy(particle=None):
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
    The application runs in a separate process in a non-blocking fashion.
    An optional particle (gyoto.std.Star or gyoto.core.Photon) can be
    provided. Changes made in the GUI affect the variable being passed.

    """
    # lazy import to not get in the way of stand-alone execution
    from .utils import recursive_value_changed_pipe_receiver

    def gtk_process(connector, obj):
        from .apps.gyotoy import GyotoyApplicationWindow
        GyotoyApplicationWindow.run(particle,
                                    blocking=True,
                                    connector=connector)

    gui_launcher(gtk_process,
                 None if particle is None else recursive_value_changed_pipe_receiver,
                 particle)

# Add edit() method to gyoto.core.Object for convenient property editing
# Note: This should ideally be achieved using SWIG's extend mechanism
def edit(self):
    """A GTK4 window for editing Gyoto object properties.

    This window provides a scrollable view of all editable properties of a
    Gyoto object, using a PropertyEditorBox as its main content.

    The application runs in a separate process in a non-blocking fashion.

    Parameters:
        obj: The Gyoto object to edit

    Note:
        The GUI runs in a separate process.

    """

    from .utils import recursive_value_changed_pipe_receiver

    def gtk_process(connector, obj):
        from .apps.gyoto_object_editor import GyotoObjectEditor
        win = GyotoObjectEditor(str(obj), blocking=True, connector=connector)
        win.present()
        win.main_loop.run()

    gui_launcher(gtk_process,
                 recursive_value_changed_pipe_receiver,
                 self)

# Monkey-patch the edit method onto gyoto.core.Object
from ..core import Object
Object.edit = edit

# Autoload submodules on demand
def __getattr__(name):
    """Autoload submodules.

    To avoid issues with loading Gtk when not needed, they should
    not be loaded by default. This will load them on demand.

    """

    # Avoid infinite recursion for standard attributes
    if name not in __all__:
        raise AttributeError(f"module '{__name__}' has no attribute '{name}'")

    fullname = f"{__name__}.{name}"

    mod = importlib.import_module(fullname)
    setattr(sys.modules[__name__], name, mod)
    return mod
