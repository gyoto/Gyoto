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

# Add edit() method to gyoto.core.Object for convenient property editing
# Note: This should ideally be achieved using SWIG's extend mechanism
def edit(self, blocking=True):
    """A GTK4 window for editing Gyoto object properties.

    This window provides a scrollable view of all editable properties of a
    Gyoto object, using a PropertyEditorBox as its main content.

    The window can run in blocking mode (manages its own GLib main loop) or
    non-blocking mode (for integration with external event loops like
    IPython).

    Parameters:
        obj: The Gyoto object to edit
        blocking (bool): If True, the window manages the GLib main
            loop.  If False, the caller must manage the event loop
            (e.g., using %gui gtk4 in IPython).

    Attributes:
        obj: The Gyoto object being edited
        main_loop: GLib.MainLoop instance (if blocking=True)
        scrolled_window: Gtk.ScrolledWindow containing the editor
        vbox: PropertyEditorBox for editing object properties

    """
    from .apps.gyoto_object_editor import GyotoObjectEditor
    GyotoObjectEditor.run(self, blocking=blocking)

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
