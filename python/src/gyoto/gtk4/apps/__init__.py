"""
Gyoto GTK4 Applications

This module provides high-level GTK4 application windows for Gyoto:

- gyotoy: Main application for simulating and visualizing geodesics
  - GyotoyApplication: GTK Application class
  - GyotoyApplicationWindow: Main window class
  - gyotoy(): Entry point for interactive use (e.g., from IPython)
  Can be called as a stand-alone application with:
    python3 -m gyoto.gtk4.apps.gyotoy

- gyoto_object_editor: Window for editing Gyoto object properties
  - GyotoObjectEditor: Window class for property editing
  Can be called as method gyoto.core.Object.edit.

"""

import importlib

__all__ = ['gyotoy', 'gyoto_object_editor']


def __getattr__(name):
    """Autoload submodules.

    To avoid issues when running submodules as scripts, they should
    not be loaded by default. This will load them on demand.

    """

    # Avoid infinite recursion for standard attributes
    if name not in __all__:
        raise AttributeError(f"module '{__name__}' has no attribute '{name}'")

    fullname = f"{__name__}.{name}"

    mod = importlib.import_module(fullname)
    setattr(sys.modules[__name__], name, mod)
    return mod
