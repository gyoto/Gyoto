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

from . import gyotoy
from . import gyoto_object_editor
