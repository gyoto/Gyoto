"""GTK4 Custom Widgets for Gyoto

This module provides a collection of GTK4 widgets specifically
designed for building Gyoto-based applications. Each widget is
tailored to handle common patterns in scientific computing and general
relativity simulations.

Available Widgets
-----------------
The following widgets are available for import:

    FilenameEditor:
        A compound widget for selecting files, combining a text entry
        with a file browser button. Useful for loading/saving XML
        configuration files.

    GyotoObjectChooser:
        A widget for choosing and editing Gyoto object kinds (metrics,
        astrobjs, spectra, etc.) with a dropdown selector and property
        editor. Supports loading objects from XML files.

    PropertyEditorBox:
        A panel for editing the properties of Gyoto
        objects. Dynamically creates appropriate editor widgets based
        on each property's type (doubles, vectors, booleans, strings,
        nested objects, etc.).

    ScientificSpin:
        A spin button with scientific notation support and optional
        unit display.  Ideal for entering physical quantities with
        large/small magnitudes.

    SimulationControls:
        A compound widget for controlling and monitoring simulations,
        including progress bar, status messages, and transport
        controls (play/pause/stop).

    VectorScientificSpin:
        A widget for editing vectors of floating-point values with
        optional unit display and dynamic resizing. Supports both
        ScientificSpin and Gtk.SpinButton as item classes.

    Viewer3D:
        A GTK4 widget embedding a Matplotlib 3D view with navigation
        toolbar.  Provides a self-contained 3D visualization component
        for trajectories and other 3D data.

Usage Examples
-------------
    # Import individual widgets
    from gyoto.gtk4.widgets import Viewer3D, ScientificSpin, SimulationControls

    # Or import all at once
    from gyoto.gtk4.widgets import *

    # Create a viewer with controls
    viewer = Viewer3D()
    controls = SimulationControls()
    controls.connect("play-pause", lambda w: viewer.draw())

See Also
--------
gyoto.gtk4.apps : High-level application windows
gyoto.gtk4.utils : GTK4 utility functions

"""

from .filename_editor import FilenameEditor
from .gyoto_object_chooser import GyotoObjectChooser
from .property_editor_box import PropertyEditorBox
from .scientific_spin import ScientificSpin
from .simulation_controls import SimulationControls
from .vector_scientific_spin import VectorScientificSpin
from .viewer_3d import Viewer3D

__all__ = [
    'FilenameEditor', 'GyotoObjectChooser', 'PropertyEditorBox',
    'ScientificSpin', 'SimulationControls', 'VectorScientificSpin', 'Viewer3D'
]
