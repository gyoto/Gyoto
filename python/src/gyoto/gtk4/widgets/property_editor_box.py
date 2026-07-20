"""PropertyEditorBox: Widget for Editing Gyoto Object Properties

This module provides a GTK4 widget for editing the properties of Gyoto
objects.

It dynamically creates appropriate editor widgets based on each
property's type.

Widget Layout
------------
    ┌───────────────────────────────────────┐
    │ (PropertyEditorBox:)                  │
    │                                       │
    │ ┌─ property1 ───────────────────────┐ │
    │ │ widget for property1              │ │
    │ └───────────────────────────────────┘ │
    │ ┌─ property2 ───────────────────────┐ │
    │ │ widget for property2              │ │
    │ └───────────────────────────────────┘ │
    │ ┌─ property3 ───────────────────────┐ │
    │ │ widget for property3              │ │
    │ └───────────────────────────────────┘ │
    │ ...                                   │
    └───────────────────────────────────────┘

Description
-----------
This widget:
- Inspects a Gyoto object's properties
- Creates appropriate GTK widgets for each property type
- Handles value changes and forwards them to the object
- Supports nested objects with recursive editing

Supported Property Types:
- double_t: ScientificSpin
- long_t/unsigned_long_t/size_t_t: Gtk.SpinButton
- bool_t: Gtk.CheckButton
- string_t: Gtk.Entry
- filename_t: FilenameEditor
- vector_double_t/vector_unsigned_long_t: VectorScientificSpin
- metric_t/screen_t/spectrum_t/astrobj_t/spectrometer_t: GyotoObjectChooser

Example:
    editor = PropertyEditorBox(obj=my_gyoto_object)
    editor.connect("value-changed",
                   lambda w, name: print(f"{name} changed"))

"""

__all__ = ['PropertyEditorBox']

import gi
gi.require_version("Gtk", "4.0")
from gi.repository import Gtk, GLib, GObject

from functools import wraps
import ctypes
import numpy
from math import isnan

from .filename_editor import FilenameEditor
from .scientific_spin import ScientificSpin
from .vector_scientific_spin import VectorScientificSpin
from .gyoto_object_chooser import GyotoObjectChooser

from ..utils import show_error_dialog

from ...core import Property, Object
from ... import core, metric, astrobj, spectrum, spectrometer

class PropertyEditorBox(Gtk.Box):
    """A compound widget for editing Gyoto object properties.

    This widget dynamically creates and manages editor widgets for
    each property of a Gyoto object, based on the property's type. It
    provides a complete UI for inspecting and modifying object
    properties.

    The widget handles:
    - Automatic widget creation based on property types
    - Value synchronization between widgets and the object
    - Error handling with user-friendly dialogs
    - Special handling for vector properties and nested objects

    Attributes:
        obj: The Gyoto object being edited
        hide: List of property names to hide from the editor
        first: List of property names to show first (before others)
        widgets: Dictionary mapping property names to their editor widgets

    Signals:
        value-changed: Emitted when a property value changes (includes
            property name)
        child-changed: Emitted when a child object is replaced
            (includes child name)
        child-mutated: Emitted when a child object's property changes
            (includes child name)
        recursive-value-changed: Emitted when a nested property changes
            (includes the full property path).

    """

    __gsignals__ = {
        "value-changed": (GObject.SignalFlags.RUN_FIRST, None, (GObject.TYPE_STRING,)),
        "child-changed": (GObject.SignalFlags.RUN_FIRST, None, (GObject.TYPE_STRING,)),
        "child-mutated": (GObject.SignalFlags.RUN_FIRST, None, (GObject.TYPE_STRING,)),
        "recursive-value-changed": (GObject.SignalFlags.RUN_FIRST, None, (GObject.TYPE_STRING,)),
    }

    @staticmethod
    def gtk_callback(method):
        """Decorator to catch Gyoto errors and show them in dialog windows.

        This decorator wraps methods that interact with Gyoto objects,
        catching any Gyoto errors and displaying them in a
        user-friendly dialog instead of printing to the terminal.

        Usage:
            @gtk_callback
            def my_method(self, widget, name, *args):
                # Method implementation that may raise Gyoto errors
                ...

        Args:
            method: The method to wrap

        Returns:
            The wrapped method with error handling

        """
        @wraps(method)
        def wrapper(self, widget=None, name='parameter', *args):
            try:
                return method(self, widget, name, *args)
            except core.Error as e:
                show_error_dialog(detail=str(e),
                                  message=f'Error setting {name}',
                                  widget=widget)
        return wrapper

    def __init__(self, obj, hide=[], first=[], connector=None,
                 *args, **kwargs):
        """Initialize the PropertyEditorBox widget.

        Args:
            obj: The Gyoto object to edit
            hide: List of property names to exclude from the editor
                (default: [])
            first: List of property names to show first, before others
                (default: [])
            connector (multiprocessing.Connection or None): If the
                editor runs in a separate process, this is used to
                send updates back to the caller.
            *args: Additional positional arguments for Gtk.Box
            **kwargs: Additional keyword arguments for Gtk.Box

        """
        if "orientation" not in kwargs:
            kwargs['orientation'] = Gtk.Orientation.VERTICAL
        if "spacing" not in kwargs:
            kwargs['spacing'] = 10
        super().__init__(*args, **kwargs)
        self.obj = obj
        self.hide = hide
        self.first = first
        self.populate_properties()
        if connector is not None:
            self.connect('recursive-value-changed',
                         self.on_recursive_value_changed_pipe_sender,
                         connector)

    def populate_properties(self):
        """Generate editor widgets for all object properties.

        This method:
        1. Gets all property names from the object
        2. Creates a frame for each property with an appropriate
           editor widget
        3. Handles special cases (like InitCoord with
           3-velocity/4-velocity)
        4. Connects signals for value changes

        """
        parameters = self.obj.getPropertyNames()
        first_and_parameters = self.first + list(parameters)

        self.widgets = dict()

        for name in first_and_parameters:
            if name not in parameters:
                continue
            if name in self.hide:
                continue
            if name in self.widgets:
                continue

            prop = self.obj.property(name)
            value = self.obj.get(prop)
            param_type = prop.type

            # Create frame for this property
            frame = Gtk.Frame()
            frame.set_label(name)
            frame.set_label_align(0.0)  # Left-align label
            frame.set_tooltip_text(prop.doc)

            hbox = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=5)
            frame.set_child(hbox)

            # Create appropriate widget based on property type
            if param_type == core.Property.double_t:
                spin = ScientificSpin(value=value,
                                      with_unit=prop.supportsUnits())
                spin.connect("value-changed", self.on_parameter_changed, name)
                spin.connect("unit-changed", self.on_unit_changed, name)
                hbox.append(spin)
                self.widgets[name] = spin

            elif param_type in (core.Property.long_t,
                                core.Property.unsigned_long_t,
                                core.Property.size_t_t):
                # Set up appropriate ranges for each integer type
                if param_type == core.Property.long_t:
                    lower = -(1 << (ctypes.sizeof(ctypes.c_long) * 8 - 1))
                    upper = (1 << (ctypes.sizeof(ctypes.c_long) * 8 - 1)) - 1
                elif param_type == core.Property.unsigned_long_t:
                    lower = 0
                    upper = (1 << (ctypes.sizeof(ctypes.c_long) * 8)) - 1
                elif param_type == core.Property.size_t_t:
                    lower = 0
                    upper = (1 << (ctypes.sizeof(ctypes.c_size_t) * 8)) - 1
                else:
                    raise RuntimeError("Bug: unhandled integer property type")

                adjustment = Gtk.Adjustment(
                    value=value,
                    lower=lower,
                    upper=upper,
                    step_increment=1,
                    page_increment=10,
                    page_size=0
                )
                spin = Gtk.SpinButton()
                spin.set_adjustment(adjustment)
                spin.set_numeric(True)
                spin.set_digits(0)
                spin.set_hexpand(True)
                hbox.append(spin)
                spin.connect("value-changed", self.on_parameter_changed, name)
                # Block scroll events to avoid accidental changes
                scroll_controller = Gtk.EventControllerScroll()
                scroll_controller.connect("scroll", lambda *args: True)
                spin.add_controller(scroll_controller)
                self.widgets[name] = spin

            elif param_type == core.Property.bool_t:
                # Create radio buttons for boolean properties
                radio_true = Gtk.CheckButton(label=prop.name)
                radio_false = Gtk.CheckButton(label=prop.name_false)
                radio_false.set_group(radio_true)
                hbox.append(radio_true)
                hbox.append(radio_false)
                radio_true.set_active(value)
                radio_false.set_active(not value)
                radio_true.connect("toggled", self.on_parameter_changed, name)
                self.widgets[name] = radio_true

            elif param_type == core.Property.string_t:
                entry = Gtk.Entry()
                entry.set_hexpand(True)
                entry.set_text(str(value))
                hbox.append(entry)
                entry.connect("activate", self.on_parameter_changed, name)
                self.widgets[name] = entry

            elif param_type == core.Property.filename_t:
                editor = FilenameEditor(str(value))
                editor.connect(
                    "value-changed",
                    self.on_parameter_changed,
                    name
                )
                hbox.append(editor)
                self.widgets[name] = editor

            elif param_type in (core.Property.vector_double_t,
                                core.Property.vector_unsigned_long_t):
                # Use VectorScientificSpin for vector properties
                itemclass = (ScientificSpin
                             if param_type == core.Property.vector_double_t
                             else Gtk.SpinButton)
                vector = VectorScientificSpin(value=value,
                                              with_unit=prop.supportsUnits(),
                                              itemclass=itemclass)
                vector.connect("value-changed", self.on_parameter_changed, name)
                vector.connect("unit-changed", self.on_unit_changed, name)
                hbox.append(vector)
                self.widgets[name] = vector

            elif param_type == core.Property.metric_t:
                # Use GyotoObjectChooser for metric properties
                chooser = GyotoObjectChooser(metric,
                                             obj=getattr(self.obj, name))
                hbox.append(chooser)
                chooser.connect("object-changed", self.on_object_changed, name)
                chooser.connect("object-mutated", self.on_object_mutated, name)
                chooser.connect("recursive-value-changed", self.on_recursive_value_changed, name)
                self.widgets[name] = chooser

            elif param_type == core.Property.screen_t:
                chooser = GyotoObjectChooser(core.Screen,
                                             obj=getattr(self.obj, name))
                hbox.append(chooser)
                chooser.connect("object-changed", self.on_object_changed, name)
                chooser.connect("object-mutated", self.on_object_mutated, name)
                chooser.connect("recursive-value-changed", self.on_recursive_value_changed, name)
                self.widgets[name] = chooser

            elif param_type == core.Property.spectrum_t:
                chooser = GyotoObjectChooser(spectrum,
                                             obj=getattr(self.obj, name))
                hbox.append(chooser)
                chooser.connect("object-changed", self.on_object_changed, name)
                chooser.connect("object-mutated", self.on_object_mutated, name)
                chooser.connect("recursive-value-changed", self.on_recursive_value_changed, name)
                self.widgets[name] = chooser

            elif param_type == core.Property.astrobj_t:
                chooser = GyotoObjectChooser(astrobj,
                                             obj=getattr(self.obj, name))
                hbox.append(chooser)
                chooser.connect("object-changed", self.on_object_changed, name)
                chooser.connect("object-mutated", self.on_object_mutated, name)
                chooser.connect("recursive-value-changed", self.on_recursive_value_changed, name)
                self.widgets[name] = chooser

            elif param_type == core.Property.spectrometer_t:
                chooser = GyotoObjectChooser(spectrometer,
                                             obj=getattr(self.obj, name))
                hbox.append(chooser)
                chooser.connect("object-changed", self.on_object_changed, name)
                chooser.connect("object-mutated", self.on_object_mutated, name)
                chooser.connect("recursive-value-changed", self.on_recursive_value_changed, name)
                self.widgets[name] = chooser

            self.append(frame)

            # Special case: InitCoord property
            if name == 'InitCoord':
                try:
                    st = astrobj.Star(self.obj)
                    isstar = True
                except Exception:
                    isstar = False
                try:
                    ph = core.Photon(self.obj)
                    isphoton = True
                except Exception:
                    isphoton = False

                if isstar or isphoton:
                    # Add 3-velocity/4-velocity radio buttons for InitCoord
                    hbox2 = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=5)
                    radio_3vel = Gtk.CheckButton(label='3-velocity')
                    radio_4vel = Gtk.CheckButton(label='4-velocity')
                    radio_4vel.set_group(radio_3vel)
                    hbox2.append(radio_3vel)
                    hbox2.append(radio_4vel)
                    radio_3vel.set_active(False)
                    radio_4vel.set_active(True)
                    radio_3vel.connect("toggled", self.on_3vel_toggled, name)
                    self.widgets[f'{name}:veltype'] = radio_3vel
                    self.widgets['InitCoord'].prepend(hbox2)

    @gtk_callback
    def on_unit_changed(self, widget, name, *args):
        """Handle unit changes from ScientificSpin or VectorScientificSpin.

        Updates the object property with the new unit.

        Args:
            widget: The spin widget that emitted the signal
            name: The name of the property being edited
            *args: Additional arguments
        """
        unit = widget.get_unit()
        value = self.obj.get(name, unit)
        widget.set_value(value)

    @gtk_callback
    def on_object_changed(self, widget, name, *args):
        """Handle object replacement from GyotoObjectChooser.

        Called when a sub-object (e.g., Metric, Spectrum) is replaced with
        a different kind.

        Args:
            widget: The GyotoObjectChooser that emitted the signal
            name: The name of the property containing the sub-object
            *args: Additional arguments
        """
        self.obj.set(name, widget.obj)
        self.emit('child-changed', name)
        self.emit('recursive-value-changed', name)

    @gtk_callback
    def on_object_mutated(self, widget, name, *args):
        """Handle property changes in sub-objects from GyotoObjectChooser.

        Called when a property of a sub-object (e.g., a Metric
        property) changes.

        Forwards the signal to the parent.

        Args:
            widget: The GyotoObjectChooser that emitted the signal
            name: The name of the property containing the sub-object
            *args: Additional arguments

        """
        self.emit('child-mutated', name)

    @gtk_callback
    def on_recursive_value_changed(self, widget, pname, name, *args):
        """Handle property changes in sub-objects from GyotoObjectChooser.

        Called when a property of a sub-object (e.g., a Metric
        property) changes. Forwards the signal to the parent with the
        full property path.

        Args:
            widget: The GyotoObjectChooser that emitted the signal.
            pname: The name of the property that changed in the sub-object.
            name: The name of the property containing the sub-object.
            *args: Additional arguments.

        """
        self.emit('recursive-value-changed', f'{name}.{pname}')

    @gtk_callback
    def on_recursive_value_changed_pipe_sender(self, wdgt, ppath,
                                               connector):
        """Send recursive value changes to another process.

        When the UI runs in a separate process, this callback sends
        property updates back to the caller via the connector.

        Args:
            wdgt: The widget that emitted the signal.
            ppath: The full path to the property that changed (e.g.,
                'Metric.Spin').
            connector: The multiprocessing.Connection to send updates to
                the caller.

        """
        descendents = ppath.split('.')
        value = self.obj
        for i in range(len(descendents)):
            obj = value
            value = obj.get(descendents[i])
        if isinstance(value, Object):
            value = str(value)
        connector.send(['update', ppath, value])
        print(f'event sent: {ppath} == {value}')

    @gtk_callback
    def on_parameter_changed(self, widget, name, *args):
        """Handle value changes from most editor widgets.

        This is the primary callback for most property changes. It extracts
        the new value from the widget and updates the object property.

        Args:
            widget: The editor widget that changed
            name: The name of the property being edited
            *args: Additional arguments
        """
        new_unit = None
        if isinstance(widget, ScientificSpin):
            new_value = widget.get_value()
            new_unit = widget.get_unit()

        elif isinstance(widget, Gtk.CheckButton):
            new_value = widget.get_active()

        elif isinstance(widget, Gtk.SpinButton):
            new_value = widget.get_value_as_int()

        elif isinstance(widget, Gtk.Entry):
            new_value = widget.get_text()

        elif isinstance(widget, FilenameEditor):
            new_value = widget.get_value()

        elif isinstance(widget, VectorScientificSpin):
            new_value = widget.get_value()
            new_unit = widget.get_unit()
            # Special case: InitCoord with 7 elements (4 pos + 3 vel)
            if len(new_value) == 7 and f'{name}:veltype' in self.widgets:
                if self.obj.Metric is not None:
                    pos = new_value[0:4]
                    vel = new_value[4:7]
                    if self.obj.kind() == 'Photon':
                        # For photons, ensure null geodesic
                        new_value = numpy.array(pos + [1.] + vel)
                        self.obj.Metric.nullifyCoord(new_value)
                    else:
                        # For stars, compute tdot and normalize
                        tdot = self.obj.Metric.SysPrimeToTdot(pos, vel)
                        if not isnan(tdot):
                            new_value = pos + [tdot] + [x*tdot for x in vel]

        # Apply the new value to the object
        if new_unit is None:
            self.obj.set(name, new_value)
        else:
            self.obj.set(name, new_value, new_unit)

        self.emit('value-changed', name)
        self.emit('recursive-value-changed', name)

    @gtk_callback
    def on_file_chooser_clicked(self, button, entry):
        """Handle file chooser button click for filename properties.

        Opens a file dialog for selecting a file.

        Args:
            button: The button that was clicked
            entry: The entry widget to update with the selected path
        """
        dialog = Gtk.FileDialog()
        dialog.open(
            self.get_root(),
            None,
            lambda dialog, result: self.on_file_selected(dialog, result, entry)
        )

    @gtk_callback
    def on_file_selected(self, dialog, result, entry):
        """Handle file selection from the dialog.

        Updates the entry widget with the selected file path.

        Args:
            dialog: The Gtk.FileDialog that completed
            result: The result of the dialog operation
            entry: The entry widget to update
        """
        try:
            file = dialog.open_finish(result)
        except GLib.Error:
            # User cancelled the dialog
            return

        if file is not None:
            entry.set_text(file.get_path())

    @gtk_callback
    def on_3vel_toggled(self, widget=None, name='InitCoord', *args):
        """Handle 3-velocity/4-velocity radio button toggles.

        This callback handles the special case of InitCoord, which can
        be represented as either 3-velocity or 4-velocity. When the
        user toggles between these representations, the values are
        converted appropriately.

        Args:
            widget: The radio button that was toggled (default: None)
            name: The property name (default: 'InitCoord')
            *args: Additional arguments

        """
        if widget is None:
            widget = self.widgets[f'{name}:veltype']
        coord = self.obj.get(name)
        if widget.get_active():
            # Switch to representing as 3-velocity
            tdot = coord[4]
            if tdot == 0.:
                tdot = 1.
            # Convert 4-velocity to 3-velocity: vel_3 = vel_4 / tdot
            pos3vel = coord[0:4] + tuple((x/tdot for x in coord[5:8]))
            self.widgets[name].set_value(pos3vel)
        else:
            # Switch back to showing 4-velocity
            self.widgets[name].set_value(coord)
