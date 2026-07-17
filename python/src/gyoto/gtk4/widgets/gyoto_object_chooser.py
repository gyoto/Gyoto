"""GyotoObjectChooser: Widget for Selecting and Editing Gyoto Objects

This module provides a GTK4 widget for choosing and editing Gyoto
object kinds (metrics, astrobjs, spectra, spectrometers, etc.) with a
dropdown selector and property editor.

Widget Layout
------------
    ┌──────────────────────────────────────────────┐
    │  ┌────────────────────────────────────────┐  │
    │  │(Gtk.DropDown:)                         │  │
    │  │  -                                   ▼ │  │
    │  │  plugin1/kind1                         │  │
    │  │  plugin1/kind2                         │  │
    │  │  ...                                   │  │
    │  │  plugin2/kind3                         │  │
    │  │  plugin2/kind4                         │  │
    │  │  ...                                   │  │
    │  │  Open...                               │  │
    │  └────────────────────────────────────────┘  │
    │  ┌────────────────────────────────────────┐  │
    │  │(PropertyEditorBox:)                    │  │
    │  │  widget for property1                  │  │
    │  │  widget for property2                  │  │
    │  │  widget for property3...               │  │
    │  └────────────────────────────────────────┘  │
    └──────────────────────────────────────────────┘

Description
-----------
This widget allows users to:
- Select a Gyoto object kind from a dropdown list
- Edit the object's properties in a PropertyEditorBox
- Load objects from XML files

The widget handles circular dependencies with PropertyEditorBox through
lazy imports.

Example:
    chooser = GyotoObjectChooser(namespace=gyoto.metric, obj=my_metric)
    chooser.connect("object-changed", lambda w: print(w.obj))

"""

__all__ = ['GyotoObjectChooser']

import gi
gi.require_version("Gtk", "4.0")
from gi.repository import Gtk, Gio, GLib

import traceback

from ..utils import show_error_dialog
from ... import core, metric, astrobj, spectrum, spectrometer

## Use lazy import to break circular dependency with property_editor_box
# from .property_editor_box import PropertyEditorBox

class GyotoObjectChooser(Gtk.Box):
    """GTK4 widget for choosing and editing Gyoto object kinds.

    This widget combines a dropdown selector with a PropertyEditorBox
    to allow users to select a Gyoto object kind and edit its
    properties. When a new kind is selected, the current object is
    replaced with a new instance of the selected kind.

    The widget supports:
    - All Gyoto namespaces: metric, astrobj, spectrum, spectrometer,
      core.Screen
    - Loading objects from XML files
    - Signal emission when objects change or are mutated

    Attributes:
        obj: The current Gyoto object being edited (None if no selection)
        namespace: The Gyoto namespace this chooser operates on
        items: List of available plugin/kind options
        dropdown: Gtk.DropDown widget for kind selection
        frame: Gtk.Frame containing the PropertyEditorBox

    Signals:
        object-changed: Emitted when a new kind is selected or object
            is loaded
        object-mutated: Emitted when a property of the current object
            changes

    Note:
        Uses lazy imports to break circular dependency with PropertyEditorBox.

    """

    _updating = False

    __gsignals__ = {
        "object-changed": (gi.repository.GObject.SignalFlags.RUN_FIRST, None, ()),
        "object-mutated": (gi.repository.GObject.SignalFlags.RUN_FIRST, None, ())
    }

    def __init__(self, namespace, obj=None):
        """Initialize the GyotoObjectChooser widget.

        Args:
            namespace: The Gyoto namespace or class to choose from.
                Valid values: gyoto.metric, gyoto.astrobj,
                gyoto.spectrum, gyoto.spectrometer, or
                gyoto.core.Screen.
            obj: Initial object to edit (default: None)

        """
        super().__init__(orientation=Gtk.Orientation.VERTICAL,
                         spacing=0)

        # Lazy import to break dependency cycle with property_editor_box
        from .property_editor_box import PropertyEditorBox

        self.obj = obj
        self.namespace = namespace

        # Build list of available kinds for the dropdown
        if self.namespace == core.Screen:
            self.items = ["-", "built-in/Screen"]
        else:
            self.items = ["-"] + [x for x in
                                  namespace.Generic.registeredPluginsSlashKinds()]
        self.items.append("Open...")

        # Create dropdown with all available kinds
        self.dropdown = Gtk.DropDown.new_from_strings(self.items)
        self.append(self.dropdown)

        # Set initial selection based on obj
        if self.obj is None:
            self.dropdown.set_selected(0)  # Select "-" (no selection)
        else:
            if self.namespace == core.Screen:
                self.dropdown.set_selected(1)  # Select "built-in/Screen"
            else:
                # Find the index matching the object's kind
                indexes = [i for i, val in enumerate(self.items)
                           if val.endswith('/'+self.obj.kind())]
                if len(indexes) == 1:
                    self.dropdown.set_selected(indexes[0])
                else:
                    # Fallback to "-" if kind not found
                    self.dropdown.set_selected(0)
                    # show_error_dialog('could not determine
                    # plugin/kind for object')

        # Connect dropdown selection changes
        self.dropdown.connect("notify::selected", self.on_dropdown_activated)

        # Frame to contain the PropertyEditorBox
        self.frame = Gtk.Frame()
        self.append(self.frame)

        # Initialize with current object if provided
        if self.obj is not None:
            # Lazy import to break dependency cycle
            from .property_editor_box import PropertyEditorBox
            box = PropertyEditorBox(self.obj, hide=['Metric'])
            self.frame.set_child(box)
            box.connect("value-changed", self.on_child_value_changed)

    def on_child_value_changed(self, widget, *args):
        """Handle value changes from the PropertyEditorBox.

        Emits 'object-mutated' signal when a property changes.

        Args:
            widget: The PropertyEditorBox that emitted the signal
            *args: Additional arguments
        """
        self.emit("object-mutated")

    def on_dropdown_activated(self, widget, *args):
        """Handle dropdown selection changes.

        Creates a new object based on the selected kind and updates the UI.

        Args:
            widget: The Gtk.DropDown that changed
            *args: Additional arguments
        """
        if self._updating:
            return

        x = widget.get_selected_item().get_string()

        if x == '-':
            # No selection
            self.obj = None
            self.frame.set_child(None)
            self.emit("object-changed")

        elif x == 'Open...':
            # Open file dialog to load from XML
            dialog = Gtk.FileDialog()

            # Set up file filters
            xml_filter = Gtk.FileFilter()
            xml_filter.set_name("XML files")
            xml_filter.add_suffix('xml')
            xml_filter.add_pattern('*.xml')

            all_filter = Gtk.FileFilter()
            all_filter.set_name("All files")
            all_filter.add_pattern('*')

            filter_list = Gio.ListStore.new(Gtk.FileFilter)
            filter_list.append(xml_filter)
            filter_list.append(all_filter)

            dialog.set_filters(filter_list)
            dialog.set_property('default-filter', xml_filter)

            dialog.open(
                self.get_root(),
                None,
                lambda dialog, result: self.on_open_file_selected(dialog, result)
            )

        else:
            # Create new object based on selected plugin/kind
            if self.namespace == core.Screen:
                self.obj = core.Screen()
            else:
                plg, knd = x.split('/')
                self.obj = self.namespace.Generic(knd, (plg,))

            # Lazy import to break dependency cycle
            from .property_editor_box import PropertyEditorBox

            # Create new PropertyEditorBox for the object
            box = PropertyEditorBox(self.obj, hide=['Metric'])
            self.frame.set_child(box)
            box.connect("value-changed", self.on_child_value_changed)
            self.emit("object-changed")

    def on_open_file_selected(self, dialog, result):
        """Handle file selection from the dialog.

        Loads the selected file and creates a new object from it.

        Args:
            dialog: The Gtk.FileDialog that completed
            result: The result of the dialog operation
        """
        try:
            file = dialog.open_finish(result)
        except GLib.Error:
            return  # Dialog was cancelled or error occurred

        factory = None
        if file is not None:
            try:
                factory = core.Factory(file.get_path())
            except core.Error as e:
                show_error_dialog(
                    message=f"Error loading XML file {file.get_path()}:",
                    detail=e.get_message(),
                    window=self.get_root()
                )

        if factory is not None:
            try:
                if self.namespace == core.Screen:
                    self.obj = factory.screen()
                elif self.namespace == astrobj:
                    self.obj = factory.astrobj()
                elif self.namespace == metric:
                    self.obj = factory.metric()
                elif self.namespace == spectrum:
                    self.obj = factory.spectrum()
                elif self.namespace == spectrometer:
                    self.obj = factory.spectrometer()
                    if self.obj is None:
                        self.dropdown.set_selected(0)
                        raise ValueError('Please select a file containing a spectrometer')
                else:
                    raise ValueError(f'Namespace unimplemented: {self.namespace}')
            except Exception:
                show_error_dialog(
                    message="Could not construct object from XML file:",
                    detail=traceback.format_exc(),
                    window=self.get_root()
                )
                return

            print(self.namespace)
            print(self.obj)

            # Update dropdown to match the loaded object's kind
            self._updating = True
            if self.namespace == core.Screen:
                self.dropdown.set_selected(1)
            else:
                indexes = [i for i, val in enumerate(self.items)
                           if val.endswith('/'+self.obj.kind())]
                if len(indexes) == 1:
                    self.dropdown.set_selected(indexes[0])
                # else: Could not find matching kind, leave as is
            self._updating = False

            # Lazy import to break dependency cycle
            from .property_editor_box import PropertyEditorBox

            # Create PropertyEditorBox for the loaded object
            box = PropertyEditorBox(self.obj, hide=['Metric'])
            self.frame.set_child(box)
            box.connect("value-changed", self.on_child_value_changed)
            self.emit("object-changed")
