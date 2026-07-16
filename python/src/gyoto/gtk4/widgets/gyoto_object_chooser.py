## Compound widget to choose/edit subobject
#
# A widget to choose a new kind for a subobject and display its own
# editor panel
#
# ┌──────────────────────────────────────────────┐
# │  ┌────────────────────────────────────────┐  │
# │  │(Gtk.DropDown:)                         │  │
# │  │  -                                   ▼ │  │
# │  │  plugin1/kind1                         │  │
# │  │  plugin1/kind2                         │  │
# │  │  ...                                   │  │
# │  │  plugin2/kind3                         │  │
# │  │  plugin2/kind4                         │  │
# │  │  ...                                   │  │
# │  └────────────────────────────────────────┘  │
# │  ┌────────────────────────────────────────┐  │
# │  │(PropertyEditorBox:)                    │  │
# │  │  widget for property1                  │  │
# │  │  widget for property2                  │  │
# │  │  widget for property3...               │  │
# │  └────────────────────────────────────────┘  │
# └──────────────────────────────────────────────┘
#

__all__ = ['GyotoObjectChooser']

import gi
gi.require_version("Gtk", "4.0")
from gi.repository import Gtk, Gio, GLib

import traceback

from ..utils import show_error_dialog
from ... import core, metric, astrobj, spectrum, spectrometer

## use lazy import due to cycle in dependencies
# from .property_editor_box import PropertyEditorBox

class GyotoObjectChooser(Gtk.Box):
    """Gtk widget for choosing and editing a Gyoto object kind

    The widget is composed of a drop-down menu from which a kind can
    be choosen, and a PropertyEditorBox for this kind. Selecting
    another kind destroys the current PropertyEditorBox and creates a
    new one.

    Parameters:

      namespace: the Gyoto namespace or class to choose from
          (gyoto.metric, gyoto.astrobj, gyoto.spectrum,
          gyoto.spectrometer or gyoto.core.Screen).

      obj (optional): original object.

    Public members:

       obj: initially, a copy of obj. Reinitialized each type a new
           kind is selected from the drop-down menu.

    Example:
      box=GyotoObjectChooser(gyoto.metric, obj=myobject.metric())

    Signals:
      This widget emits the two following signals:

      object-changed: whenever a new kind is selected, obj is
          reinitialized. The caller can connect to this signal and set
          a callback to fetch self.obj.

      object-mutated: everytime a property of obj is changed using the
          controls in the PropertyEditorBox, this signal is emitted.

    """

    _updating = False

    __gsignals__ = {
        "object-changed": (gi.repository.GObject.SignalFlags.RUN_FIRST, None, ()),
        "object-mutated": (gi.repository.GObject.SignalFlags.RUN_FIRST, None, ())
    }
    def __init__(self, namespace, obj=None):
        super().__init__(orientation=Gtk.Orientation.VERTICAL,
                         spacing=0)

        # lazy import to break dependency cycle
        from .property_editor_box import PropertyEditorBox

        self.obj=obj
        self.namespace=namespace
        if self.namespace == core.Screen:
            self.items = ["-", "built-in/Screen"]
        else:
            self.items = ["-"] + [x for x in
                                  namespace.Generic.registeredPluginsSlashKinds()]
        self.items.append("Open...")
        self.dropdown=Gtk.DropDown.new_from_strings(self.items)
        self.append(self.dropdown)

        if self.obj is None:
            self.dropdown.set_selected(0)
        else:
            if self.namespace == core.Screen:
                self.dropdown.set_selected(1)
            else:
                indexes = [i for i, val in enumerate(self.items)
                           if val.endswith('/'+self.obj.kind())]
                if len(indexes) == 1:
                    self.dropdown.set_selected(indexes[0])
                else:
                    # Let's show "-" instead of not being able to edit the object
                    self.dropdown.set_selected(0)
                    # show_error_dialog('could not determine plugin/kind for object')

        self.dropdown.connect("notify::selected", self.on_dropdown_activated)

        self.frame = Gtk.Frame()
        self.append(self.frame)
        if self.obj is not None:
            box = PropertyEditorBox(self.obj, hide=['Metric'])
            self.frame.set_child(box)
            box.connect("value-changed", self.on_child_value_changed)

    def on_child_value_changed(self, widget, *args):
        self.emit("object-mutated")
    
    def on_dropdown_activated(self, widget, *args):
        if self._updating: return

        x = widget.get_selected_item().get_string()
        if x == '-':
            self.obj=None
            self.frame.set_child(None)
            self.emit("object-changed")
        elif x == 'Open...':
            dialog = Gtk.FileDialog()

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
                lambda dialog, result:
                self.on_open_file_selected(dialog, result)
            )

        else:
            if self.namespace == core.Screen:
                self.obj = core.Screen()
            else:
                plg, knd = x.split('/')
                self.obj=self.namespace.Generic(knd, (plg,))

            # lazy import to break dependency cycle
            from .property_editor_box import PropertyEditorBox

            box = PropertyEditorBox(self.obj, hide=['Metric'])
            self.frame.set_child(box)
            box.connect("value-changed", self.on_child_value_changed)
            self.emit("object-changed")
    
    def on_open_file_selected(self, dialog, result):
        try:
            file = dialog.open_finish(result)
        except GLib.Error:
            return

        factory = None
        if file is not None:
            try:
                factory = core.Factory(file.get_path())
            except core.Error as e:
                show_error_dialog(
                    message=f"Error loading XML file{file.get_path()}:",
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
                    if self.obj == None:
                        self.dropdown.set_selected(0)
                        raise ValueError(f'Please select a file containing a spectrometer')
                else:
                    raise ValueError(f'namespace unimplemented: {self.namespace}')
            except:
                show_error_dialog(
                    message=f"Could not construct obj from XML file:",
                    detail=traceback.format_exc(),
                    window=self.get_root()
                    )
                return

            print(self.namespace)
            print(self.obj)

            self._updating = True
            if self.namespace == core.Screen:
                self.dropdown.set_selected(1)
            else:
                indexes = [i for i, val in enumerate(self.items)
                           if val.endswith('/'+self.obj.kind())]
                if len(indexes) == 1:
                    self.dropdown.set_selected(indexes[0])
                else:
                    pass
            self._updating = False


            # lazy import to break dependency cycle
            from .property_editor_box import PropertyEditorBox

            box = PropertyEditorBox(self.obj, hide=['Metric'])
            self.frame.set_child(box)
            box.connect("value-changed", self.on_child_value_changed)
            self.emit("object-changed")
