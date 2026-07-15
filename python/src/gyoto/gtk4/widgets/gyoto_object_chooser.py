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
from gi.repository import Gtk

from ... import core

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
                    show_error_dialog('could not determine plugin/kind for object')

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
        x = widget.get_selected_item().get_string()
        if x == '-':
            self.obj=None
            self.frame.set_child(None)
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
    
