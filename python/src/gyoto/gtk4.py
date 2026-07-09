from functools import wraps

import gi
gi.require_version("Gtk", "4.0")
from gi.repository import Gtk, GLib

import gyoto

class GyotoObjectChooser(Gtk.Box):
    """Gtk widget for choosing and editing a Gyoto object kind

    The widget is composed of a drop-down menu from which a kind can
    be choosen, and a PropertyEditorBox for this kind.

    Parameters:

      namespace: the Gyoto namespace to choose from (gyoto.metric,
          gyoto.astrobj, gyoto.spectrum or gyoto.spectrometer).

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
        self.obj=obj
        self.namespace=namespace
        self.items = ["-"] + [x for x in
                              namespace.Generic.registeredPluginsSlashKinds()]
        self.dropdown=Gtk.DropDown.new_from_strings(self.items)
        self.append(self.dropdown)

        if self.obj is None:
            self.dropdown.set_selected(0)
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
            box = PropertyEditorBox(self.obj)
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
            plg, knd = x.split('/')
            self.obj=self.namespace.Generic(knd, (plg,))
            box = PropertyEditorBox(self.obj)
            self.frame.set_child(box)
            box.connect("value-changed", self.on_child_value_changed)
        self.emit("object-changed")
    
    
class ScientificSpin(Gtk.Box):
    """Similar to Gtk.SpinButton in floating point notation with unit

    This widget is composed of:
        - a Gtk.Entry for the value, which accepts notations like
          1.3e-12;
        - an optional Gkt.Entry field to edit a unit;
        - a pair of buttons to increment of decrement the value.

    Parameters:
      value: the initial value;
      rel_step: the (relative) increment step (0.1);
      with_unit: whether to show the unit field (False).

    Public methods:
      Callbacks can query or set the two entry boxes using the four methods:
      get_value, set_value, get_unit and set_unit

    Signals:
      The widget emits value-changed and unit-changed when the
      respective field changes.

    """

    __gsignals__ = {
        "value-changed": (gi.repository.GObject.SignalFlags.RUN_FIRST, None, ()),
        "unit-changed": (gi.repository.GObject.SignalFlags.RUN_FIRST, None, ())
    }

    def __init__(self, value=1.0, rel_step=0.1, with_unit=False):
        super().__init__(orientation=Gtk.Orientation.HORIZONTAL,
                         spacing=0)

        self.add_css_class("scientific-spin")

        css = Gtk.CssProvider()
        css.load_from_data(b"""
        .scientific-spin button {
            padding: 0;
            margin: 0;
            min-height: 0;
            min-width: 0;
        }

        .scientific-spin button image {
            -gtk-icon-size: 10px;
        }

        .scientific-spin entry {
            padding-top: 0;
            padding-bottom: 0;
        }
        """)

        Gtk.StyleContext.add_provider_for_display(
            self.get_display(),
            css,
            Gtk.STYLE_PROVIDER_PRIORITY_APPLICATION
        )

        self.rel_step = rel_step

        self.value_entry = Gtk.Entry()
        self.value_entry.set_hexpand(True)
        self.value_entry.set_tooltip_text("value")
        self.set_value(value)

        self.append(self.value_entry)

        if with_unit:
            self.unit_entry = Gtk.Entry()
            self.unit_entry.set_hexpand(True)
            self.set_unit("")
            self.append(self.unit_entry)
            self.unit_entry.connect("activate", self.on_unit_changed)
            self.unit_entry.set_tooltip_text("unit")
        else:
            self.unit_entry=None

        buttons = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=0)

        self.up = Gtk.Button()
        self.down = Gtk.Button()

        self.up.set_icon_name("pan-up-symbolic")
        self.down.set_icon_name("pan-down-symbolic")

        self.up.add_css_class("flat")
        self.down.add_css_class("flat")
        
        self.value_entry.connect("changed", self.on_value_changed)
        self.up.connect("clicked", self.increment, +1)
        self.down.connect("clicked", self.increment, -1)

        buttons.append(self.up)
        buttons.append(self.down)

        self.append(buttons)


    def get_value(self):
        return float(self.value_entry.get_text())


    def set_value(self, value):
        self.value_entry.set_text(f"{value:.8e}")

    def get_unit(self):
        if self.unit_entry is None:
            return None
        return self.unit_entry.get_text()

    def set_unit(self, unit):
        self.unit_entry.set_text(unit)

    def increment(self, button, direction):
        try:
            value = self.get_value()
        except ValueError:
            value = 0.0

        # normal case: relative increment
        if value != 0:
            step = self.rel_step
            if direction == +1:
                value *= 1.+self.rel_step
            elif direction == -1:
                value /= 1.+self.rel_step
            else:
                raise ValueError("direction should be +1 or -1")

        # special case: start from zero
        else:
            value = direction * 1.0

        self.set_value(value)

    def on_value_changed(self, entry):
        try:
            value = float(entry.get_text())
        except ValueError:
            return
        self.emit("value-changed")

    def on_unit_changed(self, entry):
        self.emit("unit-changed")

class PropertyEditorBox(Gtk.Box):
    """
    A compound widget that presents the properties of a Gyoto object:
      - obj: the Gyoto object
      - signals: value-changed
    """

    __gsignals__ = {
        "value-changed": (gi.repository.GObject.SignalFlags.RUN_FIRST, None, ())
    }

    @staticmethod
    def gtk_callback(method):
        """Show Gyoto errors in dialog windows

        This is a deforator, use like this:

          @gtk_callback
          def method(self, widget, name *args):
              ...

        Now whenever a Gyoto error occurs in method(), it will be
        displayed in a Gtk dialog instead of the terminal.
        """
        @wraps(method)
        def wrapper(self, widget, name, *args):
            try:
                return method(self, widget, name, *args)
            except gyoto.core.Error as e:
                show_error_dialog(detail=str(e),
                                  message=f'Error setting {name}',
                                  widget=widget)
        return wrapper

    def __init__(self, obj, *args, **kwargs):
        if "orientation" not in kwargs: kwargs['orientation']=Gtk.Orientation.VERTICAL
        if "spacing" not in kwargs: kwargs['spacing']=10
        super().__init__(*args, **kwargs)
        self.obj=obj
        self.populate_properties()

    def populate_properties(self):
        """Generates widgets for object properties"""
        parameters = self.obj.getPropertyNames()

        self.widgets=dict()

        for name in parameters:
            if name in self. widgets: continue
            prop = self.obj.property(name)
            value = self.obj.get(prop)
            param_type = prop.type
            frame = Gtk.Frame()
            frame.set_label(name)
            frame.set_label_align(0.0)  # Aligner le label à gauche
            frame.set_tooltip_text(prop.doc)

            hbox = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=5)
            frame.set_child(hbox)

            if param_type == gyoto.core.Property.double_t:
                spin = ScientificSpin(value=value,
                                      with_unit=prop.supportsUnits())
                spin.connect("value-changed", self.on_parameter_changed, name)
                spin.connect("unit-changed", self.on_unit_changed, name)
                hbox.append(spin)
                self.widgets[name] = spin

            elif param_type == gyoto.core.Property.bool_t:
                radio_true = Gtk.CheckButton(label=name)
                radio_false = Gtk.CheckButton(label=prop.name_false)
                radio_false.set_group(radio_true)
                hbox.append(radio_true)
                hbox.append(radio_false)
                radio_true.set_active(value)
                radio_false.set_active(not value)
                radio_true.connect("toggled", self.on_parameter_changed, name)
                self.widgets[name] = radio_true

            elif param_type == gyoto.core.Property.metric_t:
                chooser = GyotoObjectChooser(gyoto.metric,
                                             obj=getattr(self.obj, name))
                hbox.append(chooser)
                chooser.connect("object-changed", self.on_object_changed, name)
                chooser.connect("object-mutated", self.on_object_mutated, name)
                self.widgets[name] = chooser

            elif param_type == gyoto.core.Property.spectrum_t:
                chooser = GyotoObjectChooser(gyoto.spectrum,
                                             obj=getattr(self.obj, name))
                hbox.append(chooser)
                chooser.connect("object-changed", self.on_object_changed, name)
                chooser.connect("object-mutated", self.on_object_mutated, name)
                self.widgets[name] = chooser

            elif param_type == gyoto.core.Property.astrobj_t:
                chooser = GyotoObjectChooser(gyoto.astrobj,
                                             obj=getattr(self.obj, name))
                hbox.append(chooser)
                chooser.connect("object-changed", self.on_object_changed, name)
                chooser.connect("object-mutated", self.on_object_mutated, name)
                self.widgets[name] = chooser

            elif param_type == gyoto.core.Property.spectrometer_t:
                chooser = GyotoObjectChooser(gyoto.spectrometer,
                                             obj=getattr(self.obj, name))
                hbox.append(chooser)
                chooser.connect("object-changed", self.on_object_changed, name)
                chooser.connect("object-mutated", self.on_object_mutated, name)
                self.widgets[name] = chooser

            elif param_type == "long":
                adjustment = Gtk.Adjustment(
                    value=value,
                    lower=-1e10,
                    upper=1e10,
                    step_increment=1,
                    page_increment=10,
                    page_size=0
                )
                spin = Gtk.SpinButton()
                spin.set_adjustment(adjustment)
                spin.set_numeric(True)
                spin.set_digits(20)
                hbox.append(spin)

            elif param_type == gyoto.core.Property.string_t:
                entry = Gtk.Entry()
                entry.set_hexpand(True)
                entry.set_text(str(value))
                hbox.append(entry)
                entry.connect("activate", self.on_parameter_changed, name)
                self.widgets[name] = entry

            elif param_type == "filename":
                entry = Gtk.Entry()
                entry.set_text(str(value))
                button = Gtk.Button(label="...")
                button.connect("clicked", self.on_file_chooser_clicked, entry)
                hbox.append(entry)
                hbox.append(button)

            elif param_type == "double_with_unit":
                val, unit = value
                adjustment = Gtk.Adjustment(
                    value=val,
                    lower=-1e10,
                    upper=1e10,
                    step_increment=0.1,
                    page_increment=1.0,
                    page_size=0
                )
                spin = Gtk.SpinButton()
                spin.set_adjustment(adjustment)
                spin.set_numeric(True)
                hbox.append(spin)

                unit_label = Gtk.Label(label=unit)
                hbox.append(unit_label)

                unit_entry = Gtk.Entry()
                unit_entry.set_text(unit)
                hbox.append(unit_entry)

            self.append(frame)

    @gtk_callback
    def on_unit_changed(self, widget, name, *args):
        unit = widget.get_unit()
        value = self.obj.get(name, unit)
        widget.set_value(value)

    @gtk_callback
    def on_object_changed(self, widget, name, *args):
        """Called when a sub-object is reset

        For instance if the Metric of a Star is changed to another
        kind.

        """
        self.obj.set(name, widget.obj)
        self.emit('value-changed')

    @gtk_callback
    def on_object_mutated(self, widget, name, *args):
        """Called when something changes in a sub-object

        We have nothing to do, except forward the signal.
        """
        self.emit('value-changed')

    @gtk_callback
    def on_parameter_changed(self, widget, name, *args):
        """Most widgets are connected to this callback.

        Receives:
          widget: the widget that initiated the signal
          name: the name of the Property this widget allows editing
        """
        new_unit=None
        if isinstance(widget, ScientificSpin):
            new_value = widget.get_value()
            new_unit = widget.get_unit()

        elif isinstance(widget, Gtk.CheckButton):
            new_value = widget.get_active()

        elif isinstance(widget, Gtk.SpinButton):
            new_value = widget.get_value()
            new_unit = args[0] if args else None

        elif isinstance(widget, Gtk.Entry):
            new_value = widget.get_text()
            new_unit = None

        if new_unit is None:
            self.obj.set(name, new_value)
        else:
            self.obj.set(name, new_value, new_unit)

        self.emit('value-changed')
       
class ObjectEditor(Gtk.Window):
    """A Gyoto object editor

    Displays a window giving an editable view of all the properties of
    a Gyoto object.

    Parameters:
      obj: the object to edit
      blocking: whether this object should manage the GLib loop

    """

    @staticmethod
    def run(obj, blocking=True):
        """Contruct a Gyoto ObjectEditor window and run it

        synopsis:
         ObjectEditor.run(obj, [blocking=True])
        """
        win = ObjectEditor(obj, blocking)
        win.present()
        if blocking:
            win.main_loop.run()

    def __init__(self, obj, blocking=True):
        Gtk.Window.__init__(self, title="Gyoto Object Editor")
        self.set_default_size(400, 600)
        self.obj = obj
        if blocking:
            self.main_loop=GLib.MainLoop()
        else:
            self.main_loop=None

        self.connect("close-request", self.on_close_request)

        # Main container
        self.scrolled_window = Gtk.ScrolledWindow()
        self.set_child(self.scrolled_window)

        # Vertical box for property widgets
        self.vbox = PropertyEditorBox(obj)
        self.vbox.set_margin_top(10)
        self.vbox.set_margin_bottom(10)
        self.vbox.set_margin_start(10)
        self.vbox.set_margin_end(10)
        self.scrolled_window.set_child(self.vbox)

        # We could connect to value-changed to react wheneverr obj changes 
        # self.vbox.connect("value-changed", self.on_value_changed)
        # but we don't actually need to.

    def on_close_request(self, *args):
        if self.main_loop is not None:
            GLib.idle_add(self.main_loop.quit)
        return False

    def on_file_chooser_clicked(self, button, entry):
        """Ouvre un dialogue pour choisir un fichier."""
        dialog = Gtk.FileChooserDialog(
            title="Select a file",
            parent=self,
            action=Gtk.FileChooserAction.OPEN
        )
        dialog.add_button("Cancel", Gtk.ResponseType.CANCEL)
        dialog.add_button("Open", Gtk.ResponseType.OK)
        dialog.connect("response", self.on_file_chooser_response, entry)
        dialog.show()

    def on_file_chooser_response(self, dialog, response, entry):
        """Gère la réponse du dialogue de sélection de fichier."""
        if response == Gtk.ResponseType.OK:
            file = dialog.get_file()
            if file is not None:
                entry.set_text(file.get_path())
                dialog.destroy()

def show_error_dialog(message="An error occurred", detail=None,
                      window=None, widget=None):
    """Show a pop-up error dialog

    Parameters:
      message: main message
      detail: more details about the error
      window (optional): the window this dialog should be transient of
      widget (optional): shortcut to window=widget.get_root()

    """
    if window is None and widget is not None:
        window=widget.get_root()

    dialog = Gtk.AlertDialog(
        message=message,
        buttons=["OK"]
    )

    if detail is not None: dialog.set_detail(detail)

    dialog.show(window)

# The following should be achieved using Swig's extend mechanism
def edit(self, blocking=True):
    """Display a GUI to edit this object's Properties"""
    ObjectEditor.run(self, blocking=blocking)

gyoto.core.Object.edit=edit
