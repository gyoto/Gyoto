import inspect
import gyoto
import gi
gi.require_version("Gtk", "4.0")
from gi.repository import Gtk, GLib

class GyotoObjectChooser(Gtk.Box):
    __gsignals__ = {
        "object-changed": (gi.repository.GObject.SignalFlags.RUN_FIRST, None, ()),
        "object-mutated": (gi.repository.GObject.SignalFlags.RUN_FIRST, None, ())
    }
    def __init__(self, namespace, obj=None):
        super().__init__(orientation=Gtk.Orientation.VERTICAL,
                         spacing=0)
        self.obj=obj
        self.namespace=namespace
        self.items = ["-"]
        self.classes = []
        for name in dir(namespace):
            attr = getattr(namespace, name)
            if (inspect.isclass(attr) and
                issubclass(attr, namespace.Generic) and
                attr is not namespace.Generic):
                self.classes.append(attr)
                self.items.append(name)
        self.dropdown=Gtk.DropDown.new_from_strings(self.items)
        self.append(self.dropdown)
        print("HARE", type(self.obj))
        if self.obj is None:
            self.dropdown.set_selected(0)
        elif self.obj.kind() in self.items:
            self.dropdown.set_selected(self.items.index(self.obj.kind()))
        else:
            for i in range(len(self.items)):
                x = self.items[i]
                if not hasattr(namespace, x):
                    continue
                klass = getattr(namespace, x)
                try:
                    superobj = klass(superobj)
                    self.dropdown.set_selected(i)
                    break
                except:
                    continue

        self.dropdown.connect("notify::selected", self.on_dropdown_activated)

        self.frame = Gtk.Frame()
        self.append(self.frame)
        if self.obj is not None:
            box = PropertyEditorBox(self.obj)
            self.frame.set_child(box)
            box.connect("value-changed", self.on_child_value_changed)

    def on_child_value_changed(self, widget, *args):
        print("something changed in child")
        self.emit("object-mutated")
    
    def on_dropdown_activated(self, widget, *args):
        print("child changed")
        x = widget.get_selected_item().get_string()
        print(x)
        if x == '-':
            self.obj=None
            self.frame.set_child(None)
        else:
            self.obj=getattr(self.namespace, x)()
            box = PropertyEditorBox(self.obj)
            self.frame.set_child(box)
            box.connect("value-changed", self.on_child_value_changed)
        self.emit("object-changed")
    
    
class ScientificSpin(Gtk.Box):
    """
    Scientific-number editor:
      - accepts values like -1.234e-15
      - up/down buttons
      - relative increments
    """

    __gsignals__ = {
        "value-changed": (gi.repository.GObject.SignalFlags.RUN_FIRST, None, ()),
        "unit-changed": (gi.repository.GObject.SignalFlags.RUN_FIRST, None, ())
    }

    def __init__(self, value=1.0, rel_step=0.1, exponent_step=1, with_unit=False):
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
        self.exponent_step = exponent_step

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
        # Scientific notation is usually preferable for tiny numbers
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

    def __init__(self, obj, *args, **kwargs):
        if "orientation" not in kwargs: kwargs['orientation']=Gtk.Orientation.VERTICAL
        if "spacing" not in kwargs: kwargs['spacing']=10
        super().__init__(*args, **kwargs)
        self.obj=obj
        self.populate_properties()

    def populate_properties(self):
        """Generates widgets for object properties"""
        parameters = self.obj.getPropertyNames()

        for name in parameters:
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
                print(name + " is a double_t")
                spin = ScientificSpin(value=value, with_unit=prop.supportsUnits())
                spin.connect("value-changed", self.on_parameter_changed, name)
                spin.connect("unit-changed", self.on_unit_changed, name)
                hbox.append(spin)

            elif param_type == gyoto.core.Property.bool_t:
                print(name + " is a bool_t")
                radio_true = Gtk.CheckButton(label=name)
                radio_false = Gtk.CheckButton(label=prop.name_false)
                radio_false.set_group(radio_true)
                hbox.append(radio_true)
                hbox.append(radio_false)
                radio_true.set_active(value)
                radio_false.set_active(not value)
                radio_true.connect("toggled", self.on_parameter_changed, name)

            elif param_type == gyoto.core.Property.metric_t:
                print(name + " is a metric_t")
                chooser = GyotoObjectChooser(gyoto.metric, obj=getattr(self.obj, name))
                hbox.append(chooser)
                chooser.connect("object-changed", self.on_object_changed, name)
                chooser.connect("object-mutated", self.on_object_mutated, name)

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

            elif param_type == "string":
                entry = Gtk.Entry()
                entry.set_text(str(value))
                hbox.append(entry)

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

    def on_unit_changed(self, widget, name, *args):
        unit = widget.get_unit()
        print(unit)
        value = self.obj.get(name, unit)
        widget.set_value(value)
        
    def on_object_changed(self, widget, name, *args):
        """Called when a sub-object is reset

        For instance if the Metric of a Star is changed to another
        kind.

        """
        print(type(widget.obj))
        print(widget.obj)
        setattr(self.obj, name, widget.obj)
        self.emit('value-changed')

    def on_object_mutated(self, widget, name, *args):
        """Called when something changes in a sub-object

        We have nothing to do, except forward the signal.
        """
        self.emit('value-changed')

    def on_parameter_changed(self, widget, name, *args):
        """Callback unique pour tous les paramètres.
        Reçoit :
        - widget : le widget qui a déclenché le callback
        - name : le nom du paramètre
        - *args : arguments supplémentaires (ex: unité pour double_with_unit)
        """
        print("THERE", name)
        new_unit=None
        if isinstance(widget, ScientificSpin):
            print("ScientificSpin")
            new_value = widget.get_value()
            new_unit = widget.get_unit()

        elif isinstance(widget, Gtk.CheckButton):
            print("CheckButton")
            new_value = widget.get_active()

        elif isinstance(widget, Gtk.SpinButton):
            new_value = widget.get_value()
            new_unit = args[0] if args else None

        elif isinstance(widget, Gtk.Entry):
            new_value = widget.get_text()
            new_unit = None

        print(f"Paramètre '{name}' modifié : valeur={new_value}, unité={new_unit}")

        if new_unit is None:
            self.obj.set(name, new_value)
        else:
            self.obj.set(name, new_value, new_unit)

        self.emit('value-changed')
       
class ObjectEditor(Gtk.Window):

    @staticmethod
    def run(obj, blocking=True):
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

        # Connect to 
        self.vbox.connect("value-changed", self.on_value_changed)

    def on_close_request(self, *args):
        """Callback appelé quand la fenêtre est fermée (via la croix ou Alt+F4)."""
        if self.main_loop is not None:
            GLib.idle_add(self.main_loop.quit)
        return False  # Autorise la fermeture

    def on_value_changed(self, widget, *args):
        """Callback unique pour tous les paramètres.
        Reçoit :
        - widget : le widget qui a déclenché le callback
        - name : le nom du paramètre
        - *args : arguments supplémentaires (ex: unité pour double_with_unit)
        """
        print("HERE something changed")
            
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

# eventually these will become method of gyoto.core.Object
# cleanest is to do it in the swig file gyoto.i

def edit(self, blocking=True):
    ObjectEditor.run(self, blocking=blocking)

gyoto.core.Object.edit=edit

#gg = gyoto.std.KerrBL()
#gg.edit()

#st = gyoto.std.FixedStar()
#st.Metric = gg
#st.edit()

# Exemple d'utilisation
# if __name__ == "__main__":
#     # Mockup d'un objet avec introspection
#     class MockObject:
#         def get_parameters(self):
#             return [
#                 ("Temperature", "double_with_unit", (25.0, "°C")),
#                 ("Filename", "filename", "/path/to/file.txt"),
#                 ("Count", "long", 42),
#                 ("Name", "string", "Test"),
#             ]

#     obj = MockObject()
#     app = Gtk.Application(application_id="org.example.ParameterEditor")
#     app.connect("activate", lambda app: ObjectEditor(obj).present())
#     # to be checked whether "blocking" should be True or False
#     app.run(None)
