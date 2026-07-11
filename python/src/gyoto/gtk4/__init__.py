from functools import wraps
import ctypes

import gi
gi.require_version("Gtk", "4.0")
from gi.repository import Gtk, GLib

import gyoto

GLib.set_prgname("Gyoto")

### Compond widgets suited for our needs

## A filename editor / file chooser with pop-up dialog:
# ┌───────────────────────────────┐
# │ /home/user/file.dat        📂 │
# └───────────────────────────────┘
#
class FilenameEditor(Gtk.Box):
    """A compound widget to select a file

    Composed on an entry field and a button for browsng.

    Parameters:
      value: initial value.

    Signals:
      value-changed emitted when a file is selected in the browser or
          when the user hits enter in the entry field.

    """
    __gsignals__ = {
        "value-changed": (
            gi.repository.GObject.SignalFlags.RUN_FIRST,
            None,
            ()
        )
    }

    def __init__(self, value=""):
        super().__init__(
            orientation=Gtk.Orientation.HORIZONTAL,
            spacing=4
        )

        self.entry = Gtk.Entry()
        self.entry.set_hexpand(True)
        self.entry.set_text(value)

        self.button = Gtk.Button()
        self.button.set_icon_name("folder-open-symbolic")
        self.button.set_tooltip_text("Choose a file")
        self.button.add_css_class("flat")

        self.append(self.entry)
        self.append(self.button)

        self.entry.connect(
            "activate",
            self.on_activate
        )

        self.button.connect(
            "clicked",
            self.on_browse
        )

    def get_value(self):
        return self.entry.get_text()

    def set_value(self, value):
        self.entry.set_text(value)

    def on_activate(self, entry):
        self.emit("value-changed")

    def on_browse(self, button):

        dialog = Gtk.FileDialog()

        dialog.open(
            self.get_root(),
            None,
            lambda dialog, result:
                self.on_file_selected(dialog, result)
        )

    def on_file_selected(self, dialog, result):
        try:
            file = dialog.open_finish(result)
        except GLib.Error:
            return

        if file is not None:
            self.set_value(file.get_path())
            self.emit("value-changed")

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
class GyotoObjectChooser(Gtk.Box):
    """Gtk widget for choosing and editing a Gyoto object kind

    The widget is composed of a drop-down menu from which a kind can
    be choosen, and a PropertyEditorBox for this kind. Selecting
    another kind destroys the current PropertyEditorBox and creates a
    new one.

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
        if self.namespace == gyoto.core.Screen:
            self.items = ["-", "built-in/Screen"]
        else:
            self.items = ["-"] + [x for x in
                                  namespace.Generic.registeredPluginsSlashKinds()]
        self.dropdown=Gtk.DropDown.new_from_strings(self.items)
        self.append(self.dropdown)

        if self.obj is None:
            self.dropdown.set_selected(0)
        else:
            if self.namespace == gyoto.core.Screen:
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
            if self.namespace == gyoto.core.Screen:
                self.obj = gyoto.core.Screen()
            else:
                plg, knd = x.split('/')
                self.obj=self.namespace.Generic(knd, (plg,))
            box = PropertyEditorBox(self.obj, hide=['Metric'])
            self.frame.set_child(box)
            box.connect("value-changed", self.on_child_value_changed)
        self.emit("object-changed")
    

## A SpinButton in floating-point notation with option unit field
#
# ┌────────────────────────────────────────────────┐
# │ ┌────────────────────┐ ┌───────────────┐ ┌───┐ │
# │ │ (Gtk.Entry:)       │ │ (Gtk.Entry:)  │ │ ▲ │ │
# │ │ 1.23456789e-12     │ │ m^2           │ │ ▼ │ │
# │ └────────────────────┘ └───────────────┘ └───┘ │
# └────────────────────────────────────────────────┘
#
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

        # Set this to True whenever changes should *not* emit
        # value-changed
        self._updating=False

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
        if value is not None:
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


    def set_value(self, value, *, emit=False):
        """Set the value field.

        By default, setting the value does not emit the
        ``value-changed`` signal. Set ``emit`` to ``True`` to emit the
        signal after updating the value.

        """
        was_updating=self._updating
        self._updating = was_updating or not emit
        try:
            self.value_entry.set_text(f"{value:.8e}")
        finally:
            self._updating = was_updating

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

        self.set_value(value, emit=True)

    def on_value_changed(self, entry):
        """Emit value-changed

        In some conditions, emission will be blocked by setting
        self._updating to True.

        """
        if self._updating: return
        try:
            value = float(entry.get_text())
        except ValueError:
            return

        self.emit("value-changed")

    def on_unit_changed(self, entry):
        self.emit("unit-changed")

## A collection of ScientificSpin widgets for a vector<double>
#
# ┌─────────────────────────────────────────────┐
# │ Unit: [ km ]                                │
# │                                             │
# │ 1.00000000e+00                ▲             │
# │                               ▼             │
# │                                             │
# │ 2.00000000e+00                ▲             │
# │                               ▼             │
# │                                             │
# │ 3.00000000e+00                ▲             │
# │                               ▼             │
# │                                             │
# │                         [+]  [-]            │
# └─────────────────────────────────────────────┘
#
class VectorScientificSpin(Gtk.Box):
    """Editor for std::vector<double> with optional unit."""

    __gsignals__ = {
        "value-changed": (
            gi.repository.GObject.SignalFlags.RUN_FIRST,
            None,
            ()
        ),
        "unit-changed": (
            gi.repository.GObject.SignalFlags.RUN_FIRST,
            None,
            ()
        )
    }

    def __init__(self, value=None, with_unit=False, rel_step=0.1,
                 itemclass=ScientificSpin):
        super().__init__(
            orientation=Gtk.Orientation.VERTICAL,
            spacing=6
        )
        if itemclass not in (ScientificSpin, Gtk.SpinButton):
            raise ValueError(f'itemclass must be one of ScientificSpin, Gtk.SpinButton')

        self.rel_step = rel_step
        self.spins = []
        self.itemclass=itemclass

        #
        # Optional unit field
        #
        if with_unit:
            hbox = Gtk.Box(
                orientation=Gtk.Orientation.HORIZONTAL,
                spacing=6
            )

            label = Gtk.Label(label="Unit:")
            label.set_xalign(0)

            self.unit_entry = Gtk.Entry()
            self.unit_entry.set_hexpand(True)
            self.unit_entry.connect(
                "activate",
                lambda *_: self.emit("unit-changed")
            )

            hbox.append(label)
            hbox.append(self.unit_entry)

            self.append(hbox)

        else:
            self.unit_entry = None

        #
        # Container for the ScientificSpins
        #
        self.values_box = Gtk.Box(
            orientation=Gtk.Orientation.VERTICAL,
            spacing=2
        )

        self.append(self.values_box)

        #
        # + / - buttons
        #
        buttons = Gtk.Box(
            orientation=Gtk.Orientation.HORIZONTAL,
            spacing=6
        )

        add = Gtk.Button(icon_name="list-add-symbolic")
        remove = Gtk.Button(icon_name="list-remove-symbolic")

        add.add_css_class("flat")
        remove.add_css_class("flat")

        add.connect("clicked", self.on_add)
        remove.connect("clicked", self.on_remove)

        buttons.append(add)
        buttons.append(remove)

        self.append(buttons)

        #
        # Initial values
        #
        if value is None:
            value = [0.0]

        self.set_value(value)

    #
    # Public API
    #

    def get_value(self):
        if self.itemclass == ScientificSpin:
            return [spin.get_value() for spin in self.spins]
        if self.itemclass == Gtk.SpinButton:
            return [spin.get_value_as_int() for spin in self.spins]

    def set_value(self, values):

        #
        # Remove excess widgets
        #
        while len(self.spins) > len(values):
            spin = self.spins.pop()
            self.values_box.remove(spin)

        #
        # Add missing widgets
        #
        while len(self.spins) < len(values):
            spin = self._new_spin()

            self.spins.append(spin)
            self.values_box.append(spin)

        #
        # Update values
        #
        for spin, value in zip(self.spins, values):
            spin.set_value(value)

    def get_unit(self):
        if self.unit_entry is None:
            return None
        return self.unit_entry.get_text()

    def set_unit(self, unit):
        if self.unit_entry is not None:
            self.unit_entry.set_text(unit)

    #
    # Internals
    #

    def on_spin_changed(self, spin):
        self.emit("value-changed")

    def on_add(self, button):

        if self.spins:
            value = self.spins[-1].get_value()
        else:
            value = 0.0

        spin = self._new_spin(value)
        self.spins.append(spin)
        self.values_box.append(spin)

        self.emit("value-changed")

    def on_remove(self, button):

        if len(self.spins) <= 1:
            return

        spin = self.spins.pop()

        self.values_box.remove(spin)

        self.emit("value-changed")

    def _new_spin(self, value=None):
        """Private metod to add a new ScientificSpin
        """
        if self.itemclass == ScientificSpin:
            spin = ScientificSpin(
                value=value,
                rel_step=self.rel_step,
                with_unit=False
            )
        elif self.itemclass == Gtk.SpinButton:
            lower=0
            upper=(1 << (ctypes.sizeof(ctypes.c_long) * 8))-1

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
        else:
            raise ValueError(f'itemclass must be one of ScientificSpin, Gtk.SpinButton')

        spin.connect("value-changed", self.on_spin_changed)
        return spin

## A panel to edit the properties of a gyoto.core.Object
#
# ┌───────────────────────────────────────┐
# │ (PropertyEditorBox:)                  │
# │                                       │
# │ ┌─ property1 ───────────────────────┐ │
# │ │ widget for property1              │ │
# │ └───────────────────────────────────┘ │
# │ ┌─ property2 ───────────────────────┐ │
# │ │ widget for property2              │ │
# │ └───────────────────────────────────┘ │
# │ ┌─ property3 ───────────────────────┐ │
# │ │ widget for property3              │ │
# │ └───────────────────────────────────┘ │
# │ ...                                   │
# └───────────────────────────────────────┘
#
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

    def __init__(self, obj, hide=[], *args, **kwargs):
        if "orientation" not in kwargs: kwargs['orientation']=Gtk.Orientation.VERTICAL
        if "spacing" not in kwargs: kwargs['spacing']=10
        super().__init__(*args, **kwargs)
        self.obj=obj
        self.hide=hide
        self.populate_properties()

    def populate_properties(self):
        """Generates widgets for object properties"""
        parameters = self.obj.getPropertyNames()

        self.widgets=dict()

        for name in parameters:
            if name in self.hide: continue
            if name in self.widgets: continue
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

            elif param_type in (gyoto.core.Property.long_t,
                                gyoto.core.Property.unsigned_long_t,
                                gyoto.core.Property.size_t_t):
                if param_type == gyoto.core.Property.long_t:
                    lower=-(1 << (ctypes.sizeof(ctypes.c_long) * 8 - 1))
                    upper=(1 << (ctypes.sizeof(ctypes.c_long) * 8 - 1)) - 1
                elif param_type == gyoto.core.Property.unsigned_long_t:
                    lower=0
                    upper=(1 << (ctypes.sizeof(ctypes.c_long) * 8))-1
                elif param_type == gyoto.core.Property.size_t_t :
                    lower=0
                    upper=(1 << (ctypes.sizeof(ctypes.c_size_t) * 8))-1
                else:
                    raise("bug")

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

            elif param_type == gyoto.core.Property.string_t:
                entry = Gtk.Entry()
                entry.set_hexpand(True)
                entry.set_text(str(value))
                hbox.append(entry)
                entry.connect("activate", self.on_parameter_changed, name)
                self.widgets[name] = entry

            elif param_type == gyoto.core.Property.filename_t:
                editor = FilenameEditor(str(value))
                editor.connect(
                    "value-changed",
                    self.on_parameter_changed,
                    name
                )
                hbox.append(editor)
                self.widgets[name] = editor

            elif param_type in (gyoto.core.Property.vector_double_t,
                                gyoto.core.Property.vector_unsigned_long_t):
                itemclass = (ScientificSpin
                             if param_type == gyoto.core.Property.vector_double_t
                             else Gtk.SpinButton)
                vector = VectorScientificSpin(value=value,
                                              with_unit=prop.supportsUnits(),
                                              itemclass=itemclass)
                vector.connect("value-changed", self.on_parameter_changed, name)
                vector.connect("unit-changed", self.on_unit_changed, name)
                hbox.append(vector)
                self.widgets[name] = vector

            elif param_type == gyoto.core.Property.metric_t:
                chooser = GyotoObjectChooser(gyoto.metric,
                                             obj=getattr(self.obj, name))
                hbox.append(chooser)
                chooser.connect("object-changed", self.on_object_changed, name)
                chooser.connect("object-mutated", self.on_object_mutated, name)
                self.widgets[name] = chooser

            elif param_type == gyoto.core.Property.screen_t:
                chooser = GyotoObjectChooser(gyoto.core.Screen,
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
            new_value = widget.get_value_as_int()

        elif isinstance(widget, Gtk.Entry):
            new_value = widget.get_text()

        elif isinstance(widget, FilenameEditor):
            new_value = widget.get_value()

        elif isinstance(widget, VectorScientificSpin):
            new_value = widget.get_value()
            new_unit = widget.get_unit()

        if new_unit is None:
            self.obj.set(name, new_value)
        else:
            self.obj.set(name, new_value, new_unit)

        self.emit('value-changed')

    @gtk_callback
    def on_file_chooser_clicked(self, button, entry):
        """File selection dialog for filename_t."""
        dialog = Gtk.FileDialog()
        dialog.open(
            self.get_root(),
            None,
            lambda dialog, result: self.on_file_selected(dialog, result, entry)
        )

    @gtk_callback
    def on_file_selected(self, dialog, result, entry):
        """What to do when the user selected a file using the dialog"""
        try:
            file = dialog.open_finish(result)
        except GLib.Error:
            # User cancelled
            return

        if file is not None:
            entry.set_text(file.get_path())

### Main windows for various tools

## A window to actually edit an object's properties
#
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


## This function displaysan error message in a Gtk dialog
#
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

### The following should be achieved using Swig's extend mechanism
def edit(self, blocking=True):
    """Display a GUI to edit this object's Properties"""
    ObjectEditor.run(self, blocking=blocking)

gyoto.core.Object.edit=edit
