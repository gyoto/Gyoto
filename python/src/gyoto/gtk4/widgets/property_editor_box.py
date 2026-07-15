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

__all__ = ['PropertEditorBox']

import gi
gi.require_version("Gtk", "4.0")
from gi.repository import Gtk, GLib

from functools import wraps
import ctypes

from .filename_editor import FilenameEditor
from .scientific_spin import ScientificSpin
from .vector_scientific_spin import VectorScientificSpin
from .gyoto_object_chooser import GyotoObjectChooser

from ...core import Property
from ... import core, metric, astrobj, spectrum, spectrometer

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
            except core.Error as e:
                show_error_dialog(detail=str(e),
                                  message=f'Error setting {name}',
                                  widget=widget)
        return wrapper

    def __init__(self, obj, hide=[], first=[], *args, **kwargs):
        if "orientation" not in kwargs: kwargs['orientation']=Gtk.Orientation.VERTICAL
        if "spacing" not in kwargs: kwargs['spacing']=10
        super().__init__(*args, **kwargs)
        self.obj=obj
        self.hide=hide
        self.first=first
        self.populate_properties()

    def populate_properties(self):
        """Generates widgets for object properties"""
        parameters = self.obj.getPropertyNames()
        first_and_parameters = self.first + list(parameters)

        self.widgets=dict()

        for name in first_and_parameters:
            if name not in parameters: continue
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
                if param_type == core.Property.long_t:
                    lower=-(1 << (ctypes.sizeof(ctypes.c_long) * 8 - 1))
                    upper=(1 << (ctypes.sizeof(ctypes.c_long) * 8 - 1)) - 1
                elif param_type == core.Property.unsigned_long_t:
                    lower=0
                    upper=(1 << (ctypes.sizeof(ctypes.c_long) * 8))-1
                elif param_type == core.Property.size_t_t :
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
                # block scroll events
                scroll_controller = Gtk.EventControllerScroll()
                scroll_controller.connect("scroll", lambda *args: True)
                spin.add_controller(scroll_controller)
                self.widgets[name] = spin

            elif param_type == core.Property.bool_t:
                radio_true = Gtk.CheckButton(label=name)
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
                chooser = GyotoObjectChooser(metric,
                                             obj=getattr(self.obj, name))
                hbox.append(chooser)
                chooser.connect("object-changed", self.on_object_changed, name)
                chooser.connect("object-mutated", self.on_object_mutated, name)
                self.widgets[name] = chooser

            elif param_type == core.Property.screen_t:
                chooser = GyotoObjectChooser(core.Screen,
                                             obj=getattr(self.obj, name))
                hbox.append(chooser)
                chooser.connect("object-changed", self.on_object_changed, name)
                chooser.connect("object-mutated", self.on_object_mutated, name)
                self.widgets[name] = chooser

            elif param_type == core.Property.spectrum_t:
                chooser = GyotoObjectChooser(spectrum,
                                             obj=getattr(self.obj, name))
                hbox.append(chooser)
                chooser.connect("object-changed", self.on_object_changed, name)
                chooser.connect("object-mutated", self.on_object_mutated, name)
                self.widgets[name] = chooser

            elif param_type == core.Property.astrobj_t:
                chooser = GyotoObjectChooser(astrobj,
                                             obj=getattr(self.obj, name))
                hbox.append(chooser)
                chooser.connect("object-changed", self.on_object_changed, name)
                chooser.connect("object-mutated", self.on_object_mutated, name)
                self.widgets[name] = chooser

            elif param_type == core.Property.spectrometer_t:
                chooser = GyotoObjectChooser(spectrometer,
                                             obj=getattr(self.obj, name))
                hbox.append(chooser)
                chooser.connect("object-changed", self.on_object_changed, name)
                chooser.connect("object-mutated", self.on_object_mutated, name)
                self.widgets[name] = chooser

            self.append(frame)

            # Special cases
            if name == 'InitCoord':
                try:
                    st = astrobj.Star(self.obj)
                    isstar = True
                except core.Error:
                    isstar = False
                if isstar:
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
            # special case: InitCoord
            if len(new_value) == 7 and f'{name}:veltype' in self.widgets:
                if self.obj.Metric is not None:
                    pos = new_value[0:4]
                    vel = new_value[4:7]
                    tdot = self.obj.Metric.SysPrimeToTdot(pos, vel)
                    new_value = pos + [tdot] + [x*tdot for x in vel]

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

    @gtk_callback
    def on_3vel_toggled(self, widget, name, *args):
        '''Callback handling 3-velocity/4-velocity radio buttons

        Compute the right velocity and sets the VectorScientificSpin
        for InitCoord accordingly.

        '''
        coord = self.obj.get(name)
        if widget.get_active():
            # switch to representing as 3-velocity
            tdot = coord[4]
            if tdot == 0.: tdot=1.
            pos3vel = coord[0:4] + tuple((x/tdot for x in coord[5:8]))
            self.widgets[name].set_value(pos3vel)
        else:
            # switch back to showing 4-velocity
            self.widgets[name].set_value(coord)
