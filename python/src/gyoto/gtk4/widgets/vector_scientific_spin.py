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

__all__ = ['VectorScientificSpin']

import gi
gi.require_version("Gtk", "4.0")
from gi.repository import Gtk

from .scientific_spin import ScientificSpin

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
            # Block scroll events
            scroll_controller = Gtk.EventControllerScroll()
            scroll_controller.connect("scroll", lambda *args: True)
            spin.add_controller(scroll_controller)
        else:
            raise ValueError(f'itemclass must be one of ScientificSpin, Gtk.SpinButton')

        spin.connect("value-changed", self.on_spin_changed)
        return spin

