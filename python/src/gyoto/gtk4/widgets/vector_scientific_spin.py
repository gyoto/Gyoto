"""VectorScientificSpin: GTK4 Widget for Vector<double> Editing

This module provides a custom GTK4 widget for editing vectors of
floating-point values (std::vector<double> in C++), with optional unit
display and dynamic resizing.

Widget Layout
------------
    ┌─────────────────────────────────────────────┐
    │ Unit: [ km ]                                │
    │                                             │
    │ 1.00000000e+00                ▲             │
    │                               ▼             │
    │                                             │
    │ 2.00000000e+00                ▲             │
    │                               ▼             │
    │                                             │
    │ 3.00000000e+00                ▲             │
    │                               ▼             │
    │                                             │
    │                         [+]  [-]            │
    └─────────────────────────────────────────────┘

Description
-----------
- **Unit Field**: Optional text entry for displaying measurement units
- **Value Fields**: Vertical list of ScientificSpin widgets (one per
    vector element)
- **Add/Remove Buttons**: Dynamically add or remove vector elements

Features:
- Supports both ScientificSpin and Gtk.SpinButton as item classes
- Automatic value validation
- Signal emission on changes

"""

__all__ = ['VectorScientificSpin']

import gi
gi.require_version("Gtk", "4.0")
from gi.repository import Gtk

from .scientific_spin import ScientificSpin

class VectorScientificSpin(Gtk.Box):
    """A GTK4 widget for editing vectors of real or integer values.

    This widget provides a vertical list of ScientificSpin or
    Gtk.SpinButton widgets for editing individual vector components,
    with optional unit display and buttons to add/remove elements
    dynamically.

    The widget is composed of:
        - An optional unit field at the top (if with_unit=True)
        - A vertical box containing ScientificSpin or Gtk.SpinButton
          widgets for each vector element
        - Add/Remove buttons at the bottom to modify the vector size

    Attributes:
        hold (bool): If True, don't emit signals.
        rel_step (float): Relative step size for ScientificSpin widgets
        spins (list): List of spin widgets for each vector element
        itemclass (type): Class to use for spin widgets
            (ScientificSpin or Gtk.SpinButton)
        unit_entry (Gtk.Entry or None): Unit field entry widget
        values_box (Gtk.Box): Container for the ScientificSpin widgets

    Signals:
        value-changed: Emitted when any vector component value changes
        unit-changed: Emitted when the unit changes (if unit field is
            enabled)

    Example:
        vector_spin = VectorScientificSpin(value=[1.0, 2.0, 3.0],
            with_unit=True)
        vector_spin.connect("value-changed", lambda w: print(w.get_value()))

    """

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

    hold = False

    def __init__(self, value=None, with_unit=False, rel_step=0.1,
                 itemclass=ScientificSpin):
        """Initialize the VectorScientificSpin widget.

        Args:
            value (list of float or None): Initial vector values. If
                None, defaults to [0.0].
            with_unit (bool): Whether to show a unit field (default: False).
            rel_step (float): Relative step size for ScientificSpin
                widgets (default: 0.1).
            itemclass (type): Class to use for spin widgets. Must be
                either ScientificSpin or Gtk.SpinButton (default:
                ScientificSpin).

        Raises:
            ValueError: If itemclass is not ScientificSpin or
                Gtk.SpinButton.

        """
        super().__init__(
            orientation=Gtk.Orientation.VERTICAL,
            spacing=6
        )
        if itemclass not in (ScientificSpin, Gtk.SpinButton):
            raise ValueError(f'itemclass must be one of ScientificSpin, Gtk.SpinButton')

        self.rel_step = rel_step
        self.spins = []
        self.itemclass = itemclass

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
        # + / - buttons for adding/removing vector elements
        # ■ for hold functionality
        #
        buttons = Gtk.Box(
            orientation=Gtk.Orientation.HORIZONTAL,
            spacing=6
        )

        add = Gtk.Button(icon_name="list-add-symbolic")
        remove = Gtk.Button(icon_name="list-remove-symbolic")
        stop = Gtk.ToggleButton(icon_name="media-playback-stop-symbolic")

        add.add_css_class("flat")
        remove.add_css_class("flat")
        stop.add_css_class("flat")

        add.set_tooltip_text(
            "Add item to vector"
        )
        remove.set_tooltip_text(
            "Remove last element from vector"
        )
        stop.set_tooltip_text(
            "Toggle stop mode: when active, prevents signaling upon parameter changes"
        )

        add.connect("clicked", self.on_add)
        remove.connect("clicked", self.on_remove)
        stop.connect("toggled", self.on_stop)

        buttons.append(add)
        buttons.append(remove)
        buttons.append(stop)

        self.append(buttons)

        #
        # Initialize with provided values
        #
        if value is None:
            value = [0.0]

        self.set_value(value)

    #
    # Public API
    #

    def get_value(self):
        """Get the current vector values.

        Returns:
            list: A list of numeric values from each spin widget.
                 If itemclass is ScientificSpin, returns floats.
                 If itemclass is Gtk.SpinButton, returns integers.
        """
        if self.itemclass == ScientificSpin:
            return [spin.get_value() for spin in self.spins]
        if self.itemclass == Gtk.SpinButton:
            return [spin.get_value_as_int() for spin in self.spins]

    def set_value(self, values):
        """Set the vector values.

        Dynamically adjusts the number of spin widgets to match the
        length of the input list, then updates each widget with the
        corresponding value.

        Args:
            values (list): List of numeric values to set.

        """
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
        """Get the current unit text.

        Returns:
            str or None: The current unit text, or None if unit field
                is disabled.

        """
        if self.unit_entry is None:
            return None
        return self.unit_entry.get_text()

    def set_unit(self, unit):
        """Set the unit field text.

        Args:
            unit (str): The unit text to display.
        """
        if self.unit_entry is not None:
            self.unit_entry.set_text(unit)

    #
    # Internals
    #

    def on_spin_changed(self, spin):
        """Handle value changes from individual spin widgets.

        Emits the 'value-changed' signal when any spin widget value changes.
        Args:
        
            spin (ScientificSpin or Gtk.SpinButton): The spin widget
                that changed.

        """
        if not self.hold:
            self.emit("value-changed")

    def on_add(self, button):
        """Handle click on the add button.

        Adds a new spin widget to the vector with a value equal to the last
        widget's value (or 0.0 if no widgets exist).

        Args:
            button (Gtk.Button): The add button that was clicked.
        """
        if self.spins:
            value = self.spins[-1].get_value()
        else:
            value = 0.0

        spin = self._new_spin(value)
        self.spins.append(spin)
        self.values_box.append(spin)

        if not self.hold:
            self.emit("value-changed")

    def on_remove(self, button):
        """Handle click on the remove button.

        Removes the last spin widget from the vector, but ensures at least
        one widget remains.

        Args:
            button (Gtk.Button): The remove button that was clicked.
        """
        if len(self.spins) <= 1:
            return

        spin = self.spins.pop()

        self.values_box.remove(spin)

        if not self.hold:
            self.emit("value-changed")

    def on_stop(self, button):
        """Handle click on the stop button.

        Toggles hold mode.

        Args:
            button (Gtk.ToggleButton): The stop button that was clicked.
        """
        self.hold = button.get_active()

        if not self.hold:
            self.emit("value-changed")

    def _new_spin(self, value=None):
        """Create a new spin widget.

        Creates either a ScientificSpin or Gtk.SpinButton depending on
        itemclass, connects the 'value-changed' signal, and returns
        the widget.

        Args:
            value (float or int or None): Initial value for the spin widget.

        Returns:
            ScientificSpin or Gtk.SpinButton: The new spin widget.

        Raises:
            ValueError: If itemclass is not ScientificSpin or
                Gtk.SpinButton.

        """
        if self.itemclass == ScientificSpin:
            spin = ScientificSpin(
                value=value,
                rel_step=self.rel_step,
                with_unit=False
            )
        elif self.itemclass == Gtk.SpinButton:
            # Import ctypes for spin button range
            import ctypes
            lower = 0
            upper = (1 << (ctypes.sizeof(ctypes.c_long) * 8)) - 1

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
            # Block scroll events to avoid accidental changes
            scroll_controller = Gtk.EventControllerScroll()
            scroll_controller.connect("scroll", lambda *args: True)
            spin.add_controller(scroll_controller)
        else:
            raise ValueError(f'itemclass must be one of ScientificSpin, Gtk.SpinButton')

        spin.connect("value-changed", self.on_spin_changed)
        return spin
