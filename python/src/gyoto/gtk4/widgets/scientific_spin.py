"""ScientificSpin: A GTK4 SpinButton with Scientific Notation Support

This module provides a custom GTK4 widget that combines the
functionality of a SpinButton with scientific notation support and
optional unit display.

Widget Layout
------------
    ┌────────────────────────────────────────────────┐
    │ ┌────────────────────┐ ┌───────────────┐ ┌───┐ │
    │ │ (Gtk.Entry:)       │ │ (Gtk.Entry:)  │ │ ▲ │ │
    │ │ 1.23456789e-12     │ │ m^2           │ │ ▼ │ │
    │ └────────────────────┘ └───────────────┘ └───┘ │
    └────────────────────────────────────────────────┘

Description
-----------
- **Value Entry**: Accepts scientific notation (e.g., 1.3e-12)
- **Unit Entry**: Optional field for displaying units (e.g., m^2)
- **Buttons**: Increment/decrement buttons with relative stepping
- **Keyboard Support**: Up/Down arrows and PageUp/PageDown for value
    adjustment

Public API
---------
- **Signals**: `value-changed`, `unit-changed`
- **Methods**: `get_value()`, `set_value()`, `get_unit()`, `set_unit()`

"""

__all__ = ['ScientificSpin']

import gi
gi.require_version("Gtk", "4.0")
from gi.repository import Gtk, Gdk

class ScientificSpin(Gtk.Box):
    """A custom GTK4 widget for numeric input with scientific notation.

    This widget provides a spin-button-like interface with the
    following features:
    - Scientific notation support (e.g., 1.23e-10)
    - Optional unit field for displaying measurement units
    - Increment/decrement buttons with relative stepping
    - Keyboard navigation (arrow keys, PageUp/PageDown)

    The widget is composed of:
        - A primary Gtk.Entry for the numeric value
        - An optional Gtk.Entry for the unit string
        - Up/Down buttons for value adjustment

    Attributes:
        value_entry (Gtk.Entry): The entry widget for numeric values
        unit_entry (Gtk.Entry or None): The entry widget for units (if
            enabled)
        rel_step (float): Relative step size for increment/decrement
            (default: 0.1)
        _updating (bool): Internal flag to prevent signal emission
            during updates

    Signals:
        value-changed: Emitted when the numeric value changes
        unit-changed: Emitted when the unit changes (if unit field is
            enabled)

    Example:
        spin = ScientificSpin(value=1.0, rel_step=0.1, with_unit=True)
        spin.connect("value-changed", lambda w: print(w.get_value()))

    """

    __gsignals__ = {
        "value-changed": (gi.repository.GObject.SignalFlags.RUN_FIRST, None, ()),
        "unit-changed": (gi.repository.GObject.SignalFlags.RUN_FIRST, None, ())
    }

    def __init__(self, value=1.0, rel_step=0.1, with_unit=False):
        """Initialize the ScientificSpin widget.

        Args:
            value (float): Initial numeric value (default: 1.0)
            rel_step (float): Relative step size for
              increment/decrement (default: 0.1).  E.g., 0.1 means
              each click multiplies/divides by 1.1/0.909...
            with_unit (bool): Whether to show a unit field (default: False)

        """
        super().__init__(orientation=Gtk.Orientation.HORIZONTAL,
                         spacing=0)

        # Flag to prevent signal emission during programmatic updates
        self._updating = False
        self.rel_step = rel_step

        # --- Value Entry ---
        self.value_entry = Gtk.Entry()
        self.value_entry.set_hexpand(True)
        self.value_entry.set_tooltip_text("Scientific notation value (e.g., 1.23e-10)")
        if value is not None:
            self.set_value(value)
        self.append(self.value_entry)

        # --- Unit Entry (optional) ---
        if with_unit:
            self.unit_entry = Gtk.Entry()
            self.unit_entry.set_hexpand(True)
            self.set_unit("")
            self.append(self.unit_entry)
            self.unit_entry.connect("activate", self.on_unit_changed)
            self.unit_entry.set_tooltip_text("Unit of measurement (e.g., m, s, kg)")
        else:
            self.unit_entry = None

        # --- Buttons Container ---
        buttons = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=0)

        # "-" Button
        self.down = Gtk.Button()
        self.down.set_icon_name("value-decrease-symbolic")
        self.down.set_halign(Gtk.Align.FILL)
        self.down.set_valign(Gtk.Align.FILL)
        self.down.set_has_frame(True)
        self.down.add_css_class("flat")
        self.down.set_focusable(False)

        # "+" Button
        self.up = Gtk.Button()
        self.up.set_icon_name("value-increase-symbolic")
        self.up.set_halign(Gtk.Align.FILL)
        self.up.set_valign(Gtk.Align.FILL)
        self.up.set_has_frame(True)
        self.up.add_css_class("flat")
        self.up.set_focusable(False)

        # Connect signals
        self.value_entry.connect("changed", self.on_value_changed)
        self.up.connect("clicked", self.increment, +1)
        self.down.connect("clicked", self.increment, -1)

        # Keyboard support
        key_controller = Gtk.EventControllerKey()
        key_controller.connect("key-pressed", self.on_key_press)
        self.value_entry.add_controller(key_controller)

        buttons.append(self.down)
        buttons.append(self.up)
        self.append(buttons)

        # --- Styling ---
        self.add_css_class("scientific-spin")

        # Apply CSS to match entry appearance
        css = Gtk.CssProvider()
        css.load_from_data(b"""
            /* Make buttons match entry appearance */
            .scientific-spin button {
                background: @theme_base_color;
            }
            .scientific-spin button:hover {
                background: @theme_hover_bg_color;
            }
        """)
        Gtk.StyleContext.add_provider_for_display(
            self.get_display(),
            css,
            Gtk.STYLE_PROVIDER_PRIORITY_APPLICATION
        )

    def get_value(self):
        """Get the current numeric value from the entry.

        Returns:
            float: The current value as a floating-point number

        Raises:
            ValueError: If the entry text cannot be converted to float
        """
        return float(self.value_entry.get_text())

    def set_value(self, value, *, emit=False):
        """Set the value field to a new numeric value.

        Args:
            value (float): The new value to display
            emit (bool): If True, emit the 'value-changed' signal
                after update. If False (default), suppress signal
                emission.

        Note:
            This method uses repr() to format the value, which preserves
            scientific notation for very large/small numbers.

        """
        was_updating = self._updating
        self._updating = was_updating or not emit
        try:
            self.value_entry.set_text(repr(value))
        finally:
            self._updating = was_updating

    def get_unit(self):
        """Get the current unit text.

        Returns:
            str or None: The current unit text, or None if unit field
                is disabled

        """
        if self.unit_entry is None:
            return None
        return self.unit_entry.get_text()

    def set_unit(self, unit):
        """Set the unit field to a new string.

        Args:
            unit (str): The unit text to display (e.g., "m", "s", "kg")
        """
        if self.unit_entry is not None:
            self.unit_entry.set_text(unit)

    def increment(self, button, direction):
        """Increment or decrement the value by a relative step.

        Args:
            button (Gtk.Button): The button that triggered the
                increment (unused)
            direction (int): Direction of change:
                +1: Increment by relative step
                -1: Decrement by relative step
                +10: Multiply by 10 (large increment)
                -10: Divide by 10 (large decrement)

        Raises:
            ValueError: If direction is not +1, -1, +10, or -10

        """
        try:
            value = self.get_value()
        except ValueError:
            value = 0.0

        # Normal case: relative increment (multiplicative)
        if value != 0:
            if direction == +1:
                factor = 1. + self.rel_step
            elif direction == -1:
                factor = 1. / (1. + self.rel_step)
            elif direction == +10:
                factor = 10.
            elif direction == -10:
                factor = 0.1
            else:
                raise ValueError("direction should be +1, -1, +10 or -10")
            value *= factor

        # Special case: starting from zero
        else:
            value = direction * 1.0

        self.set_value(value, emit=True)

    def on_key_press(self, controller, keyval, keycode, state):
        """Handle keyboard events for value adjustment.

        Args:
            controller (Gtk.EventControllerKey): The key controller
            keyval (Gdk.key): The key that was pressed
            keycode (int): Hardware keycode (unused)
            state (Gdk.ModifierType): Modifier state (unused)

        Returns:
            bool: True if the event was handled, False otherwise
        """
        # Check for up/down arrow keys
        if keyval == Gdk.KEY_Up:
            self.increment(None, +1)  # Increment
            return True
        elif keyval == Gdk.KEY_Down:
            self.increment(None, -1)  # Decrement
            return True
        elif keyval == Gdk.KEY_Page_Up:
            self.increment(None, +10)  # Large increment
            return True
        elif keyval == Gdk.KEY_Page_Down:
            self.increment(None, -10)  # Large decrement
            return True

        return False  # Event not handled, propagate further

    def on_value_changed(self, entry):
        """Handle value entry changes and emit value-changed signal.

        This method is connected to the 'changed' signal of
        value_entry.  It validates the input and emits the
        'value-changed' signal if valid.

        Args:
            entry (Gtk.Entry): The entry widget that changed

        """
        if self._updating:
            return
        try:
            value = float(entry.get_text())
        except ValueError:
            return  # Invalid input, don't emit

        self.emit("value-changed")

    def on_unit_changed(self, entry):
        """Handle unit entry changes and emit unit-changed signal.

        This method is connected to the 'activate' signal of unit_entry.
        It emits the 'unit-changed' signal when the unit is changed.

        Args:
            entry (Gtk.Entry): The unit entry widget that changed
        """
        self.emit("unit-changed")
