## A SpinButton in floating-point notation with option unit field
#
# ┌────────────────────────────────────────────────┐
# │ ┌────────────────────┐ ┌───────────────┐ ┌───┐ │
# │ │ (Gtk.Entry:)       │ │ (Gtk.Entry:)  │ │ ▲ │ │
# │ │ 1.23456789e-12     │ │ m^2           │ │ ▼ │ │
# │ └────────────────────┘ └───────────────┘ └───┘ │
# └────────────────────────────────────────────────┘
#

__all__ = ['ScientificSpin']

import gi
gi.require_version("Gtk", "4.0")
from gi.repository import Gtk, Gdk

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
        self.rel_step = rel_step

        # Value entry
        self.value_entry = Gtk.Entry()
        self.value_entry.set_hexpand(True)
        self.value_entry.set_tooltip_text("value")
        if value is not None:
            self.set_value(value)
        self.append(self.value_entry)

        # Unit entry
        if with_unit:
            self.unit_entry = Gtk.Entry()
            self.unit_entry.set_hexpand(True)
            self.set_unit("")
            self.append(self.unit_entry)
            self.unit_entry.connect("activate", self.on_unit_changed)
            self.unit_entry.set_tooltip_text("unit")
        else:
            self.unit_entry = None

        # Buttons container
        buttons = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=0)

        # "-" button
        self.down = Gtk.Button()
        self.down.set_icon_name("value-decrease-symbolic")
        self.down.set_halign(Gtk.Align.FILL)
        self.down.set_valign(Gtk.Align.FILL)
        self.down.set_has_frame(True)
        self.down.add_css_class("flat")
        self.down.set_focusable(False)

        # "+" button
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

        key_controller = Gtk.EventControllerKey()
        key_controller.connect("key-pressed", self.on_key_press)
        self.value_entry.add_controller(key_controller)

        buttons.append(self.down)
        buttons.append(self.up)
        self.append(buttons)

        # Minimal CSS to ensure buttons match entry height
        self.add_css_class("scientific-spin")

        css = Gtk.CssProvider()
        css.load_from_data(b"""
            /* Make buttons flat and white like entry */
            .scientific-spin button {
                /* aspect-ratio: 1; */
                background: @theme_base_color;
                /* border: none; */
                /* box-shadow: none; */
                /* min-width: 20px; */
                /* min-height: 20px; */
                /* padding: 0; */
                /* margin: 0; */
            }
            .scientific-spin button:hover {
                background: @theme_hover_bg_color;
            }
        """)
        #     /* Thin borders between widgets */
        #     .scientific-spin > box > button:first-child {
        #         border-right: 1px solid @borders;
        #     }

        #     .scientific-spin > entry + box > button:first-child {
        #         border-left: 1px solid @borders;
        #     }

        #     .scientific-spin > entry + entry + box > button:first-child {
        #         border-left: 1px solid @borders;
        #     }
        # """)
        # css.load_from_data(b"""
        #     .scientific-spin button {
        #         padding: 0;
        #         margin: 0;
        #         min-width: 20px;
        #         min-height: 20px;
        #     }
        # """)
        Gtk.StyleContext.add_provider_for_display(
            self.get_display(),
            css,
            Gtk.STYLE_PROVIDER_PRIORITY_APPLICATION
        )

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
            if direction == +1:
                factor = 1.+self.rel_step
            elif direction == -1:
                factor = 1./(1.+self.rel_step)
            elif direction == +10:
                factor=10.
            elif direction == -10:
                factor=0.1
            else:
                raise ValueError("direction should be +1, -1, +10 or -10")
            value *= factor

        # special case: start from zero
        else:
            value = direction * 1.0

        self.set_value(value, emit=True)


    def on_key_press(self, controller, keyval, keycode, state):
        """Handle up/down arrow key presses to increment/decrement the value."""

        # Check for up/down arrow keys
        if keyval == Gdk.KEY_Up:
            self.increment(None, +1)  # Increment
            return True  # Event handled
        elif keyval == Gdk.KEY_Down:
            self.increment(None, -1)  # Decrement
            return True  # Event handled
        elif keyval == Gdk.KEY_Page_Up:
            self.increment(None, +10)  # Decrement
            return True  # Event handled
        elif keyval == Gdk.KEY_Page_Down:
            self.increment(None, -10)  # Decrement
            return True  # Event handled

        return False  # Event not handled, propagate further

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

