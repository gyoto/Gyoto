import gi

gi.require_version("Gtk", "4.0")
from gi.repository import Gdk, Gtk, GObject

import matplotlib.colors


class MplColorButton(Gtk.ColorDialogButton):
    """A Gtk.ColorDialogButton using Matplotlib color representations."""

    __gsignals__ = {
        "color-changed": (
            GObject.SignalFlags.RUN_FIRST,
            None,
            (object,),
        ),
    }

    def __init__(self, color="black", **kwargs):
        dialog = Gtk.ColorDialog()
        super().__init__(dialog=dialog, **kwargs)

        self.color = color

        self.connect(
            "notify::rgba",
            self._on_rgba_changed,
        )

    @property
    def color(self):
        """The current color as a Matplotlib RGBA tuple."""
        return self._color

    @color.setter
    def color(self, color):
        rgba = matplotlib.colors.to_rgba(color)
        self._color = rgba
        self.set_rgba(self._to_gdk(rgba))

    def _on_rgba_changed(self, button, pspec):
        self._color = self._from_gdk(self.get_rgba())
        self.emit("color-changed", self._color)

    @staticmethod
    def _to_gdk(rgba):
        color = Gdk.RGBA()
        color.red = rgba[0]
        color.green = rgba[1]
        color.blue = rgba[2]
        color.alpha = rgba[3]
        return color

    @staticmethod
    def _from_gdk(rgba):
        return (
            rgba.red,
            rgba.green,
            rgba.blue,
            rgba.alpha,
        )
