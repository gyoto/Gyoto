__all__ = ['FilenameEditor']

import gi
gi.require_version("Gtk", "4.0")
from gi.repository import Gtk, Gdk, GLib

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

