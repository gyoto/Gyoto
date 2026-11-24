"""FilenameEditor: GTK4 Widget for File Selection

This module provides a compound GTK4 widget for selecting files,
combining a text entry with a file browser button.

Widget Layout
------------
    ┌───────────────────────────────┐
    │ /home/user/file.dat        📂 │
    └───────────────────────────────┘

Description
-----------
- **Entry Field**: Displays and allows editing of the file path
- **Browse Button**: Opens a file dialog to select a file
- **Signals**: Emits 'value-changed' when the path changes (via entry
    or dialog)

Example:
    editor = FilenameEditor(value="/path/to/file.txt")
    editor.connect("value-changed", lambda w: print(w.get_value()))

"""

__all__ = ['FilenameEditor']

import gi
gi.require_version("Gtk", "4.0")
from gi.repository import Gtk, Gdk, GLib

class FilenameEditor(Gtk.Box):
    """A compound widget for selecting files.

    This widget combines a Gtk.Entry for displaying/editing file paths
    with a Gtk.Button that opens a file chooser dialog. It provides a
    compact and intuitive way to select files in a GTK4 application.

    The widget is composed of:
        - A Gtk.Entry for the file path (expands horizontally)
        - A Gtk.Button with a folder icon for browsing

    Attributes:
        entry (Gtk.Entry): The text entry widget for the file path
        button (Gtk.Button): The browse button with folder icon

    Signals:
        value-changed: Emitted when the file path changes, either by:
            - User pressing Enter in the entry field
            - User selecting a file from the dialog

    """

    __gsignals__ = {
        "value-changed": (
            gi.repository.GObject.SignalFlags.RUN_FIRST,
            None,
            ()
        )
    }

    def __init__(self, value=""):
        """Initialize the FilenameEditor widget.

        Args:
            value (str): Initial file path to display (default: "")
        """
        super().__init__(
            orientation=Gtk.Orientation.HORIZONTAL,
            spacing=4
        )

        # Entry field for file path
        self.entry = Gtk.Entry()
        self.entry.set_hexpand(True)  # Expand horizontally to fill space
        self.entry.set_text(value)
        self.append(self.entry)

        # Browse button with folder icon
        self.button = Gtk.Button()
        self.button.set_icon_name("folder-open-symbolic")
        self.button.set_tooltip_text("Choose a file")
        self.button.add_css_class("flat")  # Flat styling for integration
        self.append(self.button)

        # Connect signals
        self.entry.connect(
            "activate",  # Emitted when user presses Enter
            self.on_activate
        )
        self.button.connect(
            "clicked",
            self.on_browse
        )

    def get_value(self):
        """Get the current file path from the entry.

        Returns:
            str: The current file path text
        """
        return self.entry.get_text()

    def set_value(self, value):
        """Set the file path in the entry.

        Args:
            value (str): The file path to display
        """
        self.entry.set_text(value)

    def on_activate(self, entry):
        """Handle Enter key press in the entry field.

        Emits the 'value-changed' signal when the user confirms the path.

        Args:
            entry (Gtk.Entry): The entry widget that received the
                activate signal

        """
        self.emit("value-changed")

    def on_browse(self, button):
        """Handle click on the browse button.

        Opens a Gtk.FileDialog to select a file.

        Args:
            button (Gtk.Button): The button that was clicked
        """
        dialog = Gtk.FileDialog()
        dialog.open(
            self.get_root(),  # Parent window
            None,  # Cancellation token (None = no cancellation)
            lambda dialog, result: self.on_file_selected(dialog, result)
        )

    def on_file_selected(self, dialog, result):
        """Handle file selection from the dialog.

        Updates the entry with the selected file path and emits
        'value-changed'.

        Args:
            dialog (Gtk.FileDialog): The file dialog that completed
            result (Gio.AsyncResult): The result of the dialog operation

        """
        try:
            file = dialog.open_finish(result)
        except GLib.Error:
            return  # Dialog was cancelled or error occurred

        if file is not None:
            self.set_value(file.get_path())
            self.emit("value-changed")
