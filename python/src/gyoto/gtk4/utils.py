__all__ = ['show_error_dialog']

import gi
gi.require_version("Gtk", "4.0")
from gi.repository import Gtk

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
