"""
GTK4 Utility Functions for Gyoto

This module provides utility functions for GTK4 applications in Gyoto,
including error dialog display and other common UI tasks.
"""

__all__ = ['show_error_dialog']

import gi
gi.require_version("Gtk", "4.0")
from gi.repository import Gtk

def show_error_dialog(message="An error occurred", detail=None,
                      window=None, widget=None):
    """Display an error message in a GTK4 alert dialog.

    Creates and shows a modal dialog with the specified error message
    and optional details. This is the recommended way to display
    errors to users in Gyoto GTK4 applications.

    Args:
        message (str): Main error message to display.  Default: "An
            error occurred"
        detail (str, optional): Additional details about the error.
            Displayed in an expandable section of the dialog.
        window (Gtk.Window, optional): Parent window for the dialog.
            The dialog will be transient for this window.
        widget (Gtk.Widget, optional): Alternative to window. If
            provided and window is None, uses widget's root window.

    Notes:
        If both window and widget are None, the dialog will be
        application-modal. If detail is None, no expandable section is
        shown.

    """
    # If no window specified but widget is provided, use widget's root
    if window is None and widget is not None:
        window = widget.get_root()

    # Create alert dialog with OK button
    dialog = Gtk.AlertDialog(
        message=message,
        buttons=["OK"]
    )

    # Set detail text if provided
    if detail is not None:
        dialog.set_detail(detail)

    # Show the dialog
    dialog.show(window)
