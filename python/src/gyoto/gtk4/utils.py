"""
GTK4 Utility Functions for Gyoto

This module provides utility functions for GTK4 applications in Gyoto,
including error dialog display and other common UI tasks.
"""

import warnings
import traceback
from ..core import Property, Factory, Error as GyotoError, debug

__all__ = ['show_error_dialog']

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
    # lazy import as we don't want import utils to import Gtk
    import gi
    gi.require_version("Gtk", "4.0")
    from gi.repository import Gtk

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


def gui_launcher(gui_function, handler_function, *args):
    """Launch a GTK GUI in a separate process with IPC to the caller.

    This function starts a GTK application in a subprocess and sets up
    inter-process communication (IPC) so that the GUI can send events
    back to the caller (e.g., an IPython session).

    Args:
        gui_function: Function to run in the GUI process. Receives a
            connection object for sending events back and any *args.
        handler_function: Function to handle events in the caller
            process. Receives a connection object for receiving events,
            a cleanup function, and any *args. If None, no IPC is set up.
        *args: Additional arguments passed to both gui_function and
            handler_function.

    """
    from multiprocessing import Process, Pipe
    from threading import Thread
    import atexit
    from ..core import Object

    # IPC pipe (parent = IPython, child = GTK process)
    if handler_function is None:
        snd_conn=None
    else:
        rcv_conn, snd_conn = Pipe()

    # construct p here, we need it for cleanup()
    p = Process(target=gui_function, args=(snd_conn, *args))

    def cleanup():
        try:
            if p.is_alive():
                p.terminate()
                p.join(timeout=1)
                if p.is_alive():
                    p.kill()
        except ValueError:
            pass  # Already dead

    atexit.register(cleanup)  # Clean up on main process exit

    if handler_function is not None:
        Thread(target=handler_function,
               args=(rcv_conn, cleanup, *args),
               daemon=True).start()

    # start p after the handler
    p.start()

def recursive_value_changed_pipe_receiver(connector, cleanup, obj):
    """Receive recursive-value-changed events from a GUI process.

    This function runs in the caller process and listens for events
    from a GTK GUI running in a separate process. When it receives an
    'update' event, it updates the specified property of obj.

    Args:
        connector: A multiprocessing.Connection for receiving events
            from the GUI process.
        cleanup: A function to call when the GUI process terminates or
            this receiver exits.
        obj: The Gyoto object to update when events are received.

    """
    try:
        while True:
            if connector.poll(0.01):  # non-blocking check
                try:
                    event = connector.recv()
                    print(f"event received: {event}")  # or call a callback
                    cmd = event[0]
                    if cmd == 'update':
                        ppath = event[1]
                        value = event[2]
                        print(f'updating {ppath} to {value}') 
                        descendents = ppath.split('.')
                        subobj = obj
                        pname = descendents[0]
                        print (pname)
                        for i in range(1, len(descendents)):
                            subobj = subobj.get(pname)
                            pname = descendents[i]
                            print (pname)
                        prop = subobj.property(pname)
                        if prop.type in (Property.astrobj_t, Property.metric_t,
                                         Property.screen_t, Property.spectrometer_t,
                                         Property.spectrum_t):
                            factory = Factory(value)
                            value = getattr(factory, factory.kind().lower())()
                        print(f'setting {pname} to {value}')
                        subobj.set(pname, value)
                except GyotoError as e:
                    warnings.warn(f"Gyoto object editor triggered a Gyoto error:\n{e.get_message()}")
                except EOFError:
                    # GUI process terminated, exit gracefully
                    break
    finally:
        if debug():
            warnings.warn('exiting')
        if cleanup is not None:
            cleanup()
        connector.close()
