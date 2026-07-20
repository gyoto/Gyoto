"""Gyoto Object Editor: GTK4 Window for Property Editing

This module provides a GTK4 window for editing the properties of Gyoto
objects.  It is designed for use in interactive Python sessions rather
than as a standalone application.

Note:

    The GyotoObjectEditor.run() method is wrapped in the edit() method
    of gyoto.core.Object, allowing:
        my_object.edit()

"""

__all__ = ['GyotoObjectEditor']

import gi
gi.require_version("Gtk", "4.0")
from gi.repository import Gtk, GLib

from ..widgets.property_editor_box import PropertyEditorBox
from ...core import Factory

class GyotoObjectEditor(Gtk.Window):
    """A GTK4 window for editing Gyoto object properties.

    This window provides a scrollable view of all editable properties of a
    Gyoto object, using a PropertyEditorBox as its main content.

    The window can run in blocking mode (manages its own GLib main loop) or
    non-blocking mode (for integration with external event loops like
    IPython).

    Parameters:
        obj: The Gyoto object to edit
        blocking (bool): If True, the window manages the GLib main
            loop.  If False, the caller must manage the event loop
            (e.g., using %gui gtk4 in IPython).
        connector (multiprocessing.Connection or None): If the GUI
            runs in a separate process, this is used to send updates
            back to the caller.

    Attributes:
        obj: The Gyoto object being edited
        main_loop: GLib.MainLoop instance (if blocking=True)
        scrolled_window: Gtk.ScrolledWindow containing the editor
        vbox: PropertyEditorBox for editing object properties

    """

    @staticmethod
    def run(obj, blocking=True, connector=None):
        """Construct a GyotoObjectEditor window and run it.

        This is the recommended way to use the editor. It creates a new
        window, presents it, and optionally runs the GTK main loop.

        Args:
            obj: The Gyoto object to edit
            blocking (bool): If True, run the GLib main loop.  If
                False, return immediately after presenting.
            connector (multiprocessing.Connection or None): If the GUI
                runs in a separate process, this is used to send updates
                back to the caller.

        Returns:
            GyotoObjectEditor: The created window instance

        Example:
            GyotoObjectEditor.run(my_metric, blocking=True)
            # or
            my_metric.edit()

        """
        win = GyotoObjectEditor(obj, blocking, connector)
        win.present()
        if blocking:
            win.main_loop.run()

    def __init__(self, obj, blocking=True, connector=None):
        """Initialize the GyotoObjectEditor window.

        Args:
            obj: The Gyoto object to edit
            blocking (bool): Whether to manage the GLib main loop
            connector (multiprocessing.Connection or None): If the GUI
                runs in a separate process, this is used to send updates
                back to the caller.

        """
        Gtk.Window.__init__(self, title="Gyoto Object Editor")
        self.set_default_size(400, 600)

        # handle QUIT from parent process
        self.connector = connector
        if connector is not None:
            GLib.timeout_add(50, self.check_connector)

        # obj may be the XML description of the object
        if isinstance(obj, str):
            factory = Factory(obj)
            obj = getattr(factory, factory.kind().lower())()

        self.obj = obj
        if blocking:
            self.main_loop = GLib.MainLoop()
        else:
            self.main_loop = None

        self.connect("close-request", self.on_close_request)

        # Main container: scrollable window
        self.scrolled_window = Gtk.ScrolledWindow()
        self.set_child(self.scrolled_window)

        # Vertical box for property widgets
        self.editor = PropertyEditorBox(obj, connector=connector)
        self.editor.set_margin_top(10)
        self.editor.set_margin_bottom(10)
        self.editor.set_margin_start(10)
        self.editor.set_margin_end(10)
        self.scrolled_window.set_child(self.editor)

        # Note: We could connect to value-changed to react whenever
        # obj changes self.editor.connect("value-changed",
        # self.on_value_changed) but we don't actually need to.

    def on_close_request(self, *args):
        """Handle window close request.

        Quits the main loop if the window is in blocking mode.

        Args:
            *args: GTK callback arguments

        Returns:
            bool: False to allow the window to close
        """
        if self.main_loop is not None:
            GLib.idle_add(self.main_loop.quit)
        return False

    def check_connector(self):
        """Check for QUIT commands from the parent process.

        This method is called periodically (every 50ms) via
        GLib.timeout_add to poll the inter-process communication pipe
        for a QUIT message.  When received, it quits the GTK main
        loop, allowing the process to exit gracefully.

        Returns:
            bool: False to stop the timeout (after QUIT), True to
                continue.

        """
        if self.connector is None:
            return False
        try:
            if self.connector.poll(0):
                msg = self.connector.recv()
                if msg == ('QUIT',):
                    self.close()
                    # if self.main_loop is not None:
                    #     self.main_loop.quit()
                    return False
        except:
            pass
        return True

