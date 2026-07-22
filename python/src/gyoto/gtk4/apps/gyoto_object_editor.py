"""Gyoto Object Editor: GTK4 Window for Property Editing

This module provides a GTK4 window for editing the properties of
Gyoto objects.

Note:
    The GyotoObjectEditorApplication.run_app() method is wrapped in
    the edit() method of gyoto.core.Object, allowing:
        my_object.edit()

"""

__all__ = ['GyotoObjectEditorApplication',
           'GyotoObjectEditorApplicationWindow']

import gi
gi.require_version("Gtk", "4.0")
from gi.repository import Gtk, GLib, Gio, Gdk

import argparse
import sys
import traceback

from ..widgets.property_editor_box import PropertyEditorBox
from ...core import Factory, Scenery

class GyotoObjectEditorApplication(Gtk.Application):
    """Standalone GTK application for the Gyoto Object Editor.

    This class handles the application lifecycle and window
    management.  It extends Gtk.Application to provide a proper GTK
    application structure.

    """

    def __init__(self, obj=None, connector=None, *args, **kwargs):
        """Initialize the Gyoto Object Editor GTK application.

        Args:
            obj: Initial Gyoto object to edit
            connector: Connection for inter-process communication

        """
        if 'application_id' not in kwargs:
            kwargs['application_id'] = ""
            "fr.obspm.gyoto.GyotoObjectEditor"
        if 'flags' not in kwargs:
            kwargs['flags'] = (Gio.ApplicationFlags.DEFAULT_FLAGS |
                               Gio.ApplicationFlags.NON_UNIQUE)
        super().__init__(*args, **kwargs)
        self.obj = obj
        self.windows = []

        # handle QUIT from parent process
        self.connector = connector
        if connector is not None:
            GLib.timeout_add(50, self.check_connector)

    def do_activate(self):
        """Called by GTK when the application starts.

        Creates the main window if it doesn't exist and presents it.
        """
        if not self.windows:
            window = GyotoObjectEditorApplicationWindow(
                application=self,
                obj=self.obj,
                connector=self.connector
            )
            self.windows.append(window)
        for window in self.windows:
            window.present()

    @staticmethod
    def run_app(obj=None, parsecliargs=False, *args, **kwargs):
        """Run the Gyoto Object Editor

        Parameters:
            obj: the Gyoto object to start with, or None, or the XML
                description of such an object, or the name of an XML
                file containing this description.
            parsecliargs (bool): whether to parse the command line
                arguments
            *args, **kwargs: other parameters are passed untouched
                to the GyotoApplication constructor.
        
        Returns:
            int: Application exit code

        """
        remaining = None
        if parsecliargs:
            parser = argparse.ArgumentParser(
                prog=f'{sys.argv[0]} ',
                description=__doc__,
                formatter_class=argparse.RawTextHelpFormatter
            )
            parser.add_argument(
                'xmlfile', nargs='?',
                help='XML file containing the description '
                'of a Gyoto object (optional)')
            cliargs, remaining = parser.parse_known_args()
            if 'xmlfile' in cliargs:
                obj=cliargs.xmlfile

        app = GyotoObjectEditorApplication(obj=obj, *args, **kwargs)
        return app.run(remaining)

    def remove_window(self, window):
        """Remove a window from the application's window list.

        Args:
            window: The GyotoyApplicationWindow to remove.

        """
        if window in self.windows:
            self.windows.remove(window)

    def close_all_windows(self):
        """Close all open windows and quit the application."""
        for window in self.windows[:]:
            window.close()

    def check_connector(self):
        """Check for QUIT commands from the parent process.

        This method is called periodically (every 50ms) via
        GLib.timeout_add to poll the inter-process communication
        pipe for a QUIT message.  When received, it closes all
        windows , allowing the process to exit gracefully.

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
                    self.close_all_windows()
                    return False
        except:
            traceback.print_exc()
        return True

class GyotoObjectEditorApplicationWindow(Gtk.Window):
    """A GTK4 window for editing Gyoto object properties.

    This window provides a scrollable view of all editable
    properties of a Gyoto object, using a PropertyEditorBox as its
    main content.

    Parameters:
        obj: The Gyoto object to edit
        connector (multiprocessing.Connection or None): If the GUI
            runs in a separate process, this is used to send updates
            back to the caller.

    Attributes:
        obj: The Gyoto object being edited
        scrolled_window: Gtk.ScrolledWindow containing the editor
        vbox: PropertyEditorBox for editing object properties
        filename: name of last file used or None

    """

    filename = None

    def __init__(self, application=None, obj=None, connector=None):
        """Initialize the GyotoObjectEditorApplicationWindow window.

        Args:
            application: Parent Gtk.Application instance
            obj: The Gyoto object to edit
        
            connector (multiprocessing.Connection or None): If the
                GUI runs in a separate process, this is used to send
                updates back to the caller.

        """
        super().__init__(application=application,
                         title="Gyoto Object Editor")
        self.set_default_size(400, 600)

        # obj may be the XML description of the object
        if isinstance(obj, str):
            if obj.lower().endswith('.xml'):
                self.filename = obj
            factory = Factory(obj)
            obj = getattr(factory, factory.kind().lower())()

        if obj is None:
            obj = self.default_obj()

        self.obj = obj

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

        Removes the window from the application's window list, and
        allows the window to close.

        Args:
            *args: GTK callback arguments

        Returns:
            bool: False to allow the window to close

        """
        if self.props.application is not None:
            self.props.application.remove_window(self)
        return False

    def default_obj(self):
        """Return default object
        """
        return Scenery()

# Stand-alone entry point:
if __name__ == "__main__":
    raise SystemExit(
        GyotoObjectEditorApplication.run_app(parsecliargs=True)
    )
