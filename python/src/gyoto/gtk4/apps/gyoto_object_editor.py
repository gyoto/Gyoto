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

from gettext import gettext as _

from ..widgets.property_editor_box import PropertyEditorBox
from ...core import Factory, Scenery, Error as GyotoError, Photon, Screen
from ... import core

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
            kwargs['application_id'] = "fr.obspm.gyoto.GyotoObjectEditor"
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

        # Build UI
        self.build_headerbar()
        self.build_shortcuts()

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

        # Register with application
        if application:
            application.windows.append(self)

    ####################################################################
    # UI
    ####################################################################

    def build_headerbar(self):
        """Build the window's header bar with menu button.

        Creates a header bar with a hamburger menu containing:
        - New...
          - Star
          - Photon
        - Open...
        - Save
        - Save As...
        - Help
        - Close
        - Quit

        """
        # Create title bar with hamburger button
        header = Gtk.HeaderBar()
        self.set_titlebar(header)

        # Create menu model
        menu = Gio.Menu()

        # New submenu
        new_menu = Gio.Menu()

        # Dynamically add items
        dynamic_items = {
            'core': ["Scenery", "Screen", "Photon"]
        }
        for section in 'Astrobj', 'Metric', 'Spectrum', 'Spectrometer':
            generic = getattr(core, section)
            kinds = list(generic.registeredPluginsSlashKinds())
            kinds.sort()
            dynamic_items[section] = kinds

        for section, items in dynamic_items.items():
            new_menu_section = Gio.Menu()
            for label in items:
                # Create a menu item with the label
                item = Gio.MenuItem.new(label, "win.new")
                item.set_action_and_target_value(
                    "win.new", GLib.Variant("s", section + '/' + label)
                )
                new_menu_section.append_item(item)
            new_menu.append_submenu(section, new_menu_section)

        menu_section1 = Gio.Menu()
        menu_section1.append_submenu(_("New"), new_menu)
        menu_section2 = Gio.Menu()
        menu_section2.append(_("Open…"), "win.open")
        menu_section2.append(_("Save"), "win.save")
        menu_section2.append(_("Save As…"), "win.save-as")
        menu_section3 = Gio.Menu()
        menu_section3.append(_("Help"), "win.help")
        menu_section3.append(_("Close"), "win.close")
        menu_section3.append(_("Quit"), "win.quit")

        # Main menu items
        menu.append_section(None, menu_section1)
        menu.append_section(None, menu_section2)
        menu.append_section(None, menu_section3)

        # Create menu button
        menu_button = Gtk.MenuButton(
            icon_name="open-menu-symbolic",
            menu_model=menu,
            use_underline=True
        )
        menu_button.add_css_class("flat")
        header.pack_end(menu_button)

        # Connect actions
        action_group = Gio.SimpleActionGroup()
        self.insert_action_group("win", action_group)

        action_new = Gio.SimpleAction.new("new", GLib.VariantType("s"))
        action_new.connect("activate", self.on_new_item_activated)
        action_group.add_action(action_new)

        action_group.add_action_entries([
            ("open", self.on_open, None),
            ("save", self.on_save, None),
            ("save-as", self.on_save_as, None),
            ("help", self.on_help, None),
            ("close", self.on_close, None),
            ("quit", self.on_quit, None),
        ])

    def build_shortcuts(self):
        """Create keyboard shortcuts.

        Creates keyboard shortcuts for these actions:
        - Close window: Ctrl+W,
        - Close all windows and quit: Ctrl+Q,
        - Open file: Ctrl+O,
        - Save file: Ctrl+S,
        - Save file as: Ctrl+Shift+S,
        - Help: F1,

        """
        controller = Gtk.ShortcutController()
        self.add_controller(controller)

        controller.add_shortcut(
            Gtk.Shortcut(
                trigger=Gtk.KeyvalTrigger(
                    keyval=Gdk.KEY_o,
                    modifiers=Gdk.ModifierType.CONTROL_MASK
                ),
                action=Gtk.NamedAction.new("win.open")
            )
        )

        controller.add_shortcut(
            Gtk.Shortcut(
                trigger=Gtk.KeyvalTrigger(
                    keyval=Gdk.KEY_s,
                    modifiers=Gdk.ModifierType.CONTROL_MASK
                ),
                action=Gtk.NamedAction.new("win.save")
            )
        )

        controller.add_shortcut(
            Gtk.Shortcut(
                trigger=Gtk.KeyvalTrigger(
                    keyval=Gdk.KEY_s,
                    modifiers=(
                        Gdk.ModifierType.CONTROL_MASK
                        | Gdk.ModifierType.SHIFT_MASK
                    )
                ),
                action=Gtk.NamedAction.new("win.save-as")
            )
        )

        controller.add_shortcut(
            Gtk.Shortcut(
                trigger=Gtk.KeyvalTrigger(
                    keyval=Gdk.KEY_F1,
                    modifiers=0
                ),
                action=Gtk.NamedAction.new("win.help")
            )
        )

        controller.add_shortcut(
            Gtk.Shortcut(
                trigger=Gtk.KeyvalTrigger(
                    keyval=Gdk.KEY_w,
                    modifiers=Gdk.ModifierType.CONTROL_MASK
                ),
                action=Gtk.NamedAction.new("win.close")
            )
        )

        controller.add_shortcut(
            Gtk.Shortcut(
                trigger=Gtk.KeyvalTrigger(
                    keyval=Gdk.KEY_q,
                    modifiers=Gdk.ModifierType.CONTROL_MASK
                ),
                action=Gtk.NamedAction.new("win.quit")
            )
        )

    ####################################################################
    # Callbacks
    ####################################################################

    def on_new_item_activated(self, action, parameter):
        """Create a new object and open an editor window

        Args:
            action: the Gio Action
            parameter: its parameter

        """
        obj = None
        label = parameter.get_string()  # Extract the label

        static_constructors = {
            "core/Scenery": Scenery,
            "core/Photon": Photon,
            "core/Screen": Screen,
        }

        if label in static_constructors:
            obj = static_constructors[label]()
        else:
            nspace, plg, kind = label.split('/', 2)
            obj = getattr(core, nspace)(kind, (plg,))

        if obj is None:
            raise ValueError(f'no valid constructor for {label}')

        window = GyotoObjectEditorApplicationWindow(
            application=self.props.application,
            obj=obj,
            connector=None
        )
        window.present()

    def on_open(self, *args):
        """Open a file dialog to load an object from an XML file.

        Creates a file dialog with XML and All Files filters,
        defaulting to XML. When a file is selected, it's loaded and
        the object is openned in a new editor window.

        Args:
            *args: GTK callback arguments

        """
        dialog = Gtk.FileDialog()

        # Create XML filter
        xml_filter = Gtk.FileFilter()
        xml_filter.set_name("XML files")
        xml_filter.add_suffix('xml')
        xml_filter.add_pattern('*.xml')

        # Create All Files filter
        all_filter = Gtk.FileFilter()
        all_filter.set_name("All files")
        all_filter.add_pattern('*')

        # Create list model for filters
        filter_list = Gio.ListStore.new(Gtk.FileFilter)
        filter_list.append(xml_filter)
        filter_list.append(all_filter)

        # Set filters and default
        dialog.set_filters(filter_list)
        dialog.set_property('default-filter', xml_filter)

        dialog.open(
            self,
            None,
            lambda dialog, result: self.on_open_file_selected(dialog, result)
        )

    def on_open_file_selected(self, dialog, result):
        """Handle the selection of a file to open.

        Args:
            dialog: The Gtk.FileDialog instance
            result: The result of the dialog operation

        """
        try:
            file = dialog.open_finish(result)
        except GLib.Error:
            return

        if file is not None:
            window = GyotoObjectEditorApplicationWindow(
                application=self.props.application,
                obj=file.get_path(),
                connector=None
            )
            window.present()

    def on_save(self, *args):
        """Save the object in the last XML file used.

        If the object was not read from file and not yet saved to
        file, opens a dialog.

        Args:
            *args: GTK callback arguments

        """
        if self.filename:
            try:
                Factory(self.obj).write(self.filename)
            except GyotoError as e:
                show_error_dialog(
                    message=f"Error writing XML file {self.filename}:",
                    detail=e.get_message(),
                    window=self
                )
        else:
            self.on_save_as(*args)

    def on_save_as(self, *args):
        """Open a file dialog to save the object as XML.

        Creates a file dialog with XML and All Files filters,
        defaulting to XML. When a file is selected, the object is
        saved to that file.

        Args:
            *args: GTK callback arguments

        """
        dialog = Gtk.FileDialog()

        # Create XML filter
        xml_filter = Gtk.FileFilter()
        xml_filter.set_name("XML files")
        xml_filter.add_suffix('xml')
        xml_filter.add_pattern('*.xml')

        # Create All Files filter
        all_filter = Gtk.FileFilter()
        all_filter.set_name("All files")
        all_filter.add_pattern('*')

        # Create list model for filters
        filter_list = Gio.ListStore.new(Gtk.FileFilter)
        filter_list.append(xml_filter)
        filter_list.append(all_filter)

        # Set filters and default
        dialog.set_filters(filter_list)
        dialog.set_property('default-filter', xml_filter)

        dialog.save(
            self,
            None,
            lambda dialog, result: self.on_save_file_selected(dialog, result)
        )

    def on_save_file_selected(self, dialog, result):
        """Handle the selection of a file to save to.

        Args:
            dialog: The Gtk.FileDialog instance
            result: The result of the dialog operation

        """
        try:
            file = dialog.save_finish(result)
        except GLib.Error:
            return

        if file is None:
            return

        try:
            Factory(self.obj).write(file.get_path())
            self.filename = file.get_path()
        except GyotoError as e:
            show_error_dialog(
                message=f"Error writing XML file {file.get_path()}:",
                detail=e.get_message(),
                window=self
            )

    def on_help(self, *args):
        """Display the help dialog.

        Shows a dialog with comprehensive usage information including
        keyboard shortcuts and menu button descriptions.

        Args:
            *args: GTK callback arguments

        """
        dialog = Gtk.Dialog(
            title="Help",
            transient_for=self,
            modal=False
        )
        dialog.set_default_size(600, 400)

        help_text = (
            "The Gyoto Object Editor\n"
        )

        label = Gtk.Label(
            label=help_text,
            halign=Gtk.Align.START,
            wrap=True,
            xalign=0.0
        )

        scrolled = Gtk.ScrolledWindow()
        scrolled.set_child(label)
        scrolled.set_policy(
            Gtk.PolicyType.AUTOMATIC,
            Gtk.PolicyType.AUTOMATIC
        )

        dialog.set_child(scrolled)
        dialog.present()

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

    def on_close(self, *args):
        """Handle close action from menu.

        Args:
            *args: GTK callback arguments

        """
        self.close()

    def on_quit(self, *args):
        """Handle quit all action.

        Closes all windows and quits the application.

        Args:
            *args: GTK callback arguments

        """
        if self.props.application is not None:
            self.props.application.close_all_windows()
        else:
            self.close()

    def default_obj(self):
        """Return default object
        """
        return Scenery()

# Stand-alone entry point:
if __name__ == "__main__":
    raise SystemExit(
        GyotoObjectEditorApplication.run_app(parsecliargs=True)
    )
