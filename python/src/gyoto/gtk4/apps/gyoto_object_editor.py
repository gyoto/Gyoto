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
from gi.repository import Gtk, GLib, Gio, Gdk, GObject, Pango

import os.path
import argparse
import sys
import traceback

from gettext import gettext as _

from ..widgets.property_editor_box import PropertyEditorBox
from ...core import Factory, Scenery, Error as GyotoError, Photon, Screen
from ... import core

# Custom GObject class to hold the label of a 'New >' menu item
class NewItem(GObject.Object):
    """Custom GObject class to store a label for the 'New' menu items.

    Attributes:
        label (str): The label for the menu item, typically in the format
            "section/kind" (e.g., "Astrobj/Star").
    """
    __gtype_name__ = 'NewItem'

    label = GObject.Property(type=str)

    def __init__(self, label=None):
        super().__init__()
        self.label = label


class StartupWindow(Gtk.ApplicationWindow):
    """Startup window with Open and New buttons.

    This window is shown when the application starts without an object
    to edit. It provides buttons to open an existing file or create a
    new object.

    """

    def __init__(self, application):
        """Initialize the startup window.

        Args:
            application: The parent GyotoObjectEditorApplication

        """
        super().__init__(
            application=application,
            title="Gyoto Object Editor"
        )
        self.app = application

        self.build_ui()
        self.setup_actions()
        self.setup_keyboard()

    def build_ui(self):
        """Build the startup window UI."""
        main_box = Gtk.Box(
            orientation=Gtk.Orientation.VERTICAL,
            spacing=12
        )
        main_box.set_margin_top(12)
        main_box.set_margin_bottom(12)
        main_box.set_margin_start(12)
        main_box.set_margin_end(12)
        self.set_child(main_box)

        button_box = Gtk.Box(
            orientation=Gtk.Orientation.HORIZONTAL,
            spacing=6
        )
        button_box.set_halign(Gtk.Align.CENTER)

        open_button = Gtk.Button(label="Open")
        open_button.set_action_name("win.open")
        open_button.add_css_class("flat")
        button_box.append(open_button)

        self.new_button = Gtk.Button(label="New")
        self.new_button.set_action_name("win.new")
        self.new_button.add_css_class("flat")
        button_box.append(self.new_button)

        main_box.append(button_box)

    def setup_actions(self):
        """Set up actions for this window."""
        close_action = Gio.SimpleAction.new("close", None)
        close_action.connect("activate", lambda a, p: self.close())
        self.add_action(close_action)

        open_action = Gio.SimpleAction.new("open", None)
        open_action.connect("activate", lambda a, p: self.app.on_open(self))
        self.add_action(open_action)

        new_action = Gio.SimpleAction.new("new", None)
        new_action.connect("activate",
                           lambda a, p: self.app.on_new(self.new_button))
        self.add_action(new_action)

    def setup_keyboard(self):
        """Set up keyboard handler for Escape key."""
        key_controller = Gtk.EventControllerKey()
        key_controller.connect("key-pressed", self.on_key_pressed)
        self.add_controller(key_controller)

    def on_key_pressed(self, controller, keyval, keycode, state):
        """Handle Escape key to close and quit."""
        if keyval == Gdk.KEY_Escape:
            self.close()
            self.app.quit()
            return True
        return False


# Main application class
class GyotoObjectEditorApplication(Gtk.Application):
    """Standalone GTK application for the Gyoto Object Editor.

    This class handles the application lifecycle and window
    management.  It extends Gtk.Application to provide a proper GTK
    application structure.

    """

    _next_wid = 0

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

        # App-level actions
        help_action = Gio.SimpleAction.new("help", None)
        help_action.connect("activate", lambda a, p: self.on_help())
        self.add_action(help_action)

        quit_action = Gio.SimpleAction.new("quit", None)
        quit_action.connect("activate", lambda a, p: self.on_quit())
        self.add_action(quit_action)

        # Accelerators
        self.set_accels_for_action("app.help", ["F1"])
        self.set_accels_for_action("app.quit", ["<Primary>Q"])
        self.set_accels_for_action("win.open", ["<Primary>O"])
        self.set_accels_for_action("win.new", ["<Primary>N"])
        self.set_accels_for_action("win.close", ["<Primary>W"])
        self.set_accels_for_action("win.save", ["<Primary>S"])
        self.set_accels_for_action("win.save-as", ["<Primary><Shift>S"])

    def do_activate(self):
        """Called by GTK when the application starts.

        Creates the main window if it doesn't exist and presents it.
        """
        if not self.windows:
            if self.obj is None:
                window = StartupWindow(self)
                window.present()
                return

            window = GyotoObjectEditorApplicationWindow(
                application=self,
                obj=self.obj,
                connector=self.connector
            )
        for window in self.windows:
            window.present()

    def on_open(self, window):
        """Handle win.open action.

        Opens a file dialog and creates a new editor window with the
        selected file. Closes the parent window if it's a StartupWindow.

        Args:
            action: The Gio.Action
            window: The parent for which the dialog is transient

        """
        dialog = self.create_file_dialog()
        dialog.open(
            window,
            None,
            lambda d, result: self.on_open_file_selected(
                d, result, window
            )
        )

    def on_open_file_selected(self, dialog, result, parent_window):
        """Handle file selection from open dialog.

        Args:
            dialog: The Gtk.FileDialog instance
            result: The result of the dialog operation
            parent_window: The window that triggered the open action

        """
        try:
            file = dialog.open_finish(result)
        except GLib.Error:
            return

        if file is not None:
            if isinstance(parent_window, StartupWindow):
                parent_window.close()
            window = GyotoObjectEditorApplicationWindow(
                application=self,
                obj=file.get_path(),
                connector=None
            )
            window.present()

    def on_new(self, widget):
        """Handle win.new action.

        Shows a popover for selecting object type, then creates a new
        editor window. Closes the parent window if it's a StartupWindow.

        Args:
            action: The Gio.Action
            parameter: GVariant containing the parent widget for the popover

        """
        parent_window = self.get_parent_window(widget)
        self.show_new_popover(widget, parent_window)

    def get_parent_window(self, widget):
        """Find the toplevel window for a given widget.

        Args:
            widget: A GTK widget

        Returns:
            The toplevel Gtk.Window, or None

        """
        while widget is not None:
            if isinstance(widget, Gtk.Window):
                return widget
            try:
                widget = widget.get_parent()
            except AttributeError:
                return None
        return None

    def show_new_popover(self, parent_widget, parent_window):
        """Show the New popover for object type selection.

        Args:
            parent_widget: Widget to parent the popover to
            parent_window: Window to close if it's a StartupWindow

        """
        all_new_items = self.get_all_new_items()

        # Calculate the minimum width based on the longest label
        min_width = 300
        if all_new_items:
            layout = Pango.Layout(parent_widget.get_pango_context())
            for label in all_new_items:
                layout.set_text(label)
                width, __ = layout.get_pixel_size()
                if width > min_width:
                    min_width = width
            min_width += 50

        # Create the popover
        popover = Gtk.Popover()
        popover.set_parent(parent_widget)
        popover.set_autohide(True)

        # Main box
        main_box = Gtk.Box(
            orientation=Gtk.Orientation.VERTICAL,
            spacing=6
        )
        main_box.set_size_request(min_width, 400)
        popover.set_child(main_box)

        # Title
        title = Gtk.Label(
            label=_("New File"),
            halign=Gtk.Align.START
        )
        title.add_css_class("h3")
        main_box.append(title)

        # Search entry
        search_entry = Gtk.SearchEntry()
        search_entry.set_placeholder_text(_("Search..."))
        search_entry.set_hexpand(True)
        main_box.append(search_entry)

        # Create the list model
        list_store = Gio.ListStore(item_type=NewItem)
        for item in all_new_items:
            list_store.append(NewItem(label=item))

        # Create a custom filter function
        def filter_func(item):
            search_text = search_entry.get_text().lower()
            if not search_text:
                return True
            return search_text in item.label.lower()

        # Create a custom filter
        custom_filter = Gtk.CustomFilter.new(filter_func)

        # Create filter list model
        filter_model = Gtk.FilterListModel(
            model=list_store,
            filter=custom_filter
        )

        # Create selection model
        selection = Gtk.SingleSelection(model=filter_model)

        # Create list view
        factory = Gtk.SignalListItemFactory()
        factory.connect("setup", self.on_new_list_item_setup)
        factory.connect("bind", self.on_new_list_item_bind)

        list_view = Gtk.ListView(
            model=selection,
            factory=factory
        )
        list_view.set_hexpand(True)
        list_view.set_vexpand(True)
        list_view.add_css_class("rich-list")

        # Add list view to a scrolled window
        scrolled = Gtk.ScrolledWindow()
        scrolled.set_policy(
            Gtk.PolicyType.NEVER,
            Gtk.PolicyType.AUTOMATIC
        )
        scrolled.set_child(list_view)
        scrolled.set_hexpand(True)
        scrolled.set_vexpand(True)
        scrolled.set_min_content_height(300)
        main_box.append(scrolled)

        # Buttons
        button_box = Gtk.Box(
            orientation=Gtk.Orientation.HORIZONTAL,
            spacing=6
        )
        button_box.add_css_class("linked")

        cancel_button = Gtk.Button(label=_("Cancel"))
        cancel_button.connect("clicked", lambda *_: popover.popdown())
        button_box.append(cancel_button)

        create_button = Gtk.Button(label=_("Create"))
        create_button.add_css_class("suggested-action")
        button_box.append(create_button)

        main_box.append(button_box)

        # Connect search entry
        search_entry.connect(
            "search-changed",
            lambda e: custom_filter.changed(Gtk.FilterChange.DIFFERENT)
        )

        # Define callbacks with closure over local variables
        def on_create_clicked(*args):
            if selection.get_selected_item() is None:
                return
            selected_item = selection.get_selected_item()
            label = selected_item.label

            obj = self.create_object_from_label(label)

            popover.popdown()
            if isinstance(parent_window, StartupWindow):
                parent_window.close()
            window = GyotoObjectEditorApplicationWindow(
                application=self,
                obj=obj,
                objtype=label,
                connector=None
            )
            window.present()

        def on_search_key_pressed(controller, keyval, keycode, state):
            if keyval == Gdk.KEY_Down:
                list_view.grab_focus()
                if filter_model.get_n_items() > 0:
                    selection.set_selected(0)
                return True
            return False

        def on_list_key_pressed(controller, keyval, keycode, state):
            if (
                    (hasattr(Gdk, 'KEY_Enter') and keyval == Gdk.KEY_Enter)
                    or (hasattr(Gdk, 'KEY_Return') and keyval == Gdk.KEY_Return)
            ):
                on_create_clicked()
                return True
            return False

        # Add key controller for keyboard navigation
        key_controller = Gtk.EventControllerKey()
        key_controller.connect("key-pressed", on_search_key_pressed)
        search_entry.add_controller(key_controller)

        list_key_controller = Gtk.EventControllerKey()
        list_key_controller.connect("key-pressed", on_list_key_pressed)
        list_view.add_controller(list_key_controller)

        popover.popup()


        create_button.connect("clicked", on_create_clicked)
        list_view.connect("activate", lambda l, p: on_create_clicked())

    def on_new_list_item_setup(self, factory, list_item):
        """Set up the UI for a list item in the New popover.

        Args:
            factory: The Gtk.SignalListItemFactory creating the item.
            list_item: The Gtk.ListItem to set up.

        """
        box = Gtk.Box(
            orientation=Gtk.Orientation.HORIZONTAL,
            spacing=6
        )
        list_item.set_child(box)

        label = Gtk.Label(
            halign=Gtk.Align.START,
            xalign=0
        )
        label.add_css_class("title")
        box.append(label)

    def on_new_list_item_bind(self, factory, list_item):
        """Bind the data to a list item in the New popover.

        Args:
            factory: The Gtk.SignalListItemFactory binding the item.
            list_item: The Gtk.ListItem to bind.

        """
        item = list_item.get_item()
        label = list_item.get_child().get_first_child()
        label.set_text(item.label)

    def get_all_new_items(self):
        """Return list of all available object types for New popover."""
        all_items = []
        for kind in ["Photon", "Scenery", "Screen"]:
            all_items.append(("core", kind, f"core/built-in/{kind}"))
        for section in 'Astrobj', 'Metric', 'Spectrometer', 'Spectrum':
            generic = getattr(core, section)
            kinds = list(generic.registeredPluginsSlashKinds())
            kinds.sort()
            for kind in kinds:
                all_items.append((section, kind, f"{section}/{kind}"))
        return [item[2] for item in all_items]

    def create_object_from_label(self, label):
        """Create a Gyoto object from a label string."""
        nspace, plg, kind = label.split('/', 2)
        if nspace == 'core':
            return getattr(core, kind)()
        else:
            return getattr(core, nspace)(kind, (plg,))

    def create_file_dialog(self):
        """Create and return a configured file dialog."""
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

        return dialog

    def on_help(self, *args):
        """Display the help dialog."""
        parent = self.windows[0] if self.windows else None

        dialog = Gtk.Dialog(
            title="Help",
            transient_for=parent,
            modal=False
        )
        dialog.set_default_size(600, 400)

        help_text = (
            "<span font_weight='bold' size='x-large'>"
            "Gyoto Object Editor</span>\n"
            "\n"
            "<span font_weight='bold' size='large'>"
            "Keyboard Shortcuts:</span>\n"
            "<b>•</b> Ctrl+N: Create a new object\n"
            "<b>•</b> Ctrl+O: Open an object from file\n"
            "<b>•</b> Ctrl+S: Save the current object\n"
            "<b>•</b> Ctrl+Shift+S: Save the current object as...\n"
            "<b>•</b> Ctrl+W: Close the current window\n"
            "<b>•</b> Ctrl+Q: Quit the application\n"
            "<b>•</b> F1: Show this help dialog\n"
            "\n"
            "<span font_weight='bold' size='large'>"
            "Menu Options:</span>\n"
            "<b>•</b> New: Create a new object (search and select from "
            "available types)\n"
            "<b>•</b> Open: Load an object from an XML file\n"
            "<b>•</b> Save: Save the current object to its last used file\n"
            "<b>•</b> Save As: Save the current object to a new file\n"
            "<b>•</b> Help: Show this dialog\n"
            "<b>•</b> Close: Close the current window\n"
            "<b>•</b> Quit: Close all windows and exit\n"
        )

        label = Gtk.Label(
            label=help_text,
            halign=Gtk.Align.START,
            wrap=True,
            use_markup=True,
            xalign=0.0,
            justify=Gtk.Justification.FILL
        )

        scrolled = Gtk.ScrolledWindow()
        scrolled.set_child(label)
        scrolled.set_policy(
            Gtk.PolicyType.AUTOMATIC,
            Gtk.PolicyType.AUTOMATIC
        )

        dialog.set_child(scrolled)
        dialog.present()

    def on_quit(self, *args):
        """Quit the application."""
        self.close_all_windows()
        self.quit()

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

    def get_new_wid(self):
        """Get identification number for new window

        This method should be called be new application windows to
        determine their identification number.

        Side effects:
           Increment _next_wid

        Returns:
           int: identification number for new window
        """
        wid = self._next_wid
        self._next_wid += 1
        return wid


# Main application window class
class GyotoObjectEditorApplicationWindow(Gtk.ApplicationWindow):
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
    wid = 0
    type = 'undefined'

    def __init__(self, application=None, obj=None, objtype='', connector=None):
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

        if objtype == '':
            self.type = None
        else:
            self.type = objtype

        # obj may be the XML description of the object
        if isinstance(obj, str):
            if obj.lower().endswith('.xml'):
                self.filename = obj
            factory = Factory(obj)
            obj = getattr(factory, factory.kind().lower())()
            if self.type is None:
                self.type = f'{factory.kind()}/{obj.kind()}'
        elif obj is None:
            obj = self.default_obj()
            if self.type is None:
                self.type = 'Scenery/'
        elif isinstance(obj, Photon):
            self.type = 'Photon'
        elif isinstance(obj, Screen):
            self.type = 'Screen'
        elif isinstance(obj, Scenery):
            self.type = 'Scenery'
        else:
            if self.type is None:
                for nspace, generic in {
                        'Astrobj': core.Astrobj,
                        'Metric': core.Metric,
                        'Spectrometer': core.Spectrometer,
                        'Spectrum': core.Spectrum,
                }.items():
                    if isinstance(obj, generic):
                        self.type = f'{nspace}/{obj.kind()}'
                        break

        self.obj = obj

        self.connect("close-request", self.on_close_request)

        # Build UI
        self.build_headerbar()

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
            self.wid = application.get_new_wid()

        self.set_title()

    def set_title(self, *args):
        """Set window title according to attributes

        The window title includes the ID number and filename.
        """
        if len(args):
            super().set_title(*args)
        title = ""
        if self.get_application():
            title += f"[{self.wid}] "
        title += self.type
        if self.filename:
            title += f": {os.path.basename(self.filename)}"
        super().set_title(title)

    ####################################################################
    # UI
    ####################################################################

    def build_headerbar(self):
        """Build the window's header bar with menu button.

        Creates a header bar with a hamburger menu containing:
        - New... (opens a popover with search and list)
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

        # Add "New" as a menu item that opens the popover
        new_menu_item = Gio.MenuItem.new(_("New..."), "win.new")
        menu.append_item(new_menu_item)

        menu.append(_("Open…"), "win.open")
        menu.append(_("Save"), "win.save")
        menu.append(_("Save As…"), "win.save-as")

        menu_section2 = Gio.Menu()
        menu_section2.append(_("Help"), "app.help")
        menu_section2.append(_("Close"), "win.close")
        menu_section2.append(_("Quit"), "app.quit")

        menu.append_section(None, menu_section2)

        # Create menu button
        menu_button = Gtk.MenuButton(
            icon_name="open-menu-symbolic",
            menu_model=menu,
            use_underline=True
        )
        menu_button.add_css_class("flat")
        header.pack_end(menu_button)
        self.menu_button = menu_button

        # Connect win actions
        save_action = Gio.SimpleAction.new("save", None)
        save_action.connect("activate", lambda a, p: self.on_save(a, p))
        self.add_action(save_action)

        save_as_action = Gio.SimpleAction.new("save-as", None)
        save_as_action.connect("activate",
                                lambda a, p: self.on_save_as(a, p))
        self.add_action(save_as_action)

        close_action = Gio.SimpleAction.new("close", None)
        close_action.connect("activate", lambda a, p: self.on_close(a, p))
        self.add_action(close_action)

        open_action = Gio.SimpleAction.new("open", None)
        open_action.connect(
            "activate",
            lambda a, p: self.props.application.on_open(self)
        )
        self.add_action(open_action)

        new_action = Gio.SimpleAction.new("new", None)
        new_action.connect(
            "activate",
            lambda a, p: self.props.application.on_new(self.menu_button)
        )
        self.add_action(new_action)

    ####################################################################
    # Callbacks
    ####################################################################

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

        self.set_title()

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

    def default_obj(self):
        """Return default object
        """
        return Scenery()

# Stand-alone entry point:
if __name__ == "__main__":
    raise SystemExit(
        GyotoObjectEditorApplication.run_app(parsecliargs=True)
    )
