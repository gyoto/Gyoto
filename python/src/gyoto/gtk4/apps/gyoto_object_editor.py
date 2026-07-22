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

import argparse
import sys
import traceback

from gettext import gettext as _

from ..widgets.property_editor_box import PropertyEditorBox
from ...core import Factory, Scenery, Error as GyotoError, Photon, Screen
from ... import core

# Custom GObject class to hold the label of a 'New >' menu item
class NewItem(GObject.Object):
    __gtype_name__ = 'NewItem'

    label = GObject.Property(type=str)

    def __init__(self, label=None):
        super().__init__()
        self.label = label

# Main application class
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
        menu_section2.append(_("Help"), "win.help")
        menu_section2.append(_("Close"), "win.close")
        menu_section2.append(_("Quit"), "win.quit")

        menu.append_section(None, menu_section2)

        # Create menu button
        menu_button = Gtk.MenuButton(
            icon_name="open-menu-symbolic",
            menu_model=menu,
            use_underline=True
        )
        menu_button.add_css_class("flat")
        header.pack_end(menu_button)
        self.menu_button=menu_button

        # Build the New popover (parented to the menu button)
        self.build_new_popover()

        # Connect actions
        action_group = Gio.SimpleActionGroup()
        self.insert_action_group("win", action_group)

        action_group.add_action_entries([
            ("new", self.on_new_menu_item_activated, None),
            ("open", self.on_open, None),
            ("save", self.on_save, None),
            ("save-as", self.on_save_as, None),
            ("help", self.on_help, None),
            ("close", self.on_close, None),
            ("quit", self.on_quit, None),
        ])

    def build_new_popover(self):
        """Build the New popover with search and list."""
        # Prepare all available items
        all_items = []
        for kind in ["Photon", "Scenery", "Screen"]:
            all_items.append(("core", kind, f"core/built-in/{kind}"))
        for section in 'Astrobj', 'Metric', 'Spectrometer', 'Spectrum':
            generic = getattr(core, section)
            kinds = list(generic.registeredPluginsSlashKinds())
            kinds.sort()
            for kind in kinds:
                all_items.append((section, kind, f"{section}/{kind}"))

        # Create a flat list of all items
        self.all_new_items = [item[2] for item in all_items]

        # Calculate the minimum width based on the longest label
        min_width = 300  # Default minimum width
        if self.all_new_items:
            # Create a Pango layout to measure text width
            layout = Pango.Layout(self.get_pango_context())
            for label in self.all_new_items:
                layout.set_text(label)
                width, __ = layout.get_pixel_size()
                if width > min_width:
                    min_width = width
            # Add some padding
            min_width += 50

        # Create the popover
        self.new_popover = Gtk.Popover()
        self.new_popover.set_parent(self.menu_button)
        self.new_popover.set_autohide(True)

        # Main box
        main_box = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=6)
        main_box.set_size_request(min_width, 400)  # Set minimum size for the main box
        self.new_popover.set_child(main_box)

        # Title
        title = Gtk.Label(label=_("New File"), halign=Gtk.Align.START)
        title.add_css_class("h3")
        main_box.append(title)

        # Search entry
        self.new_search_entry = Gtk.SearchEntry()
        self.new_search_entry.set_placeholder_text(_("Search..."))
        self.new_search_entry.connect("search-changed", self.on_new_search_changed)
        self.new_search_entry.set_hexpand(True)  # Allow horizontal expansion
        main_box.append(self.new_search_entry)

        # Create the list model
        self.new_list_store = Gio.ListStore(item_type=NewItem)
        for item in self.all_new_items:
            self.new_list_store.append(NewItem(label=item))

        # Create a custom filter function
        def filter_func(item):
            search_text = self.new_search_entry.get_text().lower()
            if not search_text:
                return True
            return search_text in item.label.lower()

        # Create a custom filter
        self.custom_filter = Gtk.CustomFilter.new(filter_func)

        # Create filter list model
        self.new_filter_model = Gtk.FilterListModel(
            model=self.new_list_store,
            filter=self.custom_filter
        )

        # Create selection model
        self.new_selection = Gtk.SingleSelection(model=self.new_filter_model)

        # Create list view
        factory = Gtk.SignalListItemFactory()
        factory.connect("setup", self.on_new_list_item_setup)
        factory.connect("bind", self.on_new_list_item_bind)

        self.new_list_view = Gtk.ListView(
            model=self.new_selection,
            factory=factory
        )
        self.new_list_view.set_hexpand(True)  # Expand horizontally
        self.new_list_view.set_vexpand(True)  # Expand vertically
        self.new_list_view.add_css_class("rich-list")
        self.new_list_view.connect("activate", self.on_new_list_activate)

        # Add list view to a scrolled window
        scrolled = Gtk.ScrolledWindow()
        scrolled.set_policy(Gtk.PolicyType.NEVER, Gtk.PolicyType.AUTOMATIC)  # No horizontal scroll, auto vertical
        scrolled.set_child(self.new_list_view)
        scrolled.set_hexpand(True)
        scrolled.set_vexpand(True)
        scrolled.set_min_content_height(300)  # Minimum height for the list
        main_box.append(scrolled)

        # Buttons
        button_box = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=6)
        button_box.add_css_class("linked")

        cancel_button = Gtk.Button(label=_("Cancel"))
        cancel_button.connect("clicked", lambda *_: self.new_popover.popdown())
        button_box.append(cancel_button)

        create_button = Gtk.Button(label=_("Create"))
        create_button.add_css_class("suggested-action")
        create_button.connect("clicked", self.on_new_create_clicked)
        button_box.append(create_button)

        main_box.append(button_box)

        # Add key controller to the search entry for keyboard navigation
        key_controller = Gtk.EventControllerKey()
        key_controller.connect("key-pressed", self.on_new_key_pressed)
        self.new_search_entry.add_controller(key_controller)

        # Add key controller to the list view for Up/Down navigation
        list_key_controller = Gtk.EventControllerKey()
        list_key_controller.connect("key-pressed", self.on_list_key_pressed)
        self.new_list_view.add_controller(list_key_controller)

    def filter_new_items(self, item, search_entry):
        """Filter function for the New items list."""
        search_text = search_entry.get_text().lower()
        if not search_text:
            return True
        return search_text in item.label.lower()

    def on_new_menu_item_activated(self, action, parameter, *args):
        """Show the New popover when the New menu item is activated."""
        self.new_popover.popup()

    def on_new_button_clicked(self, button, *args):
        """Show the New popover."""
        self.new_popover.popup()

    def on_new_search_changed(self, entry):
        """Update filter when search text changes."""
        self.custom_filter.changed(Gtk.FilterChange.DIFFERENT)

    def on_new_list_item_setup(self, factory, list_item):
        """Set up list items."""
        box = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=6)
        list_item.set_child(box)

        # Label for the item
        label = Gtk.Label(halign=Gtk.Align.START, xalign=0)
        label.add_css_class("title")
        box.append(label)

    def on_new_list_item_bind(self, factory, list_item):
        """Bind list items to data."""
        item = list_item.get_item()
        label = list_item.get_child().get_first_child()
        label.set_text(item.label)

    def on_new_list_activate(self, list_view, position):
        """Handle activation (Enter key) in the list."""
        self.on_new_create_clicked()

    def on_new_key_pressed(self, controller, keyval, keycode, state):
        """Handle keyboard navigation in the search entry."""
        if keyval == Gdk.KEY_Down:
            # Move focus to list view and select first item
            self.new_list_view.grab_focus()
            if self.new_filter_model.get_n_items() > 0:
                self.new_selection.set_selected(0)
            return True
        return False

    def on_new_create_clicked(self, *args):
        """Create the selected item."""
        if self.new_selection.get_selected_item() is None:
            return

        selected_item = self.new_selection.get_selected_item()
        label = selected_item.label

        nspace, plg, kind = label.split('/', 2)

        if nspace == 'core':
            obj = getattr(core, kind)()
        else:
            obj = getattr(core, nspace)(kind, (plg,))

        window = GyotoObjectEditorApplicationWindow(
            application=self.props.application,
            obj=obj,
            connector=None
        )
        window.present()

        self.new_popover.popdown()

    def on_list_key_pressed(self, controller, keyval, keycode, state):
        """Handle keyboard navigation in the list view."""
        if (
                (hasattr(Gdk, 'KEY_Enter') and keyval == Gdk.KEY_Enter)
                or (hasattr(Gdk, 'KEY_Return') and keyval == Gdk.KEY_Return)
        ):
            self.on_new_create_clicked()
            return True
        return False

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
                    keyval=Gdk.KEY_n,
                    modifiers=Gdk.ModifierType.CONTROL_MASK
                ),
                action=Gtk.CallbackAction.new(self.on_new_button_clicked)
            )
        )

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
