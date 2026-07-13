## A window to edit a gyoto object's properties

__all__ = ['ObjectEditor']

import gi
gi.require_version("Gtk", "4.0")
from gi.repository import Gtk, GLib

from ..widgets.property_editor_box import PropertyEditorBox

class ObjectEditor(Gtk.Window):
    """A Gyoto object editor

    Displays a window giving an editable view of all the properties of
    a Gyoto object.

    Parameters:
      obj: the object to edit
      blocking: whether this object should manage the GLib loop

    """

    @staticmethod
    def run(obj, blocking=True):
        """Contruct a Gyoto ObjectEditor window and run it

        synopsis:
         ObjectEditor.run(obj, [blocking=True])
        """
        win = ObjectEditor(obj, blocking)
        win.present()
        if blocking:
            win.main_loop.run()

    def __init__(self, obj, blocking=True):
        Gtk.Window.__init__(self, title="Gyoto Object Editor")
        self.set_default_size(400, 600)
        self.obj = obj
        if blocking:
            self.main_loop=GLib.MainLoop()
        else:
            self.main_loop=None

        self.connect("close-request", self.on_close_request)

        # Main container
        self.scrolled_window = Gtk.ScrolledWindow()
        self.set_child(self.scrolled_window)

        # Vertical box for property widgets
        self.vbox = PropertyEditorBox(obj)
        self.vbox.set_margin_top(10)
        self.vbox.set_margin_bottom(10)
        self.vbox.set_margin_start(10)
        self.vbox.set_margin_end(10)
        self.scrolled_window.set_child(self.vbox)

        # We could connect to value-changed to react wheneverr obj changes 
        # self.vbox.connect("value-changed", self.on_value_changed)
        # but we don't actually need to.

    def on_close_request(self, *args):
        if self.main_loop is not None:
            GLib.idle_add(self.main_loop.quit)
        return False
