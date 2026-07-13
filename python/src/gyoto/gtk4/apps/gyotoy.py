# Imports + MyApplication + MainWindow skeleton

"""
Main application window.

Layout
------

    ┌────────────────────────────────────────────────────────────────────┐
    │ MyApp                                                     ☰        │
    ├────────────────────────────────────────────────────────────────────┤
    │ ┌──────────────────────────────┬────────────────────────────────┐  │
    │ │                              │  ○ Star                        │  │
    │ │                              │  ○ Photon                      │  │
    │ │                              │                                │  │
    │ │      Matplotlib canvas       │  ┌──────────────────────────┐  │  │
    │ │                              │  │ PropertyEditorBox        │  │  │
    │ │                              │  └──────────────────────────┘  │  │
    │ └──────────────────────────────┴────────────────────────────────┘  │
    ├────────────────────────────────────────────────────────────────────┤
    │ ███████──── ☑ Interpolate Step:[1e-3▲▼] N:[100▲▼]  ⏮ ▶ ⏹           │
    └────────────────────────────────────────────────────────────────────┘
"""

from __future__ import annotations

__all__ = ['GyotoyApplication', 'MainWindow']

import gi
gi.require_version("Gtk", "4.0")
from gi.repository import Gtk, Gio, GLib

import numpy

import warnings

from gettext import gettext as _

# Custom widgets
from ..widgets.property_editor_box import PropertyEditorBox
from ..widgets.scientific_spin import ScientificSpin
from ..widgets.viewer_3d import Viewer3D
from ..widgets.simulation_controls import SimulationControls
from ..utils import show_error_dialog
from ...core import Factory, Error as GyotoError, Photon
from ...std import Star

# Main application
class GyotoyApplication(Gtk.Application):
    """Standalone GTK application."""

    def __init__(self):
        super().__init__(
            application_id="fr.obspm.gyoto.Gyotoy",
            flags=Gio.ApplicationFlags.DEFAULT_FLAGS,
        )

    def do_activate(self):
        """Called by GTK when the application starts."""

        window = self.props.active_window

        if window is None:
            window = MainWindow(application=self)

        window.present()


# Main window

class MainWindow(Gtk.ApplicationWindow):
    """Main application window."""

    blocking = True
    main_loop = None
    particle = None

    ####################################################################
    # Construction
    ####################################################################

    def __init__(self, application=None, particle=None):
        super().__init__(application=application)

        self.set_title("Gyotoy")
        self.set_default_size(1400, 900)

        # Build UI
        self.build_headerbar()
        self.build_body()
        self.connect("close-request", self.on_close_request)

        # Populate initial editor
        self.set_particle(particle)

    ####################################################################
    # UI
    ####################################################################

    def build_headerbar(self):

        # create title bar with hamberguer button
        header = Gtk.HeaderBar()
        self.set_titlebar(header)
        menu_button = Gtk.MenuButton(
            icon_name="open-menu-symbolic"
        )
        menu_button.add_css_class("flat")
        header.pack_end(menu_button)

        # attach menu to hamburger
        box = Gtk.Box(
            orientation=Gtk.Orientation.VERTICAL
        )
        popover = Gtk.Popover()
        popover.set_child(box)
        menu_button.set_popover(popover)

        # populate menu
        open_button = Gtk.Button(
            label=_("Open…")
        )
        open_button.add_css_class("flat")

        save_button = Gtk.Button(
            label=_("Save As…")
        )
        save_button.add_css_class("flat")

        quit_button = Gtk.Button(
            label=_("Quit")
        )
        quit_button.add_css_class("flat")

        box.append(open_button)
        box.append(save_button)
        box.append(Gtk.Separator())
        box.append(quit_button)

        open_button.connect("clicked", self.on_open)
        save_button.connect("clicked", self.on_save_as)
        quit_button.connect("clicked", self.on_quit)


    def build_body(self):

        #### main container: vertical box
        self.vertbox1 = Gtk.Box(
            orientation=Gtk.Orientation.VERTICAL,
            spacing=6
        )
        self.set_child(self.vertbox1)

        ### first row: paned
        self.paned = Gtk.Paned(
            orientation=Gtk.Orientation.HORIZONTAL
        )
        self.vertbox1.append(self.paned)

        ## first row, left: Viewer3d
        self.viewer = Viewer3D()
        self.paned.set_start_child(
            self.viewer
        )
        self.paned.set_resize_start_child(True)

        ## first row, right: vertical box
        right = Gtk.Box(
            orientation=Gtk.Orientation.VERTICAL,
            spacing=6
        )
        self.paned.set_end_child(right)

        # top of controil column: Star/Photon radio buttons
        frame = Gtk.Frame(
            label="Particle type"
        )
        right.append(frame)

        buttons = Gtk.Box(
            orientation=Gtk.Orientation.HORIZONTAL,
            spacing=12
        )
        frame.set_child(buttons)

        self.particle_star = Gtk.CheckButton(
            label="Star"
        )
        self.particle_photon = Gtk.CheckButton(
            label="Photon"
        )
        self.particle_photon.set_group(
            self.particle_star
        )
        buttons.append(self.particle_star)
        buttons.append(self.particle_photon)

        self.particle_star.connect(
            "toggled",
            self.on_star_selected
        )

        # I don't think we need both
        # self.photon_button.connect(
        #     "toggled",
        #     self.on_photon_selected
        # )

        # then scrollable property editor view
        scroll = Gtk.ScrolledWindow()
        scroll.set_hexpand(True)
        scroll.set_vexpand(True)
        right.append(scroll)
        self.editor_scroller = scroll

        # an additional spin button for end time
        frame = Gtk.Frame()
        frame.set_label("End time")
        frame.set_label_align(0.0)
        frame.set_tooltip_text("Run simulation from initial contition time to end time")

        hbox = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=5)
        frame.set_child(hbox)

        spin = ScientificSpin(value=0.)
        spin.connect("value-changed", self.on_parameter_changed)
        hbox.append(spin)
        right.append(frame)

        # not sure we should initialize this early
        # self.editor = PropertyEditorBox(...)
        # scroll.set_child(
        #     self.editor
        # )

        # then controls for running the integration
        self.controls = SimulationControls()
        self.vertbox1.append(
            self.controls
        )
        self.controls.connect(
            "play-pause",
            self.on_play_pause
        )

        self.controls.connect(
            "reset",
            self.on_reset
        )

        self.controls.connect(
            "stop",
            self.on_stop
        )

        self.controls.connect(
            "step-changed",
            self.on_step_changed
        )

        self.controls.connect(
            "nframes-changed",
            self.on_nframes_changed
        )

        self.controls.connect(
            "interpolate-changed",
            self.on_interpolate_changed
        )
        
    ####################################################################
    # Entry points
    ####################################################################

    @staticmethod
    def run(particle=None, blocking=True):
        """Contruct a Gyoto ObjectEditor window and run it

        synopsis:
         ObjectEditor.run(particle, [blocking=True])
        """
        win = MainWindow(particle)
        win.blocking = blocking
        win.present()
        if blocking:
            win.main_loop=GLib.MainLoop()
            win.main_loop.run()

        return win

    def on_close_request(self, *args):
        if self.main_loop is not None:
            GLib.idle_add(self.main_loop.quit)
        return False

    def on_quit(self, *args):
        self.close()

    @staticmethod
    def run_app():
        """Run as a standalone GTK application."""

        app = GyotoyApplication()
        return app.run(None)

    ####################################################################
    # UI construction
    ####################################################################

    def _create_actions(self):
        """Create Gio actions."""
        pass

    def _build_ui(self):
        """Build the complete widget hierarchy."""
        pass

    ####################################################################
    # Redraw
    ####################################################################

    def redraw(self):
        print(f'in redraw: type(self.particle):{type(self.particle)}')
        if self.particle is None:
            return
        starttime = self.particle.initCoord()[0]
        self.endtime=3000
        self.particle.xFill(self.endtime)
        interpolate = False
        if interpolate:
            interp_step=1e-2
            npoints = int(abs(self.endtime-starttime)/interp_step)
            t=numpy.linspace(starttime, self.endtime, npoints)
        else:
            npoints = self.particle.get_nelements()
            t = numpy.ndarray(npoints)
            self.particle.get_t(t)
        x=numpy.ndarray(npoints)
        y=numpy.ndarray(npoints)
        z=numpy.ndarray(npoints)
        self.particle.getCartesian(t, x, y, z)
        self.viewer.axes.plot(x, y, z)

    ####################################################################
    # Callbacks
    ####################################################################

    def __getattr__(self, name):

        # if not name.startswith("on_"):
        #     raise AttributeError(
        #         f"{type(self).__name__!s} object has no attribute {name!r}"
        #     )

        warned = False

        def callback(*args, **kwargs):
            nonlocal warned

            if not warned:
                warnings.warn(
                    f"Unimplemented method: {type(self).__name__}.{name}()",
                    RuntimeWarning,
                    stacklevel=2,
                )
                warned = True

        setattr(self, name, callback)
        return callback

    def on_open(self, *args):

        dialog = Gtk.FileDialog()

        dialog.open(
            self,
            None,
            lambda dialog, result:
                self.on_open_file_selected(dialog, result)
        )

    def on_open_file_selected(self, dialog, result):
        try:
            file = dialog.open_finish(result)
        except GLib.Error:
            return

        factory = None
        if file is not None:
            try:
                factory = Factory(file.get_path())
            except GyotoError as e:
                show_error_dialog(
                    message=f"Error loading XML file{file.get_path()}:",
                    detail=e.get_message(),
                    window=self
                )

        if factory is not None:
            kind = factory.kind()
            particle = None
            if kind in ('Astrobj', 'Scenery'):
                try:
                    ao = factory.astrobj()
                except GyotoError as e:
                    show_error_dialog(
                        message=f"Could construct Astrobj from XML file:",
                        detail=e.get_message(),
                        window=self
                    )
                try:
                    particle = Star(ao)
                except GyotoError as e:
                    show_error_dialog(
                        message=f"Could not cast Astrobj to Star:",
                        detail=e.get_message(),
                        window=self
                    )
            elif kind == 'Photon':
                try:
                    particle = factory.getPhoton()
                except GyotoError as e:
                    show_error_dialog(
                        message=f"Could construct Photon from XML file:",
                        detail=e.get_message(),
                        window=self
                    )
            else:
                show_error_dialog(
                    message=f"Could not load particle from XML file",
                    detail="XML should describe an Astrobj, Scenery or Photon",
                    window=self
                )

            if particle is not None:
                self.set_particle(particle)
        print(f'in on_open_file_selected: type(self.particle):{type(self.particle)}')

    ####################################################################
    # Setters / getters
    ####################################################################

    def set_particle(self, particle):
        if isinstance(particle, Star):
            self.particle_star.set_active(True)
        elif isinstance(particle, Photon):
            self.particle_photon.set_active(True)
        else:
            show_error_dialog(message= "Wrong type for particle:",
                              detail= repr(type(particle)),
                              window=self
            )
            return
        self.particle = particle
        print(f'in set_particle: type(self.particle):{type(self.particle)}')
        self.editor = PropertyEditorBox(particle)
        self.editor_scroller.set_child(self.editor)
        self.editor.connect('value-changed', self.on_particle_changed)
        self.redraw()

    def get_particle(self):
        return self.particle

# Widget construction (__init__)
# Header bar, menu, Matplotlib embedding, and layout
# Callbacks and helper methods
# run() / run_app()

### Stand-alone entry point:

if __name__ == "__main__":
    raise SystemExit(MainWindow.run_app())
