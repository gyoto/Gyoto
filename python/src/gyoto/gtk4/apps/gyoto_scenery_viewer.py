"""Gyoto: GTK4 Application for Gyoto Scenery Vizualization

This application provides a graphical interface for simulating and
visualizing ray-traced sceneries in spacetimes supported by the Gyoto
library.

Layout
------
    ┌────────────────────────────────────────────────────────────────┐
    │ Gyoto Scenery Viewer                                     ☰     │
    ├────────────────────────────────────────────────────────────────┤
    │ ┌──────────────────────────┬────────────────────────────────┐  │
    │ │                          │  ┌──────────────────────────┐  │  │
    │ │      Matplotlib canvas   │  │ PropertyEditorBox        │  │  │
    │ │                          │  └──────────────────────────┘  │  │
    │ └──────────────────────────┴────────────────────────────────┘  │
    ├────────────────────────────────────────────────────────────────┤
    │ ██████████████████████──────────────────────────────────────── │
    │                                                                │
    │ Status...                    N: [100  ]  ⏮        ▶        ⏹   │
    │                                                                │
    └────────────────────────────────────────────────────────────────┘

Description
-----------

- **Left Panel**: 2D Matplotlib viewer displaying the scenery
    observable quantities.
- **Right Panel**: Property editor for adjusting scenery parameters.
- **Bottom Controls**: Play/pause/stop buttons, integration
    settings, and status display.
- **Worker Process**: Heavy computations run in a background process
    to keep the UI responsive.

Usage
-----
Run as a standalone application:
    python3 -m gyoto.gtk4.apps.gyoto_scenery_viewer [-h] [filename.xml]

Or import and use programmatically:
    from gyoto.gtk4 import view_scenery
    view_scenery([scenery|'filename.xml'])

An optional scenery (gyoto.core.Scenery), or the name of an XML file
describing such a scenery, can be provided.

"""

from __future__ import annotations

__all__ = ['GyotoSceneryViewerApplication',
           'GyotoSceneryViewerApplicationWindow']

import gi
gi.require_version("Gtk", "4.0")
from gi.repository import Gtk, Gio, GLib, Gdk

import sys
import argparse
import numpy
from numpy.typing import ArrayLike
from multiprocessing import Queue, Event, Process
import time
import queue
import traceback
import warnings

from gettext import gettext as _

# Custom widgets
from ..widgets.property_editor_box import PropertyEditorBox
from ..widgets.scientific_spin import ScientificSpin
from ..widgets.viewer_2d import Viewer2D
from ..widgets.simulation_controls import SimulationControls
from ..utils import show_error_dialog
from ...utils import readScenery
from ...core import Factory, Astrobj, Error as GyotoError
from ...core import Scenery, Screen
from ...std import FixedStar, Minkowski, PowerLaw

# --- Commands for worker communication ---
RUN_SIM = 'run'
QUIT = 'quit'

def worker_func(cmd_queue, progress_queue, control_queue, pause_event,
                stop_event):
    """Persistent worker process for running simulations.

    This function runs in a separate process to avoid blocking the GTK
    main thread. It waits for commands from cmd_queue, executes
    simulations, and sends progress updates to progress_queue and
    control messages to control_queue.

    The worker can be paused via pause_event and stopped via stop_event.

    Args:
        cmd_queue: Queue for receiving commands (RUN_SIM, QUIT)
        progress_queue: Queue for sending progress updates
        control_queue: Queue for sending control messages (log, done,
            error)
        pause_event: Event to pause/resume the simulation
        stop_event: Event to stop the simulation

    """
    try:
        import numpy
        import time
        from ...core import Factory, Scenery, AstrobjProperties, array_double
        from ...core import Range, Grid, verbose

        verbose(0)

        while True:
            try:
                # Blocks until command arrives or timeout
                cmd = cmd_queue.get(timeout=1.0)
            except queue.Empty:
                continue

            if cmd[0] == QUIT:
                break

            elif cmd[0] == RUN_SIM:
                end_msg = ('done',)
                try:
                    _, scenery, ilim, jlim, nframes = cmd

                    stop_event.clear()
                    pause_event.clear()

                    # Read resolution from Scenery
                    res = scenery.Screen.Resolution

                    # Transform ilim, jlim to Gyoto pixels and sort them
                    ilim += 1
                    jlim += 1
                    ilim.sort()
                    jlim.sort()

                    # No transform that to integers, making sure we
                    # preserve the borders and don't excede the field
                    ilim[ilim < 0.5] = 0.5
                    jlim[jlim < 0.5] = 0.5
                    ilim[ilim > res+0.5] = res + 0.5
                    jlim[jlim > res+0.5] = res + 0.5

                    ilim = list(ilim)
                    jlim = list(jlim)

                    ilim[0] = int(ilim[0]+0.5)
                    upper = ilim[1]+0.5
                    if int(upper) != upper:
                        ilim[1] = int(upper)
                    else:
                        ilim[1] = int(upper)-1

                    jlim[0] = int(jlim[0]+0.5)
                    upper = jlim[1]+0.5
                    if int(upper) != upper:
                        jlim[1] = int(upper)
                    else:
                        jlim[1] = int(upper)-1

                    # Compute step between frames
                    if nframes <= 0:
                        nframes = 1

                    nrows = (jlim[1]-jlim[0]+1)
                    step = nrows // nframes

                    if step * nframes < nrows:
                        step += 1

                    if scenery.Screen and scenery.Screen.Spectrometer:
                        nsamples = scenery.Screen.Spectrometer.NSamples
                    else:
                        nsamples = None

                    aop=AstrobjProperties(scenery.Quantities,
                                          shape=(res, res),
                                          nsamples=nsamples,
                                          alloc=True)
                    for key in aop.data:
                        aop.data[key][...] = numpy.nan

                    next_update_fraction = 0.
                    last_update = time.time()

                    for k in range(nframes):
                        if stop_event.is_set():
                            end_msg = ('aborted',)
                            break
                        while pause_event.is_set() and not stop_event.is_set():
                            time.sleep(0.1)

                        first = k * step + jlim[0]
                        if first > jlim[1]:
                            break
                        last = first + step
                        if last > jlim[1]:
                            last = jlim[1]

                        ii=Range(ilim[0], ilim[1], 1)
                        jj=Range(first, last, 1)
                        grid=Grid(ii, jj, "\rj = ")
                        try:
                            scenery.rayTrace(grid, aop)
                        except GyotoError as e:
                            e.Report()

                        progress = last / nrows
                        if (progress >= next_update_fraction
                                and time.time() - last_update > 0.1):
                            progress_queue.put_nowait(
                                ('progress', progress, aop.data)
                            )
                            next_update_fraction = progress + 0.05
                            last_update = time.time()

                    progress_queue.put_nowait(
                        ('progress', progress, aop.data)
                    )
                    control_queue.put(end_msg)
                except Exception as e:
                    traceback.format_exc()
                    control_queue.put(('error', traceback.format_exc()))
    except Exception as e:
        traceback.format_exc()
        control_queue.put(('fatal_error', traceback.format_exc()))

class GyotoSceneryViewerApplication(Gtk.Application):
    """Standalone GTK application for the Gyoto Scenery Viewer.

    This class handles the application lifecycle and window management.
    It extends Gtk.Application to provide a proper GTK application
    structure and maintains a list of all open windows.

    Attributes:
        windows: List of all open GyotoSceneryViewerApplicationWindow instances.
        scenery: the initial scenery (optional)
        connector: Connection for inter-process communication


    """

    def __init__(self, scenery=None, connector=None, *args, **kwargs):
        """Initialize the Gyoto Scenery Viewer( GTK application.

        Args:
            scenery: Initial scenery to display
            connector: Connection for inter-process communication

        """
        if 'application_id' not in kwargs:
            kwargs['application_id'] = "fr.obspm.gyoto.GyotoSceneryViewer"
        if 'flags' not in kwargs:
            kwargs['flags'] = (
                Gio.ApplicationFlags.DEFAULT_FLAGS
                | Gio.ApplicationFlags.NON_UNIQUE
            )
        super().__init__(*args, **kwargs)
        self.scenery = scenery
        self.windows = []

        # Handle QUIT from parent process
        self.connector = connector
        if connector is not None:
            GLib.timeout_add(50, self.check_connector)

    def do_activate(self):
        """Called by GTK when the application starts.

        Creates the main window if it doesn't exist and presents it.

        """
        if not self.windows:
            window = GyotoSceneryViewerApplicationWindow(
                application=self,
                scenery=self.scenery,
                connector=self.connector
            )
            self.windows.append(window)
        for window in self.windows:
            window.present()

    def remove_window(self, window):
        """Remove a window from the application's window list.

        Args:
            window: The GyotoSceneryViewerApplicationWindow to remove.

        """
        if window in self.windows:
            self.windows.remove(window)

    def close_all_windows(self):
        """Close all open windows and quit the application."""
        for window in self.windows[:]:
            window.close()

    @staticmethod
    def run_app(scenery=None, parsecliargs=False, *args, **kwargs):
        """Run the Gyoto Scenery Viewer as a standalone GTK application.

        Parameters:
            scenery: the Gyoto Scenery to start with, or None, or the
                XML description of such a Scenery, or the name of an
                XML file containing this description.
            parsecliargs: whether to parse the command line arguments
            *args, **kwargs: other parameters are passed untouched to
                the GyotoSceneryViewerApplication constructor.

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
                'xmlfile',
                nargs='?',
                help='XML file containing the description of a '
                'Gyoto Scenery (optional)'
            )
            cliargs, remaining = parser.parse_known_args()
            if 'xmlfile' in cliargs:
                scenery = cliargs.xmlfile

        app = GyotoSceneryViewerApplication(scenery=scenery, *args, **kwargs)
        return app.run(remaining)

    def check_connector(self):
        """Check for QUIT commands from the parent process.

        This method is called periodically (every 50ms) via
        GLib.timeout_add to poll the inter-process communication pipe
        for a QUIT message. When received, it closes all windows.

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
        except Exception:
            pass
        return True

class GyotoSceneryViewerApplicationWindow(Gtk.ApplicationWindow):
    """Main application window for the Gyoto Scenery Viewer.

    This window contains:
    - A 2D Matplotlib viewer (left panel) for visualizing ray-traced
      quantities
        - A property editor (right panel) for adjusting the scenery
      parameters
    - Control widgets (bottom) for running/stopping simulations and
      setting parameters

    The window uses a multiprocessing worker for heavy computations to
    keep the UI responsive. Communication with the worker happens via
    multiprocessing Queues for progress updates and control messages.

    Each window instance handles exactly one scenery throughout its
    lifetime. New windows can be created for additional sceneries.

    Attributes:
        scenery: Current scenery being edited/simulated
        viewer: Viewer2D widget for 2D visualization
        editor: PropertyEditorBox for editing scenery properties
        controls: SimulationControls for play/pause/stop
        worker: Process for background computations
        simulation_running: Flag indicating if a simulation is in progress
        last_focused_widget: Last widget that had focus (for focus
            restoration)
        interpolation_step: Step size for interpolation
        connector: Connection for inter-process communication
        filename: Path to the last file used for this scenery

    """

    # Default values
    scenery = None
    viewer = None
    editor = None
    endtime = 3000
    worker = None
    simulation_running = False
    last_focused_widget = None
    connector = None
    filename = None

    # Constants
    scalar_quantities = ('Intensity', 'EmissionTime', 'MinDistance',
                         'FirstDistMin', 'Redshift', 'NbCrossEqPlane',
                         'User1', 'User2', 'User3', 'User4', 'User5')

    ####################################################################
    # Construction
    ####################################################################

    def __init__(self, application=None, scenery=None, connector=None):
        """Initialize the main window.

        Args:
            application: Parent Gtk.Application instance
            scenery: Initial scenery to display
            connector: Connection for inter-process communication

        """
        super().__init__(application=application)

        # initialize dict to store images
        self.data = {}
        # self.data will be a dict of dicts where individual images
        # are addressed as self.data[(delta, semifov)][quantity]

        # Store connector for set_scenery
        self.connector = connector

        # Handle case where scenery is a string.
        # It is either a filename or a full XML description string
        if isinstance(scenery, str):
            if scenery.lower().endswith('.xml'):
                self.filename = scenery
            scenery = Factory(scenery).scenery()

        if scenery is None:
            scenery = self.default_scenery()

        self.set_title("Gyoto Scenery Viewer")
        self.set_default_size(1024, 768)

        # Build UI
        self.build_headerbar()
        self.build_body()
        self.build_shortcuts()
        self.connect("close-request", self.on_close_request)

        # Prepare computation process
        self.cmd_queue = Queue()
        self.progress_queue = Queue()
        self.control_queue = Queue()
        self.pause_event = Event()
        self.stop_event = Event()

        # Start worker process ONCE at window creation
        self.worker = Process(
            target=worker_func,
            args=(
                self.cmd_queue,
                self.progress_queue,
                self.control_queue,
                self.pause_event,
                self.stop_event
            ),
            daemon=True
        )
        self.worker.start()

        # Poll results every 50ms
        GLib.timeout_add(50, self.process_progress)
        GLib.timeout_add(50, self.process_control)

        # Actually set self.scenery and build editor
        self.set_scenery(scenery)

        # Set the initial limits
        self.reset_limits()

        # Register window with application
        if application is not None:
            application.windows.append(self)

    ####################################################################
    # UI
    ####################################################################

    def build_headerbar(self):
        """Build the window's header bar with menu button.

        Creates a header bar with a hamburger menu containing:
        - New
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

        menu_section1 = Gio.Menu()
        menu_section1.append(_("New"), "win.new")
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

        action_group.add_action_entries([
            ("new", self.on_new, None),
            ("open", self.on_open, None),
            ("save", self.on_save, None),
            ("save-as", self.on_save_as, None),
            ("help", self.on_help, None),
            ("close", self.on_close, None),
            ("quit", self.on_quit, None),
            ("compute-and-redraw", self.compute_and_redraw, None),
        ])

    def build_shortcuts(self):
        """Create keyboard shortcuts.

        Creates keyboard shortcuts for these actions:
        - New window: Ctrl+N,
        - Close window: Ctrl+W,
        - Close all windows and quit: Ctrl+Q,
        - Open file: Ctrl+O,
        - Save file: Ctrl+S,
        - Save file as: Ctrl+Shift+S,
        - Help: F1,
        - Compute and redraw: Ctrl+R.

        """
        controller = Gtk.ShortcutController()
        self.add_controller(controller)

        controller.add_shortcut(
            Gtk.Shortcut(
                trigger=Gtk.KeyvalTrigger(
                    keyval=Gdk.KEY_n,
                    modifiers=Gdk.ModifierType.CONTROL_MASK
                ),
                action=Gtk.NamedAction.new("win.new")
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
                    keyval=Gdk.KEY_q,
                    modifiers=Gdk.ModifierType.CONTROL_MASK
                ),
                action=Gtk.NamedAction.new("win.quit")
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
                    keyval=Gdk.KEY_r,
                    modifiers=Gdk.ModifierType.CONTROL_MASK
                ),
                action=Gtk.NamedAction.new("win.compute-and-redraw")
            )
        )

    def build_body(self):
        """Build the main window body layout.

        Creates the main layout with:
        - A Paned window (horizontal) with:
          * Left: Viewer2D for 2D visualization
          * Right: Property editor and controls
        - Control widgets at the bottom

        """
        #### Main container: vertical box
        self.vertbox1 = Gtk.Box(
            orientation=Gtk.Orientation.VERTICAL,
            spacing=6
        )
        self.set_child(self.vertbox1)

        ### First row: paned window (horizontal split)
        self.paned = Gtk.Paned(
            orientation=Gtk.Orientation.HORIZONTAL
        )
        self.vertbox1.append(self.paned)

        ## First row, left: Viewer2D (2D Matplotlib canvas)
        self.left = Gtk.Box(
            orientation=Gtk.Orientation.VERTICAL,
            spacing=6
        )
        self.paned.set_start_child(self.left)

        # drop-down to select the quantioty to display
        self.quantity_dropdown = Gtk.DropDown.new_from_strings(self.scalar_quantities)
        self.quantity_dropdown.connect("notify::selected",
                                       self.on_quantity_dropdown_activated)
        self.left.append(self.quantity_dropdown)

        # actual plot widget
        self.viewer = Viewer2D()
        self.left.append(self.viewer)
        # Forbid focus to the plot, else shortcuts don't work reliably
        self.viewer.canvas.set_focusable(False)

        ## First row, right: vertical box for property editor
        self.right = Gtk.Box(
            orientation=Gtk.Orientation.VERTICAL,
            spacing=6
        )
        self.paned.set_end_child(self.right)

        # Scrollable property editor view
        scroll = Gtk.ScrolledWindow()
        scroll.set_hexpand(True)
        scroll.set_vexpand(True)
        self.right.append(scroll)
        self.editor_scroller = scroll
        self.editor_scroller.set_min_content_width(365)

        self.paned.set_resize_start_child(True)
        self.paned.set_shrink_start_child(False)
        self.paned.set_resize_end_child(False)
        self.paned.set_shrink_end_child(False)

        # Controls for running the integration (play/pause/stop)
        self.controls = SimulationControls()
        self.vertbox1.append(self.controls)
        self.controls.connect("play-pause", self.on_play_pause)
        self.controls.connect("reset", self.on_reset)
        self.controls.connect("stop", self.on_stop)

        # Set default number of frames
        self.controls.nframes.set_value(100)

    def on_close_request(self, *args):
        """Handle window close request.

        Cleanly shuts down the worker process, removes the window from
        the application's window list, and allows the window to close.

        Args:
            *args: GTK callback arguments

        Returns:
            bool: False to allow the window to close

        """
        self.cmd_queue.put((QUIT,))
        if self.worker:
            self.worker.join(timeout=2.0)
            if self.worker.is_alive():
                self.worker.terminate()
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

    def on_new(self, *args):
        """Create a new window with a default Scenery.

        Args:
            *args: GTK callback arguments

        """
        new_window = GyotoSceneryViewerApplicationWindow(
            application=self.props.application,
            scenery=self.default_scenery()
        )
        new_window.present()

    ####################################################################
    # Compute and redraw
    ####################################################################

    def compute_and_redraw(self, *args):
        """Raytrace the scenery and show the image.

        This method starts a new simulation in the worker process and
        updates the UI to reflect the running state. The actual
        computation happens in the background worker process.

        Args:
            *args: Ignored (for callback compatibility)

        """

        # Disable UI and save focus
        self.simulation_running = True
        self.controls.set_progress(0.)
        self.controls.set_running(True)
        self.controls.set_status("Integrating...")
        self.last_focused_widget = self.get_focus()
        self.right.set_sensitive(False)

        # Determine region to compute
        xlim = self.viewer.axes.get_xlim()
        ylim = self.viewer.axes.get_ylim()
        ilim, jlim = self.angles_to_pixel(xlim, ylim)

        # Send simulation command to worker
        self.cmd_queue.put((
            RUN_SIM,
            self.scenery, ilim, jlim,
            self.controls.nframes.get_value_as_int(),
        ))

    def redraw(self):
        """Draw the selected quantity from precomputed data"""
        quantity = self.quantity_dropdown.get_selected_item().get_string()
        xlim = self.viewer.axes.get_xlim()
        ylim = self.viewer.axes.get_ylim()
        self.viewer.axes.clear()

        # first determine vmin and vmax
        vmin = numpy.inf
        vmax = -numpy.inf
        for dico in self.data.values():
            temp = min(vmin, numpy.nanmin(dico[quantity]))
            if not numpy.isnan(temp):
                vmin = temp
            temp = max(vmax, numpy.nanmax(dico[quantity]))
            if not numpy.isnan(temp):
                vmax = temp

        # then show each image, starting with lower resolutions and
        # larger fields:
        for (delta, semifov), dico in sorted(self.data.items(), reverse=True):
            self.viewer.axes.imshow(
                dico[quantity],
                origin='lower',
                extent=(semifov, -semifov, -semifov, semifov),
                vmin=vmin, vmax=vmax,
            )
        self.viewer.axes.set_xlim(xlim)
        self.viewer.axes.set_ylim(ylim)
        self.viewer.canvas.draw_idle()

    def reset_limits(self, *args):
        """Reset the limits

        To the largest value needed for already computed images and
        the current Screen.

        """
        semifov = 0.

        if len(self.data):
            semifov = max(semifov for delta, semifov in self.data)

        if self.scenery.Screen:
            res = self.scenery.Screen.Resolution
            fov = self.scenery.Screen.fieldOfView("arcsec")
            delta = fov / res
            semifov2 = fov/2 + delta/2
            if semifov2 > semifov:
                semifov = semifov2

        if semifov:
            self.viewer.axes.set_xlim((semifov, -semifov))
            self.viewer.axes.set_ylim((-semifov, semifov))
            self.viewer.axes.set_aspect('equal', adjustable='box')

    def process_progress(self):
        """Process progress updates from the worker.

        This method is called periodically (every 50ms) to check for
        progress updates from the worker process. It processes all
        pending progress messages but only displays the latest one to
        avoid visual clutter.

        Returns:
            bool: True to keep the timeout source alive

        """
        msg = None
        while True:
            try:
                msg = self.progress_queue.get_nowait()
            except queue.Empty:
                break
        if msg:
            if len(msg) >= 3:
                # store computed dict is the data dict
                # with key (delta, semifov)
                res = self.scenery.Screen.Resolution
                fov = self.scenery.Screen.fieldOfView("arcsec")
                delta = fov / res
                semifov = fov/2 + delta/2
                key = delta, semifov
                if key in self.data:
                    for quantity, d in msg[2].items():
                        try:
                            # copy non-nan values into existing arrays
                            numpy.copyto(self.data[key][quantity], d,
                                         where=~numpy.isnan(d))
                        except KeyError:
                            # this quantity is not yet available
                            self.data[key] = d
                else:
                    self.data[key] = msg[2]
                #redraw
                self.redraw()
            self.controls.set_progress(msg[1])
        return True

    def process_control(self):
        """Process control messages from the worker.

        This method is called periodically (every 50ms) to check for
        control messages from the worker process (log messages,
        completion, errors). It handles different message types and
        updates the UI accordingly.

        Returns:
            bool: True to keep the timeout source alive

        """
        try:
            msg = self.control_queue.get_nowait()
            if msg[0] == 'log':
                print(msg[1])
            elif msg[0] == 'done':
                self.computation_epilogue(msg="Computation finished.")
            elif msg[0] == 'aborted':
                self.computation_epilogue(msg="Computation aborted.")
            elif msg[0] == 'error':
                self.computation_epilogue(
                    msg="Computation ended in error.",
                    error=msg[1]
                )
            elif msg[0] == 'fatal_error':
                self.computation_epilogue(
                    msg="Fatal error. Restarting worker.",
                    error=msg[1]
                )
                self.restart_worker()
        except queue.Empty:
            pass
        return True

    def computation_epilogue(self, msg="Integration finished", error=None):
        """Handle the end of a computation.

        This method is called when a simulation completes
        (successfully or not). It re-enables the UI, restores focus,
        and updates the status message.

        Args:
            msg: Status message to display
            error: Optional error message to include in tooltip

        """
        self.controls.set_running(False)
        self.controls.set_status(msg, error)
        self.right.set_sensitive(True)
        if self.last_focused_widget:
            try:
                self.last_focused_widget.grab_focus()
            except Exception:
                pass
        self.last_focused_widget = None
        self.simulation_running = False

    def restart_worker(self):
        """Restart the worker process if it died or encountered a fatal error.

        This method cleanly shuts down the existing worker (if any) and
        starts a new one. It's used for error recovery.

        """
        if self.worker and self.worker.is_alive():
            self.cmd_queue.put((QUIT,))
            self.worker.join(timeout=1.0)
            if self.worker.is_alive():
                self.worker.terminate()
        self.worker = Process(
            target=worker_func,
            args=(
                self.cmd_queue,
                self.progress_queue,
                self.control_queue,
                self.pause_event,
                self.stop_event
            ),
            daemon=True
        )
        self.worker.start()
        self.last_heartbeat = time.time()
        self.simulation_running = False

    ####################################################################
    # Callbacks
    ####################################################################

    def on_open(self, *args):
        """Open a file dialog to load an XML scenery configuration.

        Creates a file dialog with XML and All Files filters,
        defaulting to XML. When a file is selected, it's loaded and
        the scenery is updated.

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
            window = GyotoSceneryViewerApplicationWindow(
                application=self.props.application,
                scenery=file.get_path(),
                connector=None
            )
            window.present()

    def on_save(self, *args):
        """Save the current scenery in the last XML file used.

        If the scenery was not read from file and not yet saved to
        file, opens a dialog.

        Args:
            *args: GTK callback arguments

        """
        if self.filename:
            try:
                Factory(self.scenery).write(self.filename)
            except GyotoError as e:
                show_error_dialog(
                    message=f"Error writing XML file {self.filename}:",
                    detail=e.get_message(),
                    window=self
                )
        else:
            self.on_save_as(*args)

    def on_save_as(self, *args):
        """Open a file dialog to save the current scenery as XML.

        Creates a file dialog with XML and All Files filters,
        defaulting to XML. When a file is selected, the current
        scenery is saved to that file.

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
            Factory(self.scenery).write(file.get_path())
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
            "Gyotoy - Gyoto Geodesic Integration Visualizer\n"
            "\n"
            "OVERVIEW:\n"
            "Gyotoy is a GTK4 application for simulating and visualizing\n"
            "geodesics (time-like or null) in spacetimes supported by the\n"
            "Gyoto library. It provides an interactive 3D view of particle\n"
            "trajectories. Each window handles exactly one particle (Star\n"
            "or Photon).\n"
            "\n"
            "UI LAYOUT:\n"
            "- Left Panel: 3D Matplotlib viewer displaying particle\n"
            "  trajectory. Use mouse to rotate (left-click+drag), and\n"
            "  pan (middle-click+drag). Use plot toolbar for other actions\n"
            "  including saving the plot.\n"
            "- Right Panel: Property editor for adjusting metric and\n"
            "  particle parameters. Changes trigger automatic\n"
            "  recomputation.\n"
            "- Bottom: Simulation controls (play/pause/stop/reset),\n"
            "  interpolation settings, and status display.\n"
            "\n"
            "MENU BUTTONS:\n"
            "- New -> Star: Open new window with default Star (Ctrl+N)\n"
            "- New -> Photon: Open new window with default Photon\n"
            "  (Ctrl+Shift+N)\n"
            "- Open (Ctrl+O): Load an XML particle configuration file.\n"
            "- Save (Ctrl+S): Save current particle to last used file.\n"
            "- Save As (Ctrl+Shift+S): Save current particle to a new file.\n"
            "- Help (F1): Show this help dialog.\n"
            "- Close (Ctrl+W): Close current window.\n"
            "- Quit (Ctrl+Q): Close all windows and quit the application.\n"
            "\n"
            "OTHER KEYBOARD SHORTCUTS:\n"
            "- Ctrl+R: Compute and redraw trajectory.\n"
            "- Escape: Close active dialog window (error, help...).\n"
            "\n"
            "PROPERTY EDITOR:\n"
            "- Edit particle parameters: position, velocity, metric\n"
            "  properties...\n"
            "- All parameter changes immediately trigger a computation\n"
            "  and redraw unless the stop (■) at the bottom right of the\n"
            "  window is activated.\n"
            "- The InitCoord vector is by default displayed with 7 cells\n"
            "  corresponding to 4-position and 3-velocity (derivatives of\n"
            "  the space coordinates with respect to time coordinate).\n"
            "  Use the '3-velocity' and '4-velocity' radio buttons to\n"
            "  switch between this view and the 8-coordinate view\n"
            "  corresponding to 4-position and 4-velocity. By default, the\n"
            "  4-velocity is renormalized (according to mass of the\n"
            "  particle) at each keystroke. Click on the stop (■) button\n"
            "  immediately below the coordinates to temporarily inhibit\n"
            "  this behavior. Click this button again to finalize input.\n"
            "\n"
            "SIMULATION CONTROLS:\n"
            "- Reset: Reinitialize integration and clear viewer.\n"
            "- Play/Pause: Start or pause the simulation.\n"
            "- Stop: Stop the integration and inhibit/enable\n"
            "  recomputation.\n"
            "- N. frames: Number of intermediate frames to display\n"
            "  (default: 100).\n"
            "- End time: Final time for integration.\n"
            "- Interpolation step: Step size for interpolation (0 for no\n"
            "  interpolation, adaptive step used instead).\n"
            "\n"
            "WORKFLOW:\n"
            "1. If needed, open file or create new particle (Star or\n"
            "   Photon).\n"
            "2. Adjust properties in the editor.\n"
            "3. Click Play or press Ctrl+R to compute trajectory.\n"
            "4. Use 3D viewer to inspect the result.\n"
            "5. Save particle description with Save/Save As.\n"
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

    def on_quantity_dropdown_activated(self, widget, *args):
        """Handle dropdown selection changes.

        Show the selected quantity..

        Args:
            widget: The Gtk.DropDown that changed
            *args: Additional arguments
        """
        self.redraw()

    def on_recursive_value_changed(self, widget, path):
        """Handle property changes"""
        # placeholder: currently do nothing
        if path in ("Screen",
                    "Screen.Resolution",
                    "Screen.FieldOfView"):
            pass

    def on_reset(self, wdgt):
        """Handle reset button click.

        Deletes the computed images, clears the viewer and reset axes.

        Args:
            wdgt: The button that was clicked

        """
        self.data = {}
        self.viewer.axes.clear()
        self.reset_limits()
        self.viewer.canvas.draw_idle()

    def on_play_pause(self, wdgt):
        """Handle play/pause button click.

        Toggles the pause state and updates the status.

        Args:
            wdgt: The SimulationControls widget

        """
        if self.controls.running:
            self.pause_event.clear()
            self.stop_event.clear()
            wdgt.stop_button.set_active(False)
            self.controls.set_status("Integration resumed...")
            if not self.simulation_running:
                self.compute_and_redraw()
        else:
            self.pause_event.set()
            self.controls.set_status("Integration paused...")

    def on_stop(self, wdgt):
        """Handle stop button click.

        Toggles the hold state and updates the status.

        Args:
            wdgt: The SimulationControls widget

        """
        hold = wdgt.stop_button.get_active()
        if hold:
            self.stop_event.set()
            self.controls.set_status(
                "Holding integration (press play)."
            )
        else:
            self.stop_event.clear()
            self.controls.set_status("Ready for integration.")

    ####################################################################
    # Setters / getters
    ####################################################################

    def set_scenery(self, scenery):
        """Set the current scenery and update the UI.

        Args:
            scenery: The Gyoto Scenery to set.

        """
        if not isinstance(scenery, Scenery):
            show_error_dialog(
                message="Wrong type for scenery:",
                detail=repr(type(scenery)),
                window=self
            )
            return

        self.scenery = scenery

        if self.editor is not None:
            self.editor.destroy()

        self.editor = PropertyEditorBox(
            scenery,
        )

        self.editor.connect("recursive-value-changed",
                            self.on_recursive_value_changed)

        self.editor_scroller.set_child(self.editor)

    def get_scenery(self):
        """Get the current scenery.

        Returns:
            The current scenery
        """
        return self.scenery

    ####################################################################
    # Default values
    ####################################################################

    def default_scenery(self):
        """Create a default Scenery.

        Returns:
            A Scenery
        """
        sc=Scenery(
            Metric   = Minkowski(
                Mass = (4e6, "sunmass"),
                Spherical = True,
            ),
            Astrobj  = FixedStar(
                Radius = 12,
                Position = (0., 0., 0.),
                Spectrum = PowerLaw(
                    Exponent = 0.,
                    Constant = 0.01,
                ),
                Opacity = PowerLaw(
                    Exponent = 0.,
                    Constant = 0.01,
                ),
                OpticallyThin = True,
            ),
            Screen   = Screen(
                Distance    = (8, "kpc"),
                Time        = (30e3, "yr"),
                Resolution  = 64,
                Inclination = (90., "degree"),
                PALN        = 0.,
                FieldOfView = (150, "µas"),
            ),
            Delta    = 1e0,
            MinimumTime = 0.,
            Quantities = " ".join(self.scalar_quantities),
            Adaptive = True,
            NThreads = 8,
        )
        return sc

    ####################################################################
    # Helpers
    ####################################################################

    def pixel_to_angles(self,
                        i: ArrayLike,
                        j: ArrayLike,
                        unit: str='arcsec',
                        resolution: int=None,
                        fieldofview: float=None
                        ) -> tuple[float | numpy.ndarray,
                                   float | numpy.ndarray]:
        '''Transform pixel coordinates (i, j) to angle coordinates (x, y)

        Here, i and j are in Python convention (starting at
        0). Integer values of pixel coordinates correspond to pixel
        center.

        The arguments may be scalars or array-like objects. Scalar
        inputs produce scalar outputs; otherwise outputs are NumPy
        arrays.

        Keywords:
        unit: unit for angles (x, y and fieldofview)
        resolution: number of pixels in each direction of the image
            (which is always square)
        fieldofview: angular extent of the image (from center of first
            pixel to center of last pixel).

        '''
        if fieldofview is None:
            fieldofview =  self.scenery.Screen.fieldOfView(unit)
        if resolution is None:
            resolution = self.scenery.Screen.Resolution

        delta = fieldofview / resolution

        i0 = 0.5 * resolution - 0.5
        j0 = 0.5 * resolution - 0.5

        x = - (i-i0) * delta
        y = (j-j0) * delta

        return x, y

    def angles_to_pixel(self,
                        x: ArrayLike,
                        y: ArrayLike,
                        unit: str='arcsec',
                        resolution: int=None,
                        fieldofview: float=None
                        ) -> tuple[float | numpy.ndarray,
                                   float | numpy.ndarray]:
        '''Transform angle coordinates (x, y) to pixel coordinates (i, j)

        Here, i and j are in Python convention (starting at
        0). Integer values of pixel coordinates correspond to pixel
        center.

        Keywords:
        unit: unit for angles (x, y and fieldofview)
        resolution: number of pixels in each direction of the image
            (which is always square)
        fieldofview: angular extent of the image (from center of first
            pixel to center of last pixel).

        '''
        if fieldofview is None:
            fieldofview =  self.scenery.Screen.fieldOfView(unit)
        if resolution is None:
            resolution = self.scenery.Screen.Resolution

        delta = fieldofview / resolution

        i0 = 0.5 * resolution - 0.5
        j0 = 0.5 * resolution - 0.5

        i = - numpy.asarray(x) / delta + i0
        j = numpy.asarray(y) / delta + j0

        if numpy.isscalar(x): i = float(i)
        if numpy.isscalar(y): j = float(j)

        return i, j


# Stand-alone entry point:
if __name__ == "__main__":
    raise SystemExit(GyotoSceneryViewerApplication.run_app(parsecliargs=True))
