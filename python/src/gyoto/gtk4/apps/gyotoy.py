# Imports + MyApplication + GyotoyApplicationWindow skeleton

"""Gyotoy: GTK4 Application for Gyoto Geodesic Integration

This application provides a graphical interface for simulating and
visualizing geodesics (time-like or null) in spacetimes supported by
the Gyoto library.

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
    │ ████████████████████████████████████────────────────────────────── │
    │                                                                    │
    │ Status...                    N: [100  ]      ⏮        ▶        ⏹   │
    │                                                                    │
    └────────────────────────────────────────────────────────────────────┘

Description
-----------
- **Left Panel**: 3D Matplotlib viewer displaying the particle
    trajectory.
- **Right Panel**: Property editor for adjusting metric and particle
    parameters.
- **Bottom Controls**: Play/pause/stop buttons, interpolation
    settings, and status display.
- **Worker Process**: Heavy computations run in a background process
    to keep the UI responsive.

Usage
-----
Run as a standalone application:
    python3 -m gyoto.gtk4.apps.gyotoy [-h] [filename.xml]

Or import and use programmatically:
    from gyoto.gtk4.apps.gyotoy import gyotoy
    gyotoy([particle|'filename.xml'])

An optional particle (gyoto.std.Star or gyoto.core.Photon), or the
name of an XML file describing such a particle, can be provided.

"""

from __future__ import annotations

__all__ = ['GyotoyApplication', 'GyotoyApplicationWindow', 'gyotoy']

import gi
gi.require_version("Gtk", "4.0")
from gi.repository import Gtk, Gio, GLib, Gdk

import sys
import argparse
import numpy
from multiprocessing import Queue, Event, Process
import time
import queue
import traceback
import warnings

from gettext import gettext as _

# Custom widgets
from ..widgets.property_editor_box import PropertyEditorBox
from ..widgets.scientific_spin import ScientificSpin
from ..widgets.viewer_3d import Viewer3D
from ..widgets.simulation_controls import SimulationControls
from ..utils import show_error_dialog
from ...core import Factory, Astrobj, Error as GyotoError, Photon
from ...core import GYOTO_COORDKIND_SPHERICAL
from ...std import Star, KerrBL

# --- Commands for worker communication ---
RUN_SIM = 'run'
QUIT = 'quit'

def worker_func(cmd_queue, progress_queue, control_queue, pause_event, stop_event):
    """Persistent worker process for running simulations.

    This function runs in a separate process to avoid blocking the GTK
    main thread.  It waits for commands from cmd_queue, executes
    simulations, and sends progress updates to progress_queue and
    control messages to control_queue.

    The worker can be paused via pause_event and stopped via stop_event.

    Args:
        cmd_queue: Queue for receiving commands (RUN_SIM, QUIT)
        progress_queue: Queue for sending progress updates
        control_queue: Queue for sending control messages (log, done, error)
        pause_event: Event to pause/resume the simulation
        stop_event: Event to stop the simulation

    """
    try:
        import numpy
        import time
        from ...core import Factory
        from ...std import Star

        next_update_fraction = 0.
        last_update = time.time()

        while True:
            try:
                cmd = cmd_queue.get(timeout=1.0)  # Blocks until command arrives or timeout
            except queue.Empty:
                continue # No command, continue to next iteration

            if cmd[0] == QUIT:
                break  # Exit loop → process terminates

            elif cmd[0] == RUN_SIM:
                end_msg = ('done',)
                try:
                    _, particlexml, starttime, endtime, nframes, interp_step = cmd
                    stop_event.clear()
                    pause_event.clear()  # Start fresh

                    # Rebuild objects from serialized data
                    f = Factory(particlexml)
                    particle = f.photon() if f.kind() == 'Photon' else Star(f.astrobj())

                    # Ensure interp_step has the same sign as endtime-starttime
                    interp_step = numpy.sign(endtime-starttime)*abs(interp_step)

                    # Compute frame and interpolation dates
                    frametimes = numpy.linspace(starttime, endtime+interp_step, nframes + 1)

                    # If interp_step is not 0, will interpolate
                    if interp_step:
                        t = numpy.arange(starttime,
                                         endtime+interp_step,
                                         interp_step)
                        x, y, z = [numpy.full(len(t), numpy.nan, like=t)
                                   for _ in range(3)]

                    for n in range(len(frametimes) - 1):
                        if stop_event.is_set():
                            end_msg = ('aborted',)
                            break
                        while pause_event.is_set() and not stop_event.is_set():
                            time.sleep(0.1)

                        frametime = frametimes[n + 1]

                        # Actually integrate using Gyoto with built-in adaptive step
                        particle.xFill(frametime)

                        if interp_step:
                            # Extract positions at interpolated dates

                            # First check whether we could reach frametimes[n+1]:
                            # (impossible e.g. if geodesic joined the event horizon)
                            npoints = particle.get_nelements()
                            tinteg = numpy.empty(npoints)
                            particle.get_t(tinteg)

                            if interp_step > 0 :
                                frameend = numpy.min((tinteg[-1], frametime))
                                mask = (t >= frametimes[n]) * (t < frameend)
                            else:
                                frameend = numpy.max((tinteg[-1], frametime))
                                mask = (t <= frametimes[n]) * (t > frameend)
                            ttmp = t[mask]
                            if len(ttmp) == 0:
                                break
                            xtmp, ytmp, ztmp = [numpy.empty(len(ttmp), like=ttmp)
                                                for _ in range(3)]
                            particle.getCartesian(ttmp, xtmp, ytmp, ztmp)
                            x[mask], y[mask], z[mask] = xtmp, ytmp, ztmp
                            if (numpy.any(numpy.isnan(x[mask])) or
                                numpy.any(numpy.isnan(x[mask])) or
                                numpy.any(numpy.isnan(x[mask]))):
                                raise Exception('should not be NaN')
                        else:
                            # Extract positions at adaptive-step dates
                            npoints = particle.get_nelements()
                            t = numpy.empty(npoints)
                            particle.get_t(t)
                            x, y, z = [numpy.empty(npoints, like=t) for _ in range(3)]
                            particle.getCartesian(t, x, y, z)

                        progress = (frametime - starttime) / (endtime - starttime)
                        if progress >= next_update_fraction and time.time() - last_update > 0.1:
                            progress_queue.put_nowait(('progress', progress,
                                       x.tolist(), y.tolist(), z.tolist()))
                            next_update_fraction = progress + 0.05
                            last_update=time.time()

                    progress_queue.put_nowait(('progress', progress, x.tolist(), y.tolist(), z.tolist()))
                    control_queue.put(end_msg)
                except Exception as e:
                    control_queue.put(('error', traceback.format_exc()))
    except Exception as e:
        control_queue.put(('fatal_error', traceback.format_exc()))

class GyotoyApplication(Gtk.Application):
    """Standalone GTK application for Gyotoy.

    This class handles the application lifecycle and window
    management.  It extends Gtk.Application to provide a proper GTK
    application structure.

    """

    def __init__(self, particle=None, connector=None, *args, **kwargs):
        """Initialize the Gyotoy GTK application.

        Args:
            particle: Initial particle to display (Star or Photon)
            connector (multiprocessing.Connection or None): If the GUI
                runs in a separate process, this is used to send updates
                back to the caller.

        """
        if 'application_id' not in kwargs:
            kwargs['application_id'] = "fr.obspm.gyoto.Gyotoy"
        if 'flags' not in kwargs:
            kwargs['flags'] = (Gio.ApplicationFlags.DEFAULT_FLAGS |
                               Gio.ApplicationFlags.NON_UNIQUE)
        super().__init__(*args, **kwargs)
        self.particle = particle
        self.connector = connector

    def do_activate(self):
        """Called by GTK when the application starts.

        Creates the main window if it doesn't exist and presents it.
        """
        window = self.props.active_window

        if window is None:
            window = GyotoyApplicationWindow(application=self,
                                             particle=self.particle,
                                             connector=self.connector)

        window.present()

    @staticmethod
    def run_app(particle=None, parsecliargs=False, *args, **kwargs):
        """Run Gyotoy as a standalone GTK application.

        Parameters:
            particle: the particle to start with (Star or Photon), or
                None, or the XML description of such a particle, or
                the name of an XML file containing this description.
            parsecliargs (bool): whether to parse the command line
                arguments
            *args, **kwargs: other parameters are passed untouched to
                the GyotoApplication constructor.

        Returns:
            int: Application exit code

        """
        remaining = None
        if parsecliargs:
            parser = argparse.ArgumentParser(prog=f'{sys.argv[0]} ',
                                             description=__doc__,
                                             formatter_class=argparse.RawTextHelpFormatter)
            parser.add_argument('xmlfile', nargs='?',
                                help='XML file containing the description'
                                + 'of a Gyoto Star or Photon (optional)')
            cliargs, remaining = parser.parse_known_args()
            if 'xmlfile' in cliargs:
                particle=cliargs.xmlfile

        app = GyotoyApplication(particle=particle, *args, **kwargs)
        return app.run(remaining)

class GyotoyApplicationWindow(Gtk.ApplicationWindow):
    """Main application window for Gyotoy.

    This window contains:
    - A 3D Matplotlib viewer (left panel) for visualizing particle
      trajectories
    - A property editor (right panel) for adjusting particle and
      metric parameters
    - Control widgets (bottom) for running/stopping simulations and
      setting parameters

    The window uses a multiprocessing worker for heavy computations to
    keep the UI responsive. Communication with the worker happens via
    multiprocessing Queues for progress updates and control messages.

    Attributes:
        star: Default Star particle
        photon: Default Photon particle
        particle: Current particle being edited/simulated
        viewer: Viewer3D widget for 3D visualization
        editor: PropertyEditorBox for editing particle properties
        controls: SimulationControls for play/pause/stop
        worker: Process for background computations
        simulation_running: Flag indicating if a simulation is in progress
        last_focused_widget: Last widget that had focus (for focus
            restoration)
        interpolation_step: Step size for interpolation
        connector: the connector that was passed to the constructor
        connectors: dict containing either None or self.connector for
            each particle
        filenames: dict with the name of the file each particle was
            last read from or written to

    """

    # Default values
    blocking = True
    main_loop = None
    particle = None
    star = None
    photon = None
    endtime = 3000
    hold = True
    worker = None
    simulation_running = False
    last_focused_widget = None
    interpolation_step = 1.
    connector = None
    connectors = None
    filenames = None

    ####################################################################
    # Construction
    ####################################################################

    def __init__(self, application=None,
                 particle=None, star=None, photon=None,
                 connector=None):
        """Initialize the main window.

        Args:
            application: Parent Gtk.Application instance
            particle: Initial particle to display (Star or Photon)
            star: Default Star particle (created if None)
            photon: Default Photon particle (created if None)
            connector (multiprocessing.Connection or None): If the GUI
                runs in a separate process, this is used to send updates
                back to the caller.
        """
        super().__init__(application=application)

        self.connectors = {}
        self.filenames = {}
        filename = None

        # handle QUIT from parent process
        self.connector = connector
        if connector is not None:
            GLib.timeout_add(50, self.check_connector)

        # Handle case where particle is a string.
        # It is either a filename or a full XML description string
        if isinstance(particle, str):
            if particle.lower().endswith('.xml'):
                filename = particle
            factory = Factory(particle)
            particle = getattr(factory, factory.kind().lower())()

        # process particle, star and photon
        if isinstance(particle, Astrobj):
            particle = Star(particle)

        # now that particle is one of None, a Star, a Photon
        # set its name (may be None)
        self.filenames[particle] = filename

        if isinstance(star, Astrobj):
            star = Star(star)

        if star is None:
            if isinstance(particle, Star):
                self.star = particle
            else:
                self.star = self.default_star()

        if photon is None:
            if isinstance(particle, Photon):
                self.photon = particle
            else:
                self.photon = self.default_photon()

        if particle is None:
            if star is not None:
                particle = star
                self.photon.Metric = particle.Metric
            elif photon is not None:
                particle = photon
                self.star.Metric = particle.Metric
            else:
                particle = self.star

        self.set_title("Gyotoy")
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
            args=(self.cmd_queue, self.progress_queue, self.control_queue,
                  self.pause_event, self.stop_event),
            daemon=True
        )
        self.worker.start()

        # Poll results every 50ms
        GLib.timeout_add(50, self.process_progress)
        GLib.timeout_add(50, self.process_control)

        # Populate initial editor
        self.hold = False

        # one of star or photon is particle.  In the end,
        # connectors has two items, one of which may not be None.
        self.connectors[particle] = connector
        for p in self.star, self.photon:
            if p not in self.connectors:
                self.connectors[p] = None

        # also create or fill the filenames dict. One of the two may
        # have been provided.
        for p in self.star, self.photon:
            if p not in self.filenames:
                self.filenames[p] = None

        # actually set self.particle and build editor
        self.set_particle(particle)

    ####################################################################
    # UI
    ####################################################################

    def build_headerbar(self):
        """Build the window's header bar with menu button.

        Creates a header bar with a hamburger menu containing:
        - Open...
        - Save
        - Save As...
        - Help
        - Quit

        """

        # Create title bar with hamburger button
        header = Gtk.HeaderBar()
        self.set_titlebar(header)
        menu_button = Gtk.MenuButton(
            icon_name="open-menu-symbolic"
        )
        menu_button.add_css_class("flat")
        header.pack_end(menu_button)

        # Attach menu to hamburger button
        box = Gtk.Box(
            orientation=Gtk.Orientation.VERTICAL
        )
        popover = Gtk.Popover()
        popover.set_child(box)
        menu_button.set_popover(popover)

        # Populate menu with Open, Save As, and Quit buttons
        open_button = Gtk.Button(
            label=_("Open…")
        )
        open_button.add_css_class("flat")

        save_button = Gtk.Button(
            label=_("Save")
        )
        save_button.add_css_class("flat")

        save_as_button = Gtk.Button(
            label=_("Save As…")
        )
        save_as_button.add_css_class("flat")


        help_button = Gtk.Button(
            label=_("Help")
        )
        help_button.add_css_class("flat")

        quit_button = Gtk.Button(
            label=_("Quit")
        )
        quit_button.add_css_class("flat")

        box.append(open_button)
        box.append(save_button)
        box.append(save_as_button)
        box.append(help_button)
        box.append(Gtk.Separator())
        box.append(quit_button)

        open_button.connect("clicked", self.on_open)
        save_button.connect("clicked", self.on_save)
        save_as_button.connect("clicked", self.on_save_as)
        help_button.connect("clicked", self.on_help)
        quit_button.connect("clicked", self.on_quit)


    def build_shortcuts(self):
        '''Create keyboard shortcuts

        Creates keyboard shortcuts for these actions:
        - Open file: Ctrl-O,
        - Save file: Ctrl-S,
        - Save file as: Ctrl-Shift-S,
        - Help: F1,
        - Quit: Ctrl-Q,
        - Compute and redraw: Ctrl-R.

        '''
        self.add_shortcut(
            Gtk.Shortcut(
                trigger=Gtk.KeyvalTrigger(keyval=Gdk.KEY_o,
                                          modifiers=Gdk.ModifierType.CONTROL_MASK),
                action=Gtk.CallbackAction.new(self.on_open)
            )
        )

        self.add_shortcut(
            Gtk.Shortcut(
                trigger=Gtk.KeyvalTrigger(keyval=Gdk.KEY_s,
                                          modifiers=Gdk.ModifierType.CONTROL_MASK),
                action=Gtk.CallbackAction.new(self.on_save)
            )
        )

        self.add_shortcut(
            Gtk.Shortcut(
                trigger=Gtk.KeyvalTrigger(keyval=Gdk.KEY_s,
                                          modifiers=(Gdk.ModifierType.CONTROL_MASK |
                                                     Gdk.ModifierType.SHIFT_MASK)),
                action=Gtk.CallbackAction.new(self.on_save_as)
            )
        )

        self.add_shortcut(
            Gtk.Shortcut(
                trigger=Gtk.KeyvalTrigger(keyval=Gdk.KEY_q,
                                          modifiers=Gdk.ModifierType.CONTROL_MASK),
                action=Gtk.CallbackAction.new(self.on_quit)
            )
        )

        self.add_shortcut(
            Gtk.Shortcut(
                trigger=Gtk.KeyvalTrigger(keyval=Gdk.KEY_F1,
                                          modifiers=0),
                action=Gtk.CallbackAction.new(self.on_help)
            )
        )

        self.add_shortcut(
            Gtk.Shortcut(
                trigger=Gtk.KeyvalTrigger(keyval=Gdk.KEY_r,
                                          modifiers=Gdk.ModifierType.CONTROL_MASK),
                action=Gtk.CallbackAction.new(self.compute_and_redraw)
            )
        )

    def build_body(self):
        """Build the main window body layout.

        Creates the main layout with:
        - A Paned window (horizontal) with:
          * Left: Viewer3D for 3D visualization
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

        ## First row, left: Viewer3D (3D Matplotlib canvas)
        self.viewer = Viewer3D()
        self.paned.set_start_child(
            self.viewer
        )
        # forbid focus to the plot, else shortcuts don't work reliably
        self.viewer.canvas.set_focusable(False)


        ## First row, right: vertical box for property editor
        self.right = Gtk.Box(
            orientation=Gtk.Orientation.VERTICAL,
            spacing=6
        )
        self.paned.set_end_child(self.right)

        # Top of control column: Star/Photon radio buttons
        frame = Gtk.Frame(
            label="Particle type"
        )
        self.right.append(frame)

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
        # Group radio buttons so only one can be selected
        self.particle_photon.set_group(
            self.particle_star
        )
        buttons.append(self.particle_star)
        buttons.append(self.particle_photon)

        self.particle_star.connect(
            "toggled",
            self.on_star_selected
        )

        # Scrollable property editor view
        scroll = Gtk.ScrolledWindow()
        scroll.set_hexpand(True)
        scroll.set_vexpand(True)
        self.right.append(scroll)
        self.editor_scroller = scroll
        self.editor_scroller.set_min_content_width(365)

        # End time spin button
        frame = Gtk.Frame()
        frame.set_label("End time")
        frame.set_label_align(0.0)
        frame.set_tooltip_text("Run simulation from initial condition time to end time")

        hbox = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=5)
        frame.set_child(hbox)

        spin = ScientificSpin(value=0.)
        spin.set_value(self.endtime)
        spin.connect("value-changed", self.on_endtime_changed)
        hbox.append(spin)
        self.right.append(frame)

        # Interpolation step spin button
        frame = Gtk.Frame()
        frame.set_label("Interpolation step")
        frame.set_label_align(0.0)
        frame.set_tooltip_text("Interpolation step (0 for no interpolation)")

        hbox = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=5)
        frame.set_child(hbox)

        spin = ScientificSpin(value=self.interpolation_step)
        spin.connect("value-changed", self.on_interpolation_step_changed)
        hbox.append(spin)
        self.right.append(frame)

        self.paned.set_resize_start_child(True)
        self.paned.set_shrink_start_child(False)

        self.paned.set_resize_end_child(False)
        self.paned.set_shrink_end_child(False)

        # Controls for running the integration (play/pause/stop)
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

        # Set default number of frames
        self.controls.nframes.set_value(100)

    ####################################################################
    # Entry points and lifecycle control
    ####################################################################

    @staticmethod
    def run(particle=None, blocking=True, connector=None):
        """Run Gyotoy as a standalone application or embedded window.

        Args:
            particle: Initial particle to display (Star or Photon)
            blocking: If True, run the GTK main loop (for standalone use)
            connector (multiprocessing.Connection or None): If provided,
                passed to the window for inter-process communication.

        Returns:
            GyotoyApplicationWindow: The created window instance
        """
        win = GyotoyApplicationWindow(particle=particle,
                                      connector=connector)
        win.blocking = blocking
        win.present()
        if blocking:
            win.main_loop = GLib.MainLoop()
            win.main_loop.run()

        return win

    def on_close_request(self, *args):
        """Handle window close request.

        Cleanly shuts down the worker process and quits the main loop.

        Args:
            *args: GTK callback arguments

        Returns:
            bool: False to allow the window to close
        """
        self.cmd_queue.put((QUIT,))  # Signal worker to exit
        if self.worker:
            self.worker.join(timeout=2.0)  # Graceful shutdown
            if self.worker.is_alive():
                self.worker.terminate()  # Force kill if stuck
        if self.main_loop is not None:
            GLib.idle_add(self.main_loop.quit)
        return False

    def on_quit(self, *args):
        """Handle quit action from menu.

        Args:
            *args: GTK callback arguments
        """
        self.close()

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

    ####################################################################
    # Compute and redraw
    ####################################################################

    def compute_and_redraw(self, *args):
        """Compute the trajectory and redraw the plot.

        This method starts a new simulation in the worker process and
        updates the UI to reflect the running state. The actual
        computation happens in the background worker process.

        Args:
            *args: Ignored (for callback compatibility)

        """
        if self.hold or self.particle is None or self.simulation_running:
            return

        # Disable UI and save focus
        self.simulation_running = True
        self.controls.set_progress(0.)
        self.controls.set_running(True)
        self.controls.set_status("Integrating...")
        self.last_focused_widget = self.get_focus()
        self.right.set_sensitive(False)

        coord = numpy.array(self.particle.InitCoord)

        # Send simulation command to worker
        self.cmd_queue.put((
            RUN_SIM,
            str(self.particle),
            self.particle.InitCoord[0],  # starttime
            self.endtime,
            self.controls.nframes.get_value_as_int(),
            self.interpolation_step
        ))

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
            if len(msg) >= 5:
                # Convert lists back to numpy arrays and filter out NaN values
                x, y, z = numpy.array(msg[2:5])
                mask = numpy.logical_and(
                    numpy.logical_not(numpy.isnan(x)),
                    numpy.logical_not(numpy.isnan(y)),
                    numpy.logical_not(numpy.isnan(z))
                    )
                self.viewer.axes.clear()
                self.viewer.axes.plot(x[mask], y[mask], z[mask])
                self.viewer.set_equal()
                self.viewer.canvas.draw_idle()
            self.controls.set_progress(msg[1])
        return True

    def process_control(self):
        """Process control messages from the worker.

        This method is called periodically (every 50ms) to check for
        control messages from the worker process (log messages,
        completion, errors).  It handles different message types and
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
                self.computation_epilogue(msg="Computation ended in error.",
                                          error=msg[1])
            elif msg[0] == 'fatal_error':
                self.computation_epilogue(msg="Fatal error. Restarting worker.",
                                          error=msg[1])
                self.restart_worker()
        except queue.Empty:
            pass
        return True

    def computation_epilogue(self, msg="Integration finished", error=None):
        """Handle the end of a computation.

        This method is called when a simulation completes
        (successfully or not).  It re-enables the UI, restores focus,
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
                self.last_focused_widget.grab_focus()  # Restore focus
            except:
                pass  # Widget may have been destroyed
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
            args=(self.cmd_queue, self.progress_queue, self.control_queue,
                  self.pause_event, self.stop_event),
            daemon=True
        )
        self.worker.start()
        self.last_heartbeat = time.time()
        self.simulation_running = False

    ####################################################################
    # Callbacks
    ####################################################################

    def on_open(self, *args):
        """Open a file dialog to load an XML particle configuration.

        Creates a file dialog with XML and All Files filters,
        defaulting to XML.  When a file is selected, it's loaded and
        the particle is updated.

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
            lambda dialog, result:
                self.on_open_file_selected(dialog, result)
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

        factory = None
        if file is not None:
            try:
                factory = Factory(file.get_path())
            except GyotoError as e:
                show_error_dialog(
                    message=f"Error loading XML file {file.get_path()}:",
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
                        message=f"Could not construct Astrobj from XML file:",
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
                        message=f"Could not construct Photon from XML file:",
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
                self.filenames[particle] = file.get_path()
                self.set_particle(particle)

    def on_save(self, *args):
        """Save the current particle in the last XML file used.

        If the particle was not read from file and not yet saved to
        file, opens a dialog.

        Args:
            *args: GTK callback arguments

        """
        if self.filenames[self.particle]:
            try:
                Factory(self.particle).write(self.filenames[self.particle])
            except GyotoError as e:
                show_error_dialog(
                    message=f"Error writing XML file {file.get_path()}:",
                    detail=e.get_message(),
                    window=self
                )
        else:
            self.on_save_as(*args)

    def on_save_as(self, *args):
        """Open a file dialog to save the current particle as XML.

        Creates a file dialog with XML and All Files filters,
        defaulting to XML.  When a file is selected, the current
        particle is saved to that file.

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
            lambda dialog, result:
                self.on_save_file_selected(dialog, result)
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
            Factory(self.particle).write(file.get_path())
            self.filenames[self.particle] = file.get_path()
        except GyotoError as e:
            show_error_dialog(
                message=f"Error writing XML file {file.get_path()}:",
                detail=e.get_message(),
                window=self
            )

    def on_star_selected(self, wdgt):
        """Handle Star/Photon radio button toggle.

        When Star is selected, sets the current particle to the Star
        instance.  When Photon is selected, sets the current particle
        to the Photon instance.  Synchronizes the metric between
        particles.

        Args:
            wdgt: The Gtk.CheckButton that was toggled

        """
        if wdgt.get_active() and isinstance(self.star, Star):
            if self.particle is not None:
                self.star.Metric = self.particle.Metric
            self.set_particle(self.star)
        elif isinstance(self.photon, Photon):
            if self.particle is not None:
                self.photon.Metric = self.particle.Metric
            self.set_particle(self.photon)

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

        help_text = """\
Gyotoy - Gyoto Geodesic Integration Visualizer

OVERVIEW:
Gyotoy is a GTK4 application for simulating and visualizing geodesics
(time-like or null) in spacetimes supported by the Gyoto library.
It provides an interactive 3D view of particle trajectories.

UI LAYOUT:
- Left Panel: 3D Matplotlib viewer displaying particle trajectory.
  Use mouse to rotate (left-click+drag), and pan (middle-click+drag).
  Use plot toolbar for other actions including saving the plot.
- Right Panel: Property editor for adjusting metric and particle
  parameters. Changes trigger automatic recomputation.
- Bottom: Simulation controls (play/pause/stop/reset), interpolation
  settings, and status display.

MENU BUTTONS:
- Open (Ctrl+O): Load an XML particle configuration file.
- Save (Ctrl+S): Save current particle to last used file.
- Save As (Ctrl+Shift+S): Save current particle to a new file.
- Help (F1): Show this help dialog.
- Quit (Ctrl+Q): Exit the application.

OHER KEYBOARD SHORTCUTS:
- Ctrl+R: Compute and redraw trajectory.
- Escape: close active dialog window (error, help...).

PROPERTY EDITOR:
- Star/Photon radio buttons: Switch between particle types.
- Edit particle parameters: position, velocity, metric properties...
- All parameter changes immediately trigger a computation and redraw
  unless the stop (■) at the bottom right of the window is activated.
- The InitCoord vector is by default displayed with 7 cells
  corresponding to 4-position and 3-velocity (derivatives of the space
  coordinates with respect to time coordinate). Use the '3-velocity'
  and '4-velocity' radio buttons to switch between this view and the
  8-coordinate view corresponding to 4-position and 4-velocity. By
  default, the 4-velocity is renormalized (according to mass of the
  particle) at each keystroke. Click on the stop (■) button
  immediately below the coordinates to temporarily inhibit this
  behavior. Click this button again to finalize input.

SIMULATION CONTROLS:
- Reset: Reinitialize integration and clear viewer.
- Play/Pause: Start or pause the simulation.
- Stop: Stop the integration and inhibit/enable recomputation.
- N. frames: Number of intermediate frames to display (default: 100).
- End time: Final time for integration.
- Interpolation step: Step size for interpolation (0 for no
  interpolation, adaptive step used instead).

WORKFLOW:
1. Select particle type (Star or Photon).
2. Adjust properties in the editor.
3. Click Play or press Ctrl+R to compute trajectory.
4. Use 3D viewer to inspect the result.
5. Save particle description with Save/Save As.
"""
        label = Gtk.Label(
            label=help_text,
            halign=Gtk.Align.START,
            wrap=True,
            xalign=0.0
        )

        scrolled = Gtk.ScrolledWindow()
        scrolled.set_child(label)
        scrolled.set_policy(Gtk.PolicyType.AUTOMATIC,
                            Gtk.PolicyType.AUTOMATIC)

        dialog.set_child(scrolled)

        dialog.present()

    def on_value_changed(self, widget, name, *args):
        """Handle property value changes from the editor.

        When InitCoord changes, ensures the 4-velocity is normalized.
        Triggers a new simulation.

        Args:
            widget: The widget that emitted the signal
            name: The name of the property that changed
            *args: Additional arguments
        """
        if self.hold:
            return
        if name == 'InitCoord':
            coord = numpy.array(self.particle.InitCoord)
            if self.particle == self.star:
                self.particle.Metric.normalizeFourVel(coord)
            else:
                self.particle.Metric.nullifyCoord(coord)
            self.particle.InitCoord = coord
            self.editor.on_3vel_toggled(name='InitCoord')
        self.compute_and_redraw()

    def on_child_changed(self, widget, name, *args):
        """Handle metric changes from the editor.

        Synchronizes the metric between Star and Photon particles.

        Args:
            widget: The widget that emitted the signal
            name: The name of the property that changed
            *args: Additional arguments
        """
        if name == 'Metric':
            if self.particle == self.star:
                self.photon.Metric = self.star.Metric
            else:
                self.star.Metric = self.photon.Metric
        self.on_child_mutated(widget, name, *args)

    def on_child_mutated(self, widget, pname, name, *args):
        """Handle metric mutations from the editor.

        Renormalizes InitCoord for both Star and Photon when the
        metric changes.

        Args:
            widget: The widget that emitted the signal
            name: The name of the changed property in the child
            name: The name of the child
            *args: Additional arguments

        """
        # React only to events in the metric
        if name != 'Metric':
            return

        # Don't touch anything if self.hold is set
        if self.hold:
            return

        # Don't touch anything if metric is not set
        if self.particle.Metric is None:
            return

        # Normalize star init coord
        coord = numpy.array(self.star.InitCoord)
        self.star.Metric.normalizeFourVel(coord)
        self.star.InitCoord = coord

        # Nullify photon init coord
        coord = numpy.array(self.photon.InitCoord)
        self.photon.Metric.nullifyCoord(coord)
        self.photon.InitCoord = coord

        # Let the editor update InitCoord display
        self.editor.on_3vel_toggled(name='InitCoord')

        # Compute and redraw
        self.compute_and_redraw()

    def on_endtime_changed(self, wdgt):
        """Handle end time spin button value changes.

        Updates the end time and reinitializes the particle if needed.
        Triggers a new simulation.

        Args:
            wdgt: The ScientificSpin widget that changed
        """
        prev = self.endtime
        self.endtime = wdgt.get_value()
        if self.particle is not None:
            starttime = self.particle.InitCoord[0]
            if prev > starttime:
                if self.endtime < prev:
                    self.particle.reInit()
            else:
                if self.endtime > prev:
                    self.particle.reInit()
        self.compute_and_redraw()

    def on_interpolation_step_changed(self, wdgt):
        """Handle interpolation step spin button value changes.

        Updates the interpolation step and triggers a new simulation.

        Args:
            wdgt: The ScientificSpin widget that changed
        """
        self.interpolation_step = wdgt.get_value()
        self.compute_and_redraw()

    def on_reset(self, wdgt):
        """Handle reset button click.

        Reinitializes the current particle and clears the viewer.

        Args:
            wdgt: The button that was clicked
        """
        self.particle.reInit()
        self.viewer.axes.clear()
        self.viewer.set_equal()
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
        self.hold = wdgt.stop_button.get_active()
        if self.hold:
            self.stop_event.set()
            self.controls.set_status("Holding integration (press play or stop).")
        else:
            self.stop_event.clear()
            self.controls.set_status("Ready for integration.")

    ####################################################################
    # Setters / getters
    ####################################################################

    def set_particle(self, particle):
        """Set the current particle and update the UI.

        Args:
            particle: The particle to set (Star or Photon)
        """
        if isinstance(particle, Star):
            self.particle_star.set_active(True)
            self.star = particle
        elif isinstance(particle, Photon):
            self.particle_photon.set_active(True)
            self.photon = particle
        else:
            show_error_dialog(message="Wrong type for particle:",
                              detail=repr(type(particle)),
                              window=self
                              )
            return
        self.particle = particle
        # Make sure we have exactly two items in self.connectors and
        # self.filenames
        if particle not in self.connectors:
            self.connectors[particle] = None
            self.connectors = {self.star: self.connectors[self.star],
                               self.photon: self.connectors[self.photon]}
            self.filenames =  {self.star: self.filenames[self.star],
                               self.photon: self.filenames[self.photon]}
        self.editor = PropertyEditorBox(particle,
                                        first=['InitCoord', 'Metric'],
                                        connector = self.connectors[particle])
        self.editor_scroller.set_child(self.editor)
        self.editor.connect('value-changed', self.on_value_changed)
        self.editor.connect('child-changed', self.on_child_changed)
        self.editor.connect('child-mutated', self.on_child_mutated)
        # Force 3-velocity display mode for InitCoord
        self.editor.widgets['InitCoord:veltype'].set_active(True)
        self.compute_and_redraw()

    def get_particle(self):
        """Get the current particle.

        Returns:
            The current particle (Star or Photon)
        """
        return self.particle

    ####################################################################
    # Default values
    ####################################################################

    def default_metric(self):
        """Create a default KerrBL metric with high spin.

        Returns:
            KerrBL: A KerrBL metric instance with Spin=0.995
        """
        metric = KerrBL()
        metric.Spin = 0.995
        return metric

    def default_star(self):
        """Create a default Star particle.

        Returns:
            Star: A Star particle with default coordinates and metric
        """
        particle = Star()
        particle.Metric = (self.default_metric() if self.photon is None
                           else self.photon.Metric)
        particle.initCoord((0., 10.791, 1.570796326794866, 0.,
                            1.1264111886458281, 0., 0., 0.018770516047594082))
        particle.Delta = 0.01
        return particle

    def default_photon(self):
        """Create a default Photon particle.

        Returns:
            Photon: A Photon particle with default coordinates and metric
        """
        particle = Photon()
        particle.Metric = (self.default_metric() if self.star is None
                           else self.star.Metric)
        if hasattr(particle.Metric, 'Spin'):
            r = 2. * (1 + numpy.cos(2./3. * numpy.acos(-particle.Metric.Spin)))
            spherical = (particle.Metric.coordKind() == GYOTO_COORDKIND_SPHERICAL)
            coord = (0., r, 0.5*numpy.pi if spherical else 0., 0.,
                     1., 0., 0. if spherical else 1., 1./r if spherical else 0.)
            coord = numpy.array(coord)
            particle.Metric.nullifyCoord(coord)
            particle.initCoord(coord)
        return particle

# Stand-alone entry point:
if __name__ == "__main__":
    raise SystemExit(GyotoyApplication.run_app(parsecliargs=True))
