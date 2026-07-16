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
from ...core import Factory, Error as GyotoError, Photon
from ...core import GYOTO_COORDKIND_SPHERICAL
from ...std import Star, KerrBL


# --- Commands ---
RUN_SIM = 'run'
QUIT    = 'quit'

# --- Worker (runs forever) ---
def worker_func(cmd_queue, progress_queue, control_queue, pause_event, stop_event):
    """Persistent worker: waits for commands, never exits until QUIT.

    Runs in a separate process.

    Receives commands through cmd_queue and can be paused or stopped
    with pause_event and stop_event.

    Sends progress to progress_queue and control commands to control_queue.

    """
    try:
        import numpy
        import time
        from ...core import Factory
        from ...std import Star

        control_queue.put(('log', 'here'))
        next_update_fraction = 0.
        last_update = time.time()

        while True:
            try:
                cmd = cmd_queue.get(timeout=1.0)  # Blocks until command arrives
            except queue.Empty:
                continue # no command, continue

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
                        particle = f.photon() if f.kind() == 'Photon' else Star (f.astrobj())

                        # ensure interp_step as same sign as endtime.starttime
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
                            print(f'{n+1}/{len(frametimes)-1}')
                            if stop_event.is_set():
                                end_msg = ('aborted',)
                                break
                            while pause_event.is_set() and not stop_event.is_set():
                                time.sleep(0.1)

                            frametime = frametimes[n + 1]

                            # Actually integrate using Gyoto with
                            # built-in adaptive step
                            particle.xFill(frametime)

                            if interp_step:
                                # Extract positions at interpolated dates

                                # First whether we could reach
                                # frametimes[n+1]: impossible e.g. if
                                # geodesic joined the event horizon
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
                                print(tinteg)
                                print(ttmp)
                                print(frameend)
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
                        control_queue.put(('log', 'computation done'))
                        print(t)
                        control_queue.put(end_msg)
                except Exception as e:
                    control_queue.put(('error', traceback.format_exc()))
    except Exception as e:
        control_queue.put(('fatal_error', traceback.format_exc()))

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

    # Default values
    blocking  = True
    main_loop = None
    particle  = None
    star      = None
    photon    = None
    endtime   = 3000
    hold      = True
    worker    = None
    simulation_running = False
    last_focused_widget = None
    interpolation_step = 1.

    ####################################################################
    # Construction
    ####################################################################

    def __init__(self, application=None, particle=None, star=None, photon=None):
        super().__init__(application=application)

        self.star = star if star is not None else self.default_star()
        self.photon = photon if photon is not None else self.default_photon()

        self.set_title("Gyotoy")
        self.set_default_size(1024, 768)

        # Build UI
        self.build_headerbar()
        self.build_body()
        self.connect("close-request", self.on_close_request)

        # Prepare computation thread
        self.cmd_queue = Queue()
        self.progress_queue = Queue()
        self.control_queue = Queue()
        self.pause_event = Event()
        self.stop_event = Event()

        # Start worker ONCE at window creation
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
        if particle is None: particle = self.star
        self.hold=False
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

        ## first row, right: vertical box
        self.right = Gtk.Box(
            orientation=Gtk.Orientation.VERTICAL,
            spacing=6
        )
        self.paned.set_end_child(self.right)

        # top of controil column: Star/Photon radio buttons
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
        self.right.append(scroll)
        self.editor_scroller = scroll
        self.editor_scroller.set_min_content_width(365)

        # an additional spin button for end time
        frame = Gtk.Frame()
        frame.set_label("End time")
        frame.set_label_align(0.0)
        frame.set_tooltip_text("Run simulation from initial contition time to end time")

        hbox = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=5)
        frame.set_child(hbox)

        spin = ScientificSpin(value=0.)
        spin.set_value(self.endtime)
        spin.connect("value-changed", self.on_endtime_changed)
        hbox.append(spin)
        self.right.append(frame)

        # an additional spin button for interpolation step
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

        # not sure we should initialize this early
        # self.editor = PropertyEditorBox(...)
        # scroll.set_child(
        #     self.editor
        # )

        self.paned.set_resize_start_child(True)
        self.paned.set_shrink_start_child(False)

        self.paned.set_resize_end_child(False)
        self.paned.set_shrink_end_child(False)

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

        # We don't need to do anything
        # self.controls.connect(
        #     "nframes-changed",
        #     self.on_nframes_changed
        # )
        self.controls.nframes.set_value(100)

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
        self.cmd_queue.put((QUIT,))  # Signal worker to exit
        if self.worker:
            self.worker.join(timeout=2.0)  # Graceful shutdown
            if self.worker.is_alive():
                self.worker.terminate()  # Force kill if stuck
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
    # Compute and redraw
    ####################################################################

    def compute_and_redraw(self, *args):
        '''Compute the trajectory and redraw the plot

        Accepts and ignores any parameters to work as a callback.
        '''
        if self.hold or self.particle is None or self.simulation_running:
            return

        self.simulation_running = True
        self.controls.set_progress(0.)
        self.controls.set_running(True)
        self.controls.set_status("Integrating...")
        self.last_focused_widget = self.get_focus()
        self.right.set_sensitive(False)

        coord=numpy.array(self.particle.InitCoord)
        print(f'norm at start: {self.particle.Metric.norm(coord[0:4], coord[4:8])}, {coord}')
        # if self.particle == self.star:
        #     self.particle.Metric.normalizeFourVel(coord)
        # else:
        #     self.particle.Metric.nullifyCoord(coord)
        # print(f'norm at start: {self.particle.Metric.norm(coord[0:4], coord[4:8])}, {coord}')
        # self.particle.InitCoord = coord
        # self.editor.on_3vel_toggled(name='InitCoord')
        self.cmd_queue.put((
            RUN_SIM,
            str(self.particle),
            self.particle.InitCoord[0],  # starttime
            self.endtime,
            self.controls.nframes.get_value_as_int(),
            self.interpolation_step
        ))

    def process_progress(self):
        msg = None
        while True:
            try:
                msg = self.progress_queue.get_nowait()
            except queue.Empty:
                break
        if msg:
            print(msg[1])
            if len(msg) >= 5:
                x, y, z = numpy.array(msg[2:5])
                mask = numpy.logical_and(
                    numpy.logical_not(numpy.isnan(x)),
                    numpy.logical_not(numpy.isnan(y)),
                    numpy.logical_not(numpy.isnan(z))
                    )
                print(mask)
                self.viewer.axes.clear()
                self.viewer.axes.plot(x[mask], y[mask], z[mask])
                self.viewer.set_equal()
                self.viewer.canvas.draw_idle()
            self.controls.set_progress(msg[1])
        return True

    def process_control(self):
        try:
            msg = self.control_queue.get_nowait()
            if msg[0] == 'log':
                print(msg[1])
            elif msg[0] == 'done':
                self.computation_epilogue(msg = "Computation finished.")
            elif msg[0] == 'aborted':
                self.computation_epilogue(msg = "Computation aborted.")
            elif msg[0] == 'error':
                print(msg[1])
                self.computation_epilogue(msg = "Computation ended in error.",
                                          error = msg[1])
            elif msg[0] == 'fatal_error':
                self.computation_epilogue(msg = "Fatal error. Restarting worker.",
                                          error = msg[1])
                self.restart_worker()
        except queue.Empty:
            pass
        return True

    def computation_epilogue(self, msg="Integration toto finished", error=None):
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
        """Terminate old worker (if any) and start a new one."""
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

    def on_star_selected(self, wdgt):
        if wdgt.get_active() and isinstance(self.star, Star):
            if self.particle is not None:
                self.star.Metric = self.particle.Metric
            self.set_particle(self.star)
        elif isinstance(self.photon, Photon):
            if self.particle is not None:
                self.photon.Metric = self.particle.Metric
            self.set_particle(self.photon)

    def on_value_changed(self, widget, name, *args):
        '''Keep 4-velocity normalized
        '''
        if name == 'InitCoord':
            coord =  numpy.array(self.particle.InitCoord)
            if self.particle == self.star:
                self.particle.Metric.normalizeFourVel(coord)
            else:
                self.particle.Metric.nullifyCoord(coord)
            self.particle.InitCoord = coord
            self.editor.on_3vel_toggled(name='InitCoord')
        self.compute_and_redraw()

    def on_child_changed(self, widget, name, *args):
        '''Keep Metric synchronized between the two particles
        '''
        if name == 'Metric':
            if self.particle == self.star:
                self.photon.Metric = self.star.Metric
            else:
                self.star.Metric = self.photon.Metric
        self.on_child_mutated(widget, name, *args)

    def on_child_mutated(self, widget, name, *args):
        '''Renormalize InitCoord when Metric changes
        '''

        print(f'on_child_mutated: name={name}')
        if name != 'Metric': return

        # normalize star init coord
        coord = numpy.array(self.star.InitCoord)
        self.star.Metric.normalizeFourVel(coord)
        self.star.InitCoord = coord

        # nullify photon init coord
        coord = numpy.array(self.photon.InitCoord)
        self.photon.Metric.nullifyCoord(coord)
        self.photon.InitCoord = coord

        self.editor.on_3vel_toggled(name='InitCoord')
        self.compute_and_redraw()

    def on_endtime_changed(self, wdgt):
        '''Called when changing the "End time" ScientificSpin

        Effects:
        - set self.endtime
        - compute and redraw
        '''
        prev=self.endtime
        self.endtime=wdgt.get_value()
        if self.particle is not None:
            starttime=self.particle.InitCoord[0]
            if prev > starttime:
                if self.endtime < prev:
                    self.particle.reInit()
            else:
                if self.endtime > prev:
                    self.particle.reInit()
        self.compute_and_redraw()

    def on_interpolation_step_changed(self, wdgt):
        self.interpolation_step = wdgt.get_value()
        self.compute_and_redraw()

    def on_reset(self, wdgt):
        self.particle.reInit()
        self.viewer.axes.clear()
        self.viewer.set_equal()
        self.viewer.canvas.draw_idle()

    def on_play_pause(self, wdgt):
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
        self.hold = wdgt.stop_button.get_active()
        if self.hold :
            self.stop_event.set()
            self.controls.set_status("Holding integration (press play or stop).")
        else:
            self.stop_event.clear()
            self.controls.set_status("Ready for integration.")


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
        self.editor = PropertyEditorBox(particle, first=['InitCoord', 'Metric'])
        self.editor_scroller.set_child(self.editor)
        self.editor.connect('value-changed', self.on_value_changed)
        self.editor.connect('child-changed', self.on_child_changed)
        self.editor.connect('child-mutated', self.on_child_mutated)
        self.editor.widgets['InitCoord:veltype'].set_active(True)
        self.compute_and_redraw()

    def get_particle(self):
        return self.particle

    ####################################################################
    # Default values
    ####################################################################

    def default_metric(self):
        metric = KerrBL()
        metric.Spin = 0.995
        return metric

    def default_star(self):
        particle = Star()
        particle.Metric = (self.default_metric() if self.photon is None
                           else self.photon.Metric)
        particle.initCoord((0.,10.791,1.570796326794866,0.,
                            1.1264111886458281, 0.,0.,0.018770516047594082))
        particle.Delta=0.01
        return particle

    def default_photon(self):
        particle = Photon()
        particle.Metric = (self.default_metric() if self.star is None
                           else self.star.Metric)
        r = 2. * (1 + numpy.cos(2./3. * numpy.acos(-particle.Metric.Spin)))
        spherical = (particle.Metric.coordKind() == GYOTO_COORDKIND_SPHERICAL)
        coord = (0., r, 0.5*numpy.pi if spherical else 0., 0.,
                 1., 0., 0. if spherical else 1., 1./r if spherical else 0.)
        coord = numpy.array(coord)
        particle.Metric.nullifyCoord(coord)
        particle.initCoord(coord)
        return particle

# Widget construction (__init__)
# Header bar, menu, Matplotlib embedding, and layout
# Callbacks and helper methods
# run() / run_app()

### Stand-alone entry point:

if __name__ == "__main__":
    raise SystemExit(MainWindow.run_app())
