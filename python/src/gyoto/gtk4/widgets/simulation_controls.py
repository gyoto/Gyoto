# A compound widget to control/monitor a long simulation run
# ┌──────────────────────────────────────────────────────────────────────┐
# │ █████████──────  ☑ Interpolate   Step: [ 1.000e-03 ▲▼ ]   N: [100]  │
# │                                                                  │
# │                              ⏮        ▶        ⏹                  │
# └──────────────────────────────────────────────────────────────────────┘

__all__ = ['SimulationControls']

import gi
gi.require_version("Gtk", "4.0")
from gi.repository import Gtk, GObject

from .scientific_spin import ScientificSpin


class SimulationControls(Gtk.Box):
    """Controls for numerical integration / animation."""

    __gsignals__ = {
        "play-pause": (
            GObject.SignalFlags.RUN_FIRST,
            None,
            ()
        ),
        "reset": (
            GObject.SignalFlags.RUN_FIRST,
            None,
            ()
        ),
        "stop": (
            GObject.SignalFlags.RUN_FIRST,
            None,
            ()
        ),
        "step-changed": (
            GObject.SignalFlags.RUN_FIRST,
            None,
            ()
        ),
        "nframes-changed": (
            GObject.SignalFlags.RUN_FIRST,
            None,
            ()
        ),
        "interpolate-changed": (
            GObject.SignalFlags.RUN_FIRST,
            None,
            ()
        ),
    }


    def __init__(self):

        super().__init__(
            orientation=Gtk.Orientation.HORIZONTAL,
            spacing=10
        )

        self.running = False


        #
        # Progress bar
        #

        self.progress = Gtk.ProgressBar()

        self.progress.set_hexpand(True)

        self.append(self.progress)


        #
        # Interpolation
        #

        self.interpolate = Gtk.CheckButton(
            label="Interpolate"
        )

        self.interpolate.connect(
            "toggled",
            lambda *_: self.emit(
                "interpolate-changed"
            )
        )

        self.append(self.interpolate)


        #
        # Step
        #

        step_box = Gtk.Box(
            orientation=Gtk.Orientation.HORIZONTAL,
            spacing=5
        )

        step_box.append(
            Gtk.Label(label="Step:")
        )

        self.step = ScientificSpin(
            value=1e-3
        )

        self.step.connect(
            "value-changed",
            lambda *_: self.emit(
                "step-changed"
            )
        )

        step_box.append(self.step)

        self.append(step_box)


        #
        # Number of frames
        #

        frame_box = Gtk.Box(
            orientation=Gtk.Orientation.HORIZONTAL,
            spacing=5
        )

        frame_box.append(
            Gtk.Label(label="N:")
        )

        adjustment = Gtk.Adjustment(
            value=100,
            lower=1,
            upper=10**9,
            step_increment=1,
            page_increment=10
        )

        self.nframes = Gtk.SpinButton(
            adjustment=adjustment
        )

        self.nframes.set_numeric(True)

        self.nframes.connect(
            "value-changed",
            lambda *_: self.emit(
                "nframes-changed"
            )
        )

        frame_box.append(self.nframes)

        self.append(frame_box)


        #
        # Transport buttons
        #

        self.reset_button = Gtk.Button(
            icon_name="media-skip-backward-symbolic"
        )

        self.reset_button.set_tooltip_text(
            "Pause and reset integration"
        )

        self.reset_button.connect(
            "clicked",
            lambda *_: self.emit("reset")
        )


        self.play_button = Gtk.Button(
            icon_name="media-playback-start-symbolic"
        )

        self.play_button.set_tooltip_text(
            "Integrate / pause integration"
        )

        self.play_button.connect(
            "clicked",
            self.on_play_pause
        )


        self.stop_button = Gtk.ToggleButton()

        self.stop_button.set_icon_name(
            "media-playback-stop-symbolic"
        )

        self.stop_button.set_tooltip_text(
            "Push-in to inhibit redraw upon parameter changes"
        )

        self.stop_button.connect(
            "toggled",
            lambda *_: self.emit("stop")
        )


        self.append(self.reset_button)
        self.append(self.play_button)
        self.append(self.stop_button)


    ####################################################################
    # Public helpers
    ####################################################################

    def set_running(self, running):

        self.running = running

        if running:
            self.play_button.set_icon_name(
                "media-playback-pause-symbolic"
            )
        else:
            self.play_button.set_icon_name(
                "media-playback-start-symbolic"
            )


    def set_progress(self, fraction):

        self.progress.set_fraction(
            fraction
        )


    def is_stop_active(self):

        return self.stop_button.get_active()


    ####################################################################
    # Internals
    ####################################################################

    def on_play_pause(self, button):

        self.running = not self.running

        self.set_running(
            self.running
        )

        self.emit(
            "play-pause"
        )
