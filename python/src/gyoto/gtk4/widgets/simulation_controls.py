# A compound widget to control/monitor a long simulation run
# ┌──────────────────────────────────────────────────────────────────────┐
# │ █████████──────  ☑ Interpolate   Step: [ 1.000e-03 ▲▼ ]   N: [100]  │
# │                                                                  │
# │                              ⏮        ▶        ⏹                  │
# └──────────────────────────────────────────────────────────────────────┘

__all__ = ['SimulationControls']

import gi
gi.require_version("Gtk", "4.0")
from gi.repository import Gtk, GObject, Pango

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
        "nframes-changed": (
            GObject.SignalFlags.RUN_FIRST,
            None,
            ()
        ),
    }


    def __init__(self):

        super().__init__(
            orientation=Gtk.Orientation.VERTICAL,
            spacing=0
        )

        self.running = False


        #
        # Progress bar
        #

        self.progress = Gtk.ProgressBar()

        self.progress.set_hexpand(True)

        self.append(self.progress)

        #
        # Second row
        #

        hbox = Gtk.Box(
            orientation=Gtk.Orientation.HORIZONTAL,
            spacing=10
        )

        self.append(hbox)


        #
        # Status
        #

        self.status = Gtk.Label()
        self.status.set_ellipsize(Pango.EllipsizeMode.END)  # Truncate with "..."
        self.status.set_hexpand(True)  # Fixed width (300px), natural height
        self.status.set_halign(Gtk.Align.START)  # Left-align text
        self.status.set_tooltip_text("")  # Initialize empty tooltip
        hbox.append(self.status)

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

        hbox.append(frame_box)


        #
        # Transport buttons
        #

        self.reset_button = Gtk.Button(
            icon_name="media-skip-backward-symbolic"
        )

        self.reset_button.set_tooltip_text(
            "Pause and reset computation"
        )

        self.reset_button.connect(
            "clicked",
            lambda *_: self.emit("reset")
        )


        self.play_button = Gtk.Button(
            icon_name="media-playback-start-symbolic"
        )

        self.play_button.set_tooltip_text(
            "Compute / pause computation"
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


        hbox.append(self.reset_button)
        hbox.append(self.play_button)
        hbox.append(self.stop_button)


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

    def set_status(self, text, error=None):
        """Set status text (truncated) and tooltip (full text + error)."""
        self.status.set_text(text)
        if error:
            tooltip = f"{text}\n\nError: {error}"
        else:
            tooltip = text
        self.status.set_tooltip_text(tooltip)

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
