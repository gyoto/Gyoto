"""SimulationControls: GTK4 Widget for Simulation Control and Monitoring

This module provides a compound GTK4 widget for controlling and monitoring
long-running simulations, including progress display, status messages, and
computation controls.

Widget Layout
------------
    ┌──────────────────────────────────────────────────────────────────────┐
    │ ████████████████████████████──────────────────────────────────────   │
    │                                                                      │
    │ Status...            N. frames: [100  ]      ⏮        ▶        ⏹     │
    │                                                                      │
    └──────────────────────────────────────────────────────────────────────┘

Description
-----------
- **Progress Bar**: Visual indicator of simulation progress
- **Status Label**: Text display for current status (truncated with
    ellipsis)
- **N: SpinButton**: Number of frames for the simulation
- **Buttons**:
  - ⏮ (Reset): Reset the simulation
  - ▶/⏸ (Play/Pause): Start or pause the simulation
  - ⏹ (Stop): Toggle stop mode (prevents redraw on parameter changes)

Signals:
    play-pause: Emitted when play/pause button is clicked
    reset: Emitted when reset button is clicked
    stop: Emitted when stop toggle button is toggled
    nframes-changed: Emitted when number of frames changes

"""

__all__ = ['SimulationControls']

import gi
gi.require_version("Gtk", "4.0")
from gi.repository import Gtk, GObject, Pango

from .scientific_spin import ScientificSpin

class SimulationControls(Gtk.Box):
    """A compound widget for controlling and monitoring simulations.

    This widget provides a complete control panel for simulations with:
    - Progress tracking
    - Status messages
    - Play/pause/reset/stop controls
    - Frame count configuration

    The widget is designed to be used with long-running computations
    where the UI needs to remain responsive while the simulation runs
    in a background process or thread.

    Attributes:
        progress (Gtk.ProgressBar): Progress bar widget
        status (Gtk.Label): Status label with ellipsization
        nframes (Gtk.SpinButton): Number of frames spin button
        reset_button (Gtk.Button): Reset button
        play_button (Gtk.Button): Play/Pause toggle button
        stop_button (Gtk.ToggleButton): Stop toggle button
        running (bool): Current running state

    """

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
        """Initialize the SimulationControls widget.

        Creates the layout with:
        - A progress bar at the top
        - A horizontal box containing:
          * Status label (left)
          * Number of frames spin button
          * Buttons (reset, play/pause, stop)
        """
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
        # Second row: Status + N. frames + Buttons
        #
        hbox = Gtk.Box(
            orientation=Gtk.Orientation.HORIZONTAL,
            spacing=10
        )
        self.append(hbox)

        #
        # Status label
        #
        self.status = Gtk.Label()
        self.status.set_ellipsize(Pango.EllipsizeMode.END)  # Truncate with "..."
        self.status.set_hexpand(True)  # Expand horizontally
        self.status.set_halign(Gtk.Align.START)  # Left-align text
        self.status.set_tooltip_text("")  # Initialize empty tooltip
        hbox.append(self.status)

        #
        # Number of frames control
        #
        frame_box = Gtk.Box(
            orientation=Gtk.Orientation.HORIZONTAL,
            spacing=5
        )

        frame_box.append(
            Gtk.Label(label="N. frames:")
        )

        adjustment = Gtk.Adjustment(
            value=100,
            lower=1,
            upper=10**9,  # Support up to 1 billion frames
            step_increment=1,
            page_increment=10
        )

        self.nframes = Gtk.SpinButton(
            adjustment=adjustment
        )
        self.nframes.set_numeric(True)
        self.nframes.connect(
            "value-changed",
            lambda *_: self.emit("nframes-changed")
        )
        frame_box.append(self.nframes)
        hbox.append(frame_box)

        #
        # Buttons (reset, play/pause, stop)
        #
        self.reset_button = Gtk.Button(
            icon_name="media-skip-backward-symbolic"
        )
        self.reset_button.set_tooltip_text(
            "Reset the simulation (pause and clear current state)"
        )
        self.reset_button.connect(
            "clicked",
            lambda *_: self.emit("reset")
        )

        self.play_button = Gtk.Button(
            icon_name="media-playback-start-symbolic"
        )
        self.play_button.set_tooltip_text(
            "Start or pause the simulation"
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
            "Toggle stop mode: when active, prevents redraw upon parameter changes"
        )
        self.stop_button.connect(
            "toggled",
            lambda *_: self.emit("stop")
        )

        hbox.append(self.reset_button)
        hbox.append(self.play_button)
        hbox.append(self.stop_button)

    ####################################################################
    # Public API
    ####################################################################

    def set_running(self, running):
        """Update the UI to reflect the running state.

        Changes the play/pause button icon based on the running state.

        Args:
            running (bool): True if simulation is running, False if
                paused/stopped

        """
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
        """Set the progress bar fraction.

        Args:
            fraction (float): Progress value between 0.0 and 1.0
        """
        self.progress.set_fraction(fraction)

    def set_status(self, text, error=None):
        """Set status text and tooltip.

        The status text is displayed in the label (truncated with
        ellipsis if needed).  The full text (plus error details if
        provided) is shown in the tooltip.

        Args:
            text (str): Status message to display
            error (str, optional): Error message to append to tooltip

        """
        self.status.set_text(text)
        if error:
            tooltip = f"{text}\n\nError: {error}"
        else:
            tooltip = text
        self.status.set_tooltip_text(tooltip)

    def is_stop_active(self):
        """Check if stop mode is active.

        Returns:
            bool: True if stop button is toggled on, False otherwise
        """
        return self.stop_button.get_active()

    ####################################################################
    # Internals
    ####################################################################

    def on_play_pause(self, button):
        """Handle play/pause button click.

        Toggles the running state and updates the button icon.

        Args:
            button (Gtk.Button): The play/pause button that was clicked
        """
        self.running = not self.running
        self.set_running(self.running)
        self.emit("play-pause")
