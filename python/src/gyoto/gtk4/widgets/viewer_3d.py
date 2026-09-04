"""3D Matplotlib Viewer Widget for GTK4

This module provides a GTK4 widget that embeds a Matplotlib 3D axes
for visualizing 3D data. The widget is designed to be agnostic of
Gyoto-specific functionality, serving as a pure visualization
component.

Widget Layout
------------
    ┌──────────────────────────────────────────────────────────────┐
    │                                                              │
    │               matplotlib figure with 3D plot                 │
    │                                                              │
    ├──────────────────────────────────────────────────────────────┤
    │  NavigationToolbar2GTK4                                      │
    ├──────────────┬─────────────┬───────────┬────────────┬────────┤
    │ Azimuth      │ Elevation   │ Roll      │ Box Size   │Recenter│
    │ [        ]-+ │ [       ]-+ │ [     ]-+ │ [      ]-+ │  0  ⌖  │
    └──────────────┴─────────────┴───────────┴────────────┴────────┘

Description
-----------
- **Top Section**: FigureCanvasGTK4Agg displaying 3D content
- **Middle Section**: NavigationToolbar2GTK4 with standard Matplotlib
    controls
- **Bottom Section**: controls for camera orientation, box size, and
    recentering the view

The four ScientificSpin controls are kept synchronized with the
Matplotlib view. Changes made with the mouse, through the Matplotlib
API, or through the controls are reflected in the other
representation.

Example Usage
-------------
    viewer = Viewer3D()

    # Access the matplotlib axes directly
    ax = viewer.axes
    ax.plot([0, 1], [0, 1], [0, 1])  # Plot a line

    # Update the display
    viewer.draw()

    # Reset the view
    viewer.reset_view()

    # Force equal aspect ratio
    viewer.set_equal()

Signals
-------
``azim-changed``
    Emitted when the azimuth changes; the new value is passed as a
    float.
``elev-changed``
    Emitted when the elevation changes; the new value is passed as a
    float.
``roll-changed``
    Emitted when the roll changes; the new value is passed as a float.
``box-size-changed``
    Emitted when the x-axis view size changes; the new size is passed
    as a float. The size is the distance between the two values
    returned by ``get_xlim3d()``.
``center-changed``
    Emitted when the view center changes; the new center is passed as
    three floats (x, y, z).

"""

__all__ = ['Viewer3D']

import gi
gi.require_version("Gtk", "4.0")
from gi.repository import Gtk, GObject, GLib

from matplotlib.figure import Figure
from matplotlib.backends.backend_gtk4agg import FigureCanvasGTK4Agg
from matplotlib.backends.backend_gtk4 import NavigationToolbar2GTK4

from .scientific_spin import ScientificSpin


class Viewer3D(Gtk.Box):
    """GTK4 widget embedding a Matplotlib 3D view with navigation toolbar.

    This widget provides a self-contained 3D visualization component
    that can be embedded in any GTK4 application. In addition to the
    standard Matplotlib navigation toolbar, it provides controls for
    azimuth, elevation, roll, box size, and view centering.

    Attributes:
        figure (matplotlib.figure.Figure): The Matplotlib figure instance
        axes (mpl_toolkits.mplot3d.Axes3D): The 3D axes for plotting
        canvas (FigureCanvasGTK4Agg): The GTK4 canvas widget
        toolbar (NavigationToolbar2GTK4): The Matplotlib navigation toolbar
        azim_spin (ScientificSpin): Azimuth control
        elev_spin (ScientificSpin): Elevation control
        roll_spin (ScientificSpin): Roll control
        box_size_spin (ScientificSpin): View box size control

    Signals:
        azim-changed: emitted when the azimuth changes
        elev-changed: emitted when the elevation changes
        roll-changed: emitted when the roll changes
        box-size-changed: emitted when the box size changes
        center-changed: emitted when the view center changes

    Note:
        The widget uses a vertical layout with the canvas, navigation
        toolbar, and view controls from top to bottom.

    """

    __gsignals__ = {
        "azim-changed": (GObject.SignalFlags.RUN_FIRST, None, (float,)),
        "elev-changed": (GObject.SignalFlags.RUN_FIRST, None, (float,)),
        "roll-changed": (GObject.SignalFlags.RUN_FIRST, None, (float,)),
        "box-size-changed": (GObject.SignalFlags.RUN_FIRST, None, (float,)),
        "center-changed": (
            GObject.SignalFlags.RUN_FIRST, None, (float, float, float)
        ),
    }

    def __init__(self, *args, **kwargs):
        """Initialize the Viewer3D widget.

        Creates a vertical box containing a Matplotlib canvas, the standard
        NavigationToolbar2GTK4, and controls for camera orientation, box
        size, and recentering.

        Programmatic updates of the controls are guarded so that updating
        them from a Matplotlib view change cannot recursively change the
        view again.
        """
        super().__init__(
            orientation=Gtk.Orientation.VERTICAL,
            *args,
            **kwargs
        )

        self._updating_view = False
        self._view_state = None

        #
        # Matplotlib objects
        #
        self.figure = Figure(
            layout="constrained"
        )

        self.axes = self.figure.add_subplot(
            111,
            projection="3d",
            proj_type="persp"
        )

        # Set equal aspect ratio by default (cube-shaped view)
        self.axes.set_box_aspect((1, 1, 1))

        #
        # GTK canvas (embeds the Matplotlib figure)
        #
        self.canvas = FigureCanvasGTK4Agg(self.figure)
        self.canvas.set_hexpand(True)
        self.canvas.set_vexpand(True)
        self.append(self.canvas)

        #
        # Navigation toolbar (standard Matplotlib controls)
        #
        self.toolbar = NavigationToolbar2GTK4(self.canvas)
        self.append(self.toolbar)

        #
        # Additional view controls
        #
        controls = Gtk.Box(
            orientation=Gtk.Orientation.HORIZONTAL,
            spacing=4
        )
        controls.set_homogeneous(True)
        self.controls = controls
        self.append(controls)

        self.azim_spin = self._make_spin(
            "Azimuth",
            "Set the azimuth angle of the 3D view.",
            self._on_azim_changed,
            value=self.axes.azim
        )
        self.elev_spin = self._make_spin(
            "Elevation",
            "Set the elevation angle of the 3D view.",
            self._on_elev_changed,
            value=self.axes.elev
        )
        self.roll_spin = self._make_spin(
            "Roll",
            "Set the roll angle of the 3D view.",
            self._on_roll_changed,
            value=self.axes.roll
        )
        self.box_size_spin = self._make_sspin(
            "Box Size",
            "Set the size of the 3D view box.",
            self._on_box_size_changed,
            value=1.0
        )

        controls.append(self._make_frame(
            "Azimuth", self.azim_spin
        ))
        controls.append(self._make_frame(
            "Elevation", self.elev_spin
        ))
        controls.append(self._make_frame(
            "Roll", self.roll_spin
        ))
        controls.append(self._make_frame(
            "Box Size", self.box_size_spin
        ))
        controls.append(self._make_recenter_frame())

        # Keep controls synchronized with changes made through Matplotlib.
        self.canvas.mpl_connect("draw_event", self._on_draw)

        #
        # Set default camera view and initialize the controls.
        #
        self.reset_view()
        self._synchronize_controls(emit=False)

    ####################################################################
    # Public API
    ####################################################################

    def clear(self):
        """Remove all artists (plots, lines, etc.) from the 3D scene.

        This method clears the axes and resets the view to the default.
        """
        self.axes.clear()
        self.reset_view()

    def draw(self):
        """Request a redraw of the canvas.

        This schedules a redraw for the next GTK idle cycle, ensuring
        the UI remains responsive.
        """
        self.canvas.draw_idle()

    def reset_view(self, elev=30., azim=-60):
        """Restore the default camera position and orientation.

        Sets the view to a standard 3D perspective:
        - Elevation: 30 degrees
        - Azimuth: -60 degrees
        - Roll: 0 degrees
        - Box aspect: (1, 1, 1) for equal scaling on all axes
        """
        self._updating_view = True
        try:
            self.axes.view_init(
                elev=elev,
                azim=azim,
                roll=0.
            )
            self.axes.set_box_aspect((1, 1, 1))
        finally:
            self._updating_view = False

        self._synchronize_controls(emit=True)
        self.canvas.draw_idle()

    def set_equal(self):
        """Force equal scaling on all axes.

        This method adjusts the axis limits so that the 3D view has equal
        scaling in all dimensions, preventing distortion of spherical
        objects into ellipsoids when the data ranges differ.

        The view remains centered where it currently is and uses the
        maximum range across all axes as the common box size.
        """
        ax = self.axes

        xlim = ax.get_xlim3d()
        ylim = ax.get_ylim3d()
        zlim = ax.get_zlim3d()

        xmid = (xlim[0] + xlim[1]) / 2
        ymid = (ylim[0] + ylim[1]) / 2
        zmid = (zlim[0] + zlim[1]) / 2

        radius = max(
            xlim[1] - xlim[0],
            ylim[1] - ylim[0],
            zlim[1] - zlim[0],
        ) / 2

        self._updating_view = True
        try:
            ax.set_xlim3d(xmid - radius, xmid + radius)
            ax.set_ylim3d(ymid - radius, ymid + radius)
            ax.set_zlim3d(zmid - radius, zmid + radius)
            ax.set_box_aspect((1, 1, 1))
        finally:
            self._updating_view = False

        self._synchronize_controls(emit=True)
        self.canvas.draw_idle()

    ####################################################################
    # View-control construction
    ####################################################################

    def _make_spin(self, name, tooltip, callback, value):
        adjustment = Gtk.Adjustment(
            value=value,
            lower=-180,
            upper=180,
            step_increment=1,
            page_increment=10,
            page_size=0,
        )

        spin = Gtk.SpinButton(
            adjustment=adjustment,
            climb_rate=1,
            digits=0,
        )
        spin.set_wrap(True)
        spin.set_numeric(True)
        spin.set_hexpand(True)
        spin.set_tooltip_text(tooltip)
        spin.connect("value-changed", callback)

        return spin

    def _make_sspin(self, name, tooltip, callback, value):
        spin = ScientificSpin(value=value, with_unit=False)
        spin.set_hexpand(True)
        spin.set_tooltip_text(tooltip)
        spin.connect("value-changed", callback)
        return spin

    def _make_frame(self, label, child):
        frame = Gtk.Frame(label=label)
        frame.set_hexpand(True)
        frame.set_tooltip_text(child.get_tooltip_text())

        box = Gtk.Box(
            orientation=Gtk.Orientation.HORIZONTAL,
            spacing=4
        )
        box.append(child)
        frame.set_child(box)
        return frame

    def _make_recenter_frame(self):
        frame = Gtk.Frame(label="Recenter")
        frame.set_hexpand(True)
        frame.set_tooltip_text("Recenter the 3D view.")

        box = Gtk.Box(
            orientation=Gtk.Orientation.HORIZONTAL,
            spacing=2
        )
        box.set_halign(Gtk.Align.CENTER)

        origin = Gtk.Button(label="0")
        origin.add_css_class("flat")
        origin.set_tooltip_text("Center the view on the origin (0, 0, 0).")
        origin.connect("clicked", self._on_center_origin)

        data = Gtk.Button(icon_name="find-location")
        data.add_css_class("flat")
        data.set_tooltip_text("Center the view on the plotted data.")
        data.connect("clicked", self._on_center_data)

        box.append(origin)
        box.append(data)
        frame.set_child(box)

        self.origin_button = origin
        self.data_center_button = data
        return frame

    ####################################################################
    # Synchronization with Matplotlib
    ####################################################################

    def _get_view_state(self):
        xlim = self.axes.get_xlim3d()
        ylim = self.axes.get_ylim3d()
        zlim = self.axes.get_zlim3d()

        return (
            float(self.axes.azim),
            float(self.axes.elev),
            float(self.axes.roll),
            float(xlim[1] - xlim[0]),
            float((xlim[0] + xlim[1]) / 2),
            float((ylim[0] + ylim[1]) / 2),
            float((zlim[0] + zlim[1]) / 2),
        )

    def _synchronize_controls(self, emit=True):
        """Synchronize all controls with the current axes view.

        The guard prevents the controls' own ``value-changed`` handlers
        from interpreting these programmatic updates as user input.
        """
        state = self._get_view_state()
        previous = self._view_state
        self._view_state = state

        if self._updating_view:
            return

        self._updating_view = True
        try:
            for i, spin in enumerate((
                    self.azim_spin,
                    self.elev_spin,
                    self.roll_spin,
                    self.box_size_spin,
            )):
                val = spin.get_value()
                if val != state[i]:
                    spin.set_value(state[i])
        finally:
            self._updating_view = False

        if not emit or previous is None:
            return

        if state[0] != previous[0]:
            self.emit("azim-changed", state[0])
        if state[1] != previous[1]:
            self.emit("elev-changed", state[1])
        if state[2] != previous[2]:
            self.emit("roll-changed", state[2])
        if state[3] != previous[3]:
            self.emit("box-size-changed", state[3])

        old_center = previous[4:7]
        new_center = state[4:7]
        if new_center != old_center:
            self.emit("center-changed", *new_center)

    def _synchronize_controls_idle(self):
        self._synchronize_controls(emit=True)
        return GLib.SOURCE_REMOVE

    def _on_draw(self, event):
        if not self._updating_view:
            self._synchronize_controls(emit=True)

    ####################################################################
    # Control callbacks
    ####################################################################

    def _on_azim_changed(self, spin):
        if self._updating_view:
            return
        self._updating_view = True
        try:
            self.axes.view_init(
                elev=self.axes.elev,
                azim=spin.get_value(),
                roll=self.axes.roll
            )
        finally:
            self._updating_view = False
        GLib.idle_add(self._synchronize_controls_idle)
        self.canvas.draw_idle()

    def _on_elev_changed(self, spin):
        if self._updating_view:
            return
        self._updating_view = True
        try:
            self.axes.view_init(
                elev=spin.get_value(),
                azim=self.axes.azim,
                roll=self.axes.roll
            )
        finally:
            self._updating_view = False
        GLib.idle_add(self._synchronize_controls_idle)
        self.canvas.draw_idle()

    def _on_roll_changed(self, spin):
        if self._updating_view:
            return
        self._updating_view = True
        try:
            self.axes.view_init(
                elev=self.axes.elev,
                azim=self.axes.azim,
                roll=spin.get_value()
            )
        finally:
            self._updating_view = False
        GLib.idle_add(self._synchronize_controls_idle)
        self.canvas.draw_idle()

    def _on_box_size_changed(self, spin):
        if self._updating_view:
            return

        size = spin.get_value()
        if size <= 0:
            return

        xlim = self.axes.get_xlim3d()
        ylim = self.axes.get_ylim3d()
        zlim = self.axes.get_zlim3d()
        xcenter = (xlim[0] + xlim[1]) / 2
        ycenter = (ylim[0] + ylim[1]) / 2
        zcenter = (zlim[0] + zlim[1]) / 2
        half = size / 2

        self._updating_view = True
        try:
            self.axes.set_xlim3d(xcenter - half, xcenter + half)
            self.axes.set_ylim3d(ycenter - half, ycenter + half)
            self.axes.set_zlim3d(zcenter - half, zcenter + half)
        finally:
            self._updating_view = False

        GLib.idle_add(self._synchronize_controls_idle)
        self.canvas.draw_idle()

    ####################################################################
    # Recentering
    ####################################################################

    def _on_center_origin(self, button):
        self._set_center((0., 0., 0.))

    def _on_center_data(self, button):
        center = self._get_data_center()
        if center is not None:
            self._set_center(center)

    def _set_center(self, center):
        xcenter, ycenter, zcenter = center
        xlim = self.axes.get_xlim3d()
        ylim = self.axes.get_ylim3d()
        zlim = self.axes.get_zlim3d()

        xsize = xlim[1] - xlim[0]
        ysize = ylim[1] - ylim[0]
        zsize = zlim[1] - zlim[0]
        size = xsize
        half = size / 2

        self._updating_view = True
        try:
            self.axes.set_xlim3d(xcenter - half, xcenter + half)
            self.axes.set_ylim3d(ycenter - half, ycenter + half)
            self.axes.set_zlim3d(zcenter - half, zcenter + half)
        finally:
            self._updating_view = False

        self._synchronize_controls(emit=True)
        self.canvas.draw_idle()

    def _get_data_center(self):
        """Return the center of the bounding box of plotted 3D data.

        Supports the common Matplotlib 3D artist types used for lines,
        scatter plots, line collections, and polygon
        collections. Returns ``None`` when no usable 3D data are
        present.

        """
        points = []

        for line in self.axes.lines:
            try:
                x, y, z = line.get_data_3d()
                points.extend(zip(x, y, z))
            except AttributeError:
                continue

        for collection in self.axes.collections:
            offsets = getattr(collection, "_offsets3d", None)
            if offsets is not None:
                x, y, z = offsets
                points.extend(zip(x, y, z))
                continue

            segments = getattr(collection, "_segments3d", None)
            if segments is not None:
                for segment in segments:
                    points.extend(segment)
                continue

            # Poly3DCollection stores its projected 3D vertices in _vec.
            vec = getattr(collection, "_vec", None)
            if vec is not None and len(vec) >= 3:
                points.extend(zip(vec[0], vec[1], vec[2]))

        if not points:
            return None

        try:
            xmin = min(point[0] for point in points)
            xmax = max(point[0] for point in points)
            ymin = min(point[1] for point in points)
            ymax = max(point[1] for point in points)
            zmin = min(point[2] for point in points)
            zmax = max(point[2] for point in points)
        except (TypeError, ValueError):
            return None

        return (
            (xmin + xmax) / 2,
            (ymin + ymax) / 2,
            (zmin + zmax) / 2,
        )
