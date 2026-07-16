"""
3D matplotlib viewer widget.

This widget deliberately knows nothing about Gyoto.
It is only a GTK wrapper around a matplotlib Axes3D.

Example:

    viewer = Viewer3D()

    ax = viewer.axes
    ax.plot([0, 1], [0, 1], [0, 1])

    viewer.draw()
"""

__all__ = ['Viewer3d']

import gi
gi.require_version("Gtk", "4.0")
from gi.repository import Gtk

from matplotlib.figure import Figure
from matplotlib.backends.backend_gtk4agg import FigureCanvasGTK4Agg
from matplotlib.backends.backend_gtk4 import NavigationToolbar2GTK4

class Viewer3D(Gtk.Box):
    """GTK widget embedding a matplotlib 3D view."""

    def __init__(self, *args, **kwargs):

        super().__init__(
            orientation=Gtk.Orientation.VERTICAL,
            *args,
            **kwargs
        )

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

        self.axes.set_box_aspect(
            (1, 1, 1)
        )

        #
        # GTK canvas
        #

        self.canvas = FigureCanvasGTK4Agg(
            self.figure
        )

        self.canvas.set_hexpand(True)
        self.canvas.set_vexpand(True)

        self.append(self.canvas)

        #
        # Navigation toolbar
        #

        self.toolbar = NavigationToolbar2GTK4(self.canvas)
        self.append(self.toolbar)

        #
        # Default view
        #

        self.reset_view()


    ####################################################################
    # Public API
    ####################################################################

    def clear(self):
        """Remove all artists from the scene."""

        self.axes.clear()
        self.reset_view()


    def draw(self):
        """Request a redraw."""

        self.canvas.draw_idle()


    def reset_view(self):
        """Restore a reasonable default camera."""

        self.axes.view_init(
            elev=30,
            azim=-60
        )

        self.axes.set_box_aspect(
            (1, 1, 1)
        )


    # def set_equal(self):
    #     """Force equal scaling on all axes.

    #     This avoids spheres becoming ellipsoids when the ranges differ.
    #     """

    #     xlim = self.axes.get_xlim3d()
    #     ylim = self.axes.get_ylim3d()
    #     zlim = self.axes.get_zlim3d()

    #     xr = abs(xlim[1] - xlim[0])
    #     yr = abs(ylim[1] - ylim[0])
    #     zr = abs(zlim[1] - zlim[0])

    #     self.axes.set_box_aspect(
    #         (xr, yr, zr)
    #     )

    def set_equal(self):
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

        ax.set_xlim3d(xmid - radius, xmid + radius)
        ax.set_ylim3d(ymid - radius, ymid + radius)
        ax.set_zlim3d(zmid - radius, zmid + radius)

        ax.set_box_aspect((1, 1, 1))
