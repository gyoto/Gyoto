"""2D Matplotlib Viewer Widget for GTK4

This module provides a GTK4 widget that embeds a Matplotlib 2D axes
for visualizing 2D data. The widget is designed to be agnostic of
Gyoto-specific functionality, serving as a pure visualization
component.

Widget Layout
------------
    ┌──────────────────────────────────────────────────────────────┐
    │                                                              │
    │                                                              │
    │                                                              │
    │               matplotlib figure with 2D plot                 │
    │                                                              │
    │                                                              │
    │                                                              │
    │                                                              │
    ├──────────────────────────────────────────────────────────────┤
    │  ◀ Home  │  ←  →  │  ⊞  ⊟  │  🔍 Zoom  │  📏 Pan  │  💾 Save │
    └──────────────────────────────────────────────────────────────┘

Description
-----------
- **Top Section**: Matplotlib FigureCanvasGTK4Agg displaying 2D content
- **Bottom Section**: NavigationToolbar2GTK4 with standard Matplotlib
    controls
- **Features**: 2D visualization, camera control, zoom/pan tools

Example Usage
-------------
    viewer = Viewer2D()

    # Access the matplotlib axes directly
    ax = viewer.axes
    ax.imshow(data) # display image

    # Update the display
    viewer.draw()

    # Reset the view
    viewer.reset_view()

    # Force equal aspect ratio
    viewer.set_equal()

"""

__all__ = ['Viewer2D']

import gi
gi.require_version("Gtk", "4.0")
from gi.repository import Gtk

from matplotlib.figure import Figure
from matplotlib.backends.backend_gtk4agg import FigureCanvasGTK4Agg
from matplotlib.backends.backend_gtk4 import NavigationToolbar2GTK4

class Viewer2D(Gtk.Box):
    """GTK4 widget embedding a Matplotlib 2D view with navigation toolbar.

    This widget provides a self-contained 2D visualization component
    that can be embedded in any GTK4 application. It combines a
    Matplotlib FigureCanvasGTK4Agg for rendering and a
    NavigationToolbar2GTK4 for standard Matplotlib controls.

    The widget is designed to be completely independent of Gyoto,
    making it reusable for any application requiring 2D plotting.

    Attributes:
        figure (matplotlib.figure.Figure): The Matplotlib figure instance
        axes (TODO): The 2D axes for plotting
        canvas (FigureCanvasGTK4Agg): The GTK4 canvas widget
        toolbar (NavigationToolbar2GTK4): The Matplotlib navigation toolbar

    Note:
        The widget uses a vertical layout with the canvas on top and
        toolbar below.

    """

    def __init__(self, *args, **kwargs):
        """Initialize the Viewer2D widget.

        Creates a vertical box containing:
        - A Matplotlib FigureCanvasGTK4Agg for 2D rendering
        - A NavigationToolbar2GTK4 for camera controls

        Args:
            *args: Positional arguments passed to Gtk.Box
            **kwargs: Keyword arguments passed to Gtk.Box
        """
        super().__init__(
            orientation=Gtk.Orientation.VERTICAL,
            *args,
            **kwargs
        )

        #
        # Matplotlib objects
        #
        self.figure = Figure(
            layout="constrained"  # Prevents label overlap
        )

        self.axes = self.figure.add_subplot(
            111,  # Single subplot
        )

        # Set equal aspect ratio by default (cube-shaped view)
        self.axes.axis('equal')

        #
        # GTK canvas (embeds the Matplotlib figure)
        #
        self.canvas = FigureCanvasGTK4Agg(
            self.figure
        )
        self.canvas.set_hexpand(True)  # Expand horizontally
        self.canvas.set_vexpand(True)  # Expand vertically
        self.append(self.canvas)  # Add canvas to the box

        #
        # Navigation toolbar (standard Matplotlib controls)
        #
        self.toolbar = NavigationToolbar2GTK4(self.canvas)
        self.append(self.toolbar)  # Add toolbar below the canvas

    ####################################################################
    # Public API
    ####################################################################

    def clear(self):
        """Remove all artists (plots, lines, etc.) from the 2D scene.

        This method clears the axes and resets the view to the default.
        """
        self.axes.clear()

    def draw(self):
        """Request a redraw of the canvas.

        This schedules a redraw for the next GTK idle cycle, ensuring
        the UI remains responsive.
        """
        self.canvas.draw_idle()

    def set_equal(self):
        """Force equal scaling on all axes.

        This method adjusts the axis limits so that the 2D view has
        equal scaling in all dimensions, preventing distortion of
        spherical objects into ellipsoids when the data ranges differ.

        The approach centers the view on the midpoint of each axis and
        uses the maximum range across all axes to set uniform limits.

        """
        ax = self.axes

        # Get current axis limits
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()

        # Calculate midpoints for each axis
        xmid = (xlim[0] + xlim[1]) / 2
        ymid = (ylim[0] + ylim[1]) / 2

        # Find the maximum range
        radius = max(
            xlim[1] - xlim[0],
            ylim[1] - ylim[0],
        ) / 2

        # Set symmetric limits around each midpoint
        ax.set_xlim(xmid - radius, xmid + radius)
        ax.set_ylim(ymid - radius, ymid + radius)

        # Force equal aspect ratio
        ax.axis('equal')
