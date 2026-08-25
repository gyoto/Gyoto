"""VectorGyotoObjectChooser: GTK4 Widget for Complex Classes Editing

This module provides a custom GTK4 widget for editing instance of
Complex classes (Astrobj.Complex, Metric.Complex and
Spectrometer.Complex.

Widget Layout
------------
    ┌────────────────────────────────────────────┐
    │┌─────────────────────────────────────┐     │
    ││  GyotoObjectChooser                 │ - + │
    │└─────────────────────────────────────┘     │
    │┌─────────────────────────────────────┐     │
    ││  GyotoObjectChooser                 │ - + │
    │└─────────────────────────────────────┘     │
    │┌─────────────────────────────────────┐     │
    ││  GyotoObjectChooser                 │ - + │
    │└─────────────────────────────────────┘     │
    └────────────────────────────────────────────┘

Description
-----------
- **GyotoObjectChoosers**: one chooser per sub-object
- **Add/Remove Buttons**: Dynamically add or remove vector elements

Features:
- Signal emission on changes

"""

__all__ = ['VectorGyotoObjectChooser']

import gi
gi.require_version("Gtk", "4.0")
from gi.repository import Gtk

from .gyoto_object_chooser import GyotoObjectChooser
from ...core import Astrobj
from ...core import Metric
from ...core import Spectrometer
from ...core import Error as GyotoError
from ... import astrobj, metric, spectrometer

class VectorGyotoObjectChooser(Gtk.Box):
    """A GTK4 widget for editing vectors of Gyoto Objects.

    This widget provides a vertical list of GyotoObjectChosser widgets
    for editing individual vector components, with buttons to
    add/remove elements dynamically.

    Attributes:
        hold (bool): If True, don't emit signals.
        choosers (list): List of spin widgets for each vector element

    Signals:
        child-added: Emitted when a child is added.
        child-changed: Emitted when a child is replaced (includes the
            index of the child).
        child-removed: Emitted when a child is removed (includes the
            index of the child).
        recursive-value-changed: Emitted when a nested property changes
            (includes the full property path).

    """

    __gsignals__ = {
        "child-added": (gi.repository.GObject.SignalFlags.RUN_FIRST, None, ()),
        "child-changed": (gi.repository.GObject.SignalFlags.RUN_FIRST, None, (int,)),
        "child-mutated": (gi.repository.GObject.SignalFlags.RUN_FIRST, None, (int,str)),
        "child-removed": (gi.repository.GObject.SignalFlags.RUN_FIRST, None, (int,)),
        "recursive-value-changed": (gi.repository.GObject.SignalFlags.RUN_FIRST, None, (str,))
    }

    hold = False
    obj = None

    def __init__(self, value=None, dtype=None):
        """Initialize the VectorGyotoObjectChooser widget.

        Args:
            value (Gyoto Complex class): Initial vector values. If
                None, defaults to empty vector.

        """
        super().__init__(
            orientation=Gtk.Orientation.VERTICAL,
            spacing=6
        )


        if dtype is None:
            if value is None:
                raise ValueError("value and dtype can't both be None")
            if isinstance(value, Astrobj):
                dtype = astrobj.Complex
            elif isinstance(value, Metric):
                dtype = metric.Complex
            elif isinstance(value, Spectrometer):
                dtype = spectrometer.Complex
            else:
                raise ValueError(
                    "value must be an Astrobj, Metric or Spectrometer"
                )

        if dtype is astrobj.Complex:
            self.namespace = astrobj
        elif dtype is metric.Complex:
            self.namespace = metric
        elif dtype is spectrometer.Complex:
            self.namespace = spectrometer
        else:
            raise ValueError(
                'dtype should be one of CompexAstrobj, CompexMetric, ' +
                'ComplexSpectrometer'
            )

        self.dtype = dtype
        self.choosers = []
        self.items = []
        self.items_box = Gtk.Box(
            orientation=Gtk.Orientation.VERTICAL,
            spacing=2
        )

        self.append(self.items_box)
        
        # + button for adding vector elements
        buttons = Gtk.Box(
            orientation=Gtk.Orientation.HORIZONTAL,
            halign = Gtk.Align.END,
            hexpand=True,
            spacing=6
        )

        add = Gtk.Button(icon_name="list-add-symbolic")
        add.add_css_class("flat")
        add.set_tooltip_text(
            "Add item to vector"
        )
        add.connect("clicked", self.on_add)
        buttons.append(add)

        self.append(buttons)

        self._set_object(value)

    # private helpers
    def _apply_filter(self, text, force_visible=False):
        """Filter the nested PropertyEditorBox and return its match count.

        The chooser itself is kept visible by the PropertyEditorBox that owns
        it whenever either the chooser name or one of its child parameters
        matches. If the chooser name matches, its complete subtree remains
        visible.
        """
        count = 0

        for child in self.choosers:
            if child is None:
                continue

            child_count = child._apply_filter(text, force_visible=force_visible)

            self.items[self.choosers.index(child)].set_visible(
                not text or force_visible or child_count)

            count += child_count

        return count

    def _set_object(self, value):
        """Set the complext object.

        Dynamically adjusts the number of choosers to match the
        cardinal of the object, then updates each widget with the
        corresponding value.

        Args:
            value (self.dtype or None): Complex object to set

        Raises:
            ValueError: If value is not None and cannot be cast to
                self.dtype.

        """
        hold = self.hold
        self.hold = True

        if value is None:
            self.obj = None
            cardinal = 0
        else:
            try:
                self.obj = (
                    value if isinstance(value, self.dtype)
                    else self.dtype(value)
                )
            except GyotoError as e:
                raise ValueError(
                    f'object {value} cannot be cast to {self.dtype} ' +
                    f'({e.get_message()})'
                )
            cardinal = self.obj.getCardinal()
        
        #
        # Remove all widgets
        #
        while len(self.items) > 0:
            self.choosers.pop()
            item = self.items.pop()
            self.items_box.remove(item)

        #
        # Add widgets
        #
        for i in range(cardinal):
            self.on_add(obj=self.obj[i])

        self.hold=hold

    # callbacks for adding / removing child
    def on_add(self, button=None, obj=None):
        """Handle click on the add button.

        Adds a new chooser widget to the vector.

        Args:
            button (Gtk.Button): The add button that was clicked.
            obj (gyoto.core.Object): The object to represent.

        Note:
            When constructing from an existing complex object,
            self.hold must be true so that no signal is emitted and
            obj is not appended to self.obj.

        """

        item = Gtk.Box(
            orientation=Gtk.Orientation.HORIZONTAL,
            valign=Gtk.Align.START,
            spacing=2
        )

        chooser = GyotoObjectChooser(self.namespace, obj)
        chooser.connect("object-changed", self.on_child_object_changed)
        chooser.connect("object-mutated", self.on_child_object_mutated)
        chooser.connect("recursive-value-changed", self.on_child_recursive_value_changed)
        item.append(chooser)

        remove = Gtk.Button(icon_name="list-remove-symbolic")
        remove.set_vexpand(False)
        remove.add_css_class("flat")
        remove.connect("clicked", self.on_remove)
        item.append(remove)

        self.choosers.append(chooser)
        self.items.append(item)
        self.items_box.append(item)

        if not self.hold:
            self.obj.append(obj)
            self.emit("child-added")

    def on_remove(self, button):
        """Handle click on the remove button.

        Removes item from the vector.

        Args:
            button (Gtk.Button): The remove button that was clicked.
        """

        item = button.get_parent()
        index = self.items.index(item)

        self.obj.remove(index)

        self.choosers.remove(self.choosers[self.items.index(item)])
        self.items.remove(item)
        self.items_box.remove(item)

        if not self.hold:
            self.emit("child-removed", index)

    # callbacks for propagating signals
    def on_child_object_changed(self, widget, *args):
        """Handle object changes from the children.

        Emits 'object-mutated' signal when a property changes.

        Args:
            widget: The PropertyEditorBox that emitted the signal
            *args: Additional arguments
        """
        index = self.items.index(widget.get_parent())
        self.obj[index] = widget.obj
        if isinstance(widget.obj, spectrometer.Generic):
            # gyoto.spectrometer.Complex tries to hook itself to its
            # children to get notified when their properties change,
            # but the hook is not properly set when a child is
            # replaced. We fix it here.
            widget.obj.hook(self.obj)
            self.tell(widget.obj)
        elif isinstance(widget.obj, astrobj.Generic):
            # gyoto.astrobj.Complex.append sets the metric. We need to
            # emulate this behavior.
            if self.obj.Metric:
                widget.obj.Metric = self.obj.Metric
            else:
                self.obj.Metric = widget.obj.Metric
        self.emit("child-changed", index)

    def on_child_object_mutated(self, widget, name, *args):
        """Handle value changes from the children.

        Emits 'object-mutated' signal when a property changes.

        Args:
            widget: The PropertyEditorBox that emitted the signal
            *args: Additional arguments
        """
        index = self.items.index(widget.get_parent())
        self.emit("child-mutated", index, name)

    def on_child_recursive_value_changed(self, widget, ppath, *args):
        """Handle recursive value changes from the PropertyEditorBox.

        Forwards the signal to indicate that a nested property path
        has changed. Prepends ppath with [{index}] where `index` is
        the index of the mutated child.

        Args:
            widget: The PropertyEditorBox that emitted the signal.
            ppath: The full path of the property that changed.
            item: the member of self.items that contains this
                property.
            *args: Additional arguments.

        """
        index = self.items.index(widget.get_parent())
        self.emit("recursive-value-changed", f'[{index}].{ppath}')
