### Compond widgets suited for our needs
"""GTK4 custom widgets for Gyoto."""
from .filename_editor import FilenameEditor
from .gyoto_object_chooser import GyotoObjectChooser
from .property_editor_box import PropertyEditorBox
from .scientific_spin import ScientificSpin
from .simulation_controls import SimulationControls
from .vector_scientific_spin import VectorScientificSpin
from .viewer_3d import Viewer3D

__all__ = [
    'FilenameEditor', 'GyotoObjectChooser', 'PropertyEditorBox',
    'ScientificSpin', 'SimulationControls', 'VectorScientificSpin', 'Viewer3D'
]
