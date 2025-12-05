"""The General relativitY Orbit Tracer of paris Observatory

Note that importing "gyoto" is deprecated and may cease to work in a
future release. Please update your code to import gyoto.core instead.

"""
import sys
import os
import types
import importlib
import importlib.abc
import importlib.util
from pathlib import Path

# For backwards compatibility, expose gyoto.core as gyoto
from gyoto.core import *

# Provide a Pythonic wrapper around Scenery.rayTrace.
# The underlying C++-like interface remains accessible.
from gyoto import core, util
core.Scenery.rayTrace = util.rayTrace
core.Scenery.rayTrace.__doc__ += core._core.Scenery_rayTrace.__doc__
core.Scenery.__getitem__ = util.Scenery_getitem
core.Scenery.__getitem__.__doc__ += core.Scenery.rayTrace.__doc__

core.Worldline.getCartesian =  util._Worldline_getCartesian
core.Worldline.getCartesian.__doc__ = core._core.Worldline_getCartesian.__doc__
core.Worldline.getCoord =  util._Worldline_getCoord
core.Worldline.getCoord.__doc__ = core._core.Worldline_getCoord.__doc__

# import plugins as modules

class DynamicModuleLoader(importlib.abc.Loader):
    """Loader for dynamically created modules."""

    def create_module(self, spec):
        """Create a new module object."""
        fullname=spec.name
        name = fullname.split('.', 1)[1]
        try:
            # Try to load a plugin by that name
            core.requirePlugin(name)
        except core.Error:
            # If still not successful, raise the correct error
            raise ModuleNotFoundError(f"No module named '{fullname}'")

        return types.ModuleType(fullname)

    def exec_module(self, module):
        """Initialize the module dynamically.
        This is called after the module is created.
        """
        fullname=module.__name__
        name = fullname.split('.', 1)[1]

        # Implement __getattr__ for our dynamic module
        def __getattr__(clsname):
            try:
                # check whether this is a Metric
                obj = core.Metric(clsname)
                return lambda : core.Metric(clsname)
            except core.Error:
                pass

            try:
                # check whether this is an Astrobj
                obj = core.Astrobj(clsname)
                return lambda : core.Astrobj(clsname)
            except core.Error:
                pass

            try:
                # check whether this is an Spectrum
                obj = core.Spectrum(clsname)
                return lambda : core.Spectrum(clsname)
            except core.Error:
                pass

            try:
                # check whether this is an Spectrometer
                obj = core.Spectrometer(clsname)
                return lambda : core.Spectrometer(clsname)
            except core.Error:
                raise AttributeError(f"module '{fullname}' has no attribute '{clsname}'")

        module.__getattr__=__getattr__

class DynamicModuleFinder(importlib.abc.MetaPathFinder):
    """Finder for dynamic modules.
    Checks if a .py file exists for the submodule.
    If not, creates the module dynamically.
    """

    def find_spec(self, fullname, path, target=None):
        """Find the module spec for the given module name.

        Args:
            fullname: The full module name (e.g., "module.submodule").
            path: The search path (unused here).
            target: The target module (unused here).

        Returns:
            A ModuleSpec if the module should be created dynamically,
            None if the module should be loaded normally.
        """
        # Only handle submodules of "module"
        if not fullname.startswith(__name__+"."):
            return None

        # Extract the submodule name (e.g., "submodule")
        submodule_name = fullname.split('.', 1)[1]

        # Expected path for the .py file
        module_dir = Path(__file__).parent
        submodule_path = module_dir / submodule_name

        # Check for .py file
        if submodule_path.with_suffix(".py").exists():
            return None

        # Check for package (directory with __init__.py)
        if (submodule_path.is_dir()
            and (submodule_path / "__init__.py").exists()):
            return None

        # Check for compiled module (.so)
        # Example pattern: _std.cpython-313-x86_64-linux-gnu.so
        so_pattern = f"{submodule_name}.*.so"
        so_matches = list(module_dir.glob(so_pattern))
        if so_matches:
            return None

        # Otherwise, create a spec for a dynamic module
        return importlib.util.spec_from_loader(
            fullname,
            DynamicModuleLoader(),
            origin="dynamic",
        )

# Register the finder in sys.meta_path
sys.meta_path.insert(0, DynamicModuleFinder())

def __getattr__(name):
    f'''Autoload submodules

    {__name__}.__getattr__ is called when Python encounters a call to
    an attribute that does not yet exist.

    It first tries to load an existing submodule. If it fails, it
    tries to load a Gyoto plugin by that name and expose it as a
    minimal submodule.

    '''
    # __getattr__ shouldn't be called in that case, but still the
    # right answer:
    if name in sys.modules[__name__].__dict__:
        return sys.modules[__name__].__dict__[name]

    # Take care of standard attributes
    if name.startswith('__') and name.endsswith('__'):
        raise AttributeError(f"module '{__name__}' has no attribute '{name}'")

    # Now try to load a submodule
    fullname=f"{__name__}.{name}"
    mod = importlib.import_module(fullname)
    setattr(sys.modules[__name__], name, mod)
    return mod
