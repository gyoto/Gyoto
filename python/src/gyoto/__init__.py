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
from gyoto import core, util, metric, astrobj, spectrum, spectrometer
core.Scenery.rayTrace = util.rayTrace
core.Scenery.rayTrace.__doc__ += core._core.Scenery_rayTrace.__doc__
core.Scenery.__getitem__ = util.Scenery_getitem
core.Scenery.__getitem__.__doc__ += core.Scenery.rayTrace.__doc__

core.Worldline.getCartesian =  util._Worldline_getCartesian
core.Worldline.getCartesian.__doc__ = core._core.Worldline_getCartesian.__doc__
core.Worldline.getCoord =  util._Worldline_getCoord
core.Worldline.getCoord.__doc__ = core._core.Worldline_getCoord.__doc__

# import plugins as modules

class GyotoPluginLoader(importlib.abc.Loader):
    """Loader for wrapping dynamic moduls around Gyoto plugins."""

    def create_module(self, spec):
        """Create a wrapper module around a Gyoto plugin."""
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
        """Initialize a wrapper module around a Gyoto plugin.

        Populate a wrapper module around a Gyoto plugin with useful
        functions and attributes.
        """
        fullname=module.__name__
        name = fullname.split('.', 1)[1]

        # Provide generic constructors for accessing classes
        # specifically in this plugin. This is useful when a plugin
        # implements two classes with the same name in two distinct
        # registed: for instance the python plug-in implements a
        # Metric and a Spectrum class that are both called "Python".

        def Metric(clsname):
            '''obj = gyoto.<plugin>.Metric("clsname")

            Instanciate Metric of kind clsname from plugin <plugin>.

            If <plugin> only registers one class with that name,
            equivalent to: obj=gyoto.<plugin>.identifier where
            identifier is a valid Python identifier based on clsname
            (see gyoto.util.valid_identifier()).
            '''
            return core.Metric(clsname, (name,))

        def Astrobj(clsname):
            '''obj = gyoto.<plugin>.Astrobj("clsname")

            Instanciate Astrobj of kind clsname from plugin <plugin>.

            If <plugin> only registers one class with that name,
            equivalent to: obj=gyoto.<plugin>.identifier where
            identifier is a valid Python identifier based on clsname
            (see gyoto.util.valid_identifier()).
            '''
            return core.Astrobj(clsname, (name,))

        def Spectrum(clsname):
            '''obj = gyoto.<plugin>.Spectrum("clsname")

            Instanciate Spectrum of given kind from plugin <plugin>.

            If <plugin> only registers one class with that name,
            equivalent to: obj=gyoto.<plugin>.identifier where
            identifier is a valid Python identifier based on clsname
            (see gyoto.util.valid_identifier()).
            '''
            return core.Spectrum(clsname, (name,))

        def Spectrometer(clsname):
            '''obj = gyoto.<plugin>.Spectrometer("clsname")

            Instanciate Spectrometer of given kind from plugin <plugin>.

            If <plugin> only registers one class with that name,
            equivalent to: obj=gyoto.<plugin>.identifier where
            identifier is a valid Python identifier based on clsname
            (see gyoto.util.valid_identifier()).
            '''
            return core.Spectrometer(clsname, (name,))

        module.Metric = Metric
        module.Astrobj = Astrobj
        module.Spectrum = Spectrum
        module.Spectrometer = Spectrometer

        module.__file__ = __file__
        module.__qualname__ = fullname
        module.__all__ = ['Metric', 'Astrobj', 'Spectrum', 'Spectrometer',
                          '__getattr__']
        module.__doc__='''Dynamically generated module around Gyoto plugin '''+name

        for entry, namespace in ((core.getMetricRegister(), metric),
                               (core.getAstrobjRegister(), astrobj),
                               (core.getSpectrumRegister(), spectrum),
                               (core.getSpectrometerRegister(), spectrometer)):
            while entry:
                if entry.plugin() == name:
                    classname = entry.name()
                    identifier = util.valid_identifier(classname)
                    constructor = util.make_constructor(namespace, classname,
                                                        name, identifier)
                    if identifier not in module.__dict__:
                        setattr(module, identifier, constructor)
                        module.__all__.append(identifier)
                    if identifier not in namespace.__dict__:
                        setattr(namespace, identifier, constructor)
                entry = entry.next()

class GyotoPluginFinder(importlib.abc.MetaPathFinder):
    """Finder for wrapping dynamic modules around Gyoto plugins.

    Check if a static submodule by that name already exists
    (gyoto/<name>.py or gyoto/<name>/__init__.py. In this case, let
    Python handle it normally.

    Else, instruct Python to try and build it dynamically using
    GyotoPluginBuilder.
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
            GyotoPluginLoader(),
            origin="dynamic",
        )

# Register the finder in sys.meta_path
sys.meta_path.insert(0, GyotoPluginFinder())

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
