"""The General relativitY Orbit Tracer of paris Observatory

"""
import sys
import os
import types
import importlib
import importlib.abc
import importlib.util
from pathlib import Path
import warnings

# Expose gyoto.core as gyoto
from .core import *

# import gyoto.core
from gyoto import core

# avoid outputting n times the same warning about lorene while loading the rest
core.verbose(0)
from gyoto import util, animate, std, metric, astrobj, spectrum, spectrometer

# try importing the python submodule.
try:
    from gyoto import python
except:
    warnings.warn("gyoto.python is not available")

# try importing the lorene submodule
try:
    from gyoto import lorene
except:
    pass

# try importing the gtk4 submodule
try:
    from gyoto import gtk4
except:
    warnings.warn("gyoto.gtk4 is not available")

core.verbose(core.GYOTO_DEFAULT_VERBOSITY)

# Provide a Pythonic wrapper around Scenery.rayTrace.
# The underlying C++-like interface remains accessible.
core.Scenery.rayTrace = util.rayTrace
core.Scenery.rayTrace.__doc__ += core._core.Scenery_rayTrace.__doc__
core.Scenery.__getitem__ = util.Scenery_getitem
core.Scenery.__getitem__.__doc__ += core.Scenery.rayTrace.__doc__

core.Worldline.getCartesian =  util._Worldline_getCartesian
core.Worldline.getCartesian.__doc__ = core._core.Worldline_getCartesian.__doc__
core.Worldline.getCoord =  util._Worldline_getCoord
core.Worldline.getCoord.__doc__ = core._core.Worldline_getCoord.__doc__

# Generate submodules on the fly for Gyoto plugins

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

        class Metric(core.Metric):
            _plugin = name
            def __init__(self, clsname):
                super().__init__(clsname, (self._plugin,))
        Metric.__module__ = fullname
        Metric.__doc__ = '''Instanciate Metric of kind clsname from plugin '''+name+'''.

            If '''+name+''' only registers one class with that name,
            equivalent to: obj=gyoto.'''+name+'''.identifier() where
            identifier is a valid Python identifier based on clsname
            (see gyoto.util.valid_identifier()).
            '''

        class Astrobj(core.Astrobj):
            _plugin = name
            def __init__(self, clsname):
                super().__init__(clsname, (self._plugin,))
        Astrobj.__module__ = fullname
        Astrobj.__doc__ = '''Instanciate Astrobj of kind clsname from plugin '''+name+'''.

            If '''+name+''' only registers one class with that name,
            equivalent to: obj=gyoto.'''+name+'''.identifier() where
            identifier is a valid Python identifier based on clsname
            (see gyoto.util.valid_identifier()).
            '''

        class Spectrum(core.Spectrum):
            _plugin = name
            def __init__(self, clsname):
                super().__init__(clsname, (self._plugin,))
        Spectrum.__module__ = fullname
        Spectrum.__doc__ = '''Instanciate Spectrum of kind clsname from plugin '''+name+'''.

            If '''+name+''' only registers one class with that name,
            equivalent to: obj=gyoto.'''+name+'''.identifier() where
            identifier is a valid Python identifier based on clsname
            (see gyoto.util.valid_identifier()).
            '''

        class Spectrometer(core.Spectrometer):
            _plugin = name
            def __init__(self, clsname):
                super().__init__(clsname, (self._plugin,))
        Spectrometer.__module__ = fullname
        Spectrometer.__doc__ = '''Instanciate Spectrometer of kind clsname from plugin '''+name+'''.

            If '''+name+''' only registers one class with that name,
            equivalent to: obj=gyoto.'''+name+'''.identifier() where
            identifier is a valid Python identifier based on clsname
            (see gyoto.util.valid_identifier()).
            '''
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
                    klass = util.make_class(namespace, classname,
                                            name, identifier)
                    if identifier not in module.__dict__:
                        setattr(module, identifier, klass)
                        module.__all__.append(identifier)
                    if identifier not in namespace.__dict__:
                        setattr(namespace, identifier, klass)
                entry = entry.next()

class GyotoPluginFinder(importlib.abc.MetaPathFinder):
    """Finder for wrapping dynamic modules around Gyoto plugins."""

    def find_spec(self, fullname, path, target=None):
        """Find the module spec for the given module name."""
        if not fullname.startswith(__name__ + "."):
            return None

        submodule_name = fullname.split('.', 1)[1]
        module_dir = Path(__file__).parent

        # Split the submodule name into parts for nested modules
        parts = submodule_name.split('.')
        submodule_path = module_dir
        for part in parts:
            submodule_path = submodule_path / part

        # Check for .py file
        if submodule_path.with_suffix(".py").exists():
            return None

        # Check for package (directory with __init__.py)
        if (submodule_path.is_dir()
            and (submodule_path / "__init__.py").exists()):
            return None

        # Check for compiled module (.so)
        so_pattern = f"{parts[-1]}.*.so"
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
    """Autoload submodules.

    First tries to load an existing submodule. If it fails, it
    tries to load a Gyoto plugin by that name and expose it as a
    minimal submodule.
    """
    # Avoid infinite recursion for standard attributes
    if name.startswith('__') and name.endswith('__'):
        raise AttributeError(f"module '{__name__}' has no attribute '{name}'")

    fullname = f"{__name__}.{name}"

    # First, try to import the submodule normally
    try:
        mod = importlib.import_module(fullname)
        setattr(sys.modules[__name__], name, mod)
        return mod
    except ModuleNotFoundError:
        pass  # Fall back to dynamic plugin wrapper

    # If normal import fails, try to load as a Gyoto plugin
    try:
        core.requirePlugin(name)
    except core.Error:
        raise ModuleNotFoundError(f"No module named '{fullname}'")

    # Dynamically create the module
    module = types.ModuleType(fullname)
    sys.modules[fullname] = module
    setattr(sys.modules[__name__], name, module)

    # Initialize the module (same as GyotoPluginLoader.exec_module)
    GyotoPluginLoader().exec_module(module)
    return module
