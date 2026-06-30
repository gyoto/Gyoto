# Gyoto extension module for Python

## Purpose

Files in this directory (python/ in the Gyoto source tree) allow
building a Python extension module to call Gyoto routines from
Python. This allows implementing complex algorithms for instance for
computing movies or for fitting data.

This is different from, but complementary with, the files in the
plugins/python/ directory, which allow extending Gyoto with new
Metric, Astrobj and Spectrum classes implemented in Python.

## Word of caution

There are often several Python interpreters 

## Building & installing

This directory is built and installed with the rest of Gyoto. See
INSTALL.Gyoto.md in the top-most directory. Make sure the right Python
interpreter is found at the "configure" step. Beware that there is no
universal convention to install a Python module under an arbitraty
prefix. Gyoto provides two ways to specify where this module should be
installed: the arguments --with-python-install-scheme and
--with-pythondir options of the configure script. With
--with-python-install-scheme, you can select one of the schemes
provided by your system (see configure --help for a list of options
available at your site). With --with-pythondir, you can specify a
specific directory to install the module. If it is not already known
to Python, set the PYTHONPATH environment variable accordingly.

If you prefer, for instance if you want to install this module in a
virtualenv, you can also tell Gyoto to not install the module at all
and install it yourself with whatever tool you see fit (e.g. pip). At
the configure stage, use --with-python-install-scheme=none (or
--with-pythondir=nowhere). Then, running make without any argument
will still produce a wheel file under python/dist/ that you can
install with pip.

## Testing

Running 'make check' in the python/ directory will install the module
in a temporary location (the testbed/ subfolder) and run the tests in
the tests folder as well as the example*.py scripts.

One the module is installed (with 'make install' or you rown tools),
the same tests can be run from the python/ directory (assuming the
PYTHON variable is set to your python interpreter):

$PYTHON -m unittest discover -s tests -p "*.py"
$PYTHON example.py
$PYTHON example-patterndisk.py

## Usage

Once properly installed, import gyoto or its sub-modules: gyoto.core,
gyoto.std, gyoto.lorene, gyoto.utils, gyoto.metric, gyoto.astrobj,
gyoto.spectrum.

The example*.py scripts contain plenty of examples. Since the modules
expose most of the C++ API, the user and reference manuals apply, see
the doc/ directory in the Gyoto source tree and
https://gyoto.obspm.fr/
https://gyoto.obspm.fr/GyotoManual.pdf
https://gyoto.obspm.fr/namespaces.html
