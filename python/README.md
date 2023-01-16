# Gyoto extension module for Python

## Purpose

Files in this directory (python/ in the Gyoto source tree) allow
building a Python extension module to call Gyoto routines from
Python. This allows implementing complex algorithms for instance for
computing movies of for fitting data.

This is different from, but complementary with, the files in the
plugins/python/ directory, which allow extending Gyoto with new
Metric, Astrobj and Spectrum classes implemented in Python.

## Building & installing

This directory is built and installed with the rest of Gyoto. See
INSTALL.Gyoto.md in the top-most directory. Make sure the right Python
interpreter is found at the "configure" step. If installing to a
custom prefix (e.g. your home directory), set the PYTHONPATH
environment variable to where the Gyoto python module is installed at
the "make install" step.

## Usage

Once properly installed, import gyoto or its sub-modules: gyoto.core,
gyoto.std, gyoto.lorene, gyoto.utils, gyoto.metric, gyoto.astrobj,
gyoto.spectrum (but not from the python/ source tree, the build system
is currently such that this does not work).

The examples.py script contains plenty of examples. Since the modules
expose most of the C++ API, the user and reference manuals apply, see
the doc/ directory in the Gyoto source tree and
https://gyoto.obspm.fr/
https://gyoto.obspm.fr/GyotoManual.pdf
https://gyoto.obspm.fr/namespaces.html

## Known issues

It is not possible to import gyoto in Python from the python/ source
directory. Change to another directory before trying to run the
example scripts.

