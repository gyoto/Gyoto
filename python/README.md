# Python bindings for Gyoto

## Purpose

Files in this directory (python/ in the Gyoto source tree) allow
calling Gyoto routines from Python for building complex algorythms for
instance for computing movies of for fitting data. This is different
from, but complementary with, the files in the plugins/python/
directory, which allow extending Gyoto with new Metric, Astrobj and
Spectrum classes implemented in Python.

## Building & installing

This directory is built and installed with the rest of Gyoto. See
INSTALL.Gyoto.md in the top-most directory. Make sure the right Python
interpreter is found at the "configure" step. If installing to a
custom prefix (e.g. your home directory), set the PYTHONPATH
enverionment variable to where the Gyoto python module is installed at
the "make install" step.

## Known issues

It is not possible to import gyoto in Python from the python/ source
directory. Change to another directory before trying to run the
example scripts.

