# Python plug-in for Gyoto

## Purpose

Files in this directory (plugins/python/ in the Gyoto source) allow
building a Gyoto plug-in that embeds the Python interpreter to call
Python routine from within Gyoto. This allows extending Gyoto with new
Metric, Astrobj or Spectrum classes implemented in Python.

This is different from, but complementary with, the files in the
python/ directory, which allow extending Python to call Gyoto
routines from Python.

## Building and installing

This directory is built and installed with the rest of Gyoto. See
INSTALL.Gyoto.md in the top-most directory. Make sure the right Python
interpreter is found at the "configure" step.

## Caveats

Coding your new Metric, Astrobj or Spectrum class in Python is
efficient in terms of human resources, while coding it in C++ is much
more efficient in terms of computing resources.

## Usage

It is possible to compile this plug-in for the several installations
of Python that may coexist on the system. The name of the plug-in will
be the same as the name of the python executable (e.g. python3,
python3.9...) It the rest of this text, we will assume this is
`python3'.

This plug-in is not part of the plug-ins loaded by default in
Gyoto. Several possibilities exist to load it explicitly:

* set the GYOTO_PLUGINS environment variable, not forgetting to also
  list the other desired default plugins, e.g.:
  ```export GYOTO_PLUGINS=stdplug,python3,nofail:lorene```

* use the `--plugins' option of the gyoto command-line tool, e.g.:
  ```gyoto --plugins=stdplug,python3,nofail:lorene ...```

* in Python, use ```gyoto.core.requirePlugin('python3')```

* in XML, set the plugin property:
  ```<Metric kind="Python" plugin="python3">```

The content of plugins/python/doc/examples should give enough
resources to get started. The gyoto_sample_*.py modules contain
documentation about how to code a class in Python while the test*.py
scripts give an idea of how to instantiate and use such classes using
the Python bindings. One can also directly load the provided XML files
into the stand-alone gyoto executable by using the --plugins
option. Note that multi-threading is not efficient when using a Python
metric in particular, but MPI parallelization works quite well, so the
command line could look something like:

```
mpirun -np 4 gyoto --plugins=python3,stdplug -P-1 input.xml \!output.fits
```

For more information or for information on writing new classes in C++,
please refer to Sect. 9 of the Gyoto manual, available as
doc/user_guide/GyotoManual.pdf once Gyoto has been built on your
system, or at https://gyoto.obspm.fr/GyotoManual.pdf

## Copyright

This is part of Gyoto and has the same license (GPL v3+). It is
packaged in a sub-directory as an example for packaging Gyoto plug-ins.
