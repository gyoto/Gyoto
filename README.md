# GYOTO: the General relativitY Orbit Tracer of Observatoire de Paris

## What Gyoto is

[Gyoto](http://gyoto.obspm.fr/) aims at providing a framework for
computing orbits and ray-traced images in General relativity. It
consists in a library (libgyoto), utility programs, and a plug-in for
the Yorick programing language.

We request that use of Gyoto in scientific publications be properly
acknowledged. Please cite:

    GYOTO: a new general relativistic ray-tracing code, F. H. Vincent,
    T. Paumard, E. Gourgoulhon & G. Perrin, Classical and Quantum
    Gravity 28, 225011 (2011) [arXiv:1109.4769]

We also request that Gyoto modifications, extensions or plug-ins
leading to a scientific publication be made public as free software
reasonably fast (within one year after publication of the scientific
paper), for instance by contributing it directly to the Gyoto
code base. Contributors will be listed in the relevant source files as
well as in the AUTHORS file in the package.

## Copyright information

Gyoto is Copyright 2011-2019 Thibaut Paumard, Frédéric Vincent, Odele
Straub and Frédéric Lamy (To ease reading on non-UTF8 systems, French
accents are omitted in file headers).

Gyoto is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free
Software Foundation, either version 3 of the License, or (at your
option) any later version.

Gyoto is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received a copy of the [GNU General Public
License](COPYING) along with Gyoto.  If not, see
<http://www.gnu.org/licenses/>.

File python/doxyfile.py is Copyright 2008 Prabhu Ramachandran under
BSD style license. File python/numpy.i is Copyright (c) 2005-2015,
NumPy Developers under BSD 3-clause license. File bin/optionparser.h
is Copyright (C) 2012 Matthias S. Benkmann. See each file for details.

## Installation instructions

Refer to the file [INSTALL.Gyoto.md](INSTALL.Gyoto.md) for building
and installing Gyoto.

Several sample files are provided in doc/examples. You can ray-trace
those sceneries with:

  gyoto <input-file.xml> <output-file.fits>

FITS files can be read by a variety of free and proprietary
software. See http://heasarc.gsfc.nasa.gov/docs/heasarc/fits.html.

## Extending Gyoto

Custom metrics and astronomical objects can be added fairly easily by
implementing them as a Gyoto plug-in. This, of course, requires
knowledge of the C++ language. The [user
manual](http://gyoto.obspm.fr/GyotoManual.pdf) contains detailed
instructions.
                        -- Thibaut Paumard, Thu, 10 Jan 2019.
