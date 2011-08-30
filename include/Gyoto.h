/*
    Copyright 2011 Thibaut Paumard

    This file is part of Gyoto.

    Gyoto is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Gyoto is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Gyoto.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * \file Gyoto.h
 * \brief All of Gyoto
 * Gyoto.h #includes all of the other Gyoto C++ header files.
 */

#include<GyotoUtils.h>
#include<GyotoAstrobj.h>
#include<GyotoError.h>
#include<GyotoFocalPlane.h>
#include<GyotoMetric.h>
#include<GyotoSmartPointer.h>
#include<GyotoThinInfiniteDisk.h>
#include<GyotoWorldlineTime.h>
#include<GyotoCoord.h>
#include<GyotoFixedStar.h>
#include<GyotoKerr.h>
#include<GyotoPhoton.h>
#include<GyotoStar.h>
#include<GyotoWorldline.h>

/**
 * \mainpage Gyoto
<B><CENTER> General relativitY Orbit Tracer of Observatoire de Paris</CENTER></B>
 *
 * \section intro_sec Introduction
 *
 * Gyoto aims at providing a framework for computing orbits and
 * ray-traced images in General relativity. It consists of a library
 * (libgyoto), utility programs, and a plug-in for the Yorick
 * programing language.
 *
 * A graphical user interface for tracing stellar orbits is provided
 * with the Yorick plug-in (see \ref yorick_plugin_page).
 *
 * \image html gyotoy_screenshot.png
 *
 * \section downloading_sec Downloading Gyoto
 *
 * Gyoto can be downloaded at https://github.com/gyoto/Gyoto .
 *
 * \section building_sec Building Gyoto 
 *
 * \subsection complib_sec Compiling the library and the utilities 
 *
 * Just run make in the main Gyoto directory to build the library
 * (lib/libgyoto.a) and the utility programs (in bin/):
 * \code make \endcode
 *
 * \subsection compyo_sec Compiling the Yorick plug-in
 *
 * See the building section in \ref yorick_plugin_page.
 *
 * \subsection compdoc_sec Compiling (this) documentation
 *
 * From the main Gyoto directory, run
 * \code make doc \endcode
 * 
 * \section install_sec Installation
 *
 * Not provided yet.
 *  
 * \section philo_sec Philospophy and concepts
 *
 * todo...
 */

/**
 * \page yorick_plugin_page Gyoto for Yorick
 *
 * \section yoprerec_sec Prerequisites
 *
 * \subsection yo_sec Yorick
 *
 * For direct access to Gyoto, only Yorick is required. Yorick comes
 * precompiled with many distributions of linux and is present in the
 * Macports. You can also find binaries and source code at
 * http://yorick.sourceforge.net/
 *
 * One line in Yorick's yapi.h file makes the compilation fail, so
 * either use Yorick > 2.1.05 (i.e. from CVS, as of writting), or
 * change this line: PLUG_API void *ygeta_any(int iarg, long *ntot,
 * long *dims, int *typeid); to something like PLUG_API void
 * *ygeta_any(int iarg, long *ntot, long *dims, int *typeidd); As of
 * writing, a corrected yapi.h (borrowed from yorick 2.1.04) is
 * provided.
 *
 * For the graphical user interface Gyotoy, you'll need a few
 * additional pieces. First of all, I seriously doubt it will run
 * under MS Windows, but please report any news on that front.


\subsection yutils_sec Yutils

A collection of utilities for Yorick, distributed as yorick-yutils in
many Linux distributions. Actually, only pyk.i is used in Gyotoy (so
far), and as of writting, the CVS version is required.

Either get yutils > 1.4.0 from e.g. http://www.maumae.net/yorick/packages/src/
or fetch pyk.i from the CVS browser at
http://yorick.cvs.sourceforge.net/viewvc/yorick/yorick-yutils/pyk.i?view=markup


\subsection pygtk_sec PyGTK

The graphical interface is actually coded python using the pygtk
extension. So you'll need a recent version of python, PyGTK, and the
Glade library, on which PyGTK depends.


\section building_sec Building

Building the Gyoto plug-in is done the standard way for Yorick plug-ins:
\code
yorick -batch make.i
make
\endcode

You can check the package by running
make check
or by trying out the graphical interface, if you have installed all
the dependencies:
./gyotoy.i

\section installing_sec Installing

The install procedure is not yet tested.
\code
make install
\endcode
should do, mostly.

\section running_sec Running

Read gyoto.i and check.i for using Gyoto from within Yorick.

The graphical interface, Gyotoy, can be run with
yorick -i gyotoy.i

If everything goes well,
./gyotoy.i
should do, too.


 */
