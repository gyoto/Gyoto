/*
    Copyright 2011-2016, 2018-2019 Thibaut Paumard, Frédéric Vincent,
                                   Éric Gourgoulhon 

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
 *
 * Gyoto.h should include all of the other Gyoto C++ header
 * files. Note that at any given time it may be lacking some new
 * headers. It is always best to specifically include the headers that
 * you actually use.
 *
 * Beware that some symbols are resolved only in the
 * plugins (libgyoto-stdplug.so, libgyoto-lorene.so).
 */

#include<GyotoUtils.h>
#include<GyotoError.h>
#include<GyotoSmartPointer.h>
#include<GyotoWorldline.h>
#include<GyotoPhoton.h>

#include<GyotoMetric.h>
#include<GyotoKerrBL.h>
#include<GyotoKerrKS.h>
#include<GyotoChernSimons.h>
#include<GyotoRezzollaZhidenko.h>

#include<GyotoAstrobj.h>
#include<GyotoThinDisk.h>
#include<GyotoPageThorneDisk.h>
#include<GyotoDirectionalDisk.h>
#include<GyotoPatternDisk.h>
#include<GyotoPatternDiskBB.h>
#include<GyotoDynamicalDisk.h>
#include<GyotoDisk3D.h>
#include<GyotoFixedStar.h>
#include<GyotoInflateStar.h>
#include<GyotoStar.h>
#include<GyotoDeformedTorus.h>
#include<GyotoEquatorialHotSpot.h>

/**
 * \mainpage Gyoto
<B><CENTER> General relativitY Orbit Tracer of Observatoire de Paris</CENTER></B>
 *
 *
 \ref download_sec | \ref manual_sec | \ref conditions_sec
 *
 * Gyoto aims at providing a framework for computing orbits and
 * ray-traced images in General relativity. It consists of a library
 * (libgyoto), utility programs, a plug-in for the Yorick programing
 * language, and an extension for the Python 3 programing
 * language. Gyoto is known to run under Linux and Mac OS X. Please do
 * tell us if you manage to run Gyoto under a different OS. It should
 * compile and run with moderate effort on most UNIX-like systems.
 *
 * Gyoto can be expanded with plug-ins providing custom
 * Gyoto::Metric::Generic, Gyoto::Astrobj::Generic and
 * Gyoto::Spectrum::Generic classes, which describe respectively
 * analytical or numerical metrics, astronomical objects, and spectral
 * shapes for astronomical objects (see \ref
 * writing_plugins_page). Custom Gyoto::Metric::Generic,
 * Gyoto::Astrobj::Generic and Gyoto::Spectrum::Generic classes can
 * also be written in the <A HREF="https://www.python.org/">Python</A>
 * 3.x interpreted language using the \c python plug-in for Gyoto, see
 * below. Beware that a Python implementation of a custom class will
 * run significantly slower than the equivalent C++ implementation,
 * but sometimes saving on development time is better than saving on
 * computing time.
 *
 * The base distribution includes three plug-ins: the standard plug-in
 * (\c stdplug), a plug-in for using <A
 * HREF="http://www.lorene.obspm.fr/">LORENE</A>-based
 * numerical metrics (\c lorene) and a plug-in for writing custom
 * astronomical objects or metrics in the <A
 * HREF="https://www.python.org/">Python</A> 3.x interpreted
 * language (\c python).
 *
 * A graphical user interface for tracing stellar orbits is provided
 * with the Yorick plug-in (see \ref yorick_plugin_page).
 *
 * To visit the code
 * <span style="color: #ff0000"><b>Picture Gallery</b></span>:
 * click <A HREF="gallery/index.html">here</A>!
 * 
 *
 * \image html gyotoy_screenshot.png
 *
 * \section download_sec Downloading and installing
 *
 * Detailed information on installing Gyoto is available <A
 * HREF="http://github.com/gyoto/Gyoto/blob/master/INSTALL.Gyoto.md"><B>here</B></A>. The user manual below is also a valuable read.
 *
 * \section manual_sec User manual
 *
 * The user manual is available <A HREF="GyotoManual.pdf"><B>here</B></A>.
 *
 * \section conditions_sec Conditions for use
 * 
 * We request that use of Gyoto in scientific publications be properly
 * acknowledged. Please cite:
 *
 *  F. H. Vincent, T. Paumard, E. Gourgoulhon & G. Perrin: 
 *  <EM>GYOTO: a new general relativistic ray-tracing code</EM>, 
 *  Classical and Quantum Gravity <STRONG>28</STRONG>, 225011 (2011)
 *  [<A HREF="http://dx.doi.org/10.1088/0264-9381/28/22/225011">published version</A>]
 *  [<A HREF="http://arxiv.org/abs/1109.4769">preprint: arXiv:1109.4769</A>]
 *
 * We also request that Gyoto modifications, extensions or plug-ins
 * leading to a scientific publication be made public as free software
 * reasonably fast (within one year after publication of the scientific
 * paper), for instance by contributing it directly to the Gyoto
 * code base. Contributors will be listed in the relevant source files as
 * well as in the AUTHORS file in the package.
 *
 *   Gyoto is Copyright 2011-2016 Thibaut Paumard,
 *   Fr&eacute;d&eacute;ric Vincent and Odele Straub.
 *
 *  Gyoto is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 */

/**
\page yorick_plugin_page Gyoto for Yorick

\section yoprerec_sec Prerequisites

\subsection yo_sec Yorick

For direct access to Gyoto, only Yorick (>=2.2) is required. Yorick comes
precompiled with many distributions of Linux and is present in the
Macports. You can also find source code at https://github.com/dhmunro/yorick

For the graphical user interface Gyotoy, you'll need a few additional
pieces. First of all, I seriously doubt it will run under MS Windows,
but please report any news on that front.

\subsection yutils_sec Yutils

A collection of utilities for Yorick, distributed as yorick-yutils in
many Linux distributions. Actually, only pyk.i is used in Gyotoy (so
far), and as of writing, the CVS version is required.

Either get yutils > 1.4.0 from e.g. http://www.maumae.net/yorick/packages/src/
or fetch pyk.i from the CVS browser at
http://yorick.cvs.sourceforge.net/viewvc/yorick/yorick-yutils/pyk.i?view=markup


\subsection pygtk_sec PyGObject

The graphical interface gyotoy needs yorick-gy and the Gtk 3 GObject
Introspection files. Get yorick-gy from:
https://github.com/paumard/yorick-gy

Under Mac OS X (using macports), install the gtk3 package. Under a
Debian derivative GNU/Linux system, install gir1.2-gtk-3.0.  Check
yorick/README.gyotoy in the source distribution for more details.

\section building_sec Building

From the Gyoto source directory (i.e. the directory which contains the yorick/ subdirectory:
\code
make yorick
\endcode

You can check the package by running
\code
make check-yorick
\endcode

\section installing_sec Installing

\code
sudo make install-yorick
\endcode

\section running_sec Running

Read gyoto.i and check.i for using Gyoto from within Yorick.

The graphical interface, Gyotoy, can be run with
\code
yorick -i gyotoy.i
\endcode
or
\code
gyotoy
\endcode


 */

/**
 \page writing_plugins_page Writing plug-ins for Gyoto

 Al the generic Gyoto machinery for computing orbits, reading input
 files, performing ray-tracing etc. is implemented in libgyoto. On the
 other hand, all the code specific to a given metric kind or a given
 astronomical object is available in plug-ins. For instance, the
 standard plug-in libgyoto-stdplug contains the two flavors of the
 Kerr metric: Gyoto::KerrBL and Gyoto::KerrKS, as well as basic
 objects: Gyoto::FixedStar, Gyoto::Star, Gyoto::Torus,
 Gyoto::ThinInfiniteDiskBL and Gyoto::ThinInfiniteDiskKS. The
 libgyoto-lorene plug-in contains the code to access numerical metrics
 (Gyoto::LoreneMetric) as well as an example thereof:
 Gyoto::RotStar3_1. The two basic spectral shapes
 Gyoto::Spectrum::PowerLaw and Gyoto::Spectrum::BlackBody are also to
 be found in the standard plug-in.

 Gyoto can be used right away to compute stellar orbits in the Kerr
 metric or to do basic ray-tracing of accretion disks. But Gyoto is
 not limited to the basic metrics and objects we have thought of. It
 is fairly easy to add custom metrics and objects (and
 emission/absorption laws) it Gyoto, and Gyoto itself does not need to
 be modified or even re-compiled to do so: custom classes can (and
 should) be implemented as plug-ins. For an example, simply look at
 the lib/StdPlug.C file in the source distribution, and the source
 files for the objects and metrics it provides: e.g. lib/FixedStar.C and
 lib/KerrBL.C.

 To implement a new plug-in, you first need to implement a derived
 class of either the Gyoto::Astrobj, Gyoto::Metric, or
 Gyoto::Spectrum::Generic class. You don't necessarily need to
 implement everything, the Gyoto::Astrobj page explains what is
 required for an astronomical object.

 Assuming you want to be able to actually use your custom class, you
 need a way to instantiate it. This is normally the job of the
 Gyoto::Factory. You need to instruct the Gyoto::Factory how to read
 parameters for your specific class from an XML file by implementing a
 subcontractor (for a Gyoto::Astrobj, the subcontractor is a static
 method of the Gyoto::Astrobj::Subcontractor_t type). The
 subcontractor communicates with the Gyoto::Factory by means of a
 Gyoto::FactoryMessenger and basically loops calling the
 Gyoto::FactoryMessenger::getNextParameter() method (see the
 GyotoRegister.h file, unfortunately undocumented at the moment).

 You also need to register your subcontractor, so that the
 Gyoto::Factory knows it must call it when it encounters a tag of the
 form &lt;Astrobj kind=&quot;YourKind&quot;&gt; in an XML
 description. This is typically done by providing an static Init
 method in your class:
 \code
 void Gyoto::FixedStar::Init() {
   Gyoto::Astrobj::Register("FixedStar", &Gyoto::FixedStar::Subcontractor);
 }
 \endcode

 You need to make sure this Init() method is called when your plug-in
 is loaded. Assume you decide to call your plug-in MyPlug, and it
 contains a single Gyoto::Astrobj named Gyoto::MyObj. You will compile
 it under the file name libgyoto-MyPlug.so (under Linux) or
 libgyoto-MyPlug.dylib (under MacOS X). Just put this file somewhere
 where the dynamic linker can find it (any directory listed in
 $LD_LIBRARY_PATH or $DYLD_LIBRARY_PATH will be fine; /usr/local/lib/
 should also be fine). In addition to the implementation of the
 Gyoto::MyObj class, you will need to provide a function called
 __GyotoMyPlugInit() which will be exactly this:
 \code
 extern "C" void __GyotostdplugInit() {
   Gyoto::MyObj::Init();
 }
 \endcode
 This function is typically provided in a separate source file (such
 as lib/StdPlug.C in the Gyoto source) and can initialize several
 custom classes at once.

 Finally, you need to instruct Gyoto to load your plug-in at run
 time. This is done by adding the name of your plug-in to the
 GYOTO_PLUGINS environment variable. The default value for
 GYOTO_PLUGINS is "stdplug,nofail:lorene", meaning Gyoto should load
 the standard plug-in stdplug and attempt to load the lorene plug-in,
 failing only if stdplug is nowhere to be found. If you want to load
 your plug-in in addition to those, alter this variable in your shell
 (if you don't know what this means or how to do this, ask the local
 Unix guru or read the fine <A
 HREF="http://www.gnu.org/s/bash/manual/bash.html">manual</A>):
 \code
 export GYOTO_PLUGINS="stdplug,nofail:lorene,MyPlug"
 \endcode
 but if your lug-in is self-contained and your don't need the objects
 in the standard plug-ins, this will do it for you:
 \code
 export GYOTO_PLUGINS="MyPlug"
 \endcode
 This will instruct Gyoto to locate and load the file named
 libgyoto-MyPlug.(so|dylib) and to run the function named
 __GyotostdplugInit() from this library file.


 */
