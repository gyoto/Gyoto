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
 * (libgyoto), utility programs, and an extension for the Python 3
 * programing language. Gyoto is known to run under Linux and Mac OS
 * X. Please do tell us if you manage to run Gyoto under a different
 * OS. It should compile and run with moderate effort on most
 * UNIX-like systems.
 *
 * Gyoto can be expanded with plug-ins providing custom
 * Gyoto::Metric::Generic, Gyoto::Astrobj::Generic and
 * Gyoto::Spectrum::Generic classes, which describe respectively
 * analytical or numerical metrics, astronomical objects, and spectral
 * shapes for astronomical objects (see the <A
 * HREF="GyotoManual.pdf">user manual</A>). Custom
 * Gyoto::Metric::Generic, Gyoto::Astrobj::Generic and
 * Gyoto::Spectrum::Generic classes can also be written in the <A
 * HREF="https://www.python.org/">Python</A> 3.x interpreted language
 * using the \c python plug-in for Gyoto. Beware that a Python
 * implementation of a custom class will run significantly slower than
 * the equivalent C++ implementation, but sometimes saving on
 * development time is better than saving on computing time.
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
 * with the Python extension.
 *
 * To visit the code
 * <span style="color: #ff0000"><b>Picture Gallery</b></span>:
 * click <A HREF="gallery/index.html">here</A>!
 * 
 *
 * \image html gyotoy_screenshot.png width=50%
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
