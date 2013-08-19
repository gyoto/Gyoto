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

// include Metric headers
#include "GyotoKerrBL.h"
#include "GyotoKerrKS.h"
// include Astrobj headers
#include "GyotoComplexAstrobj.h"
#include "GyotoStar.h"
#include "GyotoFixedStar.h"
#include "GyotoTorus.h"
#include "GyotoThinDisk.h"
#include "GyotoPageThorneDisk.h"
#include "GyotoThinDiskPL.h"
#include "GyotoPolishDoughnut.h"
#include "GyotoThinDiskIronLine.h"
#ifdef GYOTO_USE_CFITSIO
#include "GyotoPatternDisk.h"
#include "GyotoPatternDiskBB.h"
#include "GyotoDynamicalDisk.h"
#include "GyotoDisk3D.h"
#include "GyotoDisk3D_BB.h"
#endif
// include Spectrum headers
#include "GyotoPowerLawSpectrum.h"
#include "GyotoBlackBodySpectrum.h"

using namespace Gyoto;

extern "C" void __GyotostdplugInit() {
#if defined GYOTO_USE_XERCES
  // Register Metrics
  Metric::Register("KerrBL", &(Metric::Subcontractor<Metric::KerrBL>));
  Metric::Register("KerrKS", &(Metric::Subcontractor<Metric::KerrKS>));
  // Register Astrobjs
  Astrobj::Register("Complex",   &(Astrobj::Subcontractor<Astrobj::Complex>));
  Astrobj::Register("Star",      &(Astrobj::Subcontractor<Astrobj::Star>));
  Astrobj::Register("FixedStar", &(Astrobj::Subcontractor<Astrobj::FixedStar>));
  Astrobj::Register("Torus",     &(Astrobj::Subcontractor<Astrobj::Torus>));
  Astrobj::Register("ThinDisk",  &(Astrobj::Subcontractor<Astrobj::ThinDisk>));
  Astrobj::Register("PageThorneDisk",
		    &(Astrobj::Subcontractor<Astrobj::PageThorneDisk>));
  Astrobj::Register("ThinDiskPL",  
		    &(Astrobj::Subcontractor<Astrobj::ThinDiskPL>));
  Astrobj::Register("PolishDoughnut",
		    &(Astrobj::Subcontractor<Astrobj::PolishDoughnut>));
  Astrobj::Register("ThinDiskIronLine",  
		    &(Astrobj::Subcontractor<Astrobj::ThinDiskIronLine>));
#ifdef GYOTO_USE_CFITSIO
  Astrobj::Register("PatternDisk",
		    &(Astrobj::Subcontractor<Astrobj::PatternDisk>));
  Astrobj::Register("PatternDiskBB",
		    &(Astrobj::Subcontractor<Astrobj::PatternDiskBB>));
  Astrobj::Register("DynamicalDisk",
		    &(Astrobj::Subcontractor<Astrobj::DynamicalDisk>));
  Astrobj::Register("Disk3D",
		    &(Astrobj::Subcontractor<Astrobj::Disk3D>));
  Astrobj::Register("Disk3D_BB",
		    &(Astrobj::Subcontractor<Astrobj::Disk3D_BB>));
#endif
  // Register Spectra
  Spectrum::Register("PowerLaw", 
		     &(Spectrum::Subcontractor<Spectrum::PowerLaw>));
  Spectrum::Register("BlackBody", 
		     &(Spectrum::Subcontractor<Spectrum::BlackBody>));
#endif
}
