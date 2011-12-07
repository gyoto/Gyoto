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
// include Spectrum headers
#include "GyotoPowerLawSpectrum.h"
#include "GyotoBlackBodySpectrum.h"

extern "C" void __GyotostdplugInit() {
  // Register Metrics
  Gyoto::Metric::KerrBL::Init();
  Gyoto::Metric::KerrKS::Init();
  // Register Astrobjs
  Gyoto::Astrobj::Complex::Init();
  Gyoto::Astrobj::Star::Init();
  Gyoto::Astrobj::FixedStar::Init();
  Gyoto::Astrobj::Torus::Init();
  Gyoto::Astrobj::ThinDisk::Init();
  Gyoto::Astrobj::PageThorneDisk::Init();
  // Register Spectra
  Gyoto::Spectrum::PowerLawInit();
  Gyoto::Spectrum::BlackBodyInit();
}
