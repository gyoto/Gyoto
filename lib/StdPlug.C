/*
    Copyright 2011-2018 Thibaut Paumard

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
#include "GyotoComplexMetric.h"
#include "GyotoShift.h"
#include "GyotoKerrBL.h"
#include "GyotoKerrKS.h"
#include "GyotoMinkowski.h"
#include "GyotoChernSimons.h"
#include "GyotoRezzollaZhidenko.h"
#include "GyotoHayward.h"
#include "GyotoSchwarzschildHarmonic.h"

// include Astrobj headers
#include "GyotoComplexAstrobj.h"
#include "GyotoStar.h"
#include "GyotoStarTrace.h"
#include "GyotoFixedStar.h"
#include "GyotoFreeStar.h"
#include "GyotoInflateStar.h"
#include "GyotoTorus.h"
#include "GyotoDeformedTorus.h"
#include "GyotoOscilTorus.h"
#include "GyotoThinDisk.h"
#include "GyotoPageThorneDisk.h"
#include "GyotoThinDiskPL.h"
#include "GyotoPolishDoughnut.h"
#include "GyotoThinDiskIronLine.h"
#include "GyotoEquatorialHotSpot.h"
#include "GyotoJet.h"
#include "GyotoThickDisk.h"
#include "GyotoSphericalAccretion.h"
#include "GyotoBlob.h"
#include "GyotoPlasmoid.h"
#include "GyotoFlaredDiskSynchrotron.h"
#include "GyotoThinDiskProfile.h"
#include "GyotoThinDiskGridIntensity.h"

#include "GyotoPatternDisk.h"
#include "GyotoPatternDiskBB.h"
#include "GyotoDynamicalDisk.h"
#include "GyotoDynamicalDiskBolometric.h"
#include "GyotoDisk3D.h"
#include "GyotoDynamicalDisk3D.h"
#include "GyotoDirectionalDisk.h"
#include "GyotoXillverReflection.h"

// include Spectrum headers
#include "GyotoPowerLawSpectrum.h"
#include "GyotoBlackBodySpectrum.h"
#include "GyotoThermalBremsstrahlungSpectrum.h"
#include "GyotoThermalSynchrotronSpectrum.h"
#include "GyotoPowerLawSynchrotronSpectrum.h"
#include "GyotoKappaDistributionSynchrotronSpectrum.h"

using namespace Gyoto;

extern "C" void __GyotostdplugInit() {
  // Register Metrics
  Metric::Register("Complex",   &(Metric::Subcontractor<Metric::Complex>));
  Metric::Register("Shift", &(Metric::Subcontractor<Metric::Shift>));
  Metric::Register("KerrBL", &(Metric::Subcontractor<Metric::KerrBL>));
  Metric::Register("KerrKS", &(Metric::Subcontractor<Metric::KerrKS>));
  Metric::Register("Minkowski", &(Metric::Subcontractor<Metric::Minkowski>));
  Metric::Register("ChernSimons", &(Metric::Subcontractor<Metric::ChernSimons>));
  Metric::Register("RezzollaZhidenko", &(Metric::Subcontractor<Metric::RezzollaZhidenko>));
  Metric::Register("Hayward", &(Metric::Subcontractor<Metric::Hayward>));
  Metric::Register("SchwarzschildHarmonic", &(Metric::Subcontractor<Metric::SchwarzschildHarmonic>));
  // Register Astrobjs
  Astrobj::Register("Complex",   &(Astrobj::Subcontractor<Astrobj::Complex>));
  Astrobj::Register("Star",      &(Astrobj::Subcontractor<Astrobj::Star>));
  Astrobj::Register("StarTrace", &(Astrobj::Subcontractor<Astrobj::StarTrace>));
  Astrobj::Register("FixedStar", &(Astrobj::Subcontractor<Astrobj::FixedStar>));
  Astrobj::Register("FreeStar", &(Astrobj::Subcontractor<Astrobj::FreeStar>));
  Astrobj::Register("InflateStar",      &(Astrobj::Subcontractor<Astrobj::InflateStar>));
  Astrobj::Register("Torus",     &(Astrobj::Subcontractor<Astrobj::Torus>));
  Astrobj::Register("OscilTorus",
		    &(Astrobj::Subcontractor<Astrobj::OscilTorus>));
  Astrobj::Register("DeformedTorus",
		    &(Astrobj::Subcontractor<Astrobj::DeformedTorus>));
  Astrobj::Register("ThinDisk",  &(Astrobj::Subcontractor<Astrobj::ThinDisk>));
  Astrobj::Register("PageThorneDisk",
		    &(Astrobj::Subcontractor<Astrobj::PageThorneDisk>));
  Astrobj::Register("ThinDiskPL",  
		    &(Astrobj::Subcontractor<Astrobj::ThinDiskPL>));
  Astrobj::Register("PolishDoughnut",
		    &(Astrobj::Subcontractor<Astrobj::PolishDoughnut>));
  Astrobj::Register("ThinDiskIronLine",  
		    &(Astrobj::Subcontractor<Astrobj::ThinDiskIronLine>));
  Astrobj::Register("EquatorialHotSpot",
		    &(Astrobj::Subcontractor<Astrobj::EquatorialHotSpot>));
  Astrobj::Register("PatternDisk",
		    &(Astrobj::Subcontractor<Astrobj::PatternDisk>));
  Astrobj::Register("PatternDiskBB",
		    &(Astrobj::Subcontractor<Astrobj::PatternDiskBB>));
  Astrobj::Register("DynamicalDisk",
		    &(Astrobj::Subcontractor<Astrobj::DynamicalDisk>));
  Astrobj::Register("DynamicalDiskBolometric",
		    &(Astrobj::Subcontractor<Astrobj::DynamicalDiskBolometric>));
  Astrobj::Register("Disk3D",
		    &(Astrobj::Subcontractor<Astrobj::Disk3D>));
  Astrobj::Register("DynamicalDisk3D",
		    &(Astrobj::Subcontractor<Astrobj::DynamicalDisk3D>));
  Astrobj::Register("DirectionalDisk",
		    &(Astrobj::Subcontractor<Astrobj::DirectionalDisk>));
  Astrobj::Register("Jet",
		    &(Astrobj::Subcontractor<Astrobj::Jet>));
  Astrobj::Register("ThickDisk",
		    &(Astrobj::Subcontractor<Astrobj::ThickDisk>));
  Astrobj::Register("SphericalAccretion",
		    &(Astrobj::Subcontractor<Astrobj::SphericalAccretion>));
  Astrobj::Register("ThinDiskProfile",
		    &(Astrobj::Subcontractor<Astrobj::ThinDiskProfile>));
  Astrobj::Register("Blob",
		    &(Astrobj::Subcontractor<Astrobj::Blob>));
  Astrobj::Register("Plasmoid",
		    &(Astrobj::Subcontractor<Astrobj::Plasmoid>));
  Astrobj::Register("XillverReflection",
		    &(Astrobj::Subcontractor<Astrobj::XillverReflection>));
  Astrobj::Register("FlaredDiskSynchrotron",
		    &(Astrobj::Subcontractor<Astrobj::FlaredDiskSynchrotron>));
  Astrobj::Register("ThinDiskGridIntensity",
		    &(Astrobj::Subcontractor<Astrobj::ThinDiskGridIntensity>));
  // Register Spectra
  Spectrum::Register("PowerLaw", 
		     &(Spectrum::Subcontractor<Spectrum::PowerLaw>));
  Spectrum::Register("BlackBody", 
		     &(Spectrum::Subcontractor<Spectrum::BlackBody>));
  Spectrum::Register("ThermalBremsstrahlung", 
		     &(Spectrum::Subcontractor<Spectrum::ThermalBremsstrahlung>));
  Spectrum::Register("ThermalSynchrotron", 
		     &(Spectrum::Subcontractor<Spectrum::ThermalSynchrotron>));
  Spectrum::Register("PowerLawSynchrotron", 
		     &(Spectrum::Subcontractor<Spectrum::PowerLawSynchrotron>));
  Spectrum::Register("KappaDistributionSynchrotron", 
		     &(Spectrum::Subcontractor<Spectrum::KappaDistributionSynchrotron>));
}
