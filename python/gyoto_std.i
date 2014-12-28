%module(docstring="The Gyoto standard plug-in") gyoto_std
%import gyoto.i
%{

#define GYOTO_NO_DEPRECATED
#include "GyotoConfig.h"

#include "GyotoWorldline.h"
#include "GyotoPhoton.h"
#include "GyotoUniformSphere.h"
#include "GyotoPhoton.h"
#include "GyotoScenery.h"
#include "GyotoSpectrometer.h"
#include "GyotoComplexSpectrometer.h"
#include "GyotoValue.h"

  // include Metric headers
#include "GyotoKerrBL.h"
#include "GyotoKerrKS.h"
#include "GyotoMinkowski.h"
// include Astrobj headers
#include "GyotoComplexAstrobj.h"
#include "GyotoStar.h"
#include "GyotoStarTrace.h"
#include "GyotoFixedStar.h"
#include "GyotoTorus.h"
#include "GyotoThinDisk.h"
#include "GyotoPageThorneDisk.h"
#include "GyotoThinDiskPL.h"
#include "GyotoPolishDoughnut.h"
#include "GyotoThinDiskIronLine.h"

#include "GyotoPatternDisk.h"
#include "GyotoPatternDiskBB.h"
#include "GyotoDynamicalDisk.h"
#include "GyotoDisk3D.h"
#include "GyotoDynamicalDisk3D.h"
#include "GyotoDirectionalDisk.h"

// include Spectrum headers
#include "GyotoPowerLawSpectrum.h"
#include "GyotoBlackBodySpectrum.h"
#include "GyotoThermalBremsstrahlungSpectrum.h"
%}

%array_class(double, array_double)
%array_class(double, array_unsigned_long)

%rename(__getitem__) Gyoto::Astrobj::Complex::operator[];
%rename(ComplexAstrobjPtr) ComplexPtr;
GyotoSmPtrClassDerivedPtrHdr(Astrobj, Complex, ComplexAstrobj, GyotoComplexAstrobj.h)

%extend Gyoto::SmartPointer<Gyoto::Astrobj::Complex> {
  void __setitem__(int i, Gyoto::SmartPointer<Gyoto::Astrobj::Generic> p) {
    (*$self)->operator[](i)=p;
  }
 };


GyotoSmPtrClassDerived(Astrobj, Star)
GyotoSmPtrClassDerived(Astrobj, StarTrace)
GyotoSmPtrClassDerived(Astrobj, FixedStar)
GyotoSmPtrClassDerived(Astrobj, Torus)
GyotoSmPtrClassDerived(Astrobj, PageThorneDisk)
GyotoSmPtrClassDerived(Astrobj, ThinDiskPL)
GyotoSmPtrClassDerived(Astrobj, PolishDoughnut)
GyotoSmPtrClassDerived(Astrobj, ThinDiskIronLine)
GyotoSmPtrClassDerived(Astrobj, PatternDisk)
GyotoSmPtrClassDerived(Astrobj, PatternDiskBB)
GyotoSmPtrClassDerived(Astrobj, DynamicalDisk)
GyotoSmPtrClassDerived(Astrobj, Disk3D)
GyotoSmPtrClassDerived(Astrobj, DynamicalDisk3D)
GyotoSmPtrClassDerived(Astrobj, DirectionalDisk)

GyotoSmPtrClassDerived(Metric, KerrBL)
GyotoSmPtrClassDerived(Metric, KerrKS)
GyotoSmPtrClassDerived(Metric, Minkowski)

GyotoSmPtrClassDerivedHdr(Spectrum, PowerLaw, GyotoPowerLawSpectrum.h)
GyotoSmPtrClassDerivedHdr(Spectrum, BlackBody, GyotoBlackBodySpectrum.h)
GyotoSmPtrClassDerivedHdr(Spectrum, ThermalBremsstrahlung, GyotoThermalBremsstrahlungSpectrum.h)

