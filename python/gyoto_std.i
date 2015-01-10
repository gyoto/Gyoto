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
#include "GyotoUniformSpectrometer.h"
#include "GyotoComplexSpectrometer.h"
#include "GyotoValue.h"

  // include Metric headers
#include "GyotoKerrBL.h"
#include "GyotoKerrKS.h"
#include "GyotoMinkowski.h"
// include Astrobj headers
#include "GyotoStandardAstrobj.h"
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
using namespace Gyoto;
%}

%array_class(double, array_double)
%array_class(unsigned long, array_unsigned_long)
%array_class(size_t, array_size_t)

// Typemaps to translate SmartPointer to specific classes
GyotoSmPtrTypeMapClassDerived(Astrobj, Standard)
GyotoSmPtrTypeMapClassDerived(Astrobj, UniformSphere)
GyotoSmPtrTypeMapClassDerived(Astrobj, Star)
GyotoSmPtrTypeMapClassDerived(Astrobj, StarTrace)
GyotoSmPtrTypeMapClassDerived(Astrobj, FixedStar)
GyotoSmPtrTypeMapClassDerived(Astrobj, Torus)
GyotoSmPtrTypeMapClassDerived(Astrobj, PageThorneDisk)
GyotoSmPtrTypeMapClassDerived(Astrobj, ThinDiskPL)
GyotoSmPtrTypeMapClassDerived(Astrobj, PolishDoughnut)
GyotoSmPtrTypeMapClassDerived(Astrobj, ThinDiskIronLine)
GyotoSmPtrTypeMapClassDerived(Astrobj, PatternDisk)
GyotoSmPtrTypeMapClassDerived(Astrobj, PatternDiskBB)
GyotoSmPtrTypeMapClassDerived(Astrobj, DynamicalDisk)
GyotoSmPtrTypeMapClassDerived(Astrobj, Disk3D)
GyotoSmPtrTypeMapClassDerived(Astrobj, DynamicalDisk3D)
GyotoSmPtrTypeMapClassDerived(Astrobj, DirectionalDisk)

GyotoSmPtrTypeMapClassDerived(Metric, KerrBL)
GyotoSmPtrTypeMapClassDerived(Metric, KerrKS)
GyotoSmPtrTypeMapClassDerived(Metric, Minkowski)

GyotoSmPtrTypeMapClassDerived(Spectrum, PowerLaw)
GyotoSmPtrTypeMapClassDerived(Spectrum, BlackBody)
GyotoSmPtrTypeMapClassDerived(Spectrum, ThermalBremsstrahlung)


%ignore Gyoto::Astrobj::Standard::Standard();
%ignore Gyoto::Astrobj::Standard::Standard(double radmax);
%ignore Gyoto::Astrobj::Standard::Standard(std::string kind);
%ignore Gyoto::Astrobj::Standard::Standard(const Standard& );
GyotoSmPtrClassDerivedPtrHdr(Astrobj, Standard, StandardAstrobj, GyotoStandardAstrobj.h)

%ignore Gyoto::Astrobj::UniformSphere::UniformSphere (std::string kind, SmartPointer<Metric::Generic> gg, double radius);
%ignore Gyoto::Astrobj::UniformSphere::UniformSphere (std::string kind);
%ignore Gyoto::Astrobj::UniformSphere::UniformSphere (const UniformSphere& orig);
GyotoSmPtrClassDerived(Astrobj, UniformSphere)

%extend Gyoto::Astrobj::Complex {
  Gyoto::Astrobj::Generic * __getitem__ (int i) {
    Gyoto::Astrobj::Generic * res = ($self)->operator[](i);
    res -> incRefCount();
    return res;
  }
 };
%extend Gyoto::Astrobj::Complex {
  void __setitem__(int i, Gyoto::Astrobj::Generic * p) {
    ($self)->operator[](i)=p;
  }
 };
GyotoSmPtrClassDerivedPtrHdr(Astrobj, Complex, ComplexAstrobj, GyotoComplexAstrobj.h)

// Declare specific classes
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

// Workaround cvar bug in Swig which makes help(gyoto_std) fail:
%inline {
  namespace GyotoStd {
    extern int __class__=0;
  }
}
