/*
    Copyright 2014 Thibaut Paumard

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

/*
    This is a Swig interface file. It is currently meant to provide
    Python bindings only, but it should ne be too difficult to provide
    bindings for java, Tcl or whatever other language Swig supports.
 */

// ******** SOME INITIALIZATION ********* //

// Define the module name, with a docstring
%module(docstring="The General relativitY Orbit Tracer of paris Observatory") gyoto

// Let Swig generate some documentation for the overloaded methods
%feature("autodoc", "1");

// gyoto_doc.i is generated from the doxygen comments using doxy2swig.py
%import gyoto_doc.i

// Make it possible to detect that a .h file is being processed by
// Swig rather than CPP. There should be a Swig feature for that, I
// didn't find it.
%define GYOTO_SWIGIMPORTED
%enddef

// Make sure we don't activate the deprecated method names (see
// GyotoDefs.h): it is better to fail on them to update them now.
%define GYOTO_NO_DEPRECATED
%enddef

//  ********** MACRO DEFINITIONS *********** //

// Exposing SmartPointers does not work well: it breaks automatic
// inheritance and docstring generation. We therefore provide Swig
// typemaps to convert SmartPointers to normal pointers
// automatically. We also use the ref/unref features to let Swig
// handle the reference counting for us mostly automatically (which is
// what we wouldotherwise loose by not using SmartPointers).
//
// Since SmartPointers are templates, we use macros here to provide
// typemaps for various SmartPointer specializations.

// Typemaps for Gyoto::SmartPointer<Gyoto::klass>
%define GyotoSmPtrTypeMapClass(klass)
  GyotoSmPtrTypeMap(Gyoto::klass, SWIGTYPE_p_Gyoto__ ## klass)
%enddef

// Typemaps for Gyoto::SmartPointer<Gyoto::nspace::Generic>
%define GyotoSmPtrTypeMapClassGeneric(nspace)
  GyotoSmPtrTypeMap(Gyoto::nspace::Generic, SWIGTYPE_p_Gyoto__ ## nspace ## __Generic)
%enddef

// Typemaps for Gyoto::SmartPointer<Gyoto::nspace::klass>
%define GyotoSmPtrTypeMapClassDerived(nspace, klass)
   GyotoSmPtrTypeMap(Gyoto::nspace::klass, SWIGTYPE_p_Gyoto__ ## nspace ## __ ## klass)
%enddef

// Basic macro used in the above: gtype is the Gyoto name of a type,
// e.g. Gyoto::Metric::KerrBL, while stype is the Swig name for same
// type, e.g. SWIGTYPE_p_Gyoto__Metric__KerrBL
%define GyotoSmPtrTypeMap(gtype, stype)
%typemap(in)
  Gyoto::SmartPointer<gtype>  (gtype *)
{
  int res=0;
  void *argp=0;
  res=SWIG_ConvertPtr($input, &argp, stype, 0);
  if (!SWIG_IsOK(res)) {
    SWIG_exception_fail(SWIG_ArgError(res), "argument of type '" #gtype "*'");
  }
  gtype * kp=reinterpret_cast< gtype * >(argp);
  $1 = Gyoto::SmartPointer<gtype>(kp);
}
%typemap(out)
  Gyoto::SmartPointer<gtype>  (gtype *)
{
  gtype* normal_pointer=(gtype *) (Gyoto::SmartPointer<gtype>(result));
  normal_pointer->incRefCount();
  $result = SWIG_NewPointerObj( normal_pointer, stype, SWIG_POINTER_OWN |  0 );
}
%typemap(typecheck) Gyoto::SmartPointer<gtype>, gtype * {
  void *vptr = 0;
  int res = SWIG_ConvertPtr($input, &vptr, stype, 0);
  $1 = SWIG_CheckState(res);
 }
%typemap(typecheck) Gyoto::SmartPointer< gtype > = Gyoto::SmartPointer<gtype>;
%typemap(typecheck) gtype * = Gyoto::SmartPointer<gtype>;
%typemap(typecheck) gtype const * = Gyoto::SmartPointer<gtype>;
%enddef

// Include header for a class deriving from SmartPointee, providing
// the ref and unref features
%define GyotoSmPtrClass(klass)
%feature("ref") Gyoto:: klass "$this->incRefCount();";
%feature("unref") Gyoto:: klass "$this->decRefCount(); if (!$this->getRefCount()) delete $this;";
%include Gyoto ## klass ## .h
%enddef

// Include header for a base class (e.g. GyotoMetric.h), provide
// constructor from kind string (e.g. gyoto.Metric('KerrBL')), provide
// ref/unref features
%define GyotoSmPtrClassGeneric(klass)
%rename(klass) klass ## Ptr;
%ignore Gyoto::klass::Register_;
%ignore Gyoto::klass::Register;
%ignore Gyoto::klass::initRegister;
%ignore Gyoto::klass::getSubcontractor;
%rename(klass) Gyoto::klass::Generic;
%feature("ref") Gyoto:: klass ::Generic"$this->incRefCount();";
%feature("unref") Gyoto:: klass ::Generic"$this->decRefCount(); if (!$this->getRefCount()) delete $this;";
%feature("notabstract") Gyoto::klass::Generic;
%ignore  Gyoto::klass::Generic::Generic(Gyoto::klass::Generic const &);
%ignore  Gyoto::klass::Generic::Generic(const Generic &);
%ignore  Gyoto::klass::Generic::Generic(const klass::Generic &);
%ignore  Gyoto::klass::Generic::Generic();
%ignore  Gyoto::klass::Generic::Generic(double);
%ignore  Gyoto::klass::Generic::Generic(kind_t);
%ignore  Gyoto::klass::Generic::Generic(const std::string);
%extend Gyoto::klass::Generic {
  Generic(std::string nm) {
    Gyoto::SmartPointer<Gyoto::klass::Generic> pres=
      Gyoto::klass::getSubcontractor(nm.c_str())(NULL);
    Gyoto::klass::Generic * res = (Gyoto::klass::Generic *)(pres);
    res -> incRefCount();
    return res;
  }
};
%include Gyoto ## klass ## .h
%enddef

// Include header for derived class. Parameters: nspace: namespace
// (e.g. Astrobj); klass: classname (e.g. Complex); nick: name to use
// in the flattened namespace (e.g. ComplexAstrobj); hdr: header file
// (e.g. GyotoComplexAstrobj.h). Extends class with a pseudo
// constructor for upcasting,
// e.g. cplx=gyoto_std.ComplexAstrobj(sc.strobj())
%define GyotoSmPtrClassDerivedPtrHdr(nspace, klass, nick, hdr)
%rename(nick) Gyoto::nspace::klass;
%feature("notabstract") Gyoto::nspace::klass;
%extend  Gyoto::nspace::klass {
  klass(Gyoto::nspace::Generic * base) {
    Gyoto::nspace::klass * res = dynamic_cast< Gyoto::nspace::klass * >(base);
    if (!res) Gyoto::throwError("This pointer cannot be cast to 'Gyoto::" #nspace "::" #klass "*'");
    return res;
  }
 };
%include hdr
%enddef

// Simplification of the above when nick == klass
%define GyotoSmPtrClassDerivedHdr(nspace, klass, hdr)
GyotoSmPtrClassDerivedPtrHdr(nspace, klass, klass, hdr)
%enddef

// Simplification of the above when hdr == Gyoto ## klass ## .h
%define GyotoSmPtrClassDerived(nspace, klass)
GyotoSmPtrClassDerivedHdr(nspace, klass, Gyoto ## klass ## .h)
%enddef


// ******** INCLUDES ******** //
// Include any file that is needed to compile the wrappers
%{
#define SWIG_FILE_WITH_INIT
  //#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "Python.h"
#define GYOTO_NO_DEPRECATED
#include "GyotoConfig.h"
#include "GyotoFactory.h"
#include "GyotoValue.h"
#include "GyotoProperty.h"
#include "GyotoObject.h"
#include "GyotoAstrobj.h"
#include "GyotoError.h"
#include "GyotoWorldline.h"
#include "GyotoPhoton.h"
#include "GyotoScreen.h"
#include "GyotoStandardAstrobj.h"
#include "GyotoUniformSphere.h"
#include "GyotoThinDisk.h"
#include "GyotoSpectrometer.h"
#include "GyotoComplexSpectrometer.h"
#include "GyotoUniformSpectrometer.h"
#include "GyotoRegister.h"
using namespace Gyoto;

%}

// ******** INITIALIZATION ******** //
// This will be called upon extension initialization
%init {
  Gyoto::Register::init();
  import_array();
 }

// ******** TYPEMAPS ******** //
// Actually instanciate typemaps using de macros defined above

GyotoSmPtrTypeMapClassGeneric(Metric);
GyotoSmPtrTypeMapClassGeneric(Astrobj);
GyotoSmPtrTypeMapClassGeneric(Spectrum);
GyotoSmPtrTypeMapClassGeneric(Spectrometer);

GyotoSmPtrTypeMapClass(Screen);
GyotoSmPtrTypeMapClass(Scenery);
GyotoSmPtrTypeMapClass(Photon);

GyotoSmPtrTypeMapClassDerived(Spectrometer, Complex);
GyotoSmPtrTypeMapClassDerived(Spectrometer, Uniform);

// Non-Gyoto typemaps:
// Handle std::string
%include "std_string.i";
// Handle std::vector<double> and <unsigned long int>
%include "std_vector.i";
%template(vector_double) std::vector<double>;
%template(vector_unsigned_long) std::vector<unsigned long>;
// Handle generic C arrays using a class-like interface
%include "carrays.i"
%array_class(double, array_double);
%array_class(double, array_unsigned_long);
// Handle some arrays as NumPy arrays
%include "numpy.i";
%numpy_typemaps(size_t, NPY_ULONG , size_t);
%numpy_typemaps(double, NPY_DOUBLE, size_t);

// ******** INTERFACE ******** //
// Here starts the actual pasing of the various header files

// Expose the build-time configuration variables
%include "GyotoConfig.h"

// Expose the Gyoto::Error class
%include "GyotoError.h"

// Catch all Gyoto errors and re-throw them as a Python run-time error
%exception {
	try {
	$function
	}
	catch (Gyoto::Error e) {
		PyErr_SetString(PyExc_RuntimeError, e);
		return NULL;
	}
}

// Expose Gyoto::Register::list as gyoto.listRegister
%ignore Gyoto::Register::Entry;
%ignore Gyoto::Register::init;
%rename(listRegister) Gyoto::Register::list;
%include GyotoRegister.h

%ignore Gyoto::Functor::Double_constDoubleArray;
%ignore Gyoto::Functor::Double_Double_const;
%include "GyotoFunctors.h"

%ignore Gyoto::Hook::Listener;
%ignore Gyoto::Hook::Teller;
%include "GyotoHooks.h"

%ignore Gyoto::WIP;
%include "GyotoWIP.h"

%ignore Gyoto::SmartPointer::operator();
%rename(assign) Gyoto::SmartPointer::operator=;
%include "GyotoSmartPointer.h"

%immutable Gyoto::Value::type;
%rename(assign) Gyoto::Value::operator=;
%rename(toDouble) Gyoto::Value::operator double;
%rename(toLong) Gyoto::Value::operator long;
%rename(toULong) Gyoto::Value::operator unsigned long;
%rename(toString) Gyoto::Value::operator std::string;
%rename(toVDouble) Gyoto::Value::operator std::vector<double>;
%{
  typedef unsigned long unsignedlong;
  %}
%rename(toVULong) Gyoto::Value::operator std::vector<unsigned long>;
%rename(toMetric) Gyoto::Value::operator Gyoto::SmartPointer<Gyoto::Metric::Generic>;
%rename(toAstrobj) Gyoto::Value::operator Gyoto::SmartPointer<Gyoto::Astrobj::Generic>;
%rename(toSpectrum) Gyoto::Value::operator Gyoto::SmartPointer<Gyoto::Spectrum::Generic>;
%rename(toSpectrometer) Gyoto::Value::operator Gyoto::SmartPointer<Gyoto::Spectrometer::Generic>;
%include "GyotoValue.h"
%include "GyotoObject.h"

%ignore Gyoto::Worldline::IntegState;
%ignore Gyoto::Worldline::IntegState::Generic;
%ignore Gyoto::Worldline::IntegState::Boost;
%ignore Gyoto::Worldline::IntegState::Legacy;
%include "GyotoWorldline.h"

GyotoSmPtrClass(Screen)
GyotoSmPtrClass(Scenery)
GyotoSmPtrClass(Photon)
%rename(AstrobjProperties) Gyoto::Astrobj::Properties;
%rename(increment) Gyoto::Astrobj::Properties::operator++;
GyotoSmPtrClassGeneric(Astrobj)

%ignore Gyoto::Astrobj::Standard;
%ignore Gyoto::Astrobj::UniformSphere;
%ignore Gyoto::Astrobj::ThinDisk;
%include GyotoStandardAstrobj.h
%include GyotoUniformSphere.h
%include GyotoThinDisk.h

%define _PAccessor2(member, setter)
  void setter(double *IN_ARRAY2, size_t DIM1, size_t DIM2) {
    $self->member = IN_ARRAY2;
  }
%enddef
%define _PAccessor3(member, setter)
void setter(double *IN_ARRAY3, size_t DIM1, size_t DIM2, size_t DIM3) {
    $self->member = IN_ARRAY3;
  }
%enddef
%extend Gyoto::Astrobj::Properties{
  _PAccessor2(intensity, Intensity)
  _PAccessor3(binspectrum, BinSpectrum)
  _PAccessor2(distance, MinDistance)
  _PAccessor2(first_dmin, FirstDistMin)
  _PAccessor3(impactcoords, ImpactCoords)
  _PAccessor2(redshift, Redshift)
  _PAccessor3(spectrum, Spectrum)
  _PAccessor2(time, EmissionTime)
  _PAccessor2(user1, User1)
  _PAccessor2(user2, User2)
  _PAccessor2(user3, User3)
  _PAccessor2(user4, User4)
  _PAccessor2(user5, User5)
 };

GyotoSmPtrClassGeneric(Metric)
GyotoSmPtrClassGeneric(Spectrum)
GyotoSmPtrClassGeneric(Spectrometer)

%extend Gyoto::Spectrometer::Complex {
  Gyoto::Spectrometer::Generic * __getitem__ (int i) {
    Gyoto::Spectrometer::Generic * res = ($self)->operator[](i);
    res -> incRefCount();
    return res;
  }
 };
%extend Gyoto::Spectrometer::Complex {
  void __setitem__(int i, Gyoto::Spectrometer::Generic * p) {
    ($self)->operator[](i)=p;
  }
 };
GyotoSmPtrClassDerivedPtrHdr(Spectrometer, Complex, ComplexSpectrometer, GyotoComplexSpectrometer.h)
GyotoSmPtrClassDerivedPtrHdr(Spectrometer, Uniform, UniformSpectrometer, GyotoUniformSpectrometer.h)

%include "GyotoConfig.h"
%include "GyotoDefs.h"
%include "GyotoUtils.h"
%include "GyotoFactory.h"

 // SWIG fails on nested classes. Work around this limitation:
%{
  typedef Gyoto::Screen::CoordType_e CoordType_e;
  typedef Gyoto::Screen::Coord1dSet Coord1dSet;
  typedef Gyoto::Screen::Coord2dSet Coord2dSet;
  typedef Gyoto::Screen::Grid Grid;
  typedef Gyoto::Screen::Bucket Bucket;
  typedef Gyoto::Screen::Empty Empty;
  typedef Gyoto::Screen::Range Range;
  typedef Gyoto::Screen::Indices Indices;
  typedef Gyoto::Screen::Angles Angles;
  typedef Gyoto::Screen::RepeatAngle RepeatAngle;
%}

enum CoordType_e;

%rename(increment) *::operator++;
class Coord1dSet {
public:
  const CoordType_e kind;
public:
  Coord1dSet(CoordType_e k);
  virtual void begin() =0;
  virtual bool valid() =0;
  virtual size_t size()=0;
  virtual size_t operator*() const ;
  virtual double angle() const ;
  virtual Coord1dSet& operator++()=0;
};

class Coord2dSet {
public:
  const CoordType_e kind;
  Coord2dSet(CoordType_e k);
  virtual Coord2dSet& operator++()    =0;
  virtual GYOTO_ARRAY<size_t, 2> operator*  () const;
  virtual GYOTO_ARRAY<double, 2> angles() const ;
  virtual void begin() =0;
  virtual bool valid() =0;
  virtual size_t size()=0;
};

class Grid: public Coord2dSet {
protected:
protected:
  const char * const prefix_;
  Coord1dSet &iset_;
  Coord1dSet &jset_;
public:
  Grid(Coord1dSet &iset, Coord1dSet &jset, const char * const p=NULL);
  virtual Coord2dSet& operator++();
  virtual GYOTO_ARRAY<size_t, 2> operator*  () const;
  virtual void begin();
  virtual bool valid();
  virtual size_t size();
};

class Bucket : public Coord2dSet {
protected:
  Coord1dSet &alpha_;
  Coord1dSet &delta_;
public:
  Bucket(Coord1dSet &iset, Coord1dSet &jset);
  virtual Coord2dSet& operator++();
  virtual GYOTO_ARRAY<double, 2> angles() const;
  virtual GYOTO_ARRAY<size_t, 2> operator*() const;
  virtual void begin();
  virtual bool valid();
  virtual size_t size();
};

class Empty: public Coord2dSet {
public:
  Empty();
  virtual Coord2dSet& operator++();
  virtual void begin();
  virtual bool valid();
  virtual size_t size();
};

class Range : public Coord1dSet {
protected:
  const size_t mi_, ma_, d_, sz_;
  size_t cur_;
public:
  Range(size_t mi, size_t ma, size_t d);
  void begin();
  bool valid();
  size_t size();
  Coord1dSet& operator++();
  size_t operator*() const ;
};

class Indices : public Coord1dSet {
protected:
  size_t const * const indices_;
  size_t const sz_;
  size_t i_;
public:
  Indices (size_t *IN_ARRAY1, size_t DIM1);
  void begin();
  bool valid();
  size_t size();
  Coord1dSet& operator++();
  size_t operator*() const ;
};

class Angles : public Coord1dSet {
protected:
  double const * const buf_;
  size_t const sz_;
  size_t i_;
public:
  Angles (double * IN_ARRAY1, size_t DIM1);
  void begin();
  bool valid();
  size_t size();
  Coord1dSet& operator++();
  double angle() const ;
};

class RepeatAngle : public Coord1dSet {
protected:
  double const val_;
  size_t const sz_;
  size_t i_;
public:
  RepeatAngle (double val, size_t sz);
  void begin();
  bool valid();
  size_t size();
  Coord1dSet& operator++();
  double angle() const ;
};

// Can't work out how to get swig to do this automatically
%inline {
enum {
  Property_double_t=Gyoto::Property::double_t,
  Property_long_t=Gyoto::Property::long_t,
  Property_unsigned_long_t=Gyoto::Property::unsigned_long_t,
  Property_bool_t=Gyoto::Property::bool_t,
  Property_string_t=Gyoto::Property::string_t,
  Property_filename_t=Gyoto::Property::filename_t,
  Property_vector_double_t=Gyoto::Property::vector_double_t,
  Property_vector_unsigned_long_t=Gyoto::Property::vector_unsigned_long_t,
  Property_metric_t=Gyoto::Property::metric_t,
  Property_screen_t=Gyoto::Property::screen_t,
  Property_astrobj_t=Gyoto::Property::astrobj_t,
  Property_spectrum_t=Gyoto::Property::spectrum_t,
  Property_spectrometer_t=Gyoto::Property::spectrometer_t,
  Property_empty_t=Gyoto::Property::empty_t};
}
