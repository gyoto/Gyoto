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
    Python bindings only, but it should not be too difficult to provide
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
  if (normal_pointer) normal_pointer->incRefCount();
  $result = SWIG_NewPointerObj( normal_pointer, stype, SWIG_POINTER_OWN |  0 );
}
%typemap(typecheck,precedence=SWIG_TYPECHECK_VOIDPTR)
Gyoto::SmartPointer<gtype>, gtype * {
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
%rename(klass ## Register) Gyoto::klass::Register_;
%inline {
  Gyoto::Register::Entry * get ## klass ## Register() { return Gyoto::klass::Register_; }
 }
%rename(register ## klass) Gyoto::klass::Register;
%rename(init ## klass ## Register) Gyoto::klass::initRegister;
%rename(get ## klass ## Subcontractor) Gyoto::klass::getSubcontractor;
%rename(klass) Gyoto::klass::Generic;
%feature("ref") Gyoto:: klass ::Generic"$this->incRefCount();";
%feature("unref") Gyoto:: klass ::Generic"$this->decRefCount(); if (!$this->getRefCount()) delete $this;";
// Need to mark the base classes as "notabstract" to extend them with
// a down-cast constructor
%feature("notabstract") Gyoto::klass::Generic;
// Ignore all the actual constructors are these classes are really abstract
%ignore  Gyoto::klass::Generic::Generic(Gyoto::klass::Generic const &);
%ignore  Gyoto::klass::Generic::Generic(const Generic &);
%ignore  Gyoto::klass::Generic::Generic(const klass::Generic &);
%ignore  Gyoto::klass::Generic::Generic();
%ignore  Gyoto::klass::Generic::Generic(double);
%ignore  Gyoto::klass::Generic::Generic(kind_t);
%ignore  Gyoto::klass::Generic::Generic(const std::string);
// Make a pseudo constructor for down-casting.
%extend Gyoto::klass::Generic {
  Generic(std::string nm) {
    Gyoto::SmartPointer<Gyoto::klass::Generic> pres=
      Gyoto::klass::getSubcontractor(nm.c_str())(NULL);
    Gyoto::klass::Generic * res = (Gyoto::klass::Generic *)(pres);
    if (res) res -> incRefCount();
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
#define GYOTO_NO_DEPRECATED
#include "GyotoConfig.h"
#include "GyotoDefs.h"
#include "GyotoFactory.h"
#include "GyotoFactoryMessenger.h"
#include "GyotoValue.h"
#include "GyotoProperty.h"
#include "GyotoObject.h"
#include "GyotoAstrobj.h"
#include "GyotoError.h"
#include "GyotoWorldline.h"
#include "GyotoPhoton.h"
#include "GyotoScreen.h"
#include "GyotoThinDisk.h"
#include "GyotoSpectrometer.h"
#include "GyotoComplexSpectrometer.h"
#include "GyotoUniformSpectrometer.h"
#include "GyotoRegister.h"
#include "GyotoWIP.h"
#include "GyotoConverters.h"
using namespace Gyoto;

%}

// ******** INITIALIZATION ******** //
// This will be called upon extension initialization
%init {
  Gyoto::Register::init();
#ifdef SWIGPYTHON
  import_array();
#endif
 }

// Rename operator++() -> increment() for everything
%rename(increment) *::operator++;

// ******** TYPEMAPS ******** //
// Actually instanciate typemaps using de macros defined above

GyotoSmPtrTypeMapClassGeneric(Metric);
GyotoSmPtrTypeMapClassGeneric(Astrobj);
GyotoSmPtrTypeMapClassGeneric(Spectrum);
GyotoSmPtrTypeMapClassGeneric(Spectrometer);
GyotoSmPtrTypeMapClassDerived(Astrobj, ThinDisk);

GyotoSmPtrTypeMapClass(Screen);
GyotoSmPtrTypeMapClass(Scenery);
GyotoSmPtrTypeMapClass(Photon);
GyotoSmPtrTypeMapClass(SmartPointee);

GyotoSmPtrTypeMapClassDerived(Spectrometer, Complex);
GyotoSmPtrTypeMapClassDerived(Spectrometer, Uniform);

GyotoSmPtrTypeMapClassDerived(Units, Unit);
GyotoSmPtrTypeMapClassDerived(Units, Converter);
GyotoSmPtrTypeMapClassDerived(Worldline, IntegState);
GyotoSmPtrTypeMapClassDerived(Astrobj, Properties);

// Typemaps for Gyoto::Value:

// In: cast from target language representations for all the supported
// types.
%typemap(in) Gyoto::Value {
  int res=0;
  void *argp=0;
  res=SWIG_ConvertPtr($input, &argp, SWIGTYPE_p_Gyoto__Value, 0);
  if (SWIG_IsOK(res)) {
    Gyoto::Value * temp = reinterpret_cast< Gyoto::Value * >(argp);
    $1 = *temp;
    if (SWIG_IsNewObj(res)) delete temp;
  }

  if (!SWIG_IsOK(res)) {
    res=SWIG_ConvertPtr($input, &argp, SWIGTYPE_p_Gyoto__Metric__Generic, 0);
    if (SWIG_IsOK(res)) {
      Gyoto::SmartPointer<Gyoto::Metric::Generic> temp = reinterpret_cast< Gyoto::Metric::Generic * >(argp);
      $1 = Gyoto::Value(temp);
    }
  }

  if (!SWIG_IsOK(res)) {
    res=SWIG_ConvertPtr($input, &argp, SWIGTYPE_p_Gyoto__Astrobj__Generic, 0);
    if (SWIG_IsOK(res)) {
      Gyoto::SmartPointer<Gyoto::Astrobj::Generic> temp = reinterpret_cast< Gyoto::Astrobj::Generic * >(argp);
      $1 = Gyoto::Value(temp);
    }
  }

  if (!SWIG_IsOK(res)) {
    res=SWIG_ConvertPtr($input, &argp, SWIGTYPE_p_Gyoto__Spectrum__Generic, 0);
    if (SWIG_IsOK(res)) {
      Gyoto::SmartPointer<Gyoto::Spectrum::Generic> temp = reinterpret_cast< Gyoto::Spectrum::Generic * >(argp);
      $1 = Gyoto::Value(temp);
    }
  }

  if (!SWIG_IsOK(res)) {
    res=SWIG_ConvertPtr($input, &argp, SWIGTYPE_p_Gyoto__Spectrometer__Generic, 0);
    if (SWIG_IsOK(res)) {
      Gyoto::SmartPointer<Gyoto::Spectrometer::Generic> temp = reinterpret_cast< Gyoto::Spectrometer::Generic * >(argp);
      $1 = Gyoto::Value(temp);
    }
  }

  if (!SWIG_IsOK(res)) {
    res=SWIG_ConvertPtr($input, &argp, SWIGTYPE_p_Gyoto__Screen, 0);
    if (SWIG_IsOK(res)) {
      Gyoto::SmartPointer<Gyoto::Screen> temp = reinterpret_cast< Gyoto::Screen * >(argp);
      $1 = Gyoto::Value(temp);
    }
  }

  if (!SWIG_IsOK(res)) {
    std::string *temp = 0;
    res = SWIG_AsPtr_std_string ($input, &temp) ;
    if (SWIG_IsOK(res)) $1 = Gyoto::Value(*temp);
    if (SWIG_IsNewObj(res) && temp) delete temp;
  }

  if (!SWIG_IsOK(res)) {
    std::vector<unsigned long> *temp=0;
    res = swig::traits_asptr< std::vector<unsigned long> >::asptr($input, &temp);
    if (SWIG_IsOK(res)) $1 = Gyoto::Value(*temp);
    if (SWIG_IsNewObj(res) && temp) delete temp;
  }

  if (!SWIG_IsOK(res)) {
    std::vector<double> *temp=0;
    res = swig::traits_asptr< std::vector<double> >::asptr($input, &temp);
    if (SWIG_IsOK(res)) $1 = Gyoto::Value(*temp);
    if (SWIG_IsNewObj(res) && temp) delete temp;
  }

  if (!SWIG_IsOK(res)) {
    long temp=0;
    res = SWIG_AsVal(long)($input, &temp);
    if (SWIG_IsOK(res)) $1 = Gyoto::Value(temp);
  }

  if (!SWIG_IsOK(res)) {
    double temp=0;
    res = SWIG_AsVal(double)($input, &temp);
    if (SWIG_IsOK(res)) $1 = Gyoto::Value(temp);
  }

  if (!SWIG_IsOK(res))
    SWIG_exception_fail(SWIG_ArgError(res), "argument of type 'Gyoto::Value*'");

 }
// Typecheck: should be debugged, does not seem to filter anything
%typemap(typecheck) Gyoto::Value {
  void *vptr = 0;
  int res = SWIG_ConvertPtr(argv[0], &vptr, SWIGTYPE_p_Gyoto__Value, 0);
  $1 = res;
}
// Out: cast from Gyoto::Value to language-specific representation
// for each sub-type.
%typemap(out) Gyoto::Value {
  switch ($1.type) {
  case Gyoto::Property::unsigned_long_t:
    $result = SWIG_From_unsigned_SS_long((unsigned long)($1));
    break;
  case Gyoto::Property::long_t:
    $result = SWIG_From_long(long($1));
    break;
  case Gyoto::Property::bool_t:
    $result = SWIG_From_bool(bool($1));
    break;
  case Gyoto::Property::double_t:
    $result = SWIG_From_double(double($1));
    break;
  case Gyoto::Property::filename_t:
  case Gyoto::Property::string_t:
    $result = SWIG_From_std_string(static_cast< std::string >($1));
    break;
  case Gyoto::Property::vector_double_t:
    $result = swig::from(($1).operator std::vector<double>());
    break;
  case Gyoto::Property::vector_unsigned_long_t:
    $result = swig::from(($1).operator std::vector<unsigned long>());
    break;
  case Gyoto::Property::metric_t:
    {
      Gyoto::Metric::Generic* normal_pointer=(Gyoto::Metric::Generic *) (Gyoto::SmartPointer<Gyoto::Metric::Generic>($1));
      if (normal_pointer) normal_pointer->incRefCount();
      $result = SWIG_NewPointerObj( normal_pointer, SWIGTYPE_p_Gyoto__Metric__Generic, SWIG_POINTER_OWN |  0 );
    }
    break;
  case Gyoto::Property::astrobj_t:
    {
      Gyoto::Astrobj::Generic* normal_pointer=(Gyoto::Astrobj::Generic *) (Gyoto::SmartPointer<Gyoto::Astrobj::Generic>($1));
      if (normal_pointer) normal_pointer->incRefCount();
      $result = SWIG_NewPointerObj( normal_pointer, SWIGTYPE_p_Gyoto__Astrobj__Generic, SWIG_POINTER_OWN |  0 );
    }
    break;
  case Gyoto::Property::spectrum_t:
    {
      Gyoto::Spectrum::Generic* normal_pointer=(Gyoto::Spectrum::Generic *) (Gyoto::SmartPointer<Gyoto::Spectrum::Generic>($1));
      if (normal_pointer) normal_pointer->incRefCount();
      $result = SWIG_NewPointerObj( normal_pointer, SWIGTYPE_p_Gyoto__Spectrum__Generic, SWIG_POINTER_OWN |  0 );
    }
    break;
  case Gyoto::Property::spectrometer_t:
    {
      Gyoto::Spectrometer::Generic* normal_pointer=(Gyoto::Spectrometer::Generic *) (Gyoto::SmartPointer<Gyoto::Spectrometer::Generic>($1));
      if (normal_pointer) normal_pointer->incRefCount();
      $result = SWIG_NewPointerObj( normal_pointer, SWIGTYPE_p_Gyoto__Spectrometer__Generic, SWIG_POINTER_OWN |  0 );
    }
    break;
  case Gyoto::Property::screen_t:
    {
      Gyoto::Screen* normal_pointer=(Gyoto::Screen *) (Gyoto::SmartPointer<Gyoto::Screen>($1));
      if (normal_pointer) normal_pointer->incRefCount();
      $result = SWIG_NewPointerObj( normal_pointer, SWIGTYPE_p_Gyoto__Screen, SWIG_POINTER_OWN |  0 );
    }
    break;
  default:
    $result = SWIG_NewPointerObj((new Gyoto::Value(static_cast< const Gyoto::Value& >($1))), SWIGTYPE_p_Gyoto__Value, SWIG_POINTER_OWN |  0 );
  }
}


// Non-Gyoto typemaps:
// Handle std::string
%include "std_string.i";

// Handle std::vector<double> and <unsigned long int>
%include "std_vector.i";
%template(vector_double) std::vector<double>;
%{
  typedef unsigned long unsignedlong;
  %}
%template(vector_unsigned_long) std::vector<unsigned long>;

#ifdef SWIGPYTHON
// Handle some arrays as NumPy arrays
%include "numpy.i";
%numpy_typemaps(size_t, NPY_ULONG , size_t);
%numpy_typemaps(double, NPY_DOUBLE, size_t);
#endif

// Handle generic C arrays using a class-like interface
%include "carrays.i"
%array_class(double, array_double);
%array_class(size_t, array_size_t);
%array_class(unsigned long, array_unsigned_long);
#ifdef SWIGPYTHON
// Provide conversion between generic C arrays and NumPy ndarrays
%define ExtendArrayNumPy(name, type)
%extend name {
  static name* fromnumpy1(type* IN_ARRAY1, size_t DIM1) {
    return static_cast< name * >(IN_ARRAY1);
  }
  static name* fromnumpy2(type* IN_ARRAY2, size_t DIM1, size_t DIM2) {
    return static_cast< name * >(IN_ARRAY2);
  }
  static name* fromnumpy3(type* IN_ARRAY3, size_t DIM1, size_t DIM2, size_t DIM3) {
    return static_cast< name * >(IN_ARRAY3);
  }
  static name* fromnumpy4(type* IN_ARRAY4, size_t DIM1, size_t DIM2, size_t DIM3, size_t DIM4) {
    return static_cast< name * >(IN_ARRAY4);
  }
};
%enddef

ExtendArrayNumPy(array_double, double);
ExtendArrayNumPy(array_unsigned_long, unsigned long);
ExtendArrayNumPy(array_size_t, size_t);
#endif

// ******** INTERFACE ******** //
// Here starts the actual parsing of the various header files

// Expose the build-time configuration variables
%include "GyotoConfig.h"

// Expose the global definitions and typedefs
%include "GyotoDefs.h"

// Expose the Gyoto::Error class
%include "GyotoError.h"

// Catch all Gyoto errors and re-throw them as run-time errors for the
// target language
%exception {
  try {
    $action
  }
  catch (Gyoto::Error e) {
    SWIG_Error(SWIG_RuntimeError, e);
    SWIG_fail;
  }
}

%ignore Gyoto::SmartPointer::operator();
%rename(assign) Gyoto::SmartPointer::operator=;
%include "GyotoSmartPointer.h"

// Expose Gyoto::Register::list as gyoto.listRegister
%rename(RegisterEntry) Gyoto::Register::Entry;
%rename(initRegister) Gyoto::Register::init;
%rename(listRegister) Gyoto::Register::list;
%include GyotoRegister.h

%rename(Functor__Double_constDoubleArray) Gyoto::Functor::Double_constDoubleArray;
%rename(Functor__Double_Double_const) Gyoto::Functor::Double_Double_const;
%include "GyotoFunctors.h"

%include "GyotoHooks.h"

%include "GyotoWIP.h"

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

%rename(Worldline__IntegState__Generic) Gyoto::Worldline::IntegState::Generic;
%rename(Worldline__IntegState__Boost) Gyoto::Worldline::IntegState::Boost;
%rename(Worldline__IntegState__Legacy) Gyoto::Worldline::IntegState::Legacy;
%include "GyotoWorldline.h"

GyotoSmPtrClass(Screen)
GyotoSmPtrClass(Scenery)
%ignore Gyoto::Photon::Refined;
GyotoSmPtrClass(Photon)
%rename(AstrobjProperties) Gyoto::Astrobj::Properties;
GyotoSmPtrClassGeneric(Astrobj)

GyotoSmPtrClassDerived(Astrobj, ThinDisk)

%define _PConverter(member, method)
  Gyoto::Units::Converter * method() {
  Gyoto::Units::Converter * res = $self->member;
    if (res) res -> incRefCount();
    return res;
  }
%enddef
%extend Gyoto::Astrobj::Properties{
  _PConverter(binspectrum_converter_, binSpectrumConverter)
  _PConverter(intensity_converter_, intensityConverter)
  _PConverter(spectrum_converter_, spectrumConverter)
 };

GyotoSmPtrClassGeneric(Metric)
GyotoSmPtrClassGeneric(Spectrum)
GyotoSmPtrClassGeneric(Spectrometer)

%extend Gyoto::Spectrometer::Complex {
  Gyoto::Spectrometer::Generic * __getitem__ (int i) {
    Gyoto::Spectrometer::Generic * res = ($self)->operator[](i);
    if (res) res -> incRefCount();
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
%include "GyotoUtils.h"
%include "GyotoFactory.h"
%include "GyotoFactoryMessenger.h"

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

#ifdef HAVE_BOOST_ARRAY_HPP
namespace boost {
template <typename T, size_t sz> class array {
 public:
    T& operator[](size_t c) { return buf[c] ; }
  };
}
#endif
%extend GYOTO_ARRAY {
  T __getitem__(size_t c) {
    return $self->operator[](c);
  }
};
%template(ARRAY_double_2) GYOTO_ARRAY<double, 2>;
%template(ARRAY_size_t_2) GYOTO_ARRAY<size_t, 2>;

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
  Indices (size_t *carray, size_t nel);
  void begin();
  bool valid();
  size_t size();
  Coord1dSet& operator++();
  size_t operator*() const ;
};
%extend Indices {
  Indices (size_t DIM1, size_t *IN_ARRAY1) {
    return new Indices(IN_ARRAY1, DIM1);
  }

}

class Angles : public Coord1dSet {
protected:
  double const * const buf_;
  size_t const sz_;
  size_t i_;
public:
  Angles (double * carray, size_t nel);
  void begin();
  bool valid();
  size_t size();
  Coord1dSet& operator++();
  double angle() const ;
};
%extend Angles {
  Angles (size_t DIM1, double *IN_ARRAY1) {
    return new Angles(IN_ARRAY1, DIM1);
  }

}

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

%ignore Gyoto::Property::Property;
%include "GyotoProperty.h"

%include "GyotoConverters.h"

// Workaround cvar bug in Swig which makes help(gyoto) fail:
%inline {
  namespace Gyoto {
    extern int __class__=0;
  }
}
