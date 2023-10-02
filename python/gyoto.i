/*
    Copyright 2014-2022 Thibaut Paumard

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
%module(docstring="The General relativitY Orbit Tracer of paris Observatory", package="gyoto") core

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
%feature("ref") Gyoto:: klass "$this->incRefCount();//ref";
%feature("unref") Gyoto:: klass "$this->decRefCount(); if (!$this->getRefCount()) delete $this;//unref";
// clone() method returns a new object with refCount 0
%newobject Gyoto :: klass :: clone;
%extend Gyoto::klass {
  std::string __str__() {
    return Gyoto::Factory($self).format();
  }
};
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
%feature("ref") Gyoto:: klass ::Generic"$this->incRefCount();//ref Generic";
%feature("unref") Gyoto:: klass ::Generic"$this->decRefCount(); if (!$this->getRefCount()) delete $this;//unref Generic";
%newobject Gyoto :: klass :: Generic :: clone;
// Need to mark the base classes as "notabstract" to extend them with
// a down-cast constructor
%feature("notabstract") Gyoto::klass::Generic;
// Ignore all the actual constructors as these classes are really abstract
%ignore  Gyoto::klass::Generic::Generic(Gyoto::klass::Generic const &);
%ignore  Gyoto::klass::Generic::Generic(const Generic &);
%ignore  Gyoto::klass::Generic::Generic(const klass::Generic &);
%ignore  Gyoto::klass::Generic::Generic();
%ignore  Gyoto::klass::Generic::Generic(double);
%ignore  Gyoto::klass::Generic::Generic(kind_t);
%ignore  Gyoto::klass::Generic::Generic(const std::string);
%ignore  Gyoto::klass::Generic::Generic(const int, const std::string &);
// Make a pseudo constructor for down-casting.
%extend Gyoto::klass::Generic {
  Generic(std::string nm) {
    std::vector<std::string> plugin;
    Gyoto::klass::Generic * res = NULL;
    {
      Gyoto::SmartPointer<Gyoto::klass::Generic> pres=
        Gyoto::klass::getSubcontractor(nm.c_str(), plugin)(NULL, plugin);
      res = (Gyoto::klass::Generic *)(pres);
      // We need to increment refcount, else the object is destroyed
      // when the original smartpoiter is:
      if (res) res -> incRefCount();
    }
    // Now that the original smartpointer has been detroyed, refcount is 1.
    // ref feature will increment it again, so we need to decrement is now:
    res->decRefCount();
    GYOTO_DEBUG_EXPR(res->getRefCount());
    // Special for Uniform spectrometer:
    // if 'res' can be cast to uniform spectrometer, set Kind.
    if(dynamic_cast<Gyoto::Spectrometer::Uniform*>(res))
      res->set("Kind", nm);
    // end special case
    return res;
  }
  Generic(std::string nm, std::vector<std::string> plugin) {
    GYOTO_DEBUG_EXPR(plugin.size());
    Gyoto::klass::Generic * res = NULL;
    {
      Gyoto::SmartPointer<Gyoto::klass::Generic> pres=
        Gyoto::klass::getSubcontractor(nm.c_str(), plugin)(NULL, plugin);
      res = (Gyoto::klass::Generic *)(pres);
      // We need to increment refcount, else the object is destroyed
      // when the original smartpoiter is:
      if (res) res -> incRefCount();
    }
    // Now that the original smartpointer has been detroyed, refcount is 1.
    // ref feature will increment it again, so we need to decrement is now:
    res->decRefCount();
    GYOTO_DEBUG_EXPR(res->getRefCount());
    // Special for Uniform spectrometer:
    // if 'res' can be cast to uniform spectrometer, set Kind.
    if(dynamic_cast<Gyoto::Spectrometer::Uniform*>(res))
      res->set("Kind", nm);
    // end special case
    return res;
  }
  Generic(long address) {
    Gyoto::klass::Generic * res = (Gyoto::klass::Generic *)(address);
    // Should be done by ref feature:
    // if (res) res -> incRefCount();
    return res;
  }
  Generic(Gyoto::klass::Generic *orig) {
    return orig;
  }
  std::string __str__() {
    return Gyoto::Factory($self).format();
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
%newobject Gyoto :: nspace:: klass :: clone;
%feature("notabstract") Gyoto::nspace::klass;
%extend  Gyoto::nspace::klass {
  klass(Gyoto::nspace::Generic * base) {
    Gyoto::nspace::klass * res = dynamic_cast< Gyoto::nspace::klass * >(base);
    if (!res) GYOTO_ERROR("This pointer cannot be cast to 'Gyoto::" #nspace "::" #klass "*'");
    return res;
  }
  klass(long address) {
    Gyoto::nspace::klass * res = (Gyoto::nspace::klass *)(address);
    // Should be done by ref feature:
    // if (res) res -> incRefCount();
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

// Add things for all metrics
%define GyotoSmPtrClassDerivedMetric(klass)
%extend Gyoto::Metric::klass {
  // Support this syntax:
  // vel = gg.circularVelocity(pos)
  // in addition of gg.circularVelocity(pos, vel)
  // Same for gmunu and christoffel
  void circularVelocity(double const IN_ARRAY1[4], double ARGOUT_ARRAY1[4]) {
    Gyoto::SmartPointer<Gyoto::Metric::Generic>($self)->circularVelocity(IN_ARRAY1, ARGOUT_ARRAY1);
  }
  void zamoVelocity(double const IN_ARRAY1[4], double ARGOUT_ARRAY1[4]) {
    Gyoto::SmartPointer<Gyoto::Metric::Generic>($self)->zamoVelocity(IN_ARRAY1, ARGOUT_ARRAY1);
  }
  using Gyoto::Metric::Generic::christoffel;
  void christoffel(double ARGOUT_ARRAY3[4][4][4], double const IN_ARRAY1[4]) {
    Gyoto::SmartPointer<Gyoto::Metric::Generic>($self)->christoffel(ARGOUT_ARRAY3, IN_ARRAY1);
  }
};
  GyotoSmPtrClassDerived(Metric, klass)
%enddef

// ******** INCLUDES ******** //
// Include any file that is needed to compile the wrappers
%{
#define SWIG_FILE_WITH_INIT
#define GYOTO_NO_DEPRECATED
#include "gyoto_swig.h"
using namespace Gyoto;

swig_type_info * __Gyoto_SWIGTYPE_p_Gyoto__Error() {
  return SWIGTYPE_p_Gyoto__Error;
}

%}

// ******** INITIALIZATION ******** //

// Catch all Gyoto errors and re-throw them as run-time errors for the
// target language
%exception {
  try {
    $action
  }
  catch (Gyoto::Error e) {
    //    PyErr_SetString(PyErr_NewException("gyoto.Error", NULL, NULL), e);
   SWIG_Python_Raise
     (SWIG_NewPointerObj
      ((new Gyoto::Error(static_cast<const Gyoto::Error& >(e))),
      __Gyoto_SWIGTYPE_p_Gyoto__Error(),SWIG_POINTER_NEW),
      "gyoto.Error", __Gyoto_SWIGTYPE_p_Gyoto__Error());
    SWIG_fail;
  }
}

// This will be called upon extension initialization
%init {
  if (!Gyoto::Astrobj::Register_ &&
      !Gyoto::Metric::Register_ &&
      !Gyoto::Spectrum::Register_ &&
      !Gyoto::Spectrometer::Register_)
    Gyoto::Register::init();
#ifdef SWIGPYTHON
  import_array();
#endif
 }

// Rename operator++() -> increment() for everything
%rename(increment) *::operator++;
// Rename operator=() -> assign() for everything
%rename(assign) *::operator=;
// Rename operator*() -> __ref__
// nothing to do, that's the default

// ******** TYPEMAPS ******** //
// Actually instantiate typemaps using de macros defined above

GyotoSmPtrTypeMapClassGeneric(Metric);
GyotoSmPtrTypeMapClassGeneric(Astrobj);
GyotoSmPtrTypeMapClassGeneric(Spectrum);
GyotoSmPtrTypeMapClassGeneric(Spectrometer);
GyotoSmPtrTypeMapClassDerived(Astrobj, ThinDisk);
GyotoSmPtrTypeMapClassDerived(Astrobj, Standard)

GyotoSmPtrTypeMapClass(Screen);
GyotoSmPtrTypeMapClass(Scenery);
GyotoSmPtrTypeMapClass(Photon);
%extend Gyoto::SmartPointee {
  long getPointer() {
    return long($self);
  }
 };
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
  int res=-1;
  void *argp=0;

  if (!SWIG_IsOK(res)) {
    res=SWIG_ConvertPtr($input, &argp, SWIGTYPE_p_Gyoto__Metric__Generic, 0);
    if (SWIG_IsOK(res)) {
      Gyoto::SmartPointer<Gyoto::Metric::Generic> temp = reinterpret_cast< Gyoto::Metric::Generic * >(argp);
      if(temp) $1 = Gyoto::Value(temp);
      else $1 = Gyoto::Value();
    }
  }

  if (!SWIG_IsOK(res)) {
    res=SWIG_ConvertPtr($input, &argp, SWIGTYPE_p_Gyoto__Value, 0);
    if (SWIG_IsOK(res)) {
      Gyoto::Value * temp = reinterpret_cast< Gyoto::Value * >(argp);
      $1 = *temp;
      if (SWIG_IsNewObj(res)) delete temp;
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
%template(vector_string) std::vector<std::string>;
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

%apply (double * IN_ARRAY1, size_t DIM1) {(double * dates, size_t n_dates)};
%apply (double * IN_ARRAY1, size_t DIM1) {(double const * cst, size_t const ncsts)};
%apply (double * INPLACE_ARRAY1, size_t DIM1) {(double * x1dest, size_t n1)};
%apply (double * INPLACE_ARRAY1, size_t DIM1) {(double * x2dest, size_t n2)};
%apply (double * INPLACE_ARRAY1, size_t DIM1) {(double * x3dest, size_t n3)};
%apply (double * INPLACE_ARRAY1, size_t DIM1) {(double * x0dot, size_t n0d)};
%apply (double * INPLACE_ARRAY1, size_t DIM1) {(double * x1dot, size_t n1d)};
%apply (double * INPLACE_ARRAY1, size_t DIM1) {(double * x2dot, size_t n2d)};
%apply (double * INPLACE_ARRAY1, size_t DIM1) {(double * x3dot, size_t n3d)};
// Handle all const arrays of fixed size as NumPy IN_ARRAYs
%apply (double IN_ARRAY1[ANY]) {(const double [ANY])};
%apply (double IN_ARRAY2[ANY][ANY]) {(const double [ANY][ANY])};
%apply (double IN_ARRAY3[ANY][ANY][ANY]) {(const double [ANY][ANY][ANY])};
%apply (double IN_ARRAY4[ANY][ANY][ANY][ANY]) {(const double [ANY][ANY][ANY][ANY])};
// Handle all non-const arrays of fixed size as INPLACE.
%apply (double INPLACE_ARRAY1[ANY]) {(double [ANY])};
%apply (double INPLACE_ARRAY2[ANY][ANY]) {(double [ANY][ANY])};
%apply (double INPLACE_ARRAY3[ANY][ANY][ANY]) {(double [ANY][ANY][ANY])};
%apply (double INPLACE_ARRAY4[ANY][ANY][ANY][ANY]) {(double [ANY][ANY][ANY][ANY])};
%apply (double ARGOUT_ARRAY1[ANY]) {(double ARGOUT_ARRAY1_1[ANY])}
%apply (double ARGOUT_ARRAY1[ANY]) {(double ARGOUT_ARRAY1_2[ANY])}
%apply (double ARGOUT_ARRAY1[ANY]) {(double ARGOUT_ARRAY1_3[ANY])}
%apply (double ARGOUT_ARRAY1[ANY]) {(double ARGOUT_ARRAY1_4[ANY])}


// ******** INTERFACE ******** //
// Here starts the actual parsing of the various header files

// Expose the build-time configuration variables
%include "GyotoConfig.h"

// Expose the global definitions and typedefs
%include "GyotoDefs.h"

// Expose the Gyoto::Error class
// Not a SmartPointee
%extend Gyoto::Error {
  const char *__str__() {
    return *($self);
  }
 };
%exceptionclass Gyoto::Error ;
%include "GyotoError.h"

// Expose the SmartPointer API
%ignore Gyoto::SmartPointer::operator();
%rename(assign) Gyoto::SmartPointer::operator=;
%feature("ref") Gyoto::SmartPointee "$this->incRefCount();//ref SmartPointee";
%feature("unref") Gyoto::SmartPointee "$this->decRefCount(); if (!$this->getRefCount()) delete $this;//unref SmartPointee";
%include "GyotoSmartPointer.h"

// Expose Gyoto::Register::list as gyoto.listRegister
%rename(RegisterEntry) Gyoto::Register::Entry;
%rename(initRegister) Gyoto::Register::init;
%rename(listRegister) Gyoto::Register::list;
%include GyotoRegister.h


// Not a SmartPointee
%rename(Functor__Double_constDoubleArray) Gyoto::Functor::Double_constDoubleArray;
%rename(Functor__Double_Double_const) Gyoto::Functor::Double_Double_const;
%include "GyotoFunctors.h"

// Not a SmartPointee
%include "GyotoHooks.h"

// Not a SmartPointee
%include "GyotoWIP.h"

// Worldline: not a SmartPointee
%immutable Gyoto::Value::type;
%rename(assign) Gyoto::Value::operator=;
%rename(toDouble) Gyoto::Value::operator double;
%rename(toLong) Gyoto::Value::operator long;
%rename(toULong) Gyoto::Value::operator unsigned long;
%extend Gyoto::Value {size_t toSizeT() {return size_t(*$self);}};
%rename(toString) Gyoto::Value::operator std::string;
%rename(toVDouble) Gyoto::Value::operator std::vector<double>;
%extend Gyoto::Value {
  std::vector<unsigned long> toVULong() {return ($self)->operator std::vector<unsigned long>();}
};
%rename(toMetric) Gyoto::Value::operator Gyoto::SmartPointer<Gyoto::Metric::Generic>;
%rename(toAstrobj) Gyoto::Value::operator Gyoto::SmartPointer<Gyoto::Astrobj::Generic>;
%rename(toSpectrum) Gyoto::Value::operator Gyoto::SmartPointer<Gyoto::Spectrum::Generic>;
%rename(toSpectrometer) Gyoto::Value::operator Gyoto::SmartPointer<Gyoto::Spectrometer::Generic>;
%rename(toScreen) Gyoto::Value::operator Gyoto::SmartPointer<Gyoto::Screen>;
%include "GyotoValue.h"
%include "GyotoObject.h"

%rename(Worldline__IntegState__Generic) Gyoto::Worldline::IntegState::Generic;
%rename(Worldline__IntegState__Boost) Gyoto::Worldline::IntegState::Boost;
%rename(Worldline__IntegState__Legacy) Gyoto::Worldline::IntegState::Legacy;
%extend Gyoto::Worldline {
  void get_t(double * INPLACE_ARRAY1, size_t DIM1) {
    if (DIM1 != ($self)->get_nelements()) GYOTO_ERROR("wrong output array size");
    ($self)->get_t(INPLACE_ARRAY1);
  }
  void get_xyz(
		double * x1dest, size_t n1,
		double * x2dest, size_t n2,
                double * x3dest, size_t n3) {
    if (n1 != ($self)->get_nelements() || n2 != n1 || n3 != n1)
      GYOTO_ERROR("wrong size for output array");
    ($self)->get_xyz(x1dest, x2dest, x3dest);
  }
  void getSkyPos(SmartPointer<Screen> screen,
		double * x1dest, size_t n1,
		double * x2dest, size_t n2,
                double * x3dest, size_t n3) {
    if (n1 != ($self)->get_nelements() || n2 != n1 || n3 != n1)
      GYOTO_ERROR("wrong size for output array");
    ($self)->getSkyPos(screen, x1dest, x2dest, x3dest);
  }
  void get_dot(double * x0dot, size_t n0d,
               double * x1dot, size_t n1d,
               double * x2dot, size_t n2d,
               double * x3dot, size_t n3d) {
    if (n0d != ($self)->get_nelements() || n1d != n0d || n2d != n0d || n3d != n0d)
      GYOTO_ERROR("wrong size for output array");
    ($self)->get_dot(x0dot, x1dot, x2dot, x3dot);
  }
  void get_prime(
                 double * x1dot, size_t n1d,
                 double * x2dot, size_t n2d,
                 double * x3dot, size_t n3d) {
    if (n1d != ($self)->get_nelements() || n2d != n1d || n3d != n1d)
      GYOTO_ERROR("wrong size for output array");
    ($self)->get_prime(x1dot, x2dot, x3dot);
  }
  // support this syntax:
  // vel = gg.circularVelocity(pos)
  // in addition of gg.circularVelocity(pos, vel)
  // void getInitialCoord(std::vector<double> &ARGOUT_ARRAY1) {
  //   ($self)->getInitialCoord(ARGOUT_ARRAY1);
  // }
  // void getCoord(size_t index, Gyoto::Worldline::state_type &ARGOUT_ARRAY1) {
  //   ($self)->getCoord(index, ARGOUT_ARRAY1);
  // }
  void getCartesianPos(size_t index, double ARGOUT_ARRAY1[8]) {
    ($self)->getCartesianPos(index, ARGOUT_ARRAY1);
  }
};
%include "GyotoWorldline.h"

%extend Gyoto::Screen {
  // Support this syntax:
  // pos = scr.getObserverPos()
  void getObserverPos(double ARGOUT_ARRAY1[4]) {
    ($self)->getObserverPos(ARGOUT_ARRAY1);
  }
  void getFourVel(double ARGOUT_ARRAY1[4]) {
    ($self)->getFourVel(ARGOUT_ARRAY1);
  }
 };
GyotoSmPtrClass(Screen)
GyotoSmPtrClass(Scenery)
%ignore Gyoto::Photon::Refined;
GyotoSmPtrClass(Photon)
%rename(AstrobjProperties) Gyoto::Astrobj::Properties;
GyotoSmPtrClassGeneric(Astrobj)

GyotoSmPtrClassDerived(Astrobj, ThinDisk)

%ignore Gyoto::Astrobj::Standard::Standard();
%ignore Gyoto::Astrobj::Standard::Standard(double radmax);
%ignore Gyoto::Astrobj::Standard::Standard(std::string kind);
%ignore Gyoto::Astrobj::Standard::Standard(const Standard& );
GyotoSmPtrClassDerivedPtrHdr(Astrobj, Standard, StandardAstrobj, GyotoStandardAstrobj.h)

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

%extend Gyoto::Metric::Generic {
  // Support this syntax:
  // vel = gg.circularVelocity(pos)
  // in addition of gg.circularVelocity(pos, vel)
  // Same for gmunu and christoffel
  void circularVelocity(double const IN_ARRAY1[4], double ARGOUT_ARRAY1[4]) {
    ($self)->circularVelocity(IN_ARRAY1, ARGOUT_ARRAY1);
  }
  void zamoVelocity(double const IN_ARRAY1[4], double ARGOUT_ARRAY1[4]) {
    ($self)->zamoVelocity(IN_ARRAY1, ARGOUT_ARRAY1);
  }
  void christoffel(double ARGOUT_ARRAY3[4][4][4], double const IN_ARRAY1[4]) {
    ($self)->christoffel(ARGOUT_ARRAY3, IN_ARRAY1);
  }
};
GyotoSmPtrClassGeneric(Metric)
GyotoSmPtrClassGeneric(Spectrum)
GyotoSmPtrClassGeneric(Spectrometer)


%inline {
  class myCplxSpectroIdxExcept {};
}

%exception Gyoto::Spectrometer::Complex::__getitem__ {
  try {
    $action ;
  } catch (myCplxSpectroIdxExcept e) {
    SWIG_exception_fail(SWIG_IndexError, "Index out of bounds");
  }
}

%extend Gyoto::Spectrometer::Complex {
  Gyoto::SmartPointer<Gyoto::Spectrometer::Generic> __getitem__ (size_t i) {
    if (i >= ($self)->getCardinal()) {
      throw myCplxSpectroIdxExcept();
    }
    Gyoto::SmartPointer<Gyoto::Spectrometer::Generic> res = ($self)->operator[](i);
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

// Not a class
%include "GyotoConfig.h"

// Not a class
%include "GyotoUtils.h"

// Not a SmartPointee
%include "GyotoFactory.h"

// Backwards-compatibility code introduced 2018-10-04
%extend Gyoto::Factory {
  Gyoto::SmartPointer<Gyoto::Scenery> getScenery() {
    GYOTO_WARNING << "Factory::getScenery() is deprecated, use Factory::scenery() instead\n";
    return ($self)->scenery();
  }
  Gyoto::SmartPointer<Gyoto::Photon> getPhoton() {
    GYOTO_WARNING << "Factory::getPhoton() is deprecated, use Factory::photon() instead\n";
    return ($self)->photon();
  }
 };

// Not a SmartPointee
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
%extend Coord1dSet {
  // Get value of coord set
  size_t value() const {
    return **($self);
  }
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
  size_t index() const ;
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

// Not a SmartPointee
%ignore Gyoto::Property::Property;
%include "GyotoProperty.h"

// Units and Converters are SmartPointee
%extend Gyoto::Units::Unit {
  std::string __str__() {
    return (std::string(*($self)));
  }
};
%feature("ref") Gyoto::Units::Unit "$this->incRefCount();//ref Unit";
%feature("unref") Gyoto::Units::Unit "$this->decRefCount(); if (!$this->getRefCount()) delete $this;//unref Unit";
%feature("ref") Gyoto::Units::Converter "$this->incRefCount();//ref Converter";
%feature("unref") Gyoto::Units::Converter "$this->decRefCount(); if (!$this->getRefCount()) delete $this;//unref Converter";
%include "GyotoConverters.h"

// not a SmartPointee
%include "GyotoGridData2D.h"
%include "GyotoFitsRW.h"


// Workaround cvar bug in Swig which makes help(gyoto) fail:
%inline {
  namespace Gyoto {
    extern int __class__;
  }
  int Gyoto::__class__ = 0;
}
