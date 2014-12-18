%module(docstring="The General relativitY Orbit Tracer of paris Observatory") gyoto
%feature("autodoc", "1");

%define GYOTO_SWIGIMPORTED
%enddef

%define GyotoSmPtrClass(klass)
%ignore Gyoto:: klass ;
%include Gyoto ## klass ## .h
%rename(klass) pyGyoto ## klass;
%inline {
Gyoto::SmartPointer<Gyoto::klass> pyGyoto ## klass () {
  return new Gyoto::klass();
}
}
%template(klass ## Ptr) Gyoto::SmartPointer<Gyoto::klass>;
%enddef

%define GyotoSmPtrClassGeneric(klass)
%template(klass ## Ptr) Gyoto::SmartPointer<Gyoto::klass::Generic>;
%rename(klass) klass ## Ptr;
%ignore Gyoto::klass::Register_;
%ignore Gyoto::klass::Register;
%ignore Gyoto::klass::initRegister;
%ignore Gyoto::klass::getSubcontractor;
%ignore Gyoto::klass::Generic;
%include Gyoto ## klass ## .h
%rename(klass) pyGyoto ## klass;
%inline {
Gyoto::SmartPointer<Gyoto::klass::Generic> pyGyoto ## klass (std::string const&s) {
  return Gyoto::klass::getSubcontractor(s.c_str())(NULL);
}
}
%enddef

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
using namespace Gyoto;

%}

%include "GyotoError.h"

%exception {
	try {
	$function
	}
	catch (Gyoto::Error e) {
	 	PyErr_SetString(PyExc_IndexError, e);
		return NULL;
	}
}


%include "std_string.i" 
%include "std_vector.i"
%template(vector_double) std::vector<double>;
%template(vector_unsigned_long) std::vector<unsigned long>;
%include "carrays.i"
%array_class(double, array_double)
%include "numpy.i"

%init {
  Gyoto::Register::init();
  import_array();
 }

%include "GyotoSmartPointer.h"

%include "GyotoValue.h"
%include "GyotoObject.h"

%ignore Gyoto::Worldline::IntegState;
%ignore Gyoto::Worldline::IntegState::Generic;
%ignore Gyoto::Worldline::IntegState::Boost;
%ignore Gyoto::Worldline::IntegState::Legacy;
%include "GyotoWorldline.h"

%rename (castToWorldline) pyGyotoCastToWorldline;
%inline %{
Gyoto::Worldline * pyGyotoCastToWorldline
  (Gyoto::SmartPointer<Gyoto::Astrobj::Generic> const_p) {
  Gyoto::SmartPointee * p=const_cast<Gyoto::Astrobj::Generic*>(const_p());
  Gyoto::Worldline * res = dynamic_cast<Gyoto::Worldline*>(p);
  return res;
}
%}

GyotoSmPtrClass(Screen)
GyotoSmPtrClass(Scenery)
GyotoSmPtrClass(Photon)
GyotoSmPtrClassGeneric(Astrobj)


%define _PAccessor2(member, setter)
  void setter(double *IN_ARRAY2, int DIM1, int DIM2) {
    $self->member = IN_ARRAY2;
  }
%enddef
%define _PAccessor3(member, setter)
void setter(double *IN_ARRAY3, int DIM1, int DIM2, int DIM3) {
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
  Indices (unsigned long * IN_ARRAY1, unsigned long DIM1);
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
  Angles (double const*const buf, size_t sz);
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
