#include "GyotoProperty.h"

using namespace std ;
using namespace Gyoto ;

Property::Property(Property const * const ancestors)
  : name(""), type(empty_t), parent(ancestors) {}

#define GYOTO_LOCAL(T)							\
  Property::Property(string n, set_##T##_t set, get_##T##_t get)	\
    : name(n), type(T##_t), parent(NULL) {				\
    setter.set_##T=set;							\
    getter.get_##T=get;							\
  }

GYOTO_LOCAL(long)
GYOTO_LOCAL(unsigned_long)
#if !defined(GYOTO_SIZE__T_IS_UNSIGNED_LONG)
GYOTO_LOCAL(size_t)
#endif
GYOTO_LOCAL(metric)
GYOTO_LOCAL(spectrum)
GYOTO_LOCAL(astrobj)
GYOTO_LOCAL(screen)
GYOTO_LOCAL(spectrometer)

#undef GYOTO_LOCAL

Property::Property(string n, set_double_t set, get_double_t get)
  : name(n), type(double_t), parent(NULL) {
  setter.set_double=set;
  getter.get_double=get;
  setter_unit.set_double=NULL;
  getter_unit.get_double=NULL;
}

Property::Property(string n, set_double_t set, get_double_t get,
		   set_double_unit_t setu, get_double_unit_t getu)
  : name(n), type(double_t), parent(NULL) {
  setter.set_double=set;
  getter.get_double=get;
  setter_unit.set_double=setu;
  getter_unit.get_double=getu;
}

Property::Property(string n, string nf, set_bool_t set, get_bool_t get)
  : name(n), name_false(nf), type(bool_t), parent(NULL) {
  setter.set_bool=set;
  getter.get_bool=get;
}

Property::Property(string n, set_string_t set, get_string_t get,
		   bool is_filename)
  : name(n), type(is_filename?filename_t:string_t), parent(NULL) {
  setter.set_string=set;
  getter.get_string=get;
}

Property::Property(string n,
		   set_vector_double_t set,
		   get_vector_double_t get)
  : name(n), type(vector_double_t), parent(NULL) {
  setter.set_vdouble=set;
  getter.get_vdouble=get;
  setter_unit.set_vdouble=NULL;
  getter_unit.get_vdouble=NULL;
}

Property::Property(string n,
		   set_vector_double_t set,
		   get_vector_double_t get,
		   set_vector_double_unit_t setu,
		   get_vector_double_unit_t getu)
  : name(n), type(vector_double_t), parent(NULL) {
  setter.set_vdouble=set;
  getter.get_vdouble=get;
  setter_unit.set_vdouble=setu;
  getter_unit.get_vdouble=getu;
}

Property::Property(string n,
		   set_vector_unsigned_long_t set,
		   get_vector_unsigned_long_t get)
  : name(n), type(vector_unsigned_long_t), parent(NULL) {
  setter.set_vulong=set;
  getter.get_vulong=get;
}

Property::operator bool() const { return type != empty_t; }
