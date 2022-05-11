#include "GyotoProperty.h"
#include "GyotoError.h"

using namespace std ;
using namespace Gyoto ;

Property::Property()
  : name(""), type(empty_t), parent(NULL) {}

Property::Property(std::string n, int t)
  : name(n), type(t), parent(NULL) {}

Property::Property(Property const * const ancestors)
  : name(""), type(empty_t), parent(ancestors) {}

Property::Property(std::string classname, std::string doc)
  : name(classname), type(empty_t), parent(NULL), doc(doc){}

#define GYOTO_LOCAL(T)							\
  Property::Property(string n, set_##T##_t set, get_##T##_t get, string d) \
    : name(n), type(T##_t), parent(NULL), doc(d) {			\
    setter.set_##T=set;							\
    getter.get_##T=get;							\
  }

GYOTO_LOCAL(long)
GYOTO_LOCAL(unsigned_long)
GYOTO_LOCAL(metric)
GYOTO_LOCAL(spectrum)
GYOTO_LOCAL(astrobj)
GYOTO_LOCAL(screen)
GYOTO_LOCAL(spectrometer)

#undef GYOTO_LOCAL

Property::Property(string n, set_size_t_t set, get_size_t_t get, int,
		   string d)
: name(n), type(size_t_t), parent(NULL), doc(d) {
  setter.set_size_t=set;
  getter.get_size_t=get;
}

Property::Property(string n, set_double_t set, get_double_t get, string d)
: name(n), type(double_t), parent(NULL), doc(d) {
  setter.set_double=set;
  getter.get_double=get;
  setter_unit.set_double=NULL;
  getter_unit.get_double=NULL;
}

Property::Property(string n, set_double_t set, get_double_t get,
		   set_double_unit_t setu, get_double_unit_t getu, string d)
  : name(n), type(double_t), parent(NULL), doc(d) {
  setter.set_double=set;
  getter.get_double=get;
  setter_unit.set_double=setu;
  getter_unit.get_double=getu;
}

Property::Property(string n, string nf, set_bool_t set, get_bool_t get,
		   string d)
  : name(n), name_false(nf), type(bool_t), parent(NULL), doc(d) {
  setter.set_bool=set;
  getter.get_bool=get;
}

Property::Property(string n, set_string_t set, get_string_t get,
		   bool is_filename, string d)
  : name(n), type(is_filename?filename_t:string_t), parent(NULL), doc(d) {
  setter.set_string=set;
  getter.get_string=get;
}

Property::Property(string n,
		   set_vector_double_t set,
		   get_vector_double_t get,
		   string d)
  : name(n), type(vector_double_t), parent(NULL), doc(d) {
  setter.set_vdouble=set;
  getter.get_vdouble=get;
  setter_unit.set_vdouble=NULL;
  getter_unit.get_vdouble=NULL;
}

Property::Property(string n,
		   set_vector_double_t set,
		   get_vector_double_t get,
		   set_vector_double_unit_t setu,
		   get_vector_double_unit_t getu,
		   string d)
  : name(n), type(vector_double_t), parent(NULL), doc(d) {
  setter.set_vdouble=set;
  getter.get_vdouble=get;
  setter_unit.set_vdouble=setu;
  getter_unit.get_vdouble=getu;
}

Property::Property(string n,
		   set_vector_unsigned_long_t set,
		   get_vector_unsigned_long_t get,
		   string d)
  : name(n), type(vector_unsigned_long_t), parent(NULL), doc(d) {
  setter.set_vulong=set;
  getter.get_vulong=get;
}

Property::operator bool() const { return type != empty_t || name != ""; }

Property::type_e Property::typeFromString(std::string stype) {
  if (stype=="double") {
    return Property::double_t;
  } else if (stype=="vector_double") {
    return Property::vector_double_t;
  } else if (stype=="spectrum") {
    return Property::spectrum_t;
  } else {
    GYOTO_ERROR("unimplemeted Python property type");
  }
  return Property::empty_t; // avoid warning, will never reach here
}
