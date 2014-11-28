#include "GyotoProperty.h"

using namespace std ;
using namespace Gyoto ;

Property::Property(string n, set_long_t set, get_long_t get,
		   Property const * const * ancestors)
  : name(n), type(long_t), parents(ancestors) {
  setter.set_long=set;
  getter.get_long=get;
}

Property::Property(string n,
		   set_unsigned_long_t set, get_unsigned_long_t get,
		   Property const * const * ancestors)
  : name(n), type(unsigned_long_t), parents(ancestors) {
  setter.set_unsigned_long=set;
  getter.get_unsigned_long=get;
}

Property::Property(string n, set_double_t set, get_double_t get,
		   Property const * const * ancestors)
  : name(n), type(double_t), parents(ancestors) {
  setter.set_double=set;
  getter.get_double=get;
  setter_unit.set_double=NULL;
  getter_unit.get_double=NULL;
}

Property::Property(string n, set_double_t set, get_double_t get,
		   set_double_unit_t setu, get_double_unit_t getu,
		   Property const * const * ancestors)
  : name(n), type(double_t), parents(ancestors) {
  setter.set_double=set;
  getter.get_double=get;
  setter_unit.set_double=setu;
  getter_unit.get_double=getu;
}

Property::Property(string n, string nf, set_bool_t set, get_bool_t get,
		   Property const * const * ancestors)
  : name(n), name_false(nf), type(bool_t), parents(ancestors) {
  setter.set_bool=set;
  getter.get_bool=get;
}

Property::Property(string n, set_string_t set, get_string_t get,
		   Property const * const * ancestors,
		   bool is_filename)
  : name(n), type(is_filename?filename_t:string_t), parents(ancestors) {
  setter.set_string=set;
  getter.get_string=get;
}

Property::Property(string n,
		   set_vector_double_t set,
		   get_vector_double_t get,
		   Property const * const * ancestors)
  : name(n), type(vector_double_t), parents(ancestors) {
  setter.set_vdouble=set;
  getter.get_vdouble=get;
}

Property::Property(string n,
		   set_metric_t set,
		   get_metric_t get,
		   Property const * const * ancestors)
  : name(n), type(metric_t), parents(ancestors) {
  setter.set_metric=set;
  getter.get_metric=get;
}

Property::Property(string n,
		   set_spectrum_t set,
		   get_spectrum_t get,
		   Property const * const * ancestors)
  : name(n), type(spectrum_t), parents(ancestors) {
  setter.set_spectrum=set;
  getter.get_spectrum=get;
}

Property const * Property::find(std::string n) const {
  if (this == NULL) return NULL;
  if (name == n || (type == bool_t && name_false == n)) return this;
  if (parents == NULL) return NULL;
  Property const * const * p ;
  Property const * result = NULL;
  for (p=parents;  *p != NULL; ++p) {
    result = (*p)->find(n);
    if (result) break;
  }
  return result; 
}
