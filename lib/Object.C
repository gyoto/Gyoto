#include "GyotoObject.h"
#include "GyotoError.h"
#include "GyotoFactoryMessenger.h"

#include <iostream>

using namespace std ;
using namespace Gyoto ;

GYOTO_PROPERTY_FINALIZE(Object, NULL);

/// Value

Property::Value::Value() {}
Property::Value::~Value() {}

Property::Value::Value(double val) : Double(val) {}
Property::Value::operator double() const {return Double;}

Property::Value::Value(bool val) : Bool(val) {}
Property::Value::operator bool() const {return Bool;}

Property::Value::Value(long val) : Long(val) {}
Property::Value::operator long() const {return Long;}

Property::Value::Value(std::string val) : String(val) {}
Property::Value::operator std::string() const {return String;}

Property::Value::Value(std::vector<double> val) : VDouble(val) {}
Property::Value::operator std::vector<double>() const {return VDouble;}

/// Object

Gyoto::Object::Object(std::string const &name):kind_(name) {}
Gyoto::Object::Object():kind_("") {}
Gyoto::Object::Object(Object const &o):kind_(o.kind_) {}
Gyoto::Object::~Object() {}


void Object::set(Property const &p,
		 Property::Value const &val,
		 std::string const &unit) {

  if (p.type == Property::double_t) {
    Property::set_double_unit_t setu = p.setter_unit.set_double;
    if (setu) {
      (this->*setu)(val, unit);
    } else {
      if (unit != "") throwError("Can't set this property with unit");
      set(p, val);
    }
    return;
  }

  if (unit != "")
    throwError("Can't set this property with unit (not a double)");

  set(p, val);
  return;

}

void Object::set(Property const &p, Property::Value const &val) {
# define ___local_case(type)			\
  case Property::type##_t:			\
    {						\
      Property::set_##type##_t set = p.setter.set_##type; \
      if (!set) throwError("Can't set this Property");	\
      (this->*set)(val);				\
    }							\
    break
  
  switch (p.type) {
    ___local_case(double);
    ___local_case(bool);
    ___local_case(long);
  case Property::filename_t:
    ___local_case(string);
  case Property::vector_double_t:
    {
      Property::set_vector_double_t set = p.setter.set_vdouble;
      if (!set) throwError("Can't set this Property");
      (this->*set)(val);
    }
    break;
  default:
    throwError("Unimplemented Property type in Object::set");
  }
# undef ___local_case
}

Property::Value Object::get(Property const &p,
			    std::string const &unit) const {

  if (p.type == Property::double_t) {
    Property::get_double_unit_t getu = p.getter_unit.get_double;
    if (getu) return (this->*getu)(unit);
    if (unit != "") throwError("Can't get this property with unit");
    return get(p);
  }

  if (unit != "")
    throwError("Can't set this property with unit (not a double)");

  return get(p);
}

Property::Value Object::get(Property const &p) const {
# define ___local_case(type) \
  case Property::type##_t:     \
    {			     \
    Property::get_##type##_t get = p.getter.get_##type;	\
    if (!get) throwError("Can't get this Property");	\
    val = (this->*get)();				\
    }							\
    break

  Property::Value val;
  switch (p.type) {
    ___local_case(bool);
    ___local_case(double);
    ___local_case(long);
  case Property::filename_t:
    ___local_case(string);
  case Property::vector_double_t:
    {
      Property::get_vector_double_t get = p.getter.get_vdouble;
      if (!get) throwError("Can't get this Property");
      val = (this->*get)();
    }
    break;
  default:
    throwError("Unimplemented Property type in Object::get");
  }
  return val;
# undef ___local_case
}

Property const * Object::property(std::string const pname) const {
  return getProperties() -> find(pname);
}

#ifdef GYOTO_USE_XERCES
void Object::fillProperty(Gyoto::FactoryMessenger *fmp, Property const &p) const {
  string name=p.name;
  switch (p.type) {
  case Property::bool_t:
    fmp->setParameter(get(p)?name:p.name_false);
    break;
  case Property::double_t:
    fmp->setParameter(name, double(get(p)));
    break;
  case Property::string_t:
  case Property::filename_t:
    fmp->setParameter(name, std::string(get(p)));
    break;
  case Property::vector_double_t:
    fmp->setParameter(name, get(p).VDouble);
    break;
  default:
    throwError("Property type unimplemented in Object::fillProperty()");
  }
  Property const * const * parent = p.parents;
  if (parent) {
    for ( ; *parent; ++parent) {
      fillProperty(fmp, **parent);
    } 
  }
}

void Object::fillElement(Gyoto::FactoryMessenger *fmp) const {
  fmp -> setSelfAttribute("kind", kind_);

  Property const * prop = getProperties();
  if (prop) fillProperty(fmp, *prop);
}

void Object::setParameters(Gyoto::FactoryMessenger *fmp)  {
  string name="", content="", unit="";
  if (fmp)
    while (fmp->getNextParameter(&name, &content, &unit)) {
      GYOTO_DEBUG << "Setting '" << name << "' to '" << content
		  << "' (unit='"<<unit<<"')" << endl;
      Property const * prop = property(name);
      if (!prop) {;
	// The specific setParameter() implementation may well know
	// this entity
	setParameter(name, content, unit);
      } else {
	if (prop->type == Property::filename_t)
	  content = fmp->fullPath(content);
	setParameter(*prop, name, content, unit);
      }
    }
  GYOTO_DEBUG << "Done processing parameters" << endl;
}

#endif

void Object::setParameter(Property const &p, string const &name,
			  string const & content, string const & unit) {
  Property::Value val;
  switch (p.type) {
  case Property::bool_t:
    val = (name==p.name);
    break;
  case Property::double_t:
    val = atof(content.c_str());
    set(p, val, unit);
    return;
  case Property::filename_t:
  case Property::string_t:
    val = content;
    break;
  case Property::vector_double_t:
    val = FactoryMessenger::parseArray(content);
    break;
  default:
    throwError("Property type unimplemented in Object::setParameter()");
  }
  set(p, val);
}

int Object::setParameter(string name, string content, string unit) {
  Property const * prop = property(name);
  if (!prop) return 1;
  setParameter(*prop, name, content, unit);
  return 0;
}

//// Property

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
