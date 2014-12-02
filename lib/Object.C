#include "GyotoObject.h"
#include "GyotoProperty.h"
#include "GyotoValue.h"
#include "GyotoError.h"
#include "GyotoFactoryMessenger.h"
#include "GyotoMetric.h"
#include "GyotoSpectrum.h"

#include <iostream>

using namespace std ;
using namespace Gyoto ;

// We do have a propert list. In contains a single item, wich is a
// link to the NULL pointer, meaning end-of-the list.  This is used to
// terminate the property list of our descendents in a
// forward-compatible manner, i.e., we may well add very generic
// Properties in the future.
GYOTO_PROPERTY_START(Object)
GYOTO_PROPERTY_END(Object, NULL)


Gyoto::Object::Object(std::string const &name):kind_(name) {}
Gyoto::Object::Object():kind_("") {}
Gyoto::Object::Object(Object const &o):kind_(o.kind_) {}
Gyoto::Object::~Object() {}


void Object::set(Property const &p,
		 Value val,
		 std::string const &unit) {
  GYOTO_DEBUG_EXPR(p.type);
  switch (p.type) {
  case Property::empty_t:
    throwError("Attempt to set empty_t Property");
    return;
  case Property::double_t:
    {
      Property::set_double_unit_t setu = p.setter_unit.set_double;
      if (setu) {
	GYOTO_DEBUG << "double Property which supports unit" << endl;
	(this->*setu)(val, unit);
      } else {
	GYOTO_DEBUG << "double Property which does not support unit" << endl;
	if (unit != "") throwError("Can't set this property with unit");
	set(p, val);
      }
    }
    return;
  default:
    GYOTO_DEBUG<< "Not a double_t or empty_t Property" << endl;
    if (unit != "")
      throwError("Can't set this property with unit (not a double)");
    set(p, val);
    return;
  }

}

void Object::set(Property const &p, Value val) {
# define ___local_case(type)			\
  case Property::type##_t:			\
    {						\
      GYOTO_DEBUG <<"Setting property of type " #type << endl;	\
      Property::set_##type##_t set = p.setter.set_##type;	\
      GYOTO_DEBUG_EXPR(set);					\
      if (!set) throwError("Can't set this Property");	\
      (this->*set)(val);				\
    }							\
    break
  
  switch (p.type) {
  case Property::empty_t:
    return;
    ___local_case(double);
    ___local_case(bool);
    ___local_case(long);
    ___local_case(unsigned_long);
  case Property::filename_t:
    ___local_case(string);
  case Property::vector_double_t:
    {
      Property::set_vector_double_t set = p.setter.set_vdouble;
      if (!set) throwError("Can't set this Property");
      (this->*set)(val);
    }
    break;
    ___local_case(metric);
    ___local_case(spectrum);
  default:
    throwError("Unimplemented Property type in Object::set");
  }
# undef ___local_case
}

Value Object::get(Property const &p,
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

Value Object::get(Property const &p) const {
# define ___local_case(type) \
  case Property::type##_t:     \
    {			     \
    Property::get_##type##_t get = p.getter.get_##type;	\
    if (!get) throwError("Can't get this Property");	\
    val = Value((this->*get)());			\
    }							\
    break

  Gyoto::Value val;
  switch (p.type) {
  case Property::empty_t:
    throwError("Can't get empty property");
    ___local_case(bool);
    ___local_case(double);
    ___local_case(long);
    ___local_case(unsigned_long);
  case Property::filename_t:
    ___local_case(string);
  case Property::vector_double_t:
    {
      Property::get_vector_double_t get = p.getter.get_vdouble;
      if (!get) throwError("Can't get this Property");
      val = (this->*get)();
    }
    break;
    ___local_case(metric);
    ___local_case(spectrum);
  default:
    throwError("Unimplemented Property type in Object::get");
  }
  return val;
# undef ___local_case
}

Property const * Object::property(std::string const pname) const {
  Property const * prop = getProperties(); 
  while (prop) {
    if (*prop) {
      if (prop->name == pname) return prop;
      ++prop;
    } else prop=prop->parent;
  }
  return NULL;
}

#ifdef GYOTO_USE_XERCES
void Object::fillProperty(Gyoto::FactoryMessenger *fmp, Property const &p) const {
  FactoryMessenger * childfmp=NULL;
  string name=p.name;
  switch (p.type) {
  case Property::empty_t:
    break;
  case Property::bool_t:
    fmp->setParameter(get(p)?name:p.name_false);
    break;
  case Property::long_t:
    fmp->setParameter(name, long(get(p)));
    break;
  case Property::unsigned_long_t:
    fmp->setParameter(name, (unsigned long)(get(p)));
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
  case Property::metric_t:
    fmp->metric(get(p));
    break;
  case Property::spectrum_t:
    childfmp = fmp -> makeChild ( name );
    get(p).Spectrum -> fillElement(childfmp);
    delete childfmp;
    break;
  default:
    throwError("Property type unimplemented in Object::fillProperty()");
  }
}

void Object::fillElement(Gyoto::FactoryMessenger *fmp) const {
  fmp -> setSelfAttribute("kind", kind_);
  Property const * prop = getProperties(); 
  while (prop) {
    if (*prop) {
      fillProperty(fmp, *prop);
      ++prop;
    } else prop=prop->parent;
  }
}

void Object::setParameters(Gyoto::FactoryMessenger *fmp)  {
  string name="", content="", unit="";
  FactoryMessenger * child = NULL;
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
	switch (prop->type) {
	case Property::metric_t:
	  set(*prop, fmp->metric());
	  break;
	case Property::spectrum_t:
	  content = fmp -> getAttribute("kind");
	  child = fmp -> getChild();
	  set(*prop, (*Spectrum::getSubcontractor(content))(child) );
	  delete child;
	  break;
	case Property::filename_t:
	  content = fmp->fullPath(content);
	  // no 'break;' here, we need to proceed
	default:
	  setParameter(*prop, name, content, unit);
	}
      }
    }
  GYOTO_DEBUG << "Done processing parameters" << endl;
}

#endif

void Object::setParameter(Property const &p, string const &name,
			  string const & content, string const & unit) {
  Value val;
  switch (p.type) {
  case Property::bool_t:
    val = (name==p.name);
    break;
  case Property::long_t:
    val = long(atoi(content.c_str()));
    break;
  case Property::unsigned_long_t:
    val = (unsigned long)(atoi(content.c_str()));
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
  case Property::metric_t:
    throwError("Metric can't be set using setParameter()");
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
