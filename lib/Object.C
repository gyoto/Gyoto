#include "GyotoObject.h"
#include "GyotoError.h"
#include "GyotoFactoryMessenger.h"

#include <iostream>

using namespace std ;
using namespace Gyoto ;

Property const * const Object::properties = NULL;

Gyoto::Object::Object(std::string const &name):kind_(name) {}
Gyoto::Object::Object():kind_("") {}
Gyoto::Object::Object(Object const &o):kind_(o.kind_) {}
Gyoto::Object::~Object() {}



void Object::set(Property const &p, double val) {
  if (p.type != Property::double_t)
    throwError("Property is no a double");
  Property::set_double_t set = p.setter.set_double;
  if (!set) throwError("Can't set this Property");
  (this->*set)(val);
}

void Object::set(Property const &p, double val, std::string const &unit) {
  if (p.type != Property::double_t)
    throwError("Property is no a double");
  Property::set_double_unit_t setu = p.setter_unit.set_double;
  if (setu) {
    (this->*setu)(val, unit);
  } else {
    if (unit != "") throwError("Can't set this property with unit");
    set(p, val);
  }
}

void Object::set(Property const &p, bool val) {
  if (p.type != Property::bool_t)
    throwError("Property is no a bool");
  Property::set_bool_t set = p.setter.set_bool;
  if (!set) throwError("Can't set this Property");
  (this->*set)(val);
}

void Object::set(Property const &p, string const &val) {
  throwError("called");
  if (p.type != Property::bool_t)
    throwError("Property is no a bool");
  cerr << "Object::set("<<p.name<<", "<<val<<")"<<std::endl<<std::endl;
  set(p, val==p.name);
}

// void Object::set(std::string const pname, double val) {
//   Property const * const p = property(pname);
//   if (!p) throwError ("No such Property");
//   set(*p, val);
// }

void Object::get(Property const &p, double &val) {
  if (p.type != Property::double_t)
    throwError("Property is no a double");
  Property::get_double_t get = p.getter.get_double;
  if (!get) throwError("Can't get this Property");
  val = (this->*get)();
}

void Object::get(Property const &p, double &val, std::string const &unit) {
  if (p.type != Property::double_t)
    throwError("Property is no a double");
  Property::get_double_unit_t getu = p.getter_unit.get_double;
  if (getu) {
    val = (this->*getu)(unit);
  } else {
    if (unit != "") throwError("Can't get this Property with unit");
    get(p, val);
  }
}


void Object::get(Property const &p, bool &val) {
  if (p.type != Property::bool_t)
    throwError("Property is no a bool");
  Property::get_bool_t get = p.getter.get_bool;
  if (!get) throwError("Can't get this Property");
  val = (this->*get)();
}

void Object::set(Property const &p, string &val) {
  if (p.type != Property::bool_t)
    throwError("Property is no a bool");
  bool bval;
  get(p, bval);
  if (bval) val=p.name;
  else val=p.name_false;
}

// void Object::get(std::string const pname, double &val) {
//   Property const * const p = property(pname);
//   if (!p) throwError ("No such Property");
//   get(*p, val);
// }

Property const * Object::property(std::string const pname) const {
  return getProperties() -> find(pname);
}

Property const * Object::getProperties() const {
  return properties;
}


#ifdef GYOTO_USE_XERCES
void Object::fillProperty(Gyoto::FactoryMessenger *fmp, Property const &p) {
  string name=p.name;
  switch (p.type) {
  case Property::bool_t:
    {
      bool val;
      get(p, val);
      fmp->setParameter(val?name:p.name_false);
    }
    break;
  case Property::double_t:
    {
      double val;
      get(p, val);
      fmp->setParameter(name, val);
    }
    break;
  default:
    throwError("Unimplemented");
  }
  Property const * const * parent = p.parents;
  if (parent) {
    for ( ; *parent; ++parent) {
      fillProperty(fmp, **parent);
    } 
  }
}

void Object::fillElement(Gyoto::FactoryMessenger *fmp) {
  fmp -> setSelfAttribute("kind", kind_);

  Property const * prop = getProperties();
  if (prop) fillProperty(fmp, *prop);
}

void Object::setParameters(Gyoto::FactoryMessenger *fmp)  {
  string name="", content="", unit="";
  if (fmp)
    while (fmp->getNextParameter(&name, &content, &unit))
      setParameter(name, content, unit);
}

#endif

int Object::setParameter(string name, string content, string unit) {
  Property const * prop = property(name);
  if (!prop) {
    std::string errmsg="Object::setParameter(): no such property: ";
    errmsg += name;
    throwError(errmsg);
  }
  switch (prop->type) {
  case Property::bool_t:
    Object::set(*prop, name==prop->name);
    return 0;
  case Property::double_t:
    set(*prop, atof(content.c_str()), unit);
    return 0;
  default:;
  }
  return 1;
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
