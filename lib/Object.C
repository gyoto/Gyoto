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
  if (p.type == Property::bool_t) {
    cerr << "Object::set("<<p.name<<", "<<val<<")"<<std::endl<<std::endl;
    set(p, val==p.name);
  } else if (   p.type == Property::string_t
	     || p.type == Property::filename_t) {
    Property::set_string_t set = p.setter.set_string;
    if (!set) throwError("Can't set this Property");
    (this->*set)(val);
  } else throwError("Cannot set this Property from a string");
}

void Object::set(Property const &p, vector<double> const &val) {
  if (p.type != Property::vector_double_t)
    throwError("Property is not a vector<double>");
  Property::set_vector_double_t set = p.setter.set_vdouble;
  if (!set) throwError("Can't set this Property");
  (this->*set)(val);
}

// void Object::set(std::string const pname, double val) {
//   Property const * const p = property(pname);
//   if (!p) throwError ("No such Property");
//   set(*p, val);
// }

void Object::get(Property const &p, double &val) const {
  if (p.type != Property::double_t)
    throwError("Property is no a double");
  Property::get_double_t get = p.getter.get_double;
  if (!get) throwError("Can't get this Property");
  val = (this->*get)();
}

void Object::get(Property const &p, double &val, std::string const &unit) const {
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


void Object::get(Property const &p, bool &val) const {
  if (p.type != Property::bool_t)
    throwError("Property is no a bool");
  Property::get_bool_t get = p.getter.get_bool;
  if (!get) throwError("Can't get this Property");
  val = (this->*get)();
}

void Object::get(Property const &p, string &val) const {
  if (p.type == Property::bool_t) {
    bool bval;
    get(p, bval);
    if (bval) val=p.name;
    else val=p.name_false;
  } else if (   p.type == Property::string_t
	     || p.type == Property::filename_t) {
    Property::get_string_t get = p.getter.get_string;
    if (!get) throwError("Can't get this Property as a string");
    val = (this->*get)();
  }
}

void Object::get(Property const &p, vector<double> &val) const {
  if (p.type != Property::vector_double_t)
    throwError("Property is not a vector<double>");
  Property::get_vector_double_t get = p.getter.get_vdouble;
  if (!get) throwError("Can't get this Property");
  val = (this->*get)();
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
void Object::fillProperty(Gyoto::FactoryMessenger *fmp, Property const &p) const {
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
  case Property::string_t:
  case Property::filename_t:
    {
      string val;
      get(p, val);
      fmp->setParameter(name, val);
    }
    break;
  case Property::vector_double_t:
    {
      std::vector<double> val;
      get(p, val);
      fmp->setParameter(name, val);
    }
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
      cerr << "Setting '" << name << "' to '" << content
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
  cerr << "Done processing parameters" << endl;
}

#endif

void Object::setParameter(Property const &p, string const &name,
			  string const & content, string const & unit) {
  switch (p.type) {
  case Property::bool_t:
    Object::set(p, name==p.name);
    return;
  case Property::double_t:
    set(p, atof(content.c_str()), unit);
    return;
  case Property::filename_t:
  case Property::string_t:
    set(p, content);
    return;
  case Property::vector_double_t:
    set(p, FactoryMessenger::parseArray(content));
    return;
  default:
    throwError("Property type unimplemented in Object::setParameter()");
  }
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
