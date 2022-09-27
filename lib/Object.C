/*
    Copyright 2014-2016, 2019-2020 Thibaut Paumard

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
#include "GyotoObject.h"
#include "GyotoProperty.h"
#include "GyotoValue.h"
#include "GyotoError.h"
#include "GyotoFactoryMessenger.h"
#include "GyotoMetric.h"
#include "GyotoAstrobj.h"
#include "GyotoSpectrum.h"
#include "GyotoSpectrometer.h"
#include "GyotoScreen.h"

#include <iostream>

using namespace std ;
using namespace Gyoto ;

// We do have a propert list. In contains a single item, wich is a
// link to the NULL pointer, meaning end-of-the list.  This is used to
// terminate the property list of our descendents in a
// forward-compatible manner, i.e., we may well add very generic
// Properties in the future.
GYOTO_PROPERTY_START(Gyoto::Object, "Object with properties.")
GYOTO_PROPERTY_END(Object, NULL)


Gyoto::Object::Object(std::string const &name):kind_(name), plugins_() {}
Gyoto::Object::Object():kind_(""), plugins_() {}
Gyoto::Object::Object(Object const &o):kind_(o.kind_), plugins_(o.plugins_) {}
Gyoto::Object::~Object() {}

// Output

string Gyoto::Object::kind() const {return kind_;}
void Gyoto::Object::kind(const string src) { kind_ = src;}

bool Object::isThreadSafe() const {
  /**
   * The default behaviour is to consider that everything is
   * thread-safe (for the purpose of threads in
   * Gyoto::Scenery::rayTrace()).
   *
   * For Objects that have other Object as children, we need to
   * recursively ask those children whether they are thread-safe and
   * accumulate the answer. It can be done in a generic manner as long
   * as the Object declares its children as properties, this is what
   * we do here.
   *
   * Objects that are never thread-safe must reimplement this function
   * to simply return "false", which can be done with the pair of
   * macros GYOTO_OBJECT_THREAD_SAFETY/GYOTO_PROPERTY_THREAD_UNSAFE
   * respectively in the class declaration and definition.
   *
   * Objects that need to clone children that are not declared as
   * properties in their copy constructors need to reimplement this
   * method to take care of those children.
   */
  bool safe = true;
  Property const * prop = getProperties();
  SmartPointer<SmartPointee> child=NULL;
  while (prop) {
    if (*prop) {
      switch (prop -> type) {
      case Property::metric_t:
	child=SmartPointer<Metric::Generic>(get(*prop));
	break;
      case Property::screen_t:
	child=SmartPointer<Screen>(get(prop));
	break;
      case Property::astrobj_t:
	child=SmartPointer<Astrobj::Generic>(get(*prop));
	break;
      case Property::spectrum_t:
	child=SmartPointer<Spectrum::Generic>(get(*prop));
	break;
      case Property::spectrometer_t:
	child=SmartPointer<Spectrometer::Generic>(get(*prop));
	break;
      default:
	child=NULL;
      }
      if (child) safe &= dynamic_cast<Object const*>(child()) -> isThreadSafe();
      ++prop;
    } else {
      prop=prop->parent;
    }
  }
  GYOTO_DEBUG_EXPR(safe);
  return safe;
}

void Object::set(Property const &p,
		 Value val,
		 std::string const &unit) {
  GYOTO_DEBUG_EXPR(p.type);
  switch (p.type) {
  case Property::empty_t:
    GYOTO_ERROR("Attempt to set empty_t Property");
    return;
  case Property::double_t:
    {
      Property::set_double_unit_t setu = p.setter_unit.set_double;
      if (setu) {
	GYOTO_DEBUG << "double Property which supports unit" << endl;
	(this->*setu)(val, unit);
      } else {
	GYOTO_DEBUG << "double Property which does not support unit" << endl;
	if (unit != "") GYOTO_ERROR("Can't set this property with unit");
	set(p, val);
      }
    }
    return;
  case Property::vector_double_t:
    {
      Property::set_vector_double_unit_t setu = p.setter_unit.set_vdouble;
      if (setu) {
	GYOTO_DEBUG << "vector<double> Property which supports unit" << endl;
	(this->*setu)(val, unit);
      } else {
	GYOTO_DEBUG << "vector<double> Property which does not support unit" << endl;
	if (unit != "") GYOTO_ERROR("Can't set this property with unit");
	set(p, val);
      }
    }
    return;
  default:
    GYOTO_DEBUG<< "Not a double_t, vector_double_t or empty_t Property" << endl;
    if (unit != "")
      GYOTO_ERROR("Can't set this property with unit (not a double)");
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
      if (!set) GYOTO_ERROR("Can't set this Property");	\
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
    ___local_case(size_t);
  case Property::filename_t:
    ___local_case(string);
  case Property::vector_double_t:
    {
      Property::set_vector_double_t set = p.setter.set_vdouble;
      if (!set) GYOTO_ERROR("Can't set this Property");
      (this->*set)(val);
    }
    break;
  case Property::vector_unsigned_long_t:
    {
      Property::set_vector_unsigned_long_t set = p.setter.set_vulong;
      if (!set) GYOTO_ERROR("Can't set this Property");
      (this->*set)(val);
    }
    break;
    ___local_case(metric);
    ___local_case(astrobj);
    ___local_case(spectrum);
    ___local_case(spectrometer);
    ___local_case(screen);
  default:
    GYOTO_ERROR("Unimplemented Property type in Object::set");
  }
# undef ___local_case
}

void Object::set(std::string const &pname, Value val) {
  GYOTO_DEBUG_EXPR(pname);
  Property const * p = property(pname);
  if (!p) GYOTO_ERROR("No Property by that name");
  set(*p, ((p->type == Property::bool_t && pname == p->name_false)?
	  Value(!val):val));
}

void Object::set(std::string const &pname, Value val, std::string const &unit) {
  Property const * p = property(pname);
  if (!p) GYOTO_ERROR("No Property by that name");
  set(*p, ((p->type == Property::bool_t && pname == p->name_false)?
	  Value(!val):val), unit);
}

Value Object::get(Property const &p,
			    std::string const &unit) const {

  if (p.type == Property::double_t) {
    Property::get_double_unit_t getu = p.getter_unit.get_double;
    if (getu) return (this->*getu)(unit);
    if (unit != "") GYOTO_ERROR("Can't get this property with unit");
    return get(p);
  }

  if (p.type == Property::vector_double_t) {
    Property::get_vector_double_unit_t getu = p.getter_unit.get_vdouble;
    if (getu) return (this->*getu)(unit);
    if (unit != "") GYOTO_ERROR("Can't get this property with unit");
    return get(p);
  }

  if (unit != "")
    GYOTO_ERROR("Can't set this property with unit (not a double)");

  return get(p);
}

Value Object::get(Property const &p) const {
# define ___local_case(type) \
  case Property::type##_t:     \
    {			     \
    Property::get_##type##_t get = p.getter.get_##type;	\
    if (!get) GYOTO_ERROR("Can't get this Property");	\
    val = Value((this->*get)());			\
    }							\
    break

  Gyoto::Value val;
  switch (p.type) {
  case Property::empty_t:
    GYOTO_ERROR("Can't get empty property");
    ___local_case(bool);
    ___local_case(double);
    ___local_case(long);
    ___local_case(unsigned_long);
    ___local_case(size_t);
  case Property::filename_t:
    ___local_case(string);
  case Property::vector_double_t:
    {
      Property::get_vector_double_t get = p.getter.get_vdouble;
      if (!get) GYOTO_ERROR("Can't get this Property");
      val = (this->*get)();
    }
    break;
  case Property::vector_unsigned_long_t:
    {
      Property::get_vector_unsigned_long_t get = p.getter.get_vulong;
      if (!get) GYOTO_ERROR("Can't get this Property");
      val = (this->*get)();
    }
    break;
    ___local_case(metric);
    ___local_case(astrobj);
    ___local_case(spectrum);
    ___local_case(spectrometer);
    ___local_case(screen);
  default:
    GYOTO_ERROR("Unimplemented Property type in Object::get");
  }
  return val;
# undef ___local_case
}

Value Object::get(std::string const &pname) const {
  Property const * p = property(pname);
  if (!p) GYOTO_ERROR("No Property by that name");
  Value res = get(*p);
  if (p->type == Property::bool_t && pname == p->name_false)
    return !bool(res);
  return res;
}

Value Object::get(std::string const &pname, std::string const &unit) const{
  Property const * p = property(pname);
  if (!p) GYOTO_ERROR("No Property by that name");
  Value res = get(*p, unit);
  if (p->type == Property::bool_t && pname == p->name_false)
    return !bool(res);
  return res;
}

Property const * Object::property(std::string const pname) const {
  Property const * prop = getProperties(); 
  while (prop) {
    if (*prop) {
      GYOTO_DEBUG_EXPR(prop->name);
      if (prop->name == pname ||
	  (prop->type==Property::bool_t && prop->name_false == pname))
	return prop;
      ++prop;
    } else prop=prop->parent;
  }
  return NULL;
}

#ifdef GYOTO_USE_XERCES
void Object::fillProperty(Gyoto::FactoryMessenger *fmp, Property const &p) const {
  GYOTO_DEBUG_EXPR(fmp);
  GYOTO_DEBUG_EXPR(p.name);
  GYOTO_DEBUG_EXPR(p.type);
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
  case Property::size_t_t:
    fmp->setParameter(name, (unsigned long)(size_t(get(p))));
    break;
  case Property::double_t:
    fmp->setParameter(name, double(get(p)));
    break;
  case Property::string_t:
  case Property::filename_t:
    fmp->setParameter(name, std::string(get(p)));
    break;
  case Property::vector_double_t:
    fmp->setParameter(name, get(p).operator std::vector<double>());
    break;
  case Property::vector_unsigned_long_t:
    fmp->setParameter(name, get(p).operator std::vector<unsigned long>());
    break;
  case Property::metric_t:
    fmp->metric(get(p));
    break;
  case Property::astrobj_t:
    fmp->astrobj(get(p));
    break;
  case Property::screen_t:
    fmp->screen(get(p));
    break;
  case Property::spectrum_t:
    {
      SmartPointer<Spectrum::Generic> sp=get(p);
      if (!sp) return;
      childfmp = fmp -> makeChild ( name );
      sp -> fillElement(childfmp);
      delete childfmp;
    }
    break;
  case Property::spectrometer_t:
    {
      SmartPointer<Spectrometer::Generic> spr=get(p);
      if (!spr) return;
      childfmp = fmp -> makeChild ( name );
      spr -> fillElement(childfmp);
      delete childfmp;
    }
    break;
  default:
    GYOTO_ERROR("Property type unimplemented in Object::fillProperty()");
  }
}

void Object::fillElement(Gyoto::FactoryMessenger *fmp) const {
  std::vector<std::string> const plgs=plugins();
  size_t np=plgs.size();
  if (np) {
    std::string plg(plgs[0]);
    for (size_t i=1; i<np; ++np) {
      plg += std::string(",") +plgs[i] ;
    }
    fmp -> setSelfAttribute("plugin", plg);
  }
  if (kind_ != "") fmp -> setSelfAttribute("kind", kind_);
  Property const * prop = getProperties(); 
  while (prop) {
    if (*prop) {
      if (prop->type != Property::empty_t)
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
	GYOTO_DEBUG << "'" << name << "' not found, calling setParameter()"
		    << endl;
	// The specific setParameter() implementation may well know
	// this entity
	setParameter(name, content, unit);
      } else {
	GYOTO_DEBUG << "'" << name << "' found "<< endl;
	std::vector<std::string> plugins;
	switch (prop->type) {
	case Property::metric_t:
	  set(*prop, fmp->metric());
	  break;
	case Property::astrobj_t:
	  set(*prop, fmp->astrobj());
	  break;
	case Property::screen_t:
	  set(*prop, fmp->screen());
	  break;
	case Property::spectrum_t:
	  content = fmp -> getAttribute("kind");
	  child = fmp -> getChild();
	  plugins = Gyoto::split(fmp -> getAttribute("plugin"), ",");
	  set(*prop, (*Spectrum::getSubcontractor(content, plugins))(child, plugins) );
	  delete child;
	  break;
	case Property::spectrometer_t:
	  content = fmp -> getAttribute("kind");
	  child = fmp -> getChild();
	  plugins = Gyoto::split(fmp -> getAttribute("plugin"), ",");
	  set(*prop, (*Spectrometer::getSubcontractor(content, plugins))(child, plugins) );
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
  GYOTO_DEBUG_EXPR(name);
  Value val;
  GYOTO_DEBUG_EXPR(p.type);
  GYOTO_DEBUG_EXPR(Property::double_t);
  switch (p.type) {
  case Property::bool_t:
    val = (name==p.name);
    break;
  case Property::long_t:
    val = strtol(content.c_str(), NULL, 0);
    break;
  case Property::unsigned_long_t:
    val = strtoul(content.c_str(), NULL, 0);
    break;
  case Property::size_t_t:
    val = size_t(strtoul(content.c_str(), NULL, 0));
    break;
  case Property::double_t:
    val = Gyoto::atof(content.c_str());
    GYOTO_DEBUG << "calling set(p, val, unit)" << std::endl;
    set(p, val, unit);
    return;
  case Property::filename_t:
  case Property::string_t:
    val = content;
    break;
  case Property::vector_double_t:
    val = FactoryMessenger::parseArray(content);
    set(p, val, unit);
    return;
  case Property::vector_unsigned_long_t:
    val = FactoryMessenger::parseArrayULong(content);
    break;
  case Property::metric_t:
    GYOTO_ERROR("Metric can't be set using setParameter()");
  default:
    GYOTO_ERROR("Property type unimplemented in Object::setParameter()");
  }
  GYOTO_DEBUG << "calling set" << std::endl;
  set(p, val);
}

int Object::setParameter(string name, string content, string unit) {
  Property const * prop = property(name);
  if (!prop) {
    size_t pos=name.find("::");
    if (pos != string::npos) {
      string childname = name.substr(0,pos);
      name=name.substr(pos+2);
      prop=property(childname);
      if (!prop) return 1;
      Object * obj=NULL;
      Value val=get(*prop);
      switch (prop->type) {
      case Property::screen_t:
	obj = SmartPointer<Screen>(val);
	break;
      case Property::metric_t:
	obj = SmartPointer<Metric::Generic>(val);
	break;
      case Property::astrobj_t:
	obj = SmartPointer<Astrobj::Generic>(val);
	break;
      case Property::spectrum_t:
	obj = SmartPointer<Spectrum::Generic>(val);
	break;
      case Property::spectrometer_t:
	obj = SmartPointer<Spectrometer::Generic>(val);
	break;
      default:
	GYOTO_ERROR(childname+" is not an object");
      }
      if (obj) return obj -> setParameter(name, content, unit);
      GYOTO_ERROR(childname+" not set yet");
    }
    return 1;
  }
  setParameter(*prop, name, content, unit);
  return 0;
}

std::string Object::describeProperty(Property const &p) const {
  string out=p.name+": ";
  switch (p.type) {
  case Property::empty_t:
    return "";
  case Property::bool_t:
    out = p.name + "/" + p.name_false + ": bool";
    break;
  case Property::long_t:
    out += "long";
    break;
  case Property::unsigned_long_t:
    out += "unsigned long";
    break;
  case Property::size_t_t:
    out += "size_t";
    break;
  case Property::double_t:
    out += "double";
    if (p.setter_unit.set_double) out += " with unit";
    break;
  case Property::string_t:
    out += "string";
    break;
  case Property::filename_t:
    out += "filename";
    break;
  case Property::vector_double_t:
    out += "vector<double>";
    if (p.setter_unit.set_vdouble) out += " with unit";
    break;
  case Property::vector_unsigned_long_t:
    out += "vector<unsigned> long";
    break;
  case Property::metric_t:
    out += "Gyoto::Metric::Generic";
    break;
  case Property::astrobj_t:
    out += "Gyoto::Astrobj::Generic";
    break;
  case Property::screen_t:
    out += "Gyoto::Screen";
    break;
  case Property::spectrum_t:
    out += "Gyoto::Spectrum";
    break;
  case Property::spectrometer_t:
    out += "Gyoto::Spectrometer";
    break;
  default:
    GYOTO_ERROR("Property type unimplemented in Object::fillProperty()");
  }
  return out;
}

void Object::help() const {
  Property const * prop = getProperties();
  while (prop) {
    if (*prop) {
      if (prop->type==Property::empty_t) {
	cout << prop->name << endl;
	if (prop->doc != "") cout << "  " << prop->doc << endl;
      } else {
	cout << "\t"<< describeProperty(*prop) << endl;
	if (prop->doc != "") cout << "\t  " << prop->doc << endl;
      }
      ++prop;
    } else prop=prop->parent;
  }
}
