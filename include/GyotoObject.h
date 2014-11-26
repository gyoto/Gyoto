/**
 * \file GyotoObject.h
 * \brief Introspectbale objects 
 */

/*
    Copyright 2014 Thibaut Paumard

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


#ifndef __GyotoObject_H_
#define __GyotoObject_H_

#include "GyotoConfig.h"
#include <string>
#include <vector>

namespace Gyoto {
  class Object;
  class Property;
  class FactoryMessenger;
}

/// Make an NULL-terminated array of ancestors
/**
 * Called automatically in GYOTO_PROPERTY_*
 */
#define GYOTO_PROPERTY_MAKE_ANCESTORS(name, ancestor)	\
  Property const * const name##_ancestors[] = {ancestor, NULL}

/// Define a new Property of type Bool
/*
 * Declares a static variable named "name". name and namef should not
 * be quoted.
 *
 * \param[in] class name
 * \param[in] name name of property if true;
 * \param[in] namef name of property if false;
 * \param[in] fname name of functions for setting or getting the property
 * \param[in] ancestor pointer to next Property instance
 */
#define GYOTO_PROPERTY_BOOL(class, name, namef, fname, ancestor) \
  GYOTO_PROPERTY_MAKE_ANCESTORS(name, ancestor); \
  Property const name				 \
  (#name,					 \
   #namef,								\
   (Gyoto::Property::set_bool_t)&class :: fname,			\
   (Gyoto::Property::get_bool_t)&class :: fname,			\
   name##_ancestors)

/// Define a Property of type Double
#define GYOTO_PROPERTY_DOUBLE(class, name, fname, ancestor)	\
  GYOTO_PROPERTY_MAKE_ANCESTORS(name, ancestor); \
  Property const name \
        (#name, \
	   (Gyoto::Property::set_double_t)&class::fname,	\
	   (Gyoto::Property::get_double_t)&class::fname,	\
         name##_ancestors)

/// Define a Property of type Double supporting unit
#define GYOTO_PROPERTY_DOUBLE_UNIT(class, name, fname, ancestor) \
  GYOTO_PROPERTY_MAKE_ANCESTORS(name, ancestor); \
  Property const name \
        (#name, \
	 (Gyoto::Property::set_double_t)&class::fname,	\
	 (Gyoto::Property::get_double_t)&class::fname,	\
	 (Gyoto::Property::set_double_unit_t)&class::fname,	\
	 (Gyoto::Property::get_double_unit_t)&class::fname,	\
         name##_ancestors)

/// Define a Property of type Filename
#define GYOTO_PROPERTY_FILENAME(class, name, fname, ancestor)	\
  GYOTO_PROPERTY_MAKE_ANCESTORS(name, ancestor); \
  Property const name \
        (#name, \
	   (Gyoto::Property::set_string_t)&class::fname,	\
	   (Gyoto::Property::get_string_t)&class::fname,	\
         name##_ancestors, true)

/// Define a Property of type String
#define GYOTO_PROPERTY_STRING(class, name, fname, ancestor)	\
  GYOTO_PROPERTY_MAKE_ANCESTORS(name, ancestor); \
  Property const name \
        (#name, \
	   (Gyoto::Property::set_string_t)&class::fname,	\
	   (Gyoto::Property::get_string_t)&class::fname,	\
         name##_ancestors, false)

/// Define a Property of type String
#define GYOTO_PROPERTY_VECTOR_DOUBLE(class, name, fname, ancestor)	\
  GYOTO_PROPERTY_MAKE_ANCESTORS(name, ancestor); \
  Property const name \
        (#name, \
	   (Gyoto::Property::set_vector_double_t)&class::fname,	\
	   (Gyoto::Property::get_vector_double_t)&class::fname,	\
         name##_ancestors)

/// Define class::properties and class::getProperties() 
#define GYOTO_PROPERTY_FINALIZE(class, property)		\
  Property const * const class::properties = &property;		\
  Property const * class::getProperties() const {		\
    return class::properties;					\
 }

/// Declare  class::properties and class::getProperties()
#define GYOTO_PROPERTY \
  static Property const * const properties;		\
  virtual Property const * getProperties() const



/**
 * \brief Object with properties
 */
class Gyoto::Object
{
 protected:
  std::string kind_;
 public:
  static Property const * const properties;
  Object (std::string const &kind) ;
  Object () ;
  Object (Object const &orig) ;
  virtual ~Object();
  void set(Property const &p, double val);
  void set(Property const &p, double val, std::string const &unit);
  void set(Property const &p, bool val);
  void set(Property const &p, std::string const &val);
  void set(Property const &p, std::vector<double> const &val);
  //void set(std::string const pname, double val);
  void get(Property const &p, double &val) const ;
  void get(Property const &p, double &val, std::string const &unit) const ;
  void get(Property const &p, bool &val) const ;
  void get(Property const &p, std::string &val) const ;
  void get(Property const &p, std::vector<double> &val) const;
  //void get(std::string const pname, double &val);
  Property const * property(std::string const pname) const;
  virtual Property const * getProperties() const;

#ifdef GYOTO_USE_XERCES
  void fillProperty(Gyoto::FactoryMessenger *fmp, Property const &p) ;

  /// Called from Factory
  /**
   * Object implementations should impement fillElement to save their
   * parameters to XML and call the Metric::fillElement(fmp) for the
   * shared properties
   */
  virtual void fillElement(Gyoto::FactoryMessenger *fmp) ;

  /**
   * \brief Main loop in Subcontractor_t function
   *
   * The Subcontractor_t function for each Metric kind should look
   * somewhat like this (templated as
   * Gyoto::Metric::Subcontractor<MyKind>):
\code
SmartPointer<Metric::Generic>
Gyoto::Metric::MyKind::Subcontractor(FactoryMessenger* fmp) {
  SmartPointer<MyKind> gg = new MyKind();
  gg -> setParameters(fmp);
  return gg;
}
\endcode
   *
   * Each metric kind should implement setParameter(string name,
   * string content, string unit) to interpret the individual XML
   * elements. setParameters() can be overloaded in case the specific
   * Metric class needs low level access to the FactoryMessenger. See
   * Gyoto::Astrobj::UniformSphere::setParameters().
   */
  virtual void setParameters(Gyoto::FactoryMessenger *fmp) ;
#endif
  /**
   * \brief Set parameter by name
   *
   * Assume MyKind is a subclass of Metric::Generic which has two
   * members (a string StringMember and a double DoubleMember):
\code
int MyKind::setParameter(std::string name, std::string content, std::string unit) {
 if      (name=="StringMember") setStringMember(content);
 else if (name=="DoubleMember") setDoubleMemeber(atof(content.c_str()), unit);
 else return Generic::setParameter(name, content, unit);
 return 0;
}
\endcode
   * If MyKind is not a direct subclass of Generic, it should call the
   * corresponding setParameter() implementation instead of
   * Generic::setParameter().
   *
   * \param name XML name of the parameter
   * \param content string representation of the value
   * \param unit string representation of the unit
   * \return 0 if this parameter is known, 1 if it is not.
   */
  virtual int setParameter(std::string name,
			    std::string content,
			    std::string unit);

  virtual void setParameter(Gyoto::Property const &p,
			    std::string const &name,
			    std::string const &content,
			    std::string const &unit);

};


/**
 * \brief Property
 */
class Gyoto::Property
{
 private:

 public:
  enum type_e {double_t, long_t, bool_t, string_t, filename_t,
	       vector_double_t};
  std::string name;
  std::string name_false;
  int type;
  typedef void (Object::* set_double_t)(double val);
  typedef double (Object::* get_double_t)() const;
  typedef void (Object::* set_double_unit_t)(double val,
					     std::string const &unit);
  typedef double (Object::* get_double_unit_t)(std::string const &unit) const;
  typedef void (Object::* set_long_t)(double val);
  typedef long (Object::* get_long_t)() const;
  typedef void (Object::* set_bool_t)(bool val);
  typedef bool (Object::* get_bool_t)() const;
  typedef void (Object::* set_string_t)(std::string const&);
  typedef std::string (Object::* get_string_t)() const;
  typedef void (Object::* set_fname_t)(std::string const&);
  typedef std::string (Object::* get_fname_t)() const;
  typedef void (Object::* set_vector_double_t)(std::vector<double> const&);
  typedef std::vector<double> (Object::* get_vector_double_t)() const;
  union setter_t {
    set_double_t set_double;
    set_long_t set_long;
    set_bool_t set_bool;
    set_string_t set_string;
    set_vector_double_t set_vdouble;
  };
  union getter_t {
    get_double_t get_double;
    get_long_t get_long;
    get_bool_t get_bool;
    get_string_t get_string;
    get_vector_double_t get_vdouble;
  };
  union setter_unit_t {set_double_unit_t set_double;};
  union getter_unit_t {get_double_unit_t get_double;};
  setter_t setter;
  getter_t getter;
  setter_unit_t setter_unit;
  getter_unit_t getter_unit;
  Property const * const  * const parents;
  
  Property(std::string name,
	   set_double_t set_double,
	   get_double_t get_double,
	   Property const * const * ancestors);

  Property(std::string name,
	   set_double_t set_double,
	   get_double_t get_double,
	   set_double_unit_t set_double_unit,
	   get_double_unit_t get_double_unit,
	   Property const * const * ancestors);

  Property(std::string name,
	   std::string name_false,
	   set_bool_t set_bool,
	   get_bool_t get_bool,
	   Property const * const * ancestors);

  Property(std::string name,
	   set_string_t set_string,
	   get_string_t get_string,
	   Property const * const * ancestors,
	   bool is_filename=false);

  Property(std::string name,
	   set_vector_double_t set_vdouble,
	   get_vector_double_t get_vdouble,
	   Property const * const * ancestors);

  Property const * find(std::string name) const;

};

#endif
