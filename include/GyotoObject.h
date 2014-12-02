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
  class Value;
  class FactoryMessenger;
}

/// Declare  class::properties and class::getProperties()
#define GYOTO_OBJECT \
  static Property const  properties[];		\
  virtual Property const * getProperties() const

/**
 * \brief Object with properties
 */
class Gyoto::Object
{
 protected:
  std::string kind_;
 public:
  GYOTO_OBJECT;
  Object (std::string const &kind) ;
  Object () ;
  Object (Object const &orig) ;
  virtual ~Object();
  void set(Property const &p, Value val);
  void set(Property const &p, Value val, std::string const &unit);

  Value get(Property const &p) const;
  Value get(Property const &p, std::string const &unit) const;

  Property const * property(std::string const pname) const;

#ifdef GYOTO_USE_XERCES
  virtual void fillProperty(Gyoto::FactoryMessenger *fmp, Property const &p) const ;

  /// Called from Factory
  /**
   * Object implementations should impement fillElement to save their
   * parameters to XML and call the Metric::fillElement(fmp) for the
   * shared properties
   */
  virtual void fillElement(Gyoto::FactoryMessenger *fmp) const ;

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

#endif
