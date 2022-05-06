/**
 * \file GyotoObject.h
 * \brief Introspectable objects 
 */

/*
    Copyright 2014-2016 Thibaut Paumard

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
#include "GyotoSmartPointer.h"
#include <string>
#include <vector>

namespace Gyoto {
  class Object;
  class Property;
  class Value;
  class FactoryMessenger;
}

/// Declare a pair of accessors to string member in a class declaration
/**
 * The accessors must also be defined in the .C file
 *
 * \param method name of the accessors.
 */
#define GYOTO_OBJECT_ACCESSORS_STRING(method)				\
  virtual void method(std::string const&);				\
  virtual std::string method() const;

/// Declare a pair of accessors to scalar member in a class declaration
/**
 * The accessors must also be defined in the .C file, which can be
 * done using #GYOTO_PROPERTY_ACCESSORS
 *
 * \param type data type of the memebr beeing accessed. Any scalar
 *        type (double, long, size_t, SmartPointer<Metric::Generic> etc).
 * \param method name of the accessors.
 */
#define GYOTO_OBJECT_ACCESSORS(type, method)				\
  virtual void method(type);						\
  virtual type method() const;

/// Declare a quadruplet of accessors to double member that supports unit
/**
 * The accessors must also be defined in the .C file.
 *
 * \param method name of the accessors.
 */
#define GYOTO_OBJECT_ACCESSORS_UNIT(method)				\
  GYOTO_OBJECT_ACCESSORS(double, method)				\
  virtual void method(double, std::string const &);			\
  virtual double method(std::string const &) const;

/// Declare  class::properties and class::getProperties()
/**
 * Any derived class that does define Properties (i.e. the macro
 * GYOTO_PROPERTY_START() is called somewhere in the .C file) must
 * call the #GYOTO_OBJECT macro in a "public:" section of the class
 * declaration. Else, the property list is inherited from the direct
 * parent, and calling GYOTO_PROPERTY_START() in the .C file leads to
 * a compile-time error.
 */
#define GYOTO_OBJECT \
  static Property const  properties[];			\
  virtual Property const * getProperties() const;	\
  static const std::string builtinPluginValue;		\
  virtual void plugins(std::vector<std::string> const & plugname);	\
  virtual std::vector<std::string> plugins() const

/// Declare virtual bool isThreadSafe() const
/**
 * Use this to declare that the object is not (or not always) thread
 * safe. The corresponding definition of isThreadSafe() must exist. If
 * the object is always thread unsafe (e.g. Lorene Metrics of Python
 * based objects), you can simply use GYOTO_PROPERTY_THREAD_SAFETY in
 * the corresponding .C file.
 */
#define GYOTO_OBJECT_THREAD_SAFETY		\
  virtual bool isThreadSafe() const

/// Object with properties
/**
 * The Object API allows declaring a list of Properties that can be
 * set and retrieved using a common, text-based interface. This
 * interface simplifies a lot how to read and write XML, as well as
 * writing bindings for interpreted langages (e.g. the Yorick
 * interface).
 *
 * In fact, any class member that has an interface implemented as a
 * Property can be readily read and written from/to XML as well as
 * from the Yorick plug-in, without the need for any additional code.
 *
 * To declare a Property list:
 *   1. declare (in the class declaration, .h file) and define (.C
 *      file) the pair or quadruplet of accessors for your Property
 *      (see Property class documentation;
 *   2. call the #GYOTO_OBJECT macro in in a public section of the
 *      class declaration (in the .h file):
 *      \code
 *      class A: public Object {
 *        public:
 *          GYOTO_OBJECT;
 *          A();
 *          virtual ~A();
 *        ...
 *      };
 *      \endcode
 *   3. call the various GYOTO_PROPERTY_* macros in the corresponding
 *      .C file (see the documentation of the Property class).
 *
 * It is possible to get a Property by name (Assume \a A is a class
 * deriving from Object):
 * \code
 * A myobj();
 * Property const *prop=NULL;
 * prop = myobj.property("PropertyName");
 * if (!prop) throwError("No Property by that name in this object");
 * \endcode
 * It then becomes possible to set or get the Property from or to a
 * Value:
 * \code
 * myobj.set(*prop, size_t(12));
 * size_t val = myobj.get(*prop);
 * \endcode
 * Of course the type of the Value instance and of the Property
 * instance must match. Refer to the documentation of these to classes
 * for details.
 * 
 */
class Gyoto::Object
{
 protected:
  /// The "kind" that is output in the XML entity
  /**
   * E.g. for an Astrobj, fillElement() will ensure
   * \code
   *   <Astrobj kind="kind_" ...>...</Astrobj>
   * \endcode
   * is written.
   */
  std::string kind_;

  /// The plug-ins that needs to be loaded to access this instance's class
  /**
   * E.g. for an Astrobj, fillElement() will ensure
   * \code
   *   <Astrobj ... plugin="plugins_">...</Astrobj>
   * \endcode
   * is written.
   */
  std::vector<std::string> plugins_;

 public:
  /// Whether this class is thread-safe
  /**
   * Return True if this object is thread-safe, i.e. if an instance
   * and its clone can be used in parallel threads (in the context of
   * Scenery::raytrace()). Known objects which are not thread-safe
   * include Lorene metrics and everything from the Python plug-in.
   *
   * The default implementation considers that the class itself is
   * thread safe and recurses into the declared properties to check
   * whether they are safe too. Classes that abide to the
   * Object/Property paradigm and are themselves thread-safe have
   * nothing special to do.
   *
   * Objects that clone children in their copy constructor that are
   * not declared as properties must take these children into
   * account.
   *
   * Classes that are never thread-safe must declare it. It acn be
   * easily done using GYOTO_OBJECT_THREAD_SAFETY in the class
   * declaration and GYOTO_PROPERTY_THREAD_UNSAFE in the class
   * definition.
   */
  virtual bool isThreadSafe() const;


  GYOTO_OBJECT;
  /** \fn virtual Property const * Object::getProperties() const
   *  \brief Get list of properties 
   *
   * This method is declared automatically by the #GYOTO_OBJECT macro
   * and defined automatically by the #GYOTO_PROPERTY_END macro.
   */
  /** \property static Property const  properties[]
   * \brief Property list
   * 
   * This static member is declared automatically by the #GYOTO_OBJECT
   * macro and defined automatically by the #GYOTO_PROPERTY_START,
   * #GYOTO_PROPERTY_END and GYOTO_PROPERTY_* macros.
   *
   * The list of properties is implemented as a static array of
   * Property instances. The last item in a Property of type
   * Property::empty_t, which evaluates to false, so the list can be
   * considered to be NULL-terminated (it is actually rather
   * false-terminated). This empty_t last item can be a link to
   * another Property list: for instance, the last item in
   * Gyoto::Astrobj::Standard::properties is a link to
   * Gyoto::Astrobj::Generic::properties.
   */

  /// Constructor setting kind
  Object (std::string const &kind) ;

  /// Default constructor
  Object () ;

  /// Deep copy constructor
  Object (Object const &orig) ;

  /// Virtual destructor
  virtual ~Object();

  /// Set Value of a Property
  virtual void set(Property const &p, Value val);

  /// Set Value (expressed in unit) of a Property
  virtual void set(Property const &p, Value val, std::string const &unit);

  /// Set Value of a Property
  virtual void set(std::string const &pname, Value val);

  /// Set Value (expressed in unit) of a Property
  virtual void set(std::string const &pname, Value val, std::string const &unit);

  /// Get Value of a Property
  virtual Value get(Property const &p) const;

  /// Get Value of a Property
  virtual Value get(std::string const &pname) const;

  /// Get Value of a Property, converted to unit
  virtual Value get(Property const &p, std::string const &unit) const;

  /// Get Value of a Property, converted to unit
  virtual Value get(std::string const &pname, std::string const &unit) const;

  /// Find property by name
  /**
   * Look into the Property list for a Property whose \a name (or
   * \a name_false, for a boolean Property) is \a pname. Return a const
   * pointer to the first such property found, or NULL if none is
   * found.
   */
  Property const * property(std::string const pname) const;

#ifdef GYOTO_USE_XERCES
  /// Output a single Property to XML
  /**
   * The base implementation decides what to do based on the \a
   * p.type. The format matches how setParameters() an setParameter()
   * would interpret the XML descition.
   *
   * Overriding this method should be avoided, but makes sense in some
   * cases (for instance Screen::fillProperty() selects a different
   * unit for \a Distance based on its magnitude, so that stellar
   * sizes are expressed in solar radii while smaller sizes can be
   * expressed in meters and larger sizes in parsecs).
   *
   * Overriding implementation should fall-back on calling the
   * implementation in the direct parent class:
   * \code
   * class A: public Object {};
   * class B: public A { 
   *  using B::setParameter;
   *  virtual void fillProperty(Gyoto::FactoryMessenger *fmp,
   *     			Property const &p) const ;
   * };
   * void B::fillProperty(Gyoto::FactoryMessenger *fmp,
   *     			Property const &p) const {
   *   if (name=="Duff") fmp->doSomething();
   *   else A::fillProperty(fmp, p);
   * }
   * \endcode
   */
  virtual void fillProperty(Gyoto::FactoryMessenger *fmp,
			    Property const &p) const ;

  /// Fill the XML element for this Object
  /**
   * The base implementation simply calls fillProperty() for each
   * Property defined for the Object.
   *
   * Derived classes should avoid overriding fillElement(). It may
   * make sense occasionally, e.g. to make sure that the metric is
   * output first.
   * 
   * To customize how a given Property is rendered, it is better to
   * override fillProperty().
   *
   * If this method is overridden, the implementation should in
   * general call fillElement() on the direct base.
   */
  virtual void fillElement(Gyoto::FactoryMessenger *fmp) const ;

  
  /// \brief Main loop for parsing Properties from XML description
  /**
   * This function queries the FactoryMessenger for elements to parse,
   * and tries to matche each element to a Property to set it
   * accordingly.
   *
   * Any class that tries to be buildable from XML must supply a
   * subcontractor (for base classes such as Metric, Astrobj, Spectrum
   * and Spectrometer, it is done as a template that must be
   * specialized for each class).
   *
   * This subcontractor typically looks somewhat like this:
\code
SmartPointer<Metric::Generic>
Gyoto::Metric::MyKind::Subcontractor(FactoryMessenger* fmp) {
  SmartPointer<MyKind> gg = new MyKind();
  gg -> setParameters(fmp);
  return gg;
}
\endcode
   *
   * Although this is discouraged, it is possible to override the
   * following functions to customize how XML entities are parsed:
   *   - setParameters() if low-level access to the
   *     FactoryMessenger is required;
   *   - setParameter(std::string name,
   *		      std::string content,
   *		      std::string unit)
   *     to interpret an entity that does not match a Property
   *     (e.g. alternative name);
   *   - setParameter(Gyoto::Property const &p,
   *		      std::string const &name,
   *		      std::string const &content,
   *		      std::string const &unit)
   *     to change how a Property is interpreted.
   */
  virtual void setParameters(Gyoto::FactoryMessenger *fmp) ;
#endif
  /**
   * \brief Set parameter by name
   *
   * This function is used when parsing an XML description, if no
   * Property of this \a name is found. Overriding implementation should
   * fall-back on calling the direct's parent implementation:
   *
   * \code
   * class A: public Object {};
   * class B: public A { 
   *  using B::setParameter;
   *  virtual int setParameter(std::string name,
   *			    std::string content,
   *			    std::string unit);
   * };
   * int B::setParameter(std::string name,
   *			    std::string content,
   *			    std::string unit) {
   *   if (name=="Duff") doSomething(content, unit);
   *   else return A::setParameter(name, content, unit);
   *   return 0;  // name was known
   * }
   * \endcode
   *
   * \param name XML name of the parameter (XML entity). This may have
   *        a path component, e.g. "Astrobj::Radius", in which case a
   *        property named "Astrobj" will be sought in the current
   *        object, and setParameter will be called recusrively on
   *        this Astrobj with Radius as name.
   * \param content string representation of the value
   * \param unit string representation of the unit
   * \return 0 if this parameter is known, 1 if it is not.
   */
  virtual int setParameter(std::string name,
			    std::string content,
			    std::string unit);

  /**
   * \brief Set parameter by Property (and name)
   *
   * This function is used when parsing an XML description, if
   * Property (\a p) of this \a name is found (i.e. either \a p.name or
   * \a p.name_false is equal to \a name). Implementation should fall-back
   * on calling the direct's parent implementation:
   *
   * \code
   * class A: public Object {};
   * class B: public A { 
   *  using B::setParameter;
   *  virtual void setParameter(Gyoto::Property const &p,
   *                            std::string name,
   *			        std::string content,
   *			        std::string unit);
   * };
   * void B::setParameter(Gyoto::Property const &p,
   *  			  std::string name,
   *			  std::string content,
   *			  std::string unit) {
   *   if (name=="Duff") doSomething(content, unit);
   *   else A::setParameter(p, name, content, unit);
   * }
   * \endcode
   *
   * \param p Property that matches \a name (\a p.name == \a name or
   *            \a p.name_false == \a name)
   * \param name XML name of the parameter (XML entity)
   * \param content string representation of the value
   * \param unit string representation of the unit
   */
   virtual void setParameter(Gyoto::Property const &p,
			    std::string const &name,
			    std::string const &content,
			    std::string const &unit);

   /**
    * \brief Format desrciption for a property
    *
    * Returns a string containing the name(s) and type of the
    * property, as well as whether it supports unit.
    */
   std::string describeProperty(Gyoto::Property const &p) const ;

   /**
    * \brief Print (to stdout) some help on this class
    *
    * Describe all properties that this instance supports.
    */
   void help() const ;

protected:
  /**
   * \brief Set kind_
   *
   * kind(const std::string) is protected because, for most Objects,
   * it should not be changed in runtime.
   */
  virtual void kind(const std::string); ///< Set kind_

public:
  virtual std::string kind() const; ///< Get kind_

};

#endif
