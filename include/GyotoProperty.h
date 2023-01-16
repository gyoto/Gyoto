/**
 * \file GyotoProperty.h
 * \brief Introspectable properties 
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


#ifndef __GyotoProperty_H_
#define __GyotoProperty_H_

#include "GyotoConfig.h"
#include <string>
#include <vector>

#ifndef GYOTO_PLUGIN
#define GYOTO_PLUGIN
#endif

namespace Gyoto {
  class Object;
  class Property;
  namespace Metric { class Generic; }
  namespace Astrobj { class Generic; }
  namespace Spectrum { class Generic; }
  namespace Spectrometer { class Generic; }
  class Screen;
  template <class T> class SmartPointer;
}

/// Define the class as not beeing thread-safe
/**
 * See also GYOTO_OBJECT_THREAD_SAFETY
 */
#define GYOTO_PROPERTY_THREAD_UNSAFE(class)		\
  bool class::isThreadSafe() const {return false;}

/// Define a pair of accessors to scalar member (double, long, size_t)
/**
 * Accessors must also be declared in the class declaration, which can
 * be done using #GYOTO_OBJECT_SCALAR_ACCESSORS.
 */
#define GYOTO_PROPERTY_ACCESSORS(class, type, member, method)	\
  void class::method(type v) {member=v;}				\
  type class::method() const {return member;}

/// Define a pair of accessors to scalar member (double, long, size_t)
/**
 * Accessors must also be declared in the class declaration, which can
 * be done using #GYOTO_OBJECT_SCALAR_ACCESSORS.
 *
 * This version allows performing sepcial actions in the accessors, in
 * addition to the usual stuff.
 *
 */
#define GYOTO_PROPERTY_ACCESSORS_SPECIAL(class, type, member, method, set, get) \
  void class::method(type v) {member=v; set }				\
  type class::method() const {get ; return member;}

/// Define 4 accessors to double scalar member in geometrical units
/**
 * Accessors must also be declared in the class declaration, which can
 * be done using #GYOTO_OBJECT_ACCESSORS_UNIT.
 *
 * \param class class name
 * \param member member holding the value in geometrical unit
 * \param method name for accessors
 * \metric member or expression yielding metric (which defines the
 * geometrical unit)
 */
#define GYOTO_PROPERTY_ACCESSORS_GEOMETRICAL(class, member, method, metric) \
  GYOTO_PROPERTY_ACCESSORS(class, double, member, method)		\
  void class::method(double v, std::string const &u) {			\
    member=Units::ToGeometrical(v, u, metric);				\
  }									\
  double class::method(std::string const &u) const {			\
    return Units::FromGeometrical(member, u, metric);			\
  }

/// GYOTO_PROPERTY_ACCESSORS_GEOMETRICAL + GYOTO_PROPERTY_ACCESSORS_SPECIAL
/**
 * Accessors must also be declared in the class declaration, which can
 * be done using #GYOTO_OBJECT_ACCESSORS_UNIT.
 *
 * \param class class name
 * \param member member holding the value in geometrical unit
 * \param method name for accessors
 * \metric member or expression yielding metric (which defines the
 * geometrical unit)
 */
#define GYOTO_PROPERTY_ACCESSORS_GEOMETRICAL_SPECIAL(class, member, method, metric, set, get) \
  GYOTO_PROPERTY_ACCESSORS_SPECIAL(class, double, member, method, set, get)		\
  void class::method(double v, std::string const &u) {			\
    member=Units::ToGeometrical(v, u, metric);				\
  }									\
  double class::method(std::string const &u) const {			\
    return Units::FromGeometrical(member, u, metric);			\
  }

/// Define 4 accessors to double scalar member in specified units
/**
 * Accessors must also be declared in the class declaration, which can
 * be done using #GYOTO_OBJECT_ACCESSORS_UNIT.
 *
 * \param class class name
 * \param member member holding the value in specified unit
 * \param method name for accessors
 * \unit  unit
 */
#define GYOTO_PROPERTY_ACCESSORS_UNIT(class, member, method, unit)      \
  GYOTO_PROPERTY_ACCESSORS(class, double, member, method)		\
  void class::method(double v, std::string const &u) {			\
    method(Units::Converter(u, unit)(v));				\
  }									\
  double class::method(std::string const &u) const {			\
    return Units::Converter(unit, u)(method());				\
  }

/// Start Property list
/**
 * \param class Class for which we are defining a Property list
 * \param doc Documentation for this class. Optional but recommended.
 */
#define GYOTO_PROPERTY_START(...)					\
  GYOTO_PROPERTY_CHOOSER(,##__VA_ARGS__,				\
			 "wrong number of arguments",			\
			 "wrong number of arguments",			\
			 "wrong number of arguments",			\
			 "wrong number of arguments",			\
			 GYOTO_PROPERTY_START_DOC(__VA_ARGS__),		\
			 GYOTO_PROPERTY_START_NODOC(__VA_ARGS__),	\
			 "wrong number of arguments"			\
			 )


/// Define a new Property of type bool
/*
 * Declares a Property named "name". name and namef should not
 * be quoted.
 *
 * \param[in] class Class name
 * \param[in] name Name of property if true;
 * \param[in] namef Name of property if false;
 * \param[in] fname Name of functions for setting or getting the property
 * \param[in] doc Document string (optional but recommended)
 */
#define GYOTO_PROPERTY_BOOL(...)					\
  GYOTO_PROPERTY_CHOOSER(,##__VA_ARGS__,				\
			 "wrong number of arguments",			\
			 GYOTO_PROPERTY_BOOL_DOC(__VA_ARGS__),		\
			 GYOTO_PROPERTY_BOOL_NODOC(__VA_ARGS__),	\
			 "wrong number of arguments",			\
			 "wrong number of arguments",			\
			 "wrong number of arguments",			\
			 "wrong number of arguments"			\
			 )

/// Define a Property of type double
/*
 * Declares a Property named "name". name and namef should not
 * be quoted.
 *
 * \param[in] class Class name
 * \param[in] name Name of property if true;
 * \param[in] fname Name of functions for setting or getting the property
 * \param[in] doc Document string (optional but recommended)
 */
#define GYOTO_PROPERTY_DOUBLE(...)					\
  GYOTO_PROPERTY_CHOOSER(,##__VA_ARGS__,				\
			 GYOTO_NOTHING_6,			\
			 GYOTO_NOTHING_5,			\
			 GYOTO_PROPERTY_DOUBLE_DOC(__VA_ARGS__),	\
			 GYOTO_PROPERTY_DOUBLE_NODOC(__VA_ARGS__),	\
			 GYOTO_NOTHING_2,			\
			 GYOTO_NOTHING_1,			\
			 GYOTO_NOTHING_0			\
			 )


/// Define a Property of type double with unit
/*
 * Declares a Property named "name". name and namef should not
 * be quoted.
 *
 * \param[in] class Class name
 * \param[in] name Name of property if true;
 * \param[in] fname Name of functions for setting or getting the property
 * \param[in] doc Document string (optional but recommended)
 */
#define GYOTO_PROPERTY_DOUBLE_UNIT(...)					\
  GYOTO_PROPERTY_CHOOSER(,##__VA_ARGS__,				\
			 GYOTO_NOTHING_6,			\
			 GYOTO_NOTHING_5,			\
			 GYOTO_PROPERTY_DOUBLE_UNIT_DOC(__VA_ARGS__),	\
			 GYOTO_PROPERTY_DOUBLE_UNIT_NODOC(__VA_ARGS__),	\
			 GYOTO_NOTHING_2,			\
			 GYOTO_NOTHING_1,			\
			 GYOTO_NOTHING_0			\
			 )

/// Define a Property of type vector<double>
/*
 * Declares a Property named "name". name and namef should not
 * be quoted.
 *
 * \param[in] class Class name
 * \param[in] name Name of property if true;
 * \param[in] fname Name of functions for setting or getting the property
 * \param[in] doc Document string (optional but recommended)
 */
#define GYOTO_PROPERTY_VECTOR_DOUBLE(...)				\
  GYOTO_PROPERTY_CHOOSER(,##__VA_ARGS__,				\
			 GYOTO_NOTHING_6,				\
			 GYOTO_NOTHING_5,				\
			 GYOTO_PROPERTY_VECTOR_DOUBLE_DOC(__VA_ARGS__),	\
			 GYOTO_PROPERTY_VECTOR_DOUBLE_NODOC(__VA_ARGS__), \
			 GYOTO_NOTHING_2,				\
			 GYOTO_NOTHING_1,				\
			 GYOTO_NOTHING_0				\
			 )

/// Define a Property of type vector<double> with unit
/*
 * Declares a Property named "name". name and namef should not
 * be quoted.
 *
 * \param[in] class Class name
 * \param[in] name Name of property if true;
 * \param[in] fname Name of functions for setting or getting the property
 * \param[in] doc Document string (optional but recommended)
 */
#define GYOTO_PROPERTY_VECTOR_DOUBLE_UNIT(...)				\
  GYOTO_PROPERTY_CHOOSER(,##__VA_ARGS__,				\
			 GYOTO_NOTHING_6,				\
			 GYOTO_NOTHING_5,				\
			 GYOTO_PROPERTY_VECTOR_DOUBLE_UNIT_DOC(__VA_ARGS__), \
			 GYOTO_PROPERTY_VECTOR_DOUBLE_UNIT_NODOC(__VA_ARGS__), \
			 GYOTO_NOTHING_2,				\
			 GYOTO_NOTHING_1,				\
			 GYOTO_NOTHING_0				\
			 )

/// Define a Property of type String
/*
 * Declares a Property named "name". name and namef should not
 * be quoted.
 *
 * \param[in] class Class name
 * \param[in] name Name of property if true;
 * \param[in] fname Name of functions for setting or getting the property
 * \param[in] doc Document string (optional but recommended)
 */
#define GYOTO_PROPERTY_STRING(...)					\
  GYOTO_PROPERTY_CHOOSER(,##__VA_ARGS__,				\
			 GYOTO_NOTHING_6,				\
			 GYOTO_NOTHING_5,				\
			 GYOTO_PROPERTY_STRING_DOC(__VA_ARGS__),	\
			 GYOTO_PROPERTY_STRING_NODOC(__VA_ARGS__),	\
			 GYOTO_NOTHING_2,				\
			 GYOTO_NOTHING_1,				\
			 GYOTO_NOTHING_0				\
			 )

/// Define a Property of type Filename
/*
 * Declares a Property named "name". name and namef should not
 * be quoted.
 *
 * \param[in] class Class name
 * \param[in] name Name of property if true;
 * \param[in] fname Name of functions for setting or getting the property
 * \param[in] doc Document string (optional but recommended)
 */
#define GYOTO_PROPERTY_FILENAME(...)					\
  GYOTO_PROPERTY_CHOOSER(,##__VA_ARGS__,				\
			 GYOTO_NOTHING_6,				\
			 GYOTO_NOTHING_5,				\
			 GYOTO_PROPERTY_FILENAME_DOC(__VA_ARGS__),	\
			 GYOTO_PROPERTY_FILENAME_NODOC(__VA_ARGS__),	\
			 GYOTO_NOTHING_2,				\
			 GYOTO_NOTHING_1,				\
			 GYOTO_NOTHING_0				\
			 )

/// Define a Property of type long
/*
 * Declares a Property named "name". name and namef should not
 * be quoted.
 *
 * \param[in] class Class name
 * \param[in] name Name of property if true;
 * \param[in] fname Name of functions for setting or getting the property
 * \param[in] doc Document string (optional but recommended)
 */
#define GYOTO_PROPERTY_LONG(...)				\
  GYOTO_PROPERTY_CHOOSER(,##__VA_ARGS__,				\
			 GYOTO_NOTHING_6,				\
			 GYOTO_NOTHING_5,				\
			 GYOTO_PROPERTY_LONG_DOC(__VA_ARGS__),	\
			 GYOTO_PROPERTY_LONG_NODOC(__VA_ARGS__), \
			 GYOTO_NOTHING_2,				\
			 GYOTO_NOTHING_1,				\
			 GYOTO_NOTHING_0				\
			 )

/// Define a Property of type unsigned long
/*
 * Declares a Property named "name". name and namef should not
 * be quoted.
 *
 * \param[in] class Class name
 * \param[in] name Name of property if true;
 * \param[in] fname Name of functions for setting or getting the property
 * \param[in] doc Document string (optional but recommended)
 */
#define GYOTO_PROPERTY_UNSIGNED_LONG(...)				\
  GYOTO_PROPERTY_CHOOSER(,##__VA_ARGS__,				\
			 GYOTO_NOTHING_6,				\
			 GYOTO_NOTHING_5,				\
			 GYOTO_PROPERTY_UNSIGNED_LONG_DOC(__VA_ARGS__),	\
			 GYOTO_PROPERTY_UNSIGNED_LONG_NODOC(__VA_ARGS__), \
			 GYOTO_NOTHING_2,				\
			 GYOTO_NOTHING_1,				\
			 GYOTO_NOTHING_0				\
			 )

/// Define a Property of type vector<unsigned long>
/*
 * Declares a Property named "name". name and namef should not
 * be quoted.
 *
 * \param[in] class Class name
 * \param[in] name Name of property if true;
 * \param[in] fname Name of functions for setting or getting the property
 * \param[in] doc Document string (optional but recommended)
 */
#define GYOTO_PROPERTY_VECTOR_UNSIGNED_LONG(...)			\
  GYOTO_PROPERTY_CHOOSER(,##__VA_ARGS__,				\
			 GYOTO_NOTHING_6,				\
			 GYOTO_NOTHING_5,				\
			 GYOTO_PROPERTY_VECTOR_UNSIGNED_LONG_DOC(__VA_ARGS__), \
			 GYOTO_PROPERTY_VECTOR_UNSIGNED_LONG_NODOC(__VA_ARGS__), \
			 GYOTO_NOTHING_2,				\
			 GYOTO_NOTHING_1,				\
			 GYOTO_NOTHING_0				\
			 )

/// Define a Property of type size_t
/*
 * Declares a Property named "name". name and namef should not
 * be quoted.
 *
 * \param[in] class Class name
 * \param[in] name Name of property if true;
 * \param[in] fname Name of functions for setting or getting the property
 * \param[in] doc Document string (optional but recommended)
 */
#define GYOTO_PROPERTY_SIZE_T(...)				\
  GYOTO_PROPERTY_CHOOSER(,##__VA_ARGS__,				\
			 GYOTO_NOTHING_6,				\
			 GYOTO_NOTHING_5,				\
			 GYOTO_PROPERTY_SIZE_T_DOC(__VA_ARGS__),	\
			 GYOTO_PROPERTY_SIZE_T_NODOC(__VA_ARGS__), \
			 GYOTO_NOTHING_2,				\
			 GYOTO_NOTHING_1,				\
			 GYOTO_NOTHING_0				\
			 )

/// Define a Property of type Gyoto::Metric::Generic
/*
 * Declares a Property named "name". name and namef should not
 * be quoted.
 *
 * \param[in] class Class name
 * \param[in] name Name of property if true;
 * \param[in] fname Name of functions for setting or getting the property
 * \param[in] doc Document string (optional but recommended)
 */
#define GYOTO_PROPERTY_METRIC(...)				\
  GYOTO_PROPERTY_CHOOSER(,##__VA_ARGS__,				\
			 GYOTO_NOTHING_6,				\
			 GYOTO_NOTHING_5,				\
			 GYOTO_PROPERTY_METRIC_DOC(__VA_ARGS__),	\
			 GYOTO_PROPERTY_METRIC_NODOC(__VA_ARGS__), \
			 GYOTO_NOTHING_2,				\
			 GYOTO_NOTHING_1,				\
			 GYOTO_NOTHING_0				\
			 )


/// Define a Property of type Gyoto::Spectrum::Generic
/*
 * Declares a Property named "name". name and namef should not
 * be quoted.
 *
 * \param[in] class Class name
 * \param[in] name Name of property if true;
 * \param[in] fname Name of functions for setting or getting the property
 * \param[in] doc Document string (optional but recommended)
 */
#define GYOTO_PROPERTY_SPECTRUM(...)				\
  GYOTO_PROPERTY_CHOOSER(,##__VA_ARGS__,				\
			 GYOTO_NOTHING_6,				\
			 GYOTO_NOTHING_5,				\
			 GYOTO_PROPERTY_SPECTRUM_DOC(__VA_ARGS__),	\
			 GYOTO_PROPERTY_SPECTRUM_NODOC(__VA_ARGS__), \
			 GYOTO_NOTHING_2,				\
			 GYOTO_NOTHING_1,				\
			 GYOTO_NOTHING_0				\
			 )


/// Define a Property of type Gyoto::Astrobj::Generic
/*
 * Declares a Property named "name". name and namef should not
 * be quoted.
 *
 * \param[in] class Class name
 * \param[in] name Name of property if true;
 * \param[in] fname Name of functions for setting or getting the property
 * \param[in] doc Document string (optional but recommended)
 */
#define GYOTO_PROPERTY_ASTROBJ(...)				\
  GYOTO_PROPERTY_CHOOSER(,##__VA_ARGS__,				\
			 GYOTO_NOTHING_6,				\
			 GYOTO_NOTHING_5,				\
			 GYOTO_PROPERTY_ASTROBJ_DOC(__VA_ARGS__),	\
			 GYOTO_PROPERTY_ASTROBJ_NODOC(__VA_ARGS__), \
			 GYOTO_NOTHING_2,				\
			 GYOTO_NOTHING_1,				\
			 GYOTO_NOTHING_0				\
			 )


/// Define a Property of type Gyoto::Screen
/*
 * Declares a Property named "name". name and namef should not
 * be quoted.
 *
 * \param[in] class Class name
 * \param[in] name Name of property if true;
 * \param[in] fname Name of functions for setting or getting the property
 * \param[in] doc Document string (optional but recommended)
 */
#define GYOTO_PROPERTY_SCREEN(...)				\
  GYOTO_PROPERTY_CHOOSER(,##__VA_ARGS__,				\
			 GYOTO_NOTHING_6,				\
			 GYOTO_NOTHING_5,				\
			 GYOTO_PROPERTY_SCREEN_DOC(__VA_ARGS__),	\
			 GYOTO_PROPERTY_SCREEN_NODOC(__VA_ARGS__), \
			 GYOTO_NOTHING_2,				\
			 GYOTO_NOTHING_1,				\
			 GYOTO_NOTHING_0				\
			 )


/// Define a Property of type Gyoto::Spectrometer::Generic
/*
 * Declares a Property named "name". name and namef should not
 * be quoted.
 *
 * \param[in] class Class name
 * \param[in] name Name of property if true;
 * \param[in] fname Name of functions for setting or getting the property
 * \param[in] doc Document string (optional but recommended)
 */
#define GYOTO_PROPERTY_SPECTROMETER(...)				\
  GYOTO_PROPERTY_CHOOSER(,##__VA_ARGS__,				\
			 GYOTO_NOTHING_6,				\
			 GYOTO_NOTHING_5,				\
			 GYOTO_PROPERTY_SPECTROMETER_DOC(__VA_ARGS__),	\
			 GYOTO_PROPERTY_SPECTROMETER_NODOC(__VA_ARGS__), \
			 GYOTO_NOTHING_2,				\
			 GYOTO_NOTHING_1,				\
			 GYOTO_NOTHING_0				\
			 )

/// Define class::properties and class::getProperties() 
#define GYOTO_PROPERTY_END(class, next)					\
  Property(next)};							\
  Gyoto::Property const * class::getProperties() const {		\
    return class::properties;						\
  }									\
  const std::string class::builtinPluginValue ( GYOTO_STRINGIFY(GYOTO_PLUGIN) ); \
  std::vector<std::string> class::plugins() const {			\
    if (plugins_.size() == 0) {						\
      std::vector<std::string> p;					\
      p.push_back(class::builtinPluginValue);				\
      return p;					\
    }									\
    return plugins_;							\
  }									\
  void  class::plugins(std::vector<std::string> const & plugname) {	\
    plugins_=plugname;							\
  }

/// Property that can be set and got using standard methods
/**
 * The Property API makes it easy to declare the parameters that can
 * be set in a class.
 *
 * Developpers who simply write classes (deriving from
 * Astrobj::Generic, , Metric::Generic, Spectrum::Generic) need not
 * know the inners of the Property class and interact with it only
 * using macros to declare the parameters they need to read from XML.
 *
 * To make use of the Property framework, a class
 * must derive from Gyoto::Object and use the #GYOTO_OBJECT in a
 * public section of the class declaration (i.e. in the .h
 * file). Then, in the corresponding .C file, the GYOTO_PROPERTY_*
 * macros are used as follows (note the absence of punctuation after
 * the macros):
 * \code
 * GYOTO_PROPERTY_START(MyClass)
 * GYOTO_PROPERTY_<type>(MyClass, PropertyName, accessor)
 * ...
 * GYOTO_PROPERTY_END(MyClass, ParentClass::properties)
 * \endcode
 *
 * In the above, #GYOTO_PROPERTY_START starts the definition of the
 * static member MyClass::properties. Each GYOTO_PROPERTY_<type> macro
 * declares a new property. #GYOTO_PROPERTY_END ends the definition of
 * the property list, with an optional pointer to the parent's class
 * Property list, and defines the MyClass::getProperties() method.
 *
 * The underlying accessors must always be defined, both to set and to
 * get the property. For the sake of simplicity, only a limited number
 * of data types are allowed:
 *   - double: see #GYOTO_PROPERTY_DOUBLE, #GYOTO_PROPERTY_DOUBLE_UNIT;
 *   - long: see #GYOTO_PROPERTY_LONG;
 *   - unsigned long: see #GYOTO_PROPERTY_UNSIGNED_LONG
 *     (a.k.a. size_t: see #GYOTO_PROPERTY_SIZE_T, this may break on
 *     architectures where size_t is not the same as unsigned long);
 *   - bool: see #GYOTO_PROPERTY_BOOL;
 *   - std::vector<double>: see #GYOTO_PROPERTY_VECTOR_DOUBLE
 *     and #GYOTO_PROPERTY_VECTOR_DOUBLE_UNIT;
 *   - std::vector<unsigned long>:
 *     see #GYOTO_PROPERTY_VECTOR_UNSIGNED_LONG;
 *   - Gyoto::SmartPointers to various base classes: Screen,
 *     Metric::Generic, Astrobj::Generic, Spectrum::Generic and
 *     Spectrometer::Generic. See #GYOTO_PROPERTY_METRIC,
 *     #GYOTO_PROPERTY_SCREEN, #GYOTO_PROPERTY_ASTROBJ,
 *     #GYOTO_PROPERTY_SPECTRUM and #GYOTO_PROPERTY_SPECTROMETER.
 *
 * For the floating point data-types (double and vector<double>), two
 * additional accessors supporting units can be provided. The
 * accessors must have the same name and have specific prototypes, see
 * the various function pointer typedefs, e.g. #set_double_t and
 * #get_double_t.
 *
 * The type used in these accessors may not be the same as the type of
 * the underlying class member. For instance, to read an array, it was
 * chosen to use the std::vector<type> type because it is easy to
 * read such a vector from XML and to thus determine dynamically the
 * number of elements provided. But this type is slow, so it is
 * expected that the class member will rather be a C-style array
 * (double arr[]) or something else entirely. It is not forbidden to
 * have a set of high-level accessors for the Property interface on
 * top of lower-level, more efficient accessors to be used in
 * compiled, static code:
 *
 * \code
 * void MyClass::arrayMember(double const * const tab) {
 *   for (int i=0; i<5; ++i) member_[i]=tab[i];
 * }
 * void MyClass::arrayMemberVector(std::vector<double> const &vect) {
 *   if (vect.size()!=5) throwError("Please provide 5 elements");
 *   for (int i=0; i<5; ++i) member_[i]=vect[i];
 * } 
 * double const *  MyClass::arrayMember() const {return member_;}
 * std::vector<double> MyClass::arrayMemberVector() const {
 *   std::vector<double> v(5, 0.);
 *   for (int i=0; i<5; ++i) v[i]=member_[i];
 *   return v;
 * }
 * \endcode
 *
 * In this example, assuming MyClass is based directly on Object and
 * member_ is the only parameter to read from XML, the Property list
 * may be defined as:
 * \code
 * GYOTO_PROPERTY_START(MyClass)
 * GYOTO_PROPERTY_VECTOR_DOUBLE(MyClass, ArrayMember, arrayMemberVector)
 * GYOTO_PROPERTY_END(MyClass, Object::properties)
 * \endcode
 *
 * Again, nothing more is required to read and write ArrayMember from
 * XML and from Yorick.
 */
class Gyoto::Property
{
 private:

 public:
  /// Possible type of a Property instance
  /**
   * These are all the values that Property::type may take.
   */
  enum type_e {
    /// Type is double
    double_t,
    /// Type is long
    long_t,
    /// Type is unsigned long (a.k.a. size_t)
    unsigned_long_t,
    /// Type is size_t (only if distinct from unsigned long)
    size_t_t,
    /// Type is bool
    bool_t,
    /// Type is std::string
    string_t,
    /// Type is std::string and holds a file name
    /**
     * In this case, the value may be changed to include relative or
     * absolute path. The prefix "`pwd`/" forces path to be
     * interpreted relative to the current working directory. If the
     * value starts with "/", it is interpreted as an absolute
     * path. Otherwise, when reading from XML, the path is interpreted
     * relative to the XML file.
     */
    filename_t,
    /// Type is std::vector<double>
    vector_double_t,
    /// Type is std::vector<unsigned long>
    vector_unsigned_long_t,
    /// Type is Gyoto::SmartPointer<Gyoto::Metric::Generic>
    metric_t,
    /// Type is Gyoto::SmartPointer<Gyoto::Screen::Generic>
    screen_t,
    /// Type is Gyoto::SmartPointer<Gyoto::Astrobj::Generic>
    astrobj_t,
    /// Type is Gyoto::SmartPointer<Gyoto::Spectrum::Generic>
    spectrum_t,
    /// Type is Gyoto::SmartPointer<Gyoto::Spectrometer::Generic>
    spectrometer_t,
    /// Property is empty
    /**
     * In this case (and only in this case):
     *  1. #Gyoto::Property::operator bool() const return false;
     *  2. Property::parent may be different from NULL and point to
     *     another array of Property instances, that should be
     *     interpreted as appended to this list.
     */
    empty_t};
  /// Name of this instance
  std::string name;
  /// Name if false
  /**
   * Only if #type is #empty_t
   */
  std::string name_false;
  /// Type of this instance
  int type;
  /// Prototype for an accessor to set a double
  typedef void (Object::* set_double_t)(double val);
  /// Prototype for an accessor to get a double
  typedef double (Object::* get_double_t)() const;
  /// Prototype for an accessor to set a double, with unit
  typedef void (Object::* set_double_unit_t)(double val,
					     std::string const &unit);
  /// Prototype for an accessor to get a double, with unit
  typedef double (Object::* get_double_unit_t)(std::string const &unit) const;
  /// Prototype for an accessor to set a long
  typedef void (Object::* set_long_t)(long val);
  /// Prototype for an accessor to get a long
  typedef long (Object::* get_long_t)() const;
  /// Prototype for an accessor to set an unsigned long
  typedef void (Object::* set_unsigned_long_t)(unsigned long val);
  /// Prototype for an accessor to get an unsigned long
  typedef unsigned long (Object::* get_unsigned_long_t)() const;
  /// Prototype for an accessor to set a size_t
  typedef void (Object::* set_size_t_t)(size_t val);
  /// Prototype for an accessor to get a size_t
  typedef size_t (Object::* get_size_t_t)() const;
  /// Prototype for an accessor to set a bool
  typedef void (Object::* set_bool_t)(bool val);
  /// Prototype for an accessor to get a bool
  typedef bool (Object::* get_bool_t)() const;
  /// Prototype for an accessor to set a string
  typedef void (Object::* set_string_t)(std::string const&);
  /// Prototype for an accessor to get a string
  typedef std::string (Object::* get_string_t)() const;
  /// Prototype for an accessor to set a filename
  typedef void (Object::* set_fname_t)(std::string const&);
  /// Prototype for an accessor to get a filename
  typedef std::string (Object::* get_fname_t)() const;
  /// Prototype for an accessor to set a std::vector<double>
  typedef void (Object::* set_vector_double_t)(std::vector<double> const&);
  /// Prototype for an accessor to get a std::vector<double>
  typedef std::vector<double> (Object::* get_vector_double_t)() const;
  /// Prototype for an accessor to set a std::vector<double>, with unit
  typedef void (Object::* set_vector_double_unit_t)(std::vector<double> const&, std::string const &);
  /// Prototype for an accessor to get a std::vector<double>, with unit
  typedef std::vector<double> (Object::* get_vector_double_unit_t)(std::string const &) const;
  /// Prototype for an accessor to set a std::vector<unsigned long>
  typedef void (Object::* set_vector_unsigned_long_t)(std::vector<unsigned long> const&);
  /// Prototype for an accessor to get a std::vector<unsigned long>
  typedef std::vector<unsigned long> (Object::* get_vector_unsigned_long_t)() const;

  /// Prototype for an accessor to set a Gyoto::SmartPointer<Gyoto::Metric::Generic>
  typedef void (Object::* set_metric_t)
    (Gyoto::SmartPointer<Gyoto::Metric::Generic>);
  /// Prototype for an accessor to get a Gyoto::SmartPointer<Gyoto::Metric::Generic>
  typedef Gyoto::SmartPointer<Gyoto::Metric::Generic>
    (Object::* get_metric_t)() const;

  /// Prototype for an accessor to set a Gyoto::SmartPointer<Gyoto::Screen>
  typedef void (Object::* set_screen_t)
    (Gyoto::SmartPointer<Gyoto::Screen>);
  /// Prototype for an accessor to get a Gyoto::SmartPointer<Gyoto::Screen>
  typedef Gyoto::SmartPointer<Gyoto::Screen>
    (Object::* get_screen_t)() const;

  /// Prototype for an accessor to set a Gyoto::SmartPointer<Gyoto::Astrobj::Generic>
  typedef void (Object::* set_astrobj_t)
    (Gyoto::SmartPointer<Gyoto::Astrobj::Generic>);
  /// Prototype for an accessor to get a Gyoto::SmartPointer<Gyoto::Astrobj::Generic>
  typedef Gyoto::SmartPointer<Gyoto::Astrobj::Generic>
    (Object::* get_astrobj_t)() const;

  /// Prototype for an accessor to set a Gyoto::SmartPointer<Gyoto::Spectrum::Generic>
  typedef void (Object::* set_spectrum_t)
    (Gyoto::SmartPointer<Gyoto::Spectrum::Generic>);
  /// Prototype for an accessor to get a Gyoto::SmartPointer<Gyoto::Spectrum::Generic>
  typedef Gyoto::SmartPointer<Gyoto::Spectrum::Generic>
    (Object::* get_spectrum_t)() const;

  /// Prototype for an accessor to set a Gyoto::SmartPointer<Gyoto::Spectrometer::Generic>
  typedef void (Object::* set_spectrometer_t)
    (Gyoto::SmartPointer<Gyoto::Spectrometer::Generic>);
  /// Prototype for an accessor to get a Gyoto::SmartPointer<Gyoto::Spectrometer::Generic>
  typedef Gyoto::SmartPointer<Gyoto::Spectrometer::Generic>
    (Object::* get_spectrometer_t)() const;

  /// Union holding an accessor to set any type
  /**
   * Right type is stored in #type
   */
  union setter_t {
    set_double_t set_double;
    set_long_t set_long;
    set_unsigned_long_t set_unsigned_long;
    set_size_t_t set_size_t;
    set_bool_t set_bool;
    set_string_t set_string;
    set_vector_double_t set_vdouble;
    set_vector_unsigned_long_t set_vulong;
    set_metric_t set_metric;
    set_screen_t set_screen;
    set_astrobj_t set_astrobj;
    set_spectrum_t set_spectrum;
    set_spectrometer_t set_spectrometer;
  };
  /// Union holding an accessor to get any type
  /**
   * Right type is stored in #type
   */
  union getter_t {
    get_double_t get_double;
    get_long_t get_long;
    get_unsigned_long_t get_unsigned_long;
    get_size_t_t get_size_t;
    get_bool_t get_bool;
    get_string_t get_string;
    get_vector_double_t get_vdouble;
    get_vector_unsigned_long_t get_vulong;
    get_metric_t get_metric;
    get_screen_t get_screen;
    get_astrobj_t get_astrobj;
    get_spectrum_t get_spectrum;
    get_spectrometer_t get_spectrometer;
  };
  /// Union holding an accessor to set double or vector<double> with unit
  /**
   * Right type is stored in #type
   */
  union setter_unit_t {
    set_double_unit_t set_double;
    set_vector_double_unit_t set_vdouble;
  };
  /// Union holding an accessor to get double or vector<double> with unit
  union getter_unit_t {
    get_double_unit_t get_double;
    get_vector_double_unit_t get_vdouble;
  };

  /// Pointer to the setter method
  /**
   * Right type is stored in #type
   */
  setter_t setter;
  /// Pointer to the getter method
  /**
   * Right type is stored in #type
   */
  getter_t getter;
  /// Pointer to the setter (with unit) method
  /**
   * Right type is stored in #type
   */
  setter_unit_t setter_unit;
  /// Pointer to the getter (with unit) method
  /**
   * Right type is stored in #type
   */
  getter_unit_t getter_unit;

  std::string doc;

  /// True if #Gyoto::Property::type is not #Gyoto::Property::empty_t
  operator bool() const ;

  /// If #type is #empty_t, link to another Property list
  Property const * const  parent;
  
  /// Default constructor
  Property();

  /// Constructor for #type==#empty_t
  Property(Property const * const ancestor);

  /// Constructor for class name pseudo-property
  Property(std::string classname, std::string doc="");

  /// Constructor for #type==#long_t
  Property(std::string name,
	   set_long_t set_long,
	   get_long_t get_long,
	   std::string doc);

  /// Constructor for #type==#unsigned_long_t
  Property(std::string name,
	   set_unsigned_long_t set_unsigned_long,
	   get_unsigned_long_t get_unsigned_long,
	   std::string doc);

  /// Constructor for #type==#size_t_t
  /**
   * The dummy int parameter is only there to differenciate from the
   * unsigned long constructor on platforms where size_t is a typdef
   * to unsigned long.
   */
  Property(std::string name,
	   set_size_t_t set_size_t,
	   get_size_t_t get_size_t,
	   int dummy,
	   std::string doc);

  /// Constructor for #type==#double_t, without unit support
  Property(std::string name,
	   set_double_t set_double,
	   get_double_t get_double,
	   std::string doc);

  /// Constructor for #type==#double_t, with unit support
  Property(std::string name,
	   set_double_t set_double,
	   get_double_t get_double,
	   set_double_unit_t set_double_unit,
	   get_double_unit_t get_double_unit,
	   std::string doc);

  /// Constructor for #type==#bool_t
  Property(std::string name,
	   std::string name_false,
	   set_bool_t set_bool,
	   get_bool_t get_bool,
	   std::string doc);

  /// Constructor for #type==#string_t or #filename_t
  /**
   * \param name name of the Property
   * \param set_string pointer to the setter accessor
   * \param get_string pointer to the getter accessor
   * \param is_filename true means #type=#filename_t
   */
  Property(std::string name,
	   set_string_t set_string,
	   get_string_t get_string,
	   bool is_filename,
	   std::string doc);

  /// Constructor for #type==#vector_double_t, without unit support
  Property(std::string name,
	   set_vector_double_t set_vdouble,
	   get_vector_double_t get_vdouble,
	   std::string doc);

  /// Constructor for #type==#vector_double_t, with unit support
  Property(std::string name,
	   set_vector_double_t set_vdouble,
	   get_vector_double_t get_vdouble,
	   set_vector_double_unit_t set_vdouble_unit,
	   get_vector_double_unit_t get_vdouble_unit,
	   std::string doc);

  /// Constructor for #type==#vector_unsigned_long_t
  Property(std::string name,
	   set_vector_unsigned_long_t set_vulong,
	   get_vector_unsigned_long_t get_vulong,
	   std::string doc);

  /// Constructor for #type==#metric_t
  Property(std::string name,
	   set_metric_t set_metric,
	   get_metric_t get_metric,
	   std::string doc);

  /// Constructor for #type==#screen_t
  Property(std::string name,
	   set_screen_t set_screen,
	   get_screen_t get_screen,
	   std::string doc);

  /// Constructor for #type==#astrobj_t
  Property(std::string name,
	   set_astrobj_t set_astrobj,
	   get_astrobj_t get_astrobj,
	   std::string doc);

  /// Constructor for #type==#spectrum_t
  Property(std::string name,
	   set_spectrum_t set_spectrum,
	   get_spectrum_t get_spectrum,
	   std::string doc);

  /// Constructor for #type==#spectrometer_t
  Property(std::string name,
	   set_spectrometer_t set_spectrometer,
	   get_spectrometer_t get_spectrometer,
	   std::string doc);

  /// Constructor setting only #name and #type
  /**
   * Used in the Python plug-in to provide peudo-properties
   */
  Property(std::string name, int type);

  /// Get Property::type_e value from name
  static type_e typeFromString(std::string stype);

};

/// \cond INTERNAL
#define GYOTO_PROPERTY_CHOOSER(x, A, B, C, D, E, F, FUNC, ...) FUNC

#define GYOTO_PROPERTY_START_DOC(class, doc)	\
  Property const class::properties[] = {	\
    Property (#class, doc),

#define GYOTO_PROPERTY_START_NODOC(class)		\
  GYOTO_PROPERTY_START_DOC(class, "")


#define GYOTO_PROPERTY_BOOL_DOC(class, name, namef, fname, doc)		\
  Gyoto::Property							\
  (#name,								\
   #namef,								\
   (Gyoto::Property::set_bool_t)&class :: fname,			\
   (Gyoto::Property::get_bool_t)&class :: fname,			\
   doc),
#define GYOTO_PROPERTY_BOOL_NODOC(class, name, namef, fname)	\
  GYOTO_PROPERTY_BOOL_DOC(class, name, namef, fname, "")

#define GYOTO_PROPERTY_DOUBLE_DOC(class, name, fname, doc)		\
  Gyoto::Property						\
  (#name,							\
   (Gyoto::Property::set_double_t)&class::fname,		\
   (Gyoto::Property::get_double_t)&class::fname,		\
   doc),
#define GYOTO_PROPERTY_DOUBLE_NODOC(class, name, fname)	\
  GYOTO_PROPERTY_DOUBLE_DOC(class, name, fname, "")

#define GYOTO_PROPERTY_LONG_DOC(class, name, fname, doc)	\
  Gyoto::Property						\
  (#name,							\
   (Gyoto::Property::set_long_t)&class::fname,			\
   (Gyoto::Property::get_long_t)&class::fname,			\
   doc),
#define GYOTO_PROPERTY_LONG_NODOC(class, name, fname)	\
  GYOTO_PROPERTY_LONG_DOC(class, name, fname, "")

#define GYOTO_PROPERTY_UNSIGNED_LONG_DOC(class, name, fname, doc)	\
  Gyoto::Property							\
  (#name,								\
   (Gyoto::Property::set_unsigned_long_t)&class::fname,			\
   (Gyoto::Property::get_unsigned_long_t)&class::fname,			\
   doc),
#define GYOTO_PROPERTY_UNSIGNED_LONG_NODOC(class, name, fname)	\
  GYOTO_PROPERTY_UNSIGNED_LONG_DOC(class, name, fname, "")

#define GYOTO_PROPERTY_SIZE_T_DOC(class, name, fname, doc)	\
  Gyoto::Property						\
  (#name,							\
   (Gyoto::Property::set_size_t_t)&class::fname,		\
   (Gyoto::Property::get_size_t_t)&class::fname,		\
   1,								\
   doc),
#define GYOTO_PROPERTY_SIZE_T_NODOC(class, name, fname)	\
  GYOTO_PROPERTY_SIZE_T_DOC(class, name, fname, "")

#define GYOTO_PROPERTY_DOUBLE_UNIT_DOC(class, name, fname, doc)	\
  Gyoto::Property						\
  (#name,							\
   (Gyoto::Property::set_double_t)&class::fname,		\
   (Gyoto::Property::get_double_t)&class::fname,		\
   (Gyoto::Property::set_double_unit_t)&class::fname,		\
   (Gyoto::Property::get_double_unit_t)&class::fname,		\
   doc),
#define GYOTO_PROPERTY_DOUBLE_UNIT_NODOC(class, name, fname)	\
  GYOTO_PROPERTY_DOUBLE_UNIT_DOC(class, name, fname, "")

#define GYOTO_PROPERTY_FILENAME_DOC(class, name, fname, doc)	\
  Gyoto::Property						\
  (#name,							\
   (Gyoto::Property::set_string_t)&class::fname,		\
   (Gyoto::Property::get_string_t)&class::fname,		\
   true, doc),
#define GYOTO_PROPERTY_FILENAME_NODOC(class, name, fname)	\
  GYOTO_PROPERTY_FILENAME_DOC(class, name, fname, "")

#define GYOTO_PROPERTY_STRING_DOC(class, name, fname, doc)		\
  Gyoto::Property							\
  (#name,								\
   (Gyoto::Property::set_string_t)&class::fname,			\
   (Gyoto::Property::get_string_t)&class::fname,			\
   false, doc),
#define GYOTO_PROPERTY_STRING_NODOC(class, name, fname)	\
  GYOTO_PROPERTY_STRING_DOC(class, name, fname, "")

#define GYOTO_PROPERTY_VECTOR_DOUBLE_DOC(class, name, fname, doc)	\
  Gyoto::Property							\
  (#name,								\
   (Gyoto::Property::set_vector_double_t)&class::fname,			\
   (Gyoto::Property::get_vector_double_t)&class::fname,			\
   doc),
#define GYOTO_PROPERTY_VECTOR_DOUBLE_NODOC(class, name, fname)	\
  GYOTO_PROPERTY_VECTOR_DOUBLE_DOC(class, name, fname, "")

#define GYOTO_PROPERTY_VECTOR_DOUBLE_UNIT_DOC(class, name, fname, doc)	\
  Gyoto::Property							\
  (#name,								\
   (Gyoto::Property::set_vector_double_t)&class::fname,			\
   (Gyoto::Property::get_vector_double_t)&class::fname,			\
   (Gyoto::Property::set_vector_double_unit_t)&class::fname,		\
   (Gyoto::Property::get_vector_double_unit_t)&class::fname,		\
   doc),
#define GYOTO_PROPERTY_VECTOR_DOUBLE_UNIT_NODOC(class, name, fname)	\
  GYOTO_PROPERTY_VECTOR_DOUBLE_UNIT_DOC(class, name, fname, "")

#define GYOTO_PROPERTY_VECTOR_UNSIGNED_LONG_DOC(class, name, fname, doc) \
  Gyoto::Property							\
  (#name,								\
   (Gyoto::Property::set_vector_unsigned_long_t)&class::fname,		\
   (Gyoto::Property::get_vector_unsigned_long_t)&class::fname,		\
   doc),
#define GYOTO_PROPERTY_VECTOR_UNSIGNED_LONG_NODOC(class, name, fname)	\
  GYOTO_PROPERTY_VECTOR_UNSIGNED_LONG_DOC(class, name, fname, "")

#define GYOTO_PROPERTY_METRIC_DOC(class, name, fname, doc)		\
  Gyoto::Property							\
  (#name,								\
   (Gyoto::Property::set_metric_t)&class::fname,			\
   (Gyoto::Property::get_metric_t)&class::fname,			\
   doc),
#define GYOTO_PROPERTY_METRIC_NODOC(class, name, fname)	\
  GYOTO_PROPERTY_METRIC_DOC(class, name, fname, "")

#define GYOTO_PROPERTY_SCREEN_DOC(class, name, fname, doc)		\
  Gyoto::Property							\
  (#name,								\
   (Gyoto::Property::set_screen_t)&class::fname,			\
   (Gyoto::Property::get_screen_t)&class::fname,			\
   doc),
#define GYOTO_PROPERTY_SCREEN_NODOC(class, name, fname)	\
  GYOTO_PROPERTY_SCREEN_DOC(class, name, fname, "")

#define GYOTO_PROPERTY_ASTROBJ_DOC(class, name, fname, doc)		\
  Gyoto::Property							\
  (#name,								\
   (Gyoto::Property::set_astrobj_t)&class::fname,			\
   (Gyoto::Property::get_astrobj_t)&class::fname,			\
   doc),
#define GYOTO_PROPERTY_ASTROBJ_NODOC(class, name, fname)	\
  GYOTO_PROPERTY_ASTROBJ_DOC(class, name, fname, "")

#define GYOTO_PROPERTY_SPECTRUM_DOC(class, name, fname, doc)		\
  Gyoto::Property							\
    (#name,								\
     (Gyoto::Property::set_spectrum_t)&class::fname,			\
     (Gyoto::Property::get_spectrum_t)&class::fname,			\
     doc),
#define GYOTO_PROPERTY_SPECTRUM_NODOC(class, name, fname)	\
  GYOTO_PROPERTY_SPECTRUM_DOC(class, name, fname, "")

#define GYOTO_PROPERTY_SPECTROMETER_DOC(class, name, fname, doc)	\
  Gyoto::Property							\
  (#name,								\
   (Gyoto::Property::set_spectrometer_t)&class::fname,			\
   (Gyoto::Property::get_spectrometer_t)&class::fname,			\
   doc),
#define GYOTO_PROPERTY_SPECTROMETER_NODOC(class, name, fname)	\
  GYOTO_PROPERTY_SPECTROMETER_DOC(class, name, fname, "")
// \endcond INTERNAL

#endif
