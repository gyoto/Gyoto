/**
 * \file GyotoValue.h
 * \brief Introspectable value
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


#ifndef __GyotoValue_H_
#define __GyotoValue_H_

#include "GyotoConfig.h"
#include "GyotoSmartPointer.h"
#include <string>
#include <vector>

namespace Gyoto {
  class Value;
  namespace Metric {class Generic;}
  namespace Astrobj {class Generic;}
  namespace Spectrum {class Generic;}
  namespace Spectrometer {class Generic;}
  class Screen;
}
/// Container for the value of a Property
/**
 * The Value class is very similar to the C union type (although not
 * as memory efficient): it can hold several type of values, but only
 * one at a time. Care must be taken to ensure only the member that
 * was set is retrieved. The purpose of the Value class is to be used
 * together with the Property class: code determines dynamicaly the
 * type of a Property, reads the corresponding value appropriateley
 * (e.g. from XML or from the Yorick prompt), stores the value in a
 * Value instance, and sets the Property using the Object::set()
 * method. Likewise, the Object::get() method returns a
 * Gyoto::Value. Property::type must be used to determine which member
 * of the Value is meaningful.
 *
 * Casting between Value and the various data type it can hold is
 * normally automatic, but the members can also be accessed
 * explicitly make code more easy to read and less ambiguous.
 */
class Gyoto::Value {
 public:
  /// Constructor
  Value();

  /// Destructor
  ~Value();

  /// Assignement operator
  Value& operator=(Value const&);

  /// Type of this instance
  int type;

 private:
  /// A double value
  double Double;
 public:
  /// Construct/cast from double
  Value(double);
  /// Cast to double
  operator double() const;

 private:
  /// A boolean value
  bool Bool;
 public:
  /// Construct/cast from boolean
  Value(bool);
  /// Cast to bool
  operator bool() const;

 private:
  /// A long value
  long Long;
 public:
  /// Construct/cast from long
  Value(long);
  /// Cast to long
  operator long() const;

 private:
  /// An unsigned long (a.k.a. size_t)
  unsigned long ULong;
 public:
  /// Construct/cast from unsigned long
  Value(unsigned long);
  /// Cast to unsigned long
  operator unsigned long() const;

 private:
  /// A size_t (only if distinct from unsigned long)
  size_t SizeT;
 public:
#if !defined(GYOTO_SIZE__T_IS_UNSIGNED_LONG)
  /// Construct/cast from unsigned long
  Value(size_t);
  /// Cast to unsigned long
  operator size_t() const;
#endif

 private:
  /// A string value
  std::string String;
 public:
  /// Construct/cast from string
  Value(std::string);
  /// Cast to string
  operator std::string() const;

 private:
  /// A vector of double values
  std::vector<double> VDouble;
 public:
  /// Construct/cast from vector of doubles
  Value(std::vector<double>);
  /// Cast to vector of doubles
  operator std::vector<double>() const;

 private:
  /// A vector of unsigned long values
  std::vector<unsigned long> VULong;
 public:
  /// Construct/cast from vector of unsigned long values
  Value(std::vector<unsigned long>);
  /// Cast to vector of unsigned long values
  operator std::vector<unsigned long>() const;

 private:
  /// A Metric object
  Gyoto::SmartPointer<Gyoto::Metric::Generic> Metric;
 public:
 /// Cast from Metric object
  Value(Gyoto::SmartPointer<Gyoto::Metric::Generic>);
  /// Cast to Metric object
  operator Gyoto::SmartPointer<Gyoto::Metric::Generic>() const;

 private:
  /// An Astrobj Object
  Gyoto::SmartPointer<Gyoto::Astrobj::Generic> Astrobj; 
 public:
  /// Cast from Astrobj
  Value(Gyoto::SmartPointer<Gyoto::Astrobj::Generic>);
  /// Cast to Astrobj
  operator Gyoto::SmartPointer<Gyoto::Astrobj::Generic>() const;

 private:
  /// A Spectrum object
  Gyoto::SmartPointer<Gyoto::Spectrum::Generic> Spectrum;
 public:
  /// Cast from Spectrum
  Value(Gyoto::SmartPointer<Gyoto::Spectrum::Generic>);
  /// Cast to Spectrum
  operator Gyoto::SmartPointer<Gyoto::Spectrum::Generic>() const;

 private:
  /// A Spectrometer object
  Gyoto::SmartPointer<Gyoto::Spectrometer::Generic> Spectrometer;
 public:
  /// Cast from Spectrometer
  Value(Gyoto::SmartPointer<Gyoto::Spectrometer::Generic>);
  /// Cast to Spectrometer
  operator Gyoto::SmartPointer<Gyoto::Spectrometer::Generic>() const;

 private:
  /// A Screen object
  Gyoto::SmartPointer<Gyoto::Screen> Screen;
 public:
  /// Cast from Screen
  Value(Gyoto::SmartPointer<Gyoto::Screen>);
  /// Cast to Screen
  operator Gyoto::SmartPointer<Gyoto::Screen>() const;
};

#endif
