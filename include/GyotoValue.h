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


#ifndef __GyotoValue_H_
#define __GyotoValue_H_

#include "GyotoConfig.h"
#include "GyotoSmartPointer.h"
#include <string>
#include <vector>

namespace Gyoto {
  namespace Metric {class Generic;}
  namespace Spectrum {class Generic;}
}

class Gyoto::Value {
 public:
  Value();
  ~Value();

  Value& operator=(Value const&);

  double Double;
  Value(double);
  operator double() const;

  bool Bool;
  Value(bool);
  operator bool() const;

  long Long;
  Value(long);
  operator long() const;

  std::string String;
  Value(std::string);
  operator std::string() const;

  std::vector<double> VDouble;
  Value(std::vector<double>);
  operator std::vector<double>() const;

  Gyoto::SmartPointer<Gyoto::Metric::Generic> Metric;
  Value(Gyoto::SmartPointer<Gyoto::Metric::Generic>);
  operator Gyoto::SmartPointer<Gyoto::Metric::Generic>();

  Gyoto::SmartPointer<Gyoto::Spectrum::Generic> Spectrum;
  Value(Gyoto::SmartPointer<Gyoto::Spectrum::Generic>);
  operator Gyoto::SmartPointer<Gyoto::Spectrum::Generic>();
};

#endif
