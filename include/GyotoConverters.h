/**
 * \file GyotoConverters.h
 * \brief GYOTO utilities
 *
 *  Various utilities
 */

/*
    Copyright 2011 Thibaut Paumard

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

#ifndef __GyotoConverters_H_ 
#define __GyotoConverters_H_ 

#include <GyotoUtils.h>
#include <GyotoSmartPointer.h>

#ifdef HAVE_UDUNITS
#include <udunits2.h>
#endif

#include <string>
#include <sstream>

namespace Gyoto {
  namespace Metric {
    class Generic;
  }
  namespace Units {
#ifdef HAVE_UDUNITS
    class Unit;
    class Converter;
#endif
    void Init();
    double ToMeters(double, std::string);
    double ToKilograms(double, std::string);
    double ToGeometrical(double, std::string,
			 Gyoto::SmartPointer<Gyoto::Metric::Generic>);
  }
}

#ifdef HAVE_UDUNITS
class Gyoto::Units::Unit : protected Gyoto::SmartPointee {
  friend class Gyoto::SmartPointer<Gyoto::Units::Unit>;
  friend class Gyoto::Units::Converter;
 private:
  ut_unit * unit_;
  std::string kind_;
 public:
  Unit(std::string);
  ~Unit();
  double To (double val, std::string from_unit);
  double From (double val, std::string to_unit);
  operator std::string() ;
  operator ut_unit*();
};

class Gyoto::Units::Converter : protected Gyoto::SmartPointee {
  friend class Gyoto::SmartPointer<Gyoto::Units::Converter>;
 private:
  Gyoto::SmartPointer<Gyoto::Units::Unit> from_;
  Gyoto::SmartPointer<Gyoto::Units::Unit> to_;
  cv_converter * converter_;
  void resetConverter_();
 public:
  Converter(std::string, std::string);
  Converter(Gyoto::SmartPointer<Gyoto::Units::Unit>, std::string);
  Converter(std::string, Gyoto::SmartPointer<Gyoto::Units::Unit>);
  Converter(Gyoto::SmartPointer<Gyoto::Units::Unit>,
	    Gyoto::SmartPointer<Gyoto::Units::Unit>);
  ~Converter();
  double operator()(double);
};

#endif

#endif
