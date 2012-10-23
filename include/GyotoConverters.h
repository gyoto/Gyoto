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
    ut_system * getSystem();
#endif
    void Init();
    double ToMeters(double, const std::string &,
		  const Gyoto::SmartPointer<Gyoto::Metric::Generic> &gg=NULL);
    double FromMeters(double, const std::string &,
		  const Gyoto::SmartPointer<Gyoto::Metric::Generic> &gg=NULL);
    double ToSeconds(double, const std::string &unit,
		  const Gyoto::SmartPointer<Gyoto::Metric::Generic> &gg=NULL);
    double FromSeconds(double, const std::string &unit,
		  const Gyoto::SmartPointer<Gyoto::Metric::Generic> &gg=NULL);
    double ToKilograms(double, const std::string &);
    double FromKilograms(double, const std::string &);
    double ToGeometrical(double, const std::string &,
		  const Gyoto::SmartPointer<Gyoto::Metric::Generic> &);
    double FromGeometrical(double, const std::string &,
		  const Gyoto::SmartPointer<Gyoto::Metric::Generic> &);
    bool areConvertible(const std::string &, const std::string&);
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
  Unit(const std::string &);
  ~Unit();
  double To (double val, const std::string &from_unit);
  double From (double val, const std::string &to_unit);
  operator std::string() const ;
  operator ut_unit*() const ;
};

class Gyoto::Units::Converter : protected Gyoto::SmartPointee {
  friend class Gyoto::SmartPointer<Gyoto::Units::Converter>;
 private:
  Gyoto::SmartPointer<Gyoto::Units::Unit> from_;
  Gyoto::SmartPointer<Gyoto::Units::Unit> to_;
  cv_converter * converter_;
  void resetConverter_();
 public:
  Converter(const std::string &,
	    const std::string &);
  Converter(const Gyoto::SmartPointer<Gyoto::Units::Unit>&,
	    const std::string&);
  Converter(const std::string &,
	    const Gyoto::SmartPointer<Gyoto::Units::Unit>&);
  Converter(const Gyoto::SmartPointer<Gyoto::Units::Unit>&,
	    const Gyoto::SmartPointer<Gyoto::Units::Unit>&);
  ~Converter();
  double operator()(double) const ;
};

#endif

#endif
