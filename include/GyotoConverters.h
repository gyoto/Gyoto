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

#include "GyotoUtils.h"

#ifdef HAVE_UDUNITS
#include <udunits2.h>
#endif

#include <string>
#include <sstream>

namespace Gyoto {
  namespace Units {
#ifdef HAVE_UDUNITS
    class System;
    class Unit;
#endif
    double ToMeters(double, std::string);
    double ToKilograms(double, std::string);
  }
}

#ifdef HAVE_UDUNITS
class Gyoto::Units::System {
  friend Gyoto::Units::Unit;
 protected:
  ut_system * sys_;
 public:
  System(char *);
  ~System();
};

class Gyoto::Units::Unit {
 private:
  ut_unit * unit_;
 public:
  Unit(std::string);
  ~Unit();
  double To (double val, std::string from_unit);
  double From (double val, std::string to_unit);
};
#endif

#endif
