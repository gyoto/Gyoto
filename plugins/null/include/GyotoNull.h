/*
    Copyright © 2016 Thibaut Paumard

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

/**
 * \file GyotoNull.h
 * \brief Null Astrobj, just for investigating geodesics
 *
 *  The target of ray-traced Gyoto::Photon
 */


#ifndef __GyotoNullAstrobj_H_
#define __GyotoNullAstrobj_H_

#include <iostream>
#include <fstream>
#include <iomanip>

namespace Gyoto{
  namespace Astrobj { class Null; }
  class FactoryMessenger;
  namespace Spectrum {
    class Generic;
  }
}

#include <GyotoAstrobj.h>
#include <GyotoMetric.h>

/**
 * \class Gyoto::Null. 
 * \brief Empty Astrobj
 *
 *  The target of ray-traced Gyoto::Photon
 */
class Gyoto::Astrobj::Null : public Astrobj::Generic {
  friend class Gyoto::SmartPointer<Gyoto::Astrobj::Null>;
 public:
  Null();
  Null(const Null& o);
  virtual Null * clone() const ;
  virtual ~Null() ;
  virtual int Impact(Photon *ph, size_t index, Astrobj::Properties *data=NULL);
};

#endif
