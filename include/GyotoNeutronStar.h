/**
 * \file GyotoNeutronStar.h
 * \brief Neutron star defined by its surface ; no emission
 *
 *  The target of ray-traced Gyoto::Photon
 */

/*
  Copyright (c) 2015 Frederic Vincent, 2017 Thibaut Paumard
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


#ifndef __GyotoNeutronStar_H_
#define __GyotoNeutronStar_H_

#include <GyotoStandardAstrobj.h>
#include <GyotoNumericalMetricLorene.h>

namespace Gyoto{
  namespace Astrobj { class NeutronStar; }
}

/**
 * \class Gyoto::Astrobj::NeutronStar
 * \brief Neutron star defined by its surface ; no emission
 *
 * The underlying Gyoto::Metric::Generic #gg_ instance must be a
 * Gyoto::Metric::NumericalMetricLorene describing a neutron star.
 *
 */
class Gyoto::Astrobj::NeutronStar : public Astrobj::Standard {
  friend class Gyoto::SmartPointer<Gyoto::Astrobj::NeutronStar>;

 protected:
  SmartPointer<Metric::NumericalMetricLorene> gg_; ///< Underlying metric

 public:
  GYOTO_OBJECT;
  NeutronStar(); ///< Standard constructor
  NeutronStar(std::string kin);
  NeutronStar(const NeutronStar& o); ///< Copy constructor
  virtual NeutronStar * clone() const ; ///< Cloner
  virtual ~NeutronStar() ; ///< Destructor

 public:
  virtual Gyoto::SmartPointer<Metric::Generic> metric() const; ///< Get gg_
  virtual void metric(SmartPointer<Metric::Generic> met); ///< Set gg_

  virtual double operator()(double const coord[4]);
  virtual void getVelocity(double const pos[4], double vel[4]) ;
};

#endif
