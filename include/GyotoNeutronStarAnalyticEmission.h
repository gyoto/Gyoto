/**
 * \file GyotoNeutronStarAnalyticEmission.h
 * \brief Neutron star emitting at its surface an analytic
 *  emission, typically blackbody
 *
 *  The target of ray-traced Gyoto::Photon
 */

/*
  Copyright (c) 2017-2018 Frederic Vincent, Thibaut Paumard
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


#ifndef __GyotoNeutronStarAnalyticEmission_H_
#define __GyotoNeutronStarAnalyticEmission_H_

#include <GyotoStandardAstrobj.h>
#include <GyotoNumericalMetricLorene.h>
#include <GyotoSpectrum.h>
#include <GyotoNeutronStar.h>


namespace Gyoto{
  namespace Astrobj { class NeutronStarAnalyticEmission; }
}

/**
 * \class Gyoto::Astrobj::NeutronStarAnalyticEmission
 * \brief Neutron star emitting at its surface an analytic
 *  emission, typically blackbody
 *
 * The emission law is given by #spectrum_.
 *
 */
class Gyoto::Astrobj::NeutronStarAnalyticEmission : public Astrobj::NeutronStar {
  friend class Gyoto::SmartPointer<Gyoto::Astrobj::NeutronStarAnalyticEmission>;

 protected:
  SmartPointer<Spectrum::Generic> spectrum_; ///< Emission spectrum

 public:
  GYOTO_OBJECT;
  NeutronStarAnalyticEmission(); ///< Standard constructor
  NeutronStarAnalyticEmission(const NeutronStarAnalyticEmission& o); ///< Copy constructor
  virtual NeutronStarAnalyticEmission * clone() const ; ///< Cloner
  virtual ~NeutronStarAnalyticEmission() ; ///< Destructor

 public:
  virtual void spectrum(SmartPointer<Spectrum::Generic>);
  virtual SmartPointer<Spectrum::Generic> spectrum() const;

  virtual double emission(double nu_em, double dsem,
			  state_t const &_ph, double const _obj[8]=NULL) const;

};

#endif
