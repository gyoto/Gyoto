/**
 * \file GyotoThermalBremsstrahlungSpectrum.h
 * \brief Thermal brems spectrum
 *
 */

/*
    Copyright 2014 Frederic Vincent

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

#ifndef __GyotoThermalBremsstrahlungSpectrum_H_ 
#define __GyotoThermalBremsstrahlungSpectrum_H_ 
#include "GyotoSpectrum.h"
#include <GyotoBlackBodySpectrum.h>

namespace Gyoto {
  namespace Spectrum {
    class ThermalBremsstrahlung;
  }
}

/**
 * \class Gyoto::Spectrum::ThermalBremsstrahlung
 * \brief Thermal brems spectrum
 *
 *
 *  Example XML entity:
 *  \code
 *   <Spectrum kind="ThermalBremsstrahlung">
 *   </Spectrum>
 *  \endcode
 *
 */
class Gyoto::Spectrum::ThermalBremsstrahlung : public Gyoto::Spectrum::Generic {
  friend class Gyoto::SmartPointer<Gyoto::Spectrum::ThermalBremsstrahlung>;
 protected:
  SmartPointer<Spectrum::BlackBody> spectrumBB_; ///< blackbody emission
  double cst_; ///< Scaling constant

 public:
  ThermalBremsstrahlung();

  /**
   * \brief Constructor setting T_ and cst_
   */
  virtual ThermalBremsstrahlung * clone() const; ///< Cloner

  using Gyoto::Spectrum::Generic::operator();
  virtual double operator()(double nu) const;

  double jnu(double nu, double temp,double massdensity);
  double alphanu(double nu, double temp,double massdensity);

#ifdef GYOTO_USE_XERCES
  virtual void fillElement(FactoryMessenger *fmp) const ;
#endif
};

#endif
