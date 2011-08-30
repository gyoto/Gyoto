/**
 * \file GyotoPowerLawSpectrum.h
 * \brief A constant spectrum : I_nu=cst
 *
 *  Light emitted by an astronomical object
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

#ifndef __GyotoPowerLawSpectrum_H_ 
#define __GyotoPowerLawSpectrum_H_ 
#include <GyotoSpectrum.h>

namespace Gyoto {
  namespace Spectrum {
    class PowerLaw;
  }
}


/**
 * \class Gyoto::Spectrum::PowerLaw
 * \brief I_nu=cst
 *
 *  Light emitted by e.g. a Star
 *
 */
class Gyoto::Spectrum::PowerLaw : public Gyoto::Spectrum::Generic {
  friend class Gyoto::SmartPointer<Gyoto::Spectrum::PowerLaw>;
 protected:
  double constant_;
  double exponent_;

 public:
  PowerLaw();
  PowerLaw(double exponent, double constant=1.);
  //  PowerLaw(const Spectrum &);
  virtual PowerLaw * clone() const; ///< Cloner

  double getConstant() const; ///< Get constant
  void setConstant(double);
  double getExponent() const; ///< Get exponent
  void setExponent(double);

  using Gyoto::Spectrum::Generic::operator();
  virtual double operator()(double nu) const;
                          ///< I_nu = mySpectrum(nu), nu in Hz

#ifdef GYOTO_USE_XERCES
  /**
   * Spectrum implementations should impement fillElement to save their
   * parameters to XML and call the generic implementation to save
   * generic parts.
   */

  virtual void fillElement(factoryMessenger *fmp) const ;
                                             /// < called from Factory
#endif
};

#ifdef GYOTO_USE_XERCES
namespace Gyoto {
  namespace Spectrum {
    Gyoto::SmartPointer<Gyoto::Spectrum::Generic>
      PowerLawSubcontractor(Gyoto::factoryMessenger* fmp = NULL);
    void PowerLawInit();
  }
}
#endif
#endif
