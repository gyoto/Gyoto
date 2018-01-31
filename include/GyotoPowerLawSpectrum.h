/**
 * \file GyotoPowerLawSpectrum.h
 * \brief A power law spectrum : I_nu=constant_*nu^exponent_
 *
 *  Light emitted by an astronomical object
 */

/*
    Copyright 2011-2014, 2016, 2018 Thibaut Paumard & Frederic Vincent

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
#include <string>

namespace Gyoto {
  namespace Spectrum {
    class PowerLaw;
  }
}


/**
 * \class Gyoto::Spectrum::PowerLaw
 * \brief I_nu=constant_*nu^exponent_
 *
 *  Light emitted by e.g. a Star.
 *
 *  XML stanza:
 *  \code
 *    <Spectrum kind="PowerLaw">
 *      <Exponent> 0. </Exponent>
 *      <Constant> 1. </Constant>
 *    </Spectrum>
 *  \endcode
 */
class Gyoto::Spectrum::PowerLaw : public Gyoto::Spectrum::Generic {
  friend class Gyoto::SmartPointer<Gyoto::Spectrum::PowerLaw>;
 protected:
  double constant_; ///< I_nu=constant_*nu^exponent_
  double exponent_; ///< I_nu=constant_*nu^exponent_
  double minfreq_; ///< Lower-frequency cutoff (emits at nu>=minfreq_)
  double maxfreq_; ///< Upper-frequency cutoff (emits at nu<=maxfreq_)

 public:
  GYOTO_OBJECT;

  PowerLaw();

  /**
   * \brief Constructor setting exponent_ and optionally constant_
   */
  PowerLaw(double exponent, double constant=1.);
  //  PowerLaw(const Spectrum &);
  virtual PowerLaw * clone() const; ///< Cloner

  double constant() const; ///< Get constant_
  void constant(double); ///< Set constant_
  double exponent() const; ///< Get exponent_
  void exponent(double); ///< Set exponent_
  std::vector<double> cutoff(std::string const &unit) const; ///< Get cutoffs, specifying unit
  void cutoff(std::vector<double> const &v, std::string const &unit); ///< Set cutoffs, specifying unit
  std::vector<double> cutoff() const; ///< Get cutoffs, in Hz
  void cutoff(std::vector<double> const &v); ///< Set cutoffs, in Hz

  using Gyoto::Spectrum::Generic::operator();
  virtual double operator()(double nu) const;

};

#endif
