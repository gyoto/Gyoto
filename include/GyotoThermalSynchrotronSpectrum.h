/**
 * \file GyotoThermalSynchrotronSpectrum.h
 * \brief Thermal synchrotron spectrum
 *
 */

/*
    Copyright 2018 Frederic Vincent, Thibaut Paumard

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

#ifndef __GyotoThermalSynchrotronSpectrum_H_ 
#define __GyotoThermalSynchrotronSpectrum_H_ 
#include "GyotoSpectrum.h"
#include <GyotoBlackBodySpectrum.h>

namespace Gyoto {
  namespace Spectrum {
    class ThermalSynchrotron;
  }
}

/**
 * \class Gyoto::Spectrum::ThermalSynchrotron
 * \brief Thermal synchrotron spectrum
 *
 *
 *  Example XML entity:
 *  \code
 *   <Spectrum kind="ThermalSynchrotron">
 *   </Spectrum>
 *  \endcode
 *
 */
class Gyoto::Spectrum::ThermalSynchrotron : public Gyoto::Spectrum::Generic {
  friend class Gyoto::SmartPointer<Gyoto::Spectrum::ThermalSynchrotron>;
 protected:
  SmartPointer<Spectrum::BlackBody> spectrumBB_; ///< blackbody emission
  double T_; ///< Temperature
  double numberdensityCGS_; ///< Number density in CGS UNITS (careful)
  double angle_B_pem_; ///< Angle between Bfield and emission direction (rad)
  double cyclotron_freq_; ///< Cyclotron frequency (e*B / 2*pi*me*c)
  bool angle_averaged_; ///< Boolean for angle averaging
  double bessel_K2_; ///< Bessel K2 function

 public:
  GYOTO_OBJECT;

  ThermalSynchrotron();
  ThermalSynchrotron(const ThermalSynchrotron &);
  virtual ThermalSynchrotron * clone() const; ///< Cloner

  using Gyoto::Spectrum::Generic::operator();
 /**
   * This function returns the optically thick Inu
   * which is not defined here, returns an error
   *
   * \param nu frequency in Hz
   */
  virtual double operator()(double nu) const;
 /**
   * This function returns the optically thin increment
   * to intensity dI_nu = j_nu*ds*exp(-alpha_nu*ds)
   * in SI units
   *
   * \param nu frequency in Hz
   * \param ds length element in SI (careful to this)
   */
#ifndef GYOTO_SWIGIMPORTED
  virtual double operator()(double nu,double ,double ds) const;
#endif
  // NB: the second argument, opacity in the Spectrum API
  // is useless here

  double temperature() const;
  void temperature(double tt);
  double numberdensityCGS() const;
  void numberdensityCGS(double rho);
  double angle_B_pem() const;
  void angle_B_pem(double rho);
  double cyclotron_freq() const;
  void cyclotron_freq(double rho);
  bool angle_averaged() const;
  void angle_averaged(bool ang);
  double besselK2() const;
  void besselK2(double bessel);
  
 /**
   * Returns the emission coefficient j_nu in cgs units
   * i.e. erg cm^-3 s^-1 ster^-1 Hz^-1
   *
   * \param nu frequency in Hz
   */
  double jnuCGS(double nu) const;
  /**
   * Returns the Stokes Q emission coefficient j_nu in cgs units
   * i.e. erg cm^-3 s^-1 ster^-1 Hz^-1
   *
   * \param nu frequency in Hz
   */
  double jQnuCGS(double nu) const;
  /**
   * Returns the Stokes U emission coefficient j_nu in cgs units
   * i.e. erg cm^-3 s^-1 ster^-1 Hz^-1
   *
   * \param nu frequency in Hz
   */
  double jUnuCGS(double nu) const;
  /**
   * Returns the Stokes V emission coefficient j_nu in cgs units
   * i.e. erg cm^-3 s^-1 ster^-1 Hz^-1
   *
   * \param nu frequency in Hz
   */
  double jVnuCGS(double nu) const;

 /**
   * Returns the absorption coefficient alpha_nu in cgs units [cm^-1]
   *
   * \param nu frequency in Hz
   */
  double alphanuCGS(double nu) const;
  /**
   * Returns the Stokes Q absorption coefficient alpha_nu in cgs units [cm^-1]
   *
   * \param nu frequency in Hz
   */
  double alphaQnuCGS(double nu) const;
  /**
   * Returns the Stokes U absorption coefficient alpha_nu in cgs units [cm^-1]
   *
   * \param nu frequency in Hz
   */
  double alphaUnuCGS(double nu) const;
  /**
   * Returns the Stokes V absorption coefficient alpha_nu in cgs units [cm^-1]
   *
   * \param nu frequency in Hz
   */
  double alphaVnuCGS(double nu) const;

  /**
   * Returns the Stokes Q Faraday rotation coefficient in cgs units [cm^-1]
   * 
   *  \param nu frequency in Hz
   */
  double rQnuCGS(double nu) const;
  /**
   * Returns the Stokes U Faraday rotation coefficient in cgs units [cm^-1]
   * 
   *  \param nu frequency in Hz
   */
  double rUnuCGS(double nu) const;
  /**
   * Returns the Stokes V Faraday rotation coefficient in cgs units [cm^-1]
   * 
   *  \param nu frequency in Hz
   */
  double rVnuCGS(double nu) const;


   /**
   * Returns the emission and absorption coef in SI
   *
   */
  void radiativeQ(double jnu[], // output
		    double anu[], // output
		    double const nu_ems[],
		    size_t nbnu
		    ) ;
  /**
   * Returns the emission, absorption and Fraday rotation coef in SI for the 4 Stokes parameters
   *
   */
  void radiativeQ(double jInu[], double jQnu[], double jUnu[], double jVnu[],
        double aInu[], double aQnu[], double aUnu[], double aVnu[],
        double rQnu[], double rUnu[], double rVnu[],
        double const nu_ems[],
        size_t nbnu );

};

#endif
