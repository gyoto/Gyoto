/**
 * \file GyotoSphericalAccretion.h
 * \brief A spherically-symmetric accretion flow radially
 * falling onto the central object.
 * 
 * Density is assumed to follow a r^{-2} law while temperature
 * is a power law with a specified slope.
 *
 * This astrobj emits thermal synchrotron.
 */

/*
    Copyright 2021 Frederic Vincent, Thibaut Paumard

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

#ifndef __GyotoSphericalAccretion_H_ 
#define __GyotoSphericalAccretion_H_ 

#include <iostream>
#include <fstream>
#include <iomanip>

namespace Gyoto{
  namespace Astrobj { class SphericalAccretion; }
}

#include <GyotoStandardAstrobj.h>
#include <GyotoKappaDistributionSynchrotronSpectrum.h>
#include <GyotoThermalSynchrotronSpectrum.h>

/**
 * \class Gyoto::Astrobj::SphericalAccretion
 * \brief A spherically-symmetric accretion flow radially
 * falling onto the central object.
 * 
 * Density is assumed to follow a r^{-2} law while temperature
 * is a power law with a specified slope.
 *
 * This astrobj emits thermal synchrotron.
 */

class Gyoto::Astrobj::SphericalAccretion
: public Astrobj::Standard,
  public Hook::Listener
{
  friend class Gyoto::SmartPointer<Gyoto::Astrobj::SphericalAccretion>;
 private:
  SmartPointer<Spectrum::ThermalSynchrotron> spectrumThermalSynch_;
  bool use_selfabsorption_; ///< True if selfabs is used in radiative transfer
  double sphericalAccretionInnerRadius_; ///< Inner radius of flow in M units
  double numberDensityAtInnerRadius_cgs_; ///< electron nb density at inner radius (cgs)
  double densitySlope_; ///< electron density \propto r^{-densitySlope_}
  double temperatureAtInnerRadius_; ///< electron temperature at inner radius (K)
  double temperatureSlope_; ///< electron temperature \propto z^temperatureSlope_
  double magnetizationParameter_; ///< P<SUB>magn</SUB>/(n<SUB>e</SUB> m<SUB>p</SUB> c<SUP>2</SUP>)

  // Constructors - Destructor
  // -------------------------
 public:
  GYOTO_OBJECT;
  GYOTO_OBJECT_THREAD_SAFETY;
  
  SphericalAccretion(); ///< Standard constructor
  
  SphericalAccretion(const SphericalAccretion& ) ;///< Copy constructor
  virtual SphericalAccretion* clone () const; ///< Cloner
  
  virtual ~SphericalAccretion() ;                        ///< Destructor
  
  // Accessors
  // ---------
 public:
  void useSelfAbsorption(bool abs) ;
  bool useSelfAbsorption() const;
  void sphericalAccretionInnerRadius(double hh);
  double sphericalAccretionInnerRadius() const;
  double numberDensityAtInnerRadius() const;
  double numberDensityAtInnerRadius(std::string const &unit) const;
  void numberDensityAtInnerRadius(double ne);
  void numberDensityAtInnerRadius(double dens, std::string const &unit);
  void densitySlope(double ss);
  double densitySlope()const;
  void temperatureAtInnerRadius(double tt);
  double temperatureAtInnerRadius()const;
  void temperatureSlope(double ss);
  double temperatureSlope()const;
  void magnetizationParameter(double rr);
  double magnetizationParameter()const;
  
 public:
  using Generic::metric;
  virtual void metric(SmartPointer<Metric::Generic>);
    
  virtual double operator()(double const coord[4]) ;

  virtual void radiativeQ(double Inu[], double Taunu[], 
			  double const nu_em[], size_t nbnu,
			  double dsem, state_t const &coord_ph,
			  double const coord_obj[8]=NULL) const ;
  virtual void getVelocity(double const pos[4], double vel[4]) ;

};

#endif
