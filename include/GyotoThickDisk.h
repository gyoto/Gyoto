/**
 * \file GyotoThickDisk.h
 * \brief A thick accretion disk described by its opening
 * angle between the BH spin axis and the disk surface,
 * and its inner radius. 
 * 
 * Density is assumed to follow a r^{-2} law while temperature
 * is a power law with a specified slope.
 *
 * This astrobj emits thermal synchrotron.
 */

/*
    Copyright 2019 Frederic Vincent, Thibaut Paumard

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

#ifndef __GyotoThickDisk_H_ 
#define __GyotoThickDisk_H_ 

#include <iostream>
#include <fstream>
#include <iomanip>

namespace Gyoto{
  namespace Astrobj { class ThickDisk; }
}

#include <GyotoStandardAstrobj.h>
#include <GyotoKappaDistributionSynchrotronSpectrum.h>
#include <GyotoThermalSynchrotronSpectrum.h>

/**
 * \class Gyoto::Astrobj::ThickDisk
 * \brief A thick accretion disk described by its opening
 * angle between the BH spin axis and the disk surface,
 * and its inner radius. 
 * 
 * Density is assumed to follow a r^{-2} law while temperature
 * is a power law with a specified slope.
 *
 * This astrobj emits thermal synchrotron.
 */

class Gyoto::Astrobj::ThickDisk
: public Astrobj::Standard,
  public Hook::Listener
{
  friend class Gyoto::SmartPointer<Gyoto::Astrobj::ThickDisk>;
 private:
  SmartPointer<Spectrum::ThermalSynchrotron> spectrumThermalSynch_;
  double thickDiskOpeningAngle_; ///< ThickDisk opening angle (rad)
  double thickDiskInnerRadius_; ///< Inner disk radius in M units
  double numberDensityAtInnerRadius_cgs_; ///< electron nb density at inner radius (cgs)
  double temperatureAtInnerRadius_; ///< electron temperature at inner radius (K)
  double temperatureSlope_; ///< electron temperature \propto z^temperatureSlope_
  double magnetizationParameter_; ///< P<SUB>magn</SUB>/(n<SUB>e</SUB> m<SUB>p</SUB> c<SUP>2</SUP>)
  double veloZAMONorm_; ///< ZAMO-observed velocity norm below ISCO
  double Vphi_over_V_; ///< Vphi/V where V is the ZAMO-observed velocity below ISCO expressed in a unit-vector basis
  double radius_interpol_vphi_; ///< radius at which Vphi/V reaches zero (typically the horizon); Vphi/V is interpolated between 1 at ISCO and 0 at this radius

  // Constructors - Destructor
  // -------------------------
 public:
  GYOTO_OBJECT;
  GYOTO_OBJECT_THREAD_SAFETY;
  
  ThickDisk(); ///< Standard constructor
  
  ThickDisk(const ThickDisk& ) ;///< Copy constructor
  virtual ThickDisk* clone () const; ///< Cloner
  
  virtual ~ThickDisk() ;                        ///< Destructor
  
  // Accessors
  // ---------
 public:
  void thickDiskOpeningAngle(double ang);
  double thickDiskOpeningAngle() const;
  void thickDiskInnerRadius(double hh);
  double thickDiskInnerRadius() const;
  double numberDensityAtInnerRadius() const;
  double numberDensityAtInnerRadius(std::string const &unit) const;
  void numberDensityAtInnerRadius(double ne);
  void numberDensityAtInnerRadius(double dens, std::string const &unit);
  void temperatureAtInnerRadius(double tt);
  double temperatureAtInnerRadius()const;
  void temperatureSlope(double ss);
  double temperatureSlope()const;
  void magnetizationParameter(double rr);
  double magnetizationParameter()const;
  void velocityBelowIsco(std::vector<double> const &v);
  std::vector<double> velocityBelowIsco() const;
  void velocityBelowIscoInterpol(std::vector<double> const &v);
  std::vector<double> velocityBelowIscoInterpol() const;
  
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
