/**
 * \file GyotoThickDisk.h
 * \brief A thick accretion disk described by its 
 * inner radius and the fwhm of the Gaussian factor affecting
 * the density out of the equatorial plane. 
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
 * \brief A thick accretion disk described by its 
 * inner radius and the fwhm of the Gaussian factor affecting
 * the density out of the equatorial plane. 
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
  double thickDiskInnerRadius_; ///< Inner disk radius in M units
  double thickDiskZGaussianSigma_; ///< Stdev of the Gaussian modulating the density along z in units of cylindrical radius
  bool use_selfabsorption_; ///< True if selfabs is used in radiative transfer
  double alpha_veloparam_; ///< alpha such that u^r = u^r_circ + (1-alpha)*(u^r_rad - u^r_circ)
  double beta_veloparam_; ///< beta such that Omega = Omega_circ + (1-beta)*(Omega_rad - Omega_circ)
  double numberDensityAtInnerRadius_cgs_; ///< electron nb density at inner radius (cgs)
  double temperatureAtInnerRadius_; ///< electron temperature at inner radius (K)
  double temperatureSlope_; ///< electron temperature \propto r^{-temperatureSlope_}
  double densitySlope_; ///< electron density \propto r^{-densitySlope_}
  double magnetizationParameter_; ///< P<SUB>magn</SUB>/(n<SUB>e</SUB> m<SUB>p</SUB> c<SUP>2</SUP>)
  std::string magneticConfig_; ///< Specify the magnetic field configuration for polarisation

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
  void thickDiskInnerRadius(double hh);
  double thickDiskInnerRadius() const;
  void thickDiskZGaussianSigma(double sig);
  double thickDiskZGaussianSigma() const;
  void useSelfAbsorption(bool abs) ;
  bool useSelfAbsorption() const;
  void veloParam(std::vector<double> const &v);
  std::vector<double> veloParam() const;
  double numberDensityAtInnerRadius() const;
  double numberDensityAtInnerRadius(std::string const &unit) const;
  void numberDensityAtInnerRadius(double ne);
  void numberDensityAtInnerRadius(double dens, std::string const &unit);
  void temperatureAtInnerRadius(double tt);
  double temperatureAtInnerRadius()const;
  void temperatureSlope(double ss);
  double temperatureSlope()const;
  void densitySlope(double ss);
  double densitySlope()const;
  void magnetizationParameter(double rr);
  double magnetizationParameter()const;
  void magneticConfiguration(std::string config);
  std::string magneticConfiguration() const;
  
 public:
  using Generic::metric;
  virtual void metric(SmartPointer<Metric::Generic>);
    
  virtual double operator()(double const coord[4]) ;

  virtual void radiativeQ(double Inu[], double Taunu[], 
			  double const nu_em[], size_t nbnu,
			  double dsem, state_t const &coord_ph,
			  double const coord_obj[8]=NULL) const ;
  virtual void getVelocity(double const pos[4], double vel[4]) ;

  virtual void radiativeQ(double *Inu, double *Qnu, double *Unu,
      double *Vnu,
      Eigen::Matrix4d *Onu,
      double const *nuem , size_t nbnu,
      double dsem,
      state_t const &cph,
      double const *co) const;

};

#endif
