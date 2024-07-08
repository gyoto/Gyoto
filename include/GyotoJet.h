/**
 * \file GyotoJet.h
 * \brief Simple jet model with thermal or kappa-distribution
 * synchrotron emission from Pandya et al. (2016)
 *
 * This class implements a jet defined as the volume contained
 * between the two conical/parabolic surfaces defined either 
 * by the angles jetInnerOpeningAngle_
 * and jetOuterOpeningAngle_, or by the parabola parameter such that
 * z = param * rcyl^2. The apex of the cone/parabola is at (0,0).
 * The region below jetInnerRadius_ is removed, so this quantity is
 * the "base of the jet".
 *
 * The Lorentz factor is assumed constant at gammaJet_.
 * The electron number density at the base of the jet is baseNumberDensity_cgs_,
 * its r-evolution is dedictated by mass conservation.
 * The electron temperature is baseTemperature_, its r-evolution is assumed
 * to follow a power law r^temperatureSlope_. The magnetic field
 * amplitude is defined by the magnetization parameter,
 * magnetizationParameter_.
 * 
 * The jet emits synchrotron radiation, assuming that the electrons
 * follow a thermal or kappa distribution, ie the smooth gluing of a thermal
 * distribution at low electron Lorentz factor, to a power-law distribution
 * at high electron Lorentz factor. This distribution, as well as the
 * resulting emission and absorption coefficients are taken from:
 * Pandya et al., ApJ, 822, 34 (2016), section 5.3.3
 */

/*
    Copyright 2017-2024 Frederic Vincent, Thibaut Paumard, Paloma Thévenet

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

#ifndef __GyotoJet_H_ 
#define __GyotoJet_H_ 

#include <iostream>
#include <fstream>
#include <iomanip>

namespace Gyoto{
  namespace Astrobj { class Jet; }
}

#include <GyotoStandardAstrobj.h>
#include <GyotoKappaDistributionSynchrotronSpectrum.h>
#include <GyotoThermalSynchrotronSpectrum.h>

/**
 * \class Gyoto::Astrobj::Jet
 * \brief Simple jet model with thermal or kappa-distribution
 * synchrotron emission from Pandya et al. (2016)
 *
 * This class implements a jet defined as the volume contained
 * between the two conical surfaces defined by angles jetInnerOpeningAngle_
 * and jetOuterOpeningAngle_, with apex located on the black hole rotation
 * axis at altitude jetBaseHeight_ in units of M.
 * 
 * The Lorentz factor is assumed constant at gammaJet_.
 * The electron number density at the base of the jet is baseNumberDensity_cgs_,
 * its z-evolution is dedictated by mass conservation.
 * The electron temperature is baseTemperature_, its z-evolution is assumed
 * to follow a power law z^temperatureSlope_. The magnetic field
 * amplitude is defined by the magnetization parameter,
 * magnetizationParameter_.
 *
 * The jet emits synchrotron radiation, assuming that the electrons
 * follow a thermal or kappa distribution, ie the smooth gluing of a thermal
 * distribution at low electron Lorentz factor, to a power-law distribution
 * at high electron Lorentz factor. This distribution, as well as the
 * resulting emission and absorption coefficients are taken from:
 * Pandya et al., ApJ, 822, 34 (2016), section 5.3.3
 */

class Gyoto::Astrobj::Jet
: public Astrobj::Standard,
  public Hook::Listener
{
  friend class Gyoto::SmartPointer<Gyoto::Astrobj::Jet>;
 private:
  SmartPointer<Spectrum::KappaDistributionSynchrotron> spectrumKappaSynch_;
  SmartPointer<Spectrum::ThermalSynchrotron> spectrumThermalSynch_;
  
  bool parabolic_; ///< True when the jet sheath has a parabolic shape; if false the shape will be conical (following Vincent+19 torus+jet paper)

  bool outflowing_; ///< True when the jet is outflowing. Else, inflowing.
  
  double jetShapeInnerParabolaParam_; ///< The jet shape inner boundary follows z = jetShapeInnerParabolaParam_ * rcyl^2, where rcyl is the cylindrical radius; used when parabolic_=True
  double jetShapeOuterParabolaParam_; ///< The jet shape outer boundary follows z = jetShapeOuterParabolaParam_ * rcyl^2, where rcyl is the cylindrical radius; used when parabolic_=True

  double jetOuterOpeningAngle_; ///< Jet outer opening angle (rad); used when parabolic_=False
  double jetInnerOpeningAngle_; ///< Jet inner opening angle (rad); used when parabolic_=False

  double jetVphiOverVr_; ///< ratio V^(phi)/V^(r) in orthonormal basis where V is ZAMO-measured jet velocity

  double jetStagnationRadius_; ///< Jet outflowing above, inflowing below
  double jetInnerRadius_; ///< Jet inner radius, or "base of the jet", used for scaling the thermo quantities.
  double gammaJet_; ///< Constant Lorentz factor in jet (same def for parabolic and conical jet)
  double baseNumberDensity_cgs_; ///< electron nb density at jet base (cgs)
  double baseTemperature_; ///< electron temperature at jet base (K)
  double temperatureSlope_; ///< electron temperature \propto z^temperatureSlope_
  double magnetizationParameter_; ///< P<SUB>magn</SUB>/(n<SUB>e</SUB> m<SUB>p</SUB> c<SUP>2</SUP>)
  std::string magneticConfig_; ///< Magnetic field configuration

  // Constructors - Destructor
  // -------------------------
 public:
  GYOTO_OBJECT;
  GYOTO_OBJECT_THREAD_SAFETY;
  
  Jet(); ///< Standard constructor
  
  Jet(const Jet& ) ;///< Copy constructor
  virtual Jet* clone () const; ///< Cloner
  
  virtual ~Jet() ;                        ///< Destructor
  
  // Accessors
  // ---------
 public:
  void parabolic(bool parabol);
  bool parabolic() const ;

  void outflowing(bool out);
  bool outflowing() const ;
  
  void jetShapeInnerParabolaParam(double param);
  double jetShapeInnerParabolaParam() const;
  void jetShapeOuterParabolaParam(double param);
  double jetShapeOuterParabolaParam() const;

  void jetOuterOpeningAngle(double ang);
  double jetOuterOpeningAngle() const;
  void jetInnerOpeningAngle(double ang);
  double jetInnerOpeningAngle() const;

  void jetStagnationRadius(double param);
  double jetStagnationRadius() const;
  void jetVphiOverVr(double alpha);
  double jetVphiOverVr()const;
  void jetInnerRadius(double hh);
  double jetInnerRadius() const;
  void gammaJet(double gam);
  double gammaJet() const;
  double baseNumberDensity() const;
  double baseNumberDensity(std::string const &unit) const;
  void baseNumberDensity(double ne);
  void baseNumberDensity(double dens, std::string const &unit);
  void baseTemperature(double tt);
  double baseTemperature()const;
  void temperatureSlope(double ss);
  double temperatureSlope()const;
  void magnetizationParameter(double rr);
  double magnetizationParameter()const;
  void kappaIndex(double index);
  double kappaIndex()const;
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

  virtual void radiativeQ(double Inu[], double Qnu[], double Unu[], 
              double Vnu[], Eigen::Matrix4d Onu[],
              double const nu_ems[], size_t nbnu, double dsem,
              state_t const &coord_ph, double const coord_obj[8]) const;

};

#endif
