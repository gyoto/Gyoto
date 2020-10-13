/**
 * \file GyotoPlasmoid.h
 * \brief Plasmoid sphere formed by magnetic reconnection following a Star orbit, emitting synchrotron,
 * with two distributions of electrons:
 * one thermal at "low" temperature and one kappa at "high" temperature
 *
 */

/*
    Copyright 2019 Frederic Vincent, Thibaut Paumard, Nicolas Aimar

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


#ifndef __GyotoPlasmoid_H_ 
#define __GyotoPlasmoid_H_ 

namespace Gyoto{
  namespace Astrobj { class Plasmoid; }
}

#include <GyotoMetric.h>
#include <GyotoStar.h>
#include <GyotoPowerLawSynchrotronSpectrum.h>
#include <GyotoThermalSynchrotronSpectrum.h>

#ifdef GYOTO_USE_XERCES
#include <GyotoRegister.h>
#endif

#include <string>

/**
 * \class Gyoto::Astrobj::Plasmoid
 * \brief Plasmoid sphere formed by magnetic reconnection following a Star orbit, emitting synchrotron,
 * with two distributions of electrons:
 * one thermal at "low" temperature and one kappa at "high" temperature
 *
 */
class Gyoto::Astrobj::Plasmoid :
  public Gyoto::Astrobj::Star {
  friend class Gyoto::SmartPointer<Gyoto::Astrobj::Plasmoid>;
  
  // Data : 
  // -----
 private:
  double numberDensity_cgs_; ///< cgs-unit number density of plasmoid
  double temperature_ini_; ///< initial temperature of plasmoid before reconnection
  double temperature_rec_; ///< temperature of plasmoid after reconnection
  double magnetizationParameter_; ///< magnetization parameter
  double powerIndex_; ///< plasmoid synchrotron power law distribution index
  SmartPointer<Spectrum::ThermalSynchrotronSpectrum> spectrumLowThermalSynch_; // thermal-distribution synchrotron spectrum at low Temperature
  SmartPointer<Spectrum::ThermalSynchrotronSpectrum> spectrumHighThermalSynch_; // thermal-distribution synchrotron spectrum at high Temperature
  SmartPointer<Spectrum::PowerLawSynchrotronSpectrum> spectrumPowerLawSynch_; // power law distribution synchrotron spectrum

  // Constructors - Destructor
  // -------------------------
 public:
  GYOTO_OBJECT; // This object has a (non-inherited) Property list

 /**
  * Create Plasmoid object with undefined initial conditions. One needs to
  * set the coordinate system, the metric, and the initial position
  * and velocity before integrating the orbit. setInititialCondition()
  * can be used for that.
  */
  Plasmoid(); ///< Default constructor
  
  Plasmoid(const Plasmoid& orig); ///< Copy constructor
  virtual Plasmoid * clone() const ;

  virtual ~Plasmoid() ;                        ///< Destructor
  
 public:
  virtual std::string className() const ; ///< "Plasmoid"
  virtual std::string className_l() const ; ///< "inflate_star"

 public:
  double numberDensity() const;
  double numberDensity(std::string const &unit) const;
  void numberDensity(double ne);
  void numberDensity(double dens, std::string const &unit);
  double temperature_ini() const;
  void temperature_ini(double tt);
  double temperature_rec() const;
  void temperature_rec(double tt);
  void magnetizationParameter(double rr);
  double magnetizationParameter() const;
  double kappaIndex() const;
  void kappaIndex(double);
  
  virtual void radiativeQ(double Inu[], double Taunu[], 
			  double const nu_em[], size_t nbnu,
			  double dsem, state_t const &coord_ph,
			  double const coord_obj[8]=NULL) const ;

};


#endif