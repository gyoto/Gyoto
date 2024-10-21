/**
 * \file GyotoBlob.h
 * \brief Spherical blob of plasma following some (geodesic or not) orbit, 
 * emitting synchrotron radiation that can be modulated in time
 * or in space.
 *
 * The blob's orbit can be circular equatorial (defined by the user-provided constant Omega=dphi/dt),
 * or helical (constant v^r, v^phi found from Newtonian ang mom assumed constant).
 *
 * The blob's volume is defined by that of the background UniformSphere, see the UniformSphere::radius_
 * keyword. 
 *
 * The emission of the blob can be modulated by a temporal
 * Gaussian defined by its peak time timeRef_ and temporal sigma timeSigma_
 * (careful with the units of these). It can also be modulated by a
 * spatial Gaussian such that the UniformSphere::radius_ keyword is interpreted
 * as the 3-sigma spatial extension of the blob. 
 *
 * The electron distribution function can be thermal, kappa or power-law.
 *
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


#ifndef __GyotoBlob_H_ 
#define __GyotoBlob_H_ 

namespace Gyoto{
  namespace Astrobj { class Blob; }
}

#include <GyotoMetric.h>
#include <GyotoUniformSphere.h>
#include <GyotoKappaDistributionSynchrotronSpectrum.h>
#include <GyotoPowerLawSynchrotronSpectrum.h>
#include <GyotoThermalSynchrotronSpectrum.h>

#ifdef GYOTO_USE_XERCES
#include <GyotoRegister.h>
#endif

#include <string>

/**
 * \class Gyoto::Astrobj::Blob
 * \brief Blob of plasma following a given orbit, emitting synchrotron,
 * with Gaussian density and temperature
 *
 */
class Gyoto::Astrobj::Blob :
  public Gyoto::Astrobj::UniformSphere {
  friend class Gyoto::SmartPointer<Gyoto::Astrobj::Blob>;
  
  // Data : 
  // -----
 private:
  bool time_gauss_modulated_; ///< True if blob emission time-modulated by a Gaussian with parameters timeRef_M_ and timeSigma_M_
  bool space_gauss_modulated_; ///< True if blob emitting volume is a Gaussian with 3-sigma extension coinciding with the blob's radius_
  double* init4Coord_; ///< Initial 4-coordinate of the Blob, eg (t,r,theta,phi)
  double* init3Velo_; ///< Initial 3-velocity of the Blob, eg (dr/dt, dtheta/dt, dphi/dt)
  std::string blobMotionType_; ///< Type of motion of the Blob, "Equatorial" is circular constant Omega in equatorial plane; "Helical" if at constant theta, with v^r=cst, and v^phi \propto 1/r^2 (constancy of Newtonian ang mom)
  double numberDensity_cgs_; ///< cgs-unit number density of hotspot
  double temperature_; ///< temperature of hotspot
  double timeRef_M_; ///< M-unit reference time for Gaussian hotspot evolution
  double timeSigma_M_; ///< M-unit temporal sigma for Gaussian hotspot evolution
  double magnetizationParameter_; ///< magnetization parameter
  double kappaIndex_; ///< hotspot synchrotron kappa-distribution index
  SmartPointer<Spectrum::KappaDistributionSynchrotron> spectrumKappaSynch_; // kappa-distribution synchrotron spectrum
  SmartPointer<Spectrum::PowerLawSynchrotron> spectrumPLSynch_; // PL-distribution synchrotron spectrum
  SmartPointer<Spectrum::ThermalSynchrotron> spectrumThermalSynch_; // Thermal distribution synchrotron spectrum
  std::string magneticConfig_; // Magnetic field geometry (toroidal, vertical)
  std::string electronDistrib_; // Electron distribution (thermal, kappa)

  // Constructors - Destructor
  // -------------------------
 public:
  GYOTO_OBJECT; // This object has a (non-inherited) Property list

 /**
  * Create Blob object with undefined initial conditions. One needs to
  * set the coordinate system, the metric, and the initial position
  * and velocity before integrating the orbit. setInititialCondition()
  * can be used for that.
  */
  Blob(); ///< Default constructor
  
  Blob(const Blob& orig); ///< Copy constructor
  virtual Blob * clone() const ;

  virtual ~Blob() ;                        ///< Destructor
  
 public:
  void electronDistribution(const std::string &kind);
  std::string electronDistribution() const;

  void blobMotionType(const std::string &kind) ;
  std::string blobMotionType() const ;
  
  virtual std::string className() const ; ///< "Blob"
  virtual std::string className_l() const ; ///< "inflate_star"

 public:

  bool timeGaussianModulated() const;
  void timeGaussianModulated(bool timemod) ;
  bool spaceGaussianModulated() const;
  void spaceGaussianModulated(bool spacemod) ;
  void init4Coord(std::vector<double> const &v) ;
  std::vector<double> init4Coord() const ;
  void init3Velo(std::vector<double> const &v);
  std::vector<double> init3Velo() const ;
  double numberDensity() const;
  double numberDensity(std::string const &unit) const;
  void numberDensity(double ne);
  void numberDensity(double dens, std::string const &unit);
  double temperature() const;
  void temperature(double tt);
  double timeRef() const;
  double timeRef(std::string const &unit) const;
  void timeRef(double tt);
  void timeRef(double tt, std::string const &unit);
  double timeSigma() const;
  double timeSigma(std::string const &unit) const;
  void timeSigma(double tt);
  void timeSigma(double tt, std::string const &unit);
  void magnetizationParameter(double rr);
  double magnetizationParameter() const;
  double kappaIndex() const;
  void kappaIndex(double);
  void magneticConfiguration(std::string config);
  std::string magneticConfiguration() const;
  
  virtual void radiativeQ(double Inu[], double Taunu[], 
			  double const nu_em[], size_t nbnu,
			  double dsem, state_t const &coord_ph,
			  double const coord_obj[8]=NULL) const ;

  virtual void radiativeQ(double Inu[], double Qnu[], double Unu[], 
              double Vnu[], Eigen::Matrix4d Onu[],
              double const nu_ems[], size_t nbnu, double dsem,
              state_t const &coord_ph, double const coord_obj[8]) const;
  
  void getVelocity(double const pos[4], double vel[4]);
  
  void getCartesian(double const * const dates, size_t const n_dates,
		    double * const x, double * const y, double * const z, 
		    double * const xprime=NULL, double * const yprime=NULL, double * const zprime=NULL);
  

};


#endif
