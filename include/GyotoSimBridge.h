/**
 * \file GyotoSimBridge.h
 * \brief Object which computes radiative quantities from numerical simulations files
 *
 *  The accretion-ejection is described by a set of FITS files for a set of different times
 */

/*
    Copyright 2024 Aimar Nicolas

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

#ifndef __GyotoSimBridge_H_ 
#define __GyotoSimBridge_H_ 

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstring>
#include <string>

#include <GyotoBlackBodySpectrum.h>
#include <GyotoKappaDistributionSynchrotronSpectrum.h>
#include <GyotoPowerLawSynchrotronSpectrum.h>
#include <GyotoThermalSynchrotronSpectrum.h>
#include <GyotoFitsRW.h>
#include <GyotoStandardAstrobj.h>
#include <GyotoMetric.h>

namespace Gyoto{
  namespace Astrobj { class SimBridge; }
}

/**
 * \class Gyoto::Astrobj::SimBridge
 * \brief Object that read physical quantities from a set of FITS files comming from numerical simulations.
 * 
 * The coordinate system could be spherical or cartesian (thus the names of quantities which depends on the coordinate system, like the velocity, are labeled with numbers from 0 to 3 
 * instead of {t, r, \theta, \phi} or {t, x, y z}).
 * 
 * To make the interpolation, the user have to provide the arrays of the grid for each dimensions with the following labels in the FITS files :
 * - 'X0' for the time
 * - 'X1' for radius or X (in spherical or cartesian coordinate system)
 * - 'X2' for \theta or Y (in spherical or cartesian coordinate system)
 * - 'X3' for \phi or Z (in spherical or cartesian coordinate system)
 * - 'FREQ' for the frequency array (in case the quantities are radiative coefficients like for (GR)PIC simulations) (Optional)
 * 
 * The grid DOES NOT need to be regularly spaced.
 * Each FITS file should contain in the header of the primary HDU the length of each grid arrays and the time of the FITS file.
 * 
 * 
 * If the numerical simulation is not 3D, the grid in the FITS files still must be 3D.
 * The dimension not simulated has only one point, i.e. the point of the simulation.
 * For exemple, one wants to make images of a 2D disk in the equatorial plane from GRMHD simulation.
 * If the coordinate system is spherical, then the 'X2' array which correspond to the polar angle \theta constains only one value: PI/2.
 * If the coordinate system is cartezian, then the 'X3' array which correspond to the Z axis constains only one value: 0.
 * 
 * 
 * For GRMHD simulations, the required physical quantities (in the FITS files) are 'NUMBERDENSITY', 'TEMPERATURE' and the velocity (see below).
 * In this case, the user have to set to true the flag 'temperature(bool)'.
 * 
 * The code will try to read the magnetic field components 'B0', 'B1", 'B2', 'B3'.
 * If not provided, the user can specify the a magnetic field configuration which will be defined everywhere thanks to the function magneticConfiguration(std::string).
 * The possible arguments are "None", "Toroidal", "Radial", "Vertical".
 * If both the magnetic field is not in the FITS files and the magneticConfiguration is to "None" (by default), the radiative coefficients will be averaged over pitch angle.
 * Although, if the magnetic field is not provided in FITS files, the user MUST set the magnetization parameter to define the magnetic field strength with the function magnetization(double).
 * 
 * 
 * For (GR)PIC simulations (or equivalent), the required quantities are 'J_I', the total emission coefficient and the velocity (see below).
 * All the other radiative (polarized) quantities ('J_Q', 'J_U', 'J_V', 'ALPHA_I', 'ALPHA_Q', 'ALPHA_U', 'ALPHA_V', 'R_Q', 'R_U', 'R_V') are optional.
 * In that case, the frequency array MUST be set (even with only one value).
 * 
 * 
 * Either from GRMHD and (GR)PIC simulations, the velocity could be the 3D velocity (in the emitter frame) or the 4-velocity.
 * The names for the velocity components are 'VELOCITY0', 'VELOCITY1', VELOCITY2', 'VELOCITY3'.
 * If 'VELOCITY0' is not provided, the code assume that the 3 other components (which are mendatory) are the 3D velocity in the emitter frame 
 * and will compute the 4-velocity from the metric.
 * 
 *
 */
class Gyoto::Astrobj::SimBridge : public Gyoto::Astrobj::Standard, public FitsRW {
  friend class Gyoto::SmartPointer<Gyoto::Astrobj::SimBridge>;
 
 protected:
  SmartPointer<Spectrum::BlackBody> spectrumBB_; ///< Black Body
  SmartPointer<Spectrum::KappaDistributionSynchrotron> spectrumKappaSynch_; // kappa-distribution synchrotron spectrum
  SmartPointer<Spectrum::PowerLawSynchrotron> spectrumPLSynch_; // PL-distribution synchrotron spectrum
  SmartPointer<Spectrum::ThermalSynchrotron> spectrumThermalSynch_; // Thermal distribution synchrotron spectrum
  ///< emission law
 private:
  std::string dirname_; ///< FITS files directory
  std::string fname_; ///< FITS files prefix (without the number neither the extension, i.e. '.fits')
  bool temperature_; ///< 1 if temperature is given in fits data file, 0 if emission coef is directly given
  std::string emission_; // Type of emission : Black Body or synchrotron from electron distribution (thermal, PL, kappa)
  double PLindex_; ///< power law index such that density_elec(E) &prop; E<SUP>-p</SUP>
  double gammaMin_; ///< minimum value of gamma for power law energy density
  double gammaMax_; ///< maximum value of gamma for power law energy density
  
  std::string magneticConfig_; ///< Magnetic field geometry (toroidal, vertical) if magnetic field not present in FITS files
  double magnetizationParameter_; ///< magnetization parameter if magnetic field not present in FITS files

  double openingAngle_; ///< angle that define the geometry of the object FOR testing purpose only.
  
  double floortemperature_; ///< 

  //Number of dimensions (maximum 4D : [t, r, theta, phi] or [t, x, y, z]). When the number of dimension is lower than 4, the array of the complement dimension(s) is of length 1 (in the FITS files headers).
  double* time_array_;
  double* x1_array_;
  double* x2_array_;
  double* x3_array_;
  double* nu_array_;
  int ntime_;
  int nx1_;
  int nx2_;
  int nx3_;
  int nnu_; // optional

  std::string* boundCond_;

  // Constructors - Destructor
  // -------------------------
 public:
  GYOTO_OBJECT;

  SimBridge(); ///< Standard constructor
  
  SimBridge(const SimBridge& ) ;///< Copy constructor
  virtual SimBridge* clone () const; ///< Cloner
  
  virtual ~SimBridge(); ///< Destructor

  virtual std::string className() const ; ///< "SimBridge"
  virtual std::string className_l() const ; ///< "simbridge"
  
  // Accessors
  // ---------
 public:
  void directory(std::string const &d);
  std::string directory() const;
  void filename(std::string d);
  std::string filename() const;
  void PLindex(double pl);
  double PLindex()const;
  void gammaMin(double gmin);
  double gammaMin() const;
  void gammaMax(double gmax);
  double gammaMax() const;
  void temperature(bool t);
  bool temperature() const;
  void floorTemperature(double t);
  double floorTemperature()const;
  void magneticConfiguration(std::string config);
  std::string magneticConfiguration() const;
  void magnetization(double ss);
  double magnetization() const;
  void boundaryConditions(std::string x0BC, std::string x1BC, std::string x2BC, std::string x3BC, std::string freqBC="None"); ///< the array of string must contain 4 strings even if the number of dimensions is lower than 4.
  std::string* boundaryConditions() const;
  void emissionType(std::string const &kind);
  std::string emissionType() const;

  virtual void radiativeQ(double Inu[], double Taunu[], double const nu_em[], size_t nbnu,
			  double dsem, state_t const &coord_ph, double const coord_obj[8]=NULL) const ;

  virtual void radiativeQ(double *Inu, double *Qnu, double *Unu, double *Vnu,
			  Eigen::Matrix4d *Onu, double const *nuem , size_t nbnu, double dsem,
			  state_t const &coord_ph, double const *coord_obj) const ;

  virtual void getVelocity(double const pos[4], double vel[4]);

  virtual double operator()(double const coord[4]);

  void openingAngle(double angle);
  double openingAngle() const;

};

#endif
