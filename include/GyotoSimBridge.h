/**
 * \file GyotoSimBridge.h
 * \brief Object which computes radiative quantities from numerical simulations files
 *
 *  The accretion-ejection is described by a set of FITS files for a set of different times
 */

/*
    Copyright 2024 Aimar Nicolas, Irene Urso

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
 * To make the interpolation, the user has to provide the arrays of the grid for each dimensions with the following labels in the FITS files :
 * - 'X0' for the time
 * - 'X1' for radius or X (in spherical or cartesian coordinate system)
 * - 'X2' for \theta or Y (in spherical or cartesian coordinate system)
 * - 'X3' for \phi or Z (in spherical or cartesian coordinate system)
 * - 'FREQ' for the frequency array (in case the quantities are radiative coefficients like for (GR)PIC simulations) (Optional)
 * 
 * The grid DOES NOT need to be regularly spaced, but MUST NOT chanfe from time to time.
 * Each FITS file should contain in the header of the primary HDU the length of each grid arrays and the time of the FITS file.
 * 
 * 
 * If the numerical simulation is not 3D, the grid in the FITS files still must be 3D.
 * The dimension not simulated has only one point, i.e. the point of the simulation.
 * For exeample, one wants to make images of a 2D disk in the equatorial plane from GRMHD simulation.
 * If the coordinate system is spherical, then the 'X2' array which correspond to the polar angle \theta constains only one value: PI/2.
 * If the coordinate system is cartezian, then the 'X3' array which correspond to the Z axis constains only one value: 0.
 * 
 * 
 * For GRMHD simulations, the required physical quantities (in the FITS files) are 'NUMBERDENSITY', 'TEMPERATURE' and the velocity (see below).
 * In this case, the user has to set to true the flag 'temperature(bool)'.
 * 
 * The code will try to read the magnetic field components 'B0', 'B1", 'B2', 'B3'.
 * If not provided, the user can specify the magnetic field configuration which will be defined everywhere thanks to the function magneticConfiguration(std::string).
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
 * Either from GRMHD and (GR)PIC simulations, the velocity MUST be the 3D velocity.
 * The names for the velocity components are 'VELOCITY1', VELOCITY2', 'VELOCITY3'.
 * The code will NOT search for 'VELOCITY0'! The 4-velocity will be computed with the metric from the metric.
 * 
 *
 * The geometry of the emitting region is a sphere with the radius being either the maximum radius of the simulation or the maximum radius of integration set by rMax(double).
 * It is highly recommended to NOT use this object as it is and create a new Astrobj that inherits from SimBridge which redefines the operator() function
 * to specify a custom (more appropriate) geometry of the emitting region.
 * Another way to reduce the computation time is to properly set the rmax_.
 * 
 */
class Gyoto::Astrobj::SimBridge : public Gyoto::Astrobj::Standard, public FitsRW {
  friend class Gyoto::SmartPointer<Gyoto::Astrobj::SimBridge>;
 
 protected:
  SmartPointer<Spectrum::BlackBody> spectrumBB_; ///< Black Body
  SmartPointer<Spectrum::KappaDistributionSynchrotron> spectrumKappaSynch_; // kappa-distribution synchrotron spectrum
  SmartPointer<Spectrum::PowerLawSynchrotron> spectrumPLSynch_; // PL-distribution synchrotron spectrum
  SmartPointer<Spectrum::ThermalSynchrotron> spectrumThermalSynch_; // Thermal distribution synchrotron spectrum

 private:
  std::string dirname_; ///< FITS files directory
  std::string fname_; ///< FITS files prefix (without the number neither the extension, i.e. '.fits')
  bool temperature_; ///< 1 if temperature is given in fits data file, 0 if emission coef is directly given
  bool circularmotion_; ///< 1 if velocity is given in fits data file, 0 if circularmotion is directly given
  bool BinFile_; ///< Define if the magnetic field is saved in FITS file or not
  std::string emission_; // Type of emission : Black Body or synchrotron from electron distribution (thermal, PL, kappa)
  double PLindex_; ///< power law index such that density_elec(E) &prop; E<SUP>-p</SUP>
  double gammaMin_; ///< minimum value of gamma for power law energy density
  double gammaMax_; ///< maximum value of gamma for power law energy density
  
  std::string magneticConfig_; ///< Magnetic field geometry (toroidal, vertical) if magnetic field not present in FITS files
  double magnetizationParameter_; ///< magnetization parameter if magnetic field not present in FITS files
  
  double floortemperature_; ///< Set minimum of temperature

  int emisInFile_[4]; ///< Flag that specify which emission coefficients are in the FITS files
  int absInFile_[4]; ///< Flag that specify which absorption coefficients are in the FITS files
  int rotInFile_[3]; ///< Flag that specify which rotation coefficients are in the FITS files

  protected:
  //Number of dimensions (maximum 4D : [t, r, theta, phi] or [t, x, y, z]). When the number of dimension is lower than 4, the array of the complement dimension(s) is of length 1 (in the FITS files headers).
  double* time_array_; ///< array containing the time evolution of each FITS files
  double* x1_array_; ///< First spatial dimension array (radius in spherical, X in Cartesian)
  double* x2_array_; ///< Second spatial dimension array (\theta in spherical, Y in Cartesian)
  double* x3_array_; ///< Third spatial dimension array (\phi in spherical, Z in Cartesian)
  double* nu_array_; ///< frequency array if quantities in FITS files are radiative coefficients

  int ntime_; ///< length of time_array_
  int nx1_; ///< length of x1_array_
  int nx2_; ///< length of x2_array_
  int nx3_; ///< length of x3_array_
  int nnu_; ///< ///< length of nu_array_, optional

  std::string* boundCond_; ///< Table of string which store the boundary conditions of all dimensions.

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
  void directory(std::string const &d); ///< Set the directory where the FITS files are stored
  std::string directory() const; ///< Get the directory where the FITS files will be searched
  void filename(std::string const &d); ///< Set the prefix of the FITS filenames, it should be followed by a 4-digit number (filled by zeros)
  std::string filename() const; ///< Get the FITS filename prefix
  void PLindex(double pl); ///< Set the power law index when the emission is set to "Power-Law" or "Kappa"
  double PLindex()const; ///< Get the power law index of the electron distribution function
  void gammaMin(double gmin); ///< Set the minimum gamma factor when the emission is set to "Power-Law" or "Kappa"
  double gammaMin() const; ///< Get the minimum gamma factor of the electron distribution function
  void gammaMax(double gmax); ///< Set the maximum gamma factor when the emission is set to "Power-Law" or "Kappa"
  double gammaMax() const; ///< Get the maximum gamma factor of the electron distribution function
  void circularMotion(bool t); ///< Set the flag to select the type of quantities in FITS files (GRMHD or GRPIC outputs)
  bool circularMotion() const; ///< Get the flag which select the type of quantities in FITS files
  void floorTemperature(double t); ///< Set the minimum of temperature (GRMHD case)
  double floorTemperature()const; ///< get the minimum of temperature (GRMHD case)
  void magneticConfiguration(std::string const &config); ///< Set the magnetic field configuration wanted if it is not founded in the FITS file (GRMHD case) 
  std::string magneticConfiguration() const; ///< Get the magnetic field configuration set by user.
  void magnetization(double ss); ///< Set the magnetization parameter to define magnetic field strength if not founded in FITS files (GRMHD case)
  double magnetization() const; ///< Get magnetization paramater
  void boundaryConditions(std::string const &sbc); ///< Set the boundary condition of all axis.
  std::string boundaryConditions() const; ///< Get the boundary condition array
  void emissionType(std::string const &kind); ///< Set the emission type ("BlackBody", or synchrotron "Thermal", "Power-Law" or "Kappa").
  std::string emissionType() const; ///< Get emission type

  virtual void radiativeQ(double Inu[], double Taunu[], double const nu_em[], size_t nbnu,
			  double dsem, state_t const &coord_ph, double const coord_obj[8]=NULL) const ;

  virtual void radiativeQ(double *Inu, double *Qnu, double *Unu, double *Vnu,
			  Eigen::Matrix4d *Onu, double const *nuem , size_t nbnu, double dsem,
			  state_t const &coord_ph, double const *coord_obj) const ;

  virtual void getVelocity(double const pos[4], double vel[4]);

  virtual double operator()(double const coord[4]);

  private:
  int getIndex(double const tcur) const;

};

#endif
