/**
 * \file GyotoDynamicalDisk3D.h
 * \brief A geometrically thick, optically thin disk, evolving dynamically,
 *  with black body emission. 
 *
 *  The disk is described by a set of FITS files for a set of different times
 */

/*
    Copyright 2011-2016, 2018 Frederic Vincent, Thibaut Paumard

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

#ifndef __GyotoDynamicalDisk3D_H_ 
#define __GyotoDynamicalDisk3D_H_ 

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstring>

namespace Gyoto{
  namespace Astrobj { class DynamicalDisk3D; }
}

#include <GyotoDisk3D.h>
#include <GyotoBlackBodySpectrum.h>

/**
 * \class Gyoto::Astrobj::DynamicalDisk3D
 * \brief Geometrically thick optically thin disk 
 *  read from a set of FITS files. 
 * 
 *   This class describes a PatternDiskBB that evolves dynamically. 
 *   It is described by a set of FITS files for different times.
 *   Its emission is blackbody.
 *
 *   The disk is assumed to be described by a regular, 
 *   constant in time, grid.
 *  
 *   The metric must be Kerr in BL coordinates.
 *
 */
class Gyoto::Astrobj::DynamicalDisk3D : public Astrobj::Disk3D {
  friend class Gyoto::SmartPointer<Gyoto::Astrobj::DynamicalDisk3D>;
 protected:
  SmartPointer<Spectrum::BlackBody> spectrumBB_; ///< disk black body
  ///< emission law
 private:
  bool temperature_; ///< 1 if temperature is given in fits data file, 0 if emission coef is directly given
  char* dirname_; ///< FITS files directory
  double tinit_; ///< Time of the first FITS file
  double dt_; ///< Time increment between two FITS (assumed constant)
  int nb_times_; ///< Number of times
  double PLindex_; ///< power law index such that density_elec(E) &prop; E<SUP>-p</SUP>
  bool novel_; ///< put to true if velocity of emitting particle is not provided
  double floortemperature_; ///< if non-zero, emission and absorption are 0 for temperatures below this floor, emission=blackbody and absorption is infty for temperatures above (this is a kind of fake optically thick case, when the emitting surface is inside the grid, not at the boundary of the grid)

  /**
   * An array of arrays of dimensionality double[nr_][nz_][nphi_][nnu_]. 
   * In FITS format, the first dimension is nu, the second phi, the third
   * z, the last r. It contains emission coef, or a closely related quantity.
   */
  double ** emission_array_;

  /**
   * An array of arrays of dimensionality double[nr_][nz_][nphi_][nnu_]. 
   * In FITS format, the first dimension is nu, the second phi, the third
   * z, the last r. It contains absorption coef, or a closely related quantity.
   */
  double ** absorption_array_;

  /**
   * An array of arrays of dimensionality double[nr_][nz_][nphi_][3].
   * In FITS format, the second dimension is phi, and the third r. 
   * The first plane in the first FITS dimention is dphi/dt, 
   * the second dz/dt, the third dr/dt.
   */
  double ** velocity_array_; ///< velocity(r, z, phi)

  // Constructors - Destructor
  // -------------------------
 public:
  GYOTO_OBJECT;
  GYOTO_OBJECT_THREAD_SAFETY;

  DynamicalDisk3D(); ///< Standard constructor
  
  DynamicalDisk3D(const DynamicalDisk3D& ) ;///< Copy constructor
  virtual DynamicalDisk3D* clone () const; ///< Cloner
  
  virtual ~DynamicalDisk3D() ;                        ///< Destructor
  
  // Accessors
  // ---------
 public:
  using Generic::metric;
  void metric(SmartPointer<Metric::Generic> gg);
  void file(std::string const &f);
  std::string file() const;
  void tinit(double t);
  double tinit()const;
  void dt(double t);
  double dt()const;
  void PLindex(double t);
  double PLindex()const;
  void floorTemperature(double t);
  double floorTemperature()const;
  void temperature(bool t);
  bool temperature() const;
  void withVelocity(bool t);
  bool withVelocity() const;

  // Stuff
  // -----
  /// Compute emission at one grid date.
  double emission1date(double nu_em, double dsem,
		  state_t const &c_ph, double const c_obj[8]) const;

  using Disk3D::emission;
  /// Interpolate emission between grid dates.
  virtual double emission(double nu_em, double dsem,
			  state_t const &c_ph, double const c_obj[8]=NULL) const;

  /// Compute transmission at one grid date.
  double transmission1date(double nu_em, double dsem,
		  state_t const &c_ph, double const c_obj[8]) const;

  /// Interpolate transmission between grid dates.
  double transmission(double nu_em, double dsem,
			  state_t const &cph, double const co[8]) const;

  void getVelocity(double const pos[4], double vel[4]);
  double const * getVelocity() const;
  
 protected:

  /// Set underlying Disk3D pointers to a specific date slice.
  /**
   * \param iq Index of the date slice.
   */
  void copyQuantities(int iq) ;

};

#endif
