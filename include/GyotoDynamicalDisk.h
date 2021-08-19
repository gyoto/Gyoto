/**
 * \file GyotoDynamicalDisk.h
 * \brief A geometrically thin, optically thick disk, evolving dynamically
 *
 *  The disk is described by a set of FITS files for a set of different times
 */

/*
    Copyright 2011-2015, 2018 Frederic Vincent, Thibaut Paumard

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

#ifndef __GyotoDynamicalDisk_H_ 
#define __GyotoDynamicalDisk_H_ 

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstring>

namespace Gyoto{
  namespace Astrobj { class DynamicalDisk; }
}

//#include <GyotoMetric.h>
#include <GyotoPatternDiskBB.h>

/**
 * \class Gyoto::Astrobj::DynamicalDisk
 * \brief Geometrically thin disk read from a set of FITS files
 * 
 *   This class describes a PatternDiskBB that evolves dynamically. 
 *   It is described by a set of FITS files.
 *
 */
class Gyoto::Astrobj::DynamicalDisk : public Astrobj::PatternDiskBB {
  friend class Gyoto::SmartPointer<Gyoto::Astrobj::DynamicalDisk>;
 private:
  char* dirname_; ///< FITS files directory
  double tinit_; ///< date of the first FITS file
  double dt_; ///< Time increment between two FITS (assumed constant)
  int nb_times_; ///< Number of dates
  int nnu_, nphi_, nr_; ///< Grid dimensions (assumed constant)

  /// Array of PatternDisk::emission_ arrays
  double ** emission_array_;

  /// Array of PatternDisk::velocity_ arrays
  double ** velocity_array_; ///< 

  /// Array of PatternDisk::radius_ arrays
  double ** radius_array_; ///< radius vector

  // Constructors - Destructor
  // -------------------------
 public:
  GYOTO_OBJECT;

  DynamicalDisk(); ///< Standard constructor
  
  DynamicalDisk(const DynamicalDisk& ) ;///< Copy constructor
  virtual DynamicalDisk* clone () const; ///< Cloner
  
  virtual ~DynamicalDisk() ;                        ///< Destructor
  
  // Accessors
  // ---------
 public:

  std::string file() const;
  void file(std::string const &fname);
  void tinit(double t);
  double tinit()const;
  void dt(double t);
  double dt()const;

  // fillProperty is overridden to set only the directory 
  void fillProperty(Gyoto::FactoryMessenger *fmp, Property const &p) const;

  using PatternDiskBB::emission;
  virtual double emission(double nu_em, double dsem,
			  state_t const &c_ph, double const c_obj[8]=NULL) const;

  void getVelocity(double const pos[4], double vel[4]);
  double const * getVelocity() const;
  
 protected:

  /// Set underlying PatternDisk pointers to a specific date slice.
  /**
   * \param iq Index of the date slice.
   */
  void copyQuantities(int iq) ;

  void nullifyQuantities() ;


};

#endif
