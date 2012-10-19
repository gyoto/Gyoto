/**
 * \file GyotoDynamicalDisk.h
 * \brief A geometrically thin, optically thick disk, evolving dynamically
 *
 *  The disk is described by a set of FITS files for a set of different times
 */

/*
    Copyright 2011 Frederic Vincent, Thibaut Paumard

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
  double tinit_; ///< Time of the first FITS file
  double dt_; ///< Time increment between two FITS (assumed constant)
  int nb_times_; ///< Number of times

  /**
   * An array of dimensionality double[nr_][nphi_][nnu_]. In FITS
   * format, the first dimension is nu, the second phi, and the third
   * r.
   */
  double ** emission_array_;

  double ** opacity_array_; ///< same dimenstions as emission, or NULL

  /**
   * An array of dimensionality double[nr_][nphi_][2]. In FITS format,
   * the second dimension is phi, and the third r. The first plane in
   * the first FITS dimention is dphi/dt, the second dr/dt.
   */
  double ** velocity_array_; ///< velocity(, r, phi)

  /**
   * In case of adaptive grid.
   */
  double ** radius_array_; ///< radius vector

  double * dnu_array_;
  double * nu0_array_;
  size_t * nnu_array_;

  double * dphi_array_;
  size_t * nphi_array_;

  double * dr_array_;
  size_t * nr_array_;
  //  double r0_; // this is rin_

  // Constructors - Destructor
  // -------------------------
 public:
  DynamicalDisk(); ///< Standard constructor
  
  DynamicalDisk(const DynamicalDisk& ) ;///< Copy constructor
  virtual DynamicalDisk* clone () const; ///< Cloner
  
  virtual ~DynamicalDisk() ;                        ///< Destructor
  
  // Accessors
  // ---------
 public:

  virtual int setParameter(std::string name,
			   std::string content,
			   std::string unit);

  virtual double emission(double nu_em, double dsem,
			  double c_ph[8], double c_obj[8]) const;

  void getVelocity(double const pos[4], double vel[4]);
  double const * const getVelocity() const;
  
 protected:

  void copyQuantities(int iq) ;

 public:
#ifdef GYOTO_USE_XERCES
  virtual void fillElement(FactoryMessenger *fmp) const ;
#endif

};

#endif
