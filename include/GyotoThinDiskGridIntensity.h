/**
 * \file GyotoThinDiskGridIntensity.h
 * \brief A disk defined from a 2D grid in the equatorial plane
 * and extrapolated in the vertical direction with H/r<<1
 */

/*
    Copyright 2019-2021 Frederic Vincent, Thibaut Paumard

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

#ifndef __GyotoThinDiskGridIntensity_H_ 
#define __GyotoThinDiskGridIntensity_H_ 

#include <iostream>
#include <fstream>
#include <iomanip>

namespace Gyoto{
  namespace Astrobj { class ThinDiskGridIntensity; }
  class GridData2D;
}

#include <GyotoThinDisk.h>
#include <GyotoGridData2D.h>

/**
 * \class Gyoto::Astrobj::ThinDiskGridIntensity
 */

class Gyoto::Astrobj::ThinDiskGridIntensity
: public Astrobj::ThinDisk,
  public GridData2D,
  public Hook::Listener
{
  friend class Gyoto::SmartPointer<Gyoto::Astrobj::ThinDiskGridIntensity>;
 private:
  std::string filename_; ///< Optional FITS file name containing the arrays
  /**
   * An array of dimensionality double[nr_][nphi_][nt_]. In FITS
   * format, the first dimension is t, the second phi, and the third
   * r.
   */
  double * intensity_; ///< Intensity (&nu;, r, &phi;)
  double * time_array_; /// 1D Vector containing the times values of each time steps (dt not constant)
  double deltat_;///< time translation

 public:
  GYOTO_OBJECT;
  GYOTO_OBJECT_THREAD_SAFETY;
  
  // Constructors - Destructor
  // -------------------------
  ThinDiskGridIntensity(); ///< Standard constructor
  
  ThinDiskGridIntensity(const ThinDiskGridIntensity& ) ;///< Copy constructor
  virtual ThinDiskGridIntensity* clone () const; ///< Cloner
  
  virtual ~ThinDiskGridIntensity() ;                        ///< Destructor
  
  // Accessors
  // ---------
 public:
  void file(std::string const &f) ;
  std::string file() const;
  /*
    timeTranslation shifts the value of GridData2D::tmin_ and tmax_,
    allowing to scan the full simulation without having to change
    the value of the Screen observation time (which is typically
    not provided in M unit in the XML). 
    Choosing a negative timeTranslation, i.e. performing tmin_,tmax_-=dt, 
    amounts to increasing the Screen observation time by the same value, 
    tobs+=dt.

   */
  void timeTranslation_inMunit(double const dt) ;
  double timeTranslation_inMunit() const ;
  //void dt(double dd);
  void copyIntensity(double const *const intensity,
		     size_t const naxes[3]);
  double const * getIntensity() const;
  void copyTimeArray(double const *const time_array, size_t const ntimes);
  double const * getTimeArray() const;
 public:
  using Generic::metric;
  std::vector<size_t> fitsRead(std::string filename) ;
  virtual double emission(double nu_em, double dsem,
			  state_t const &c_ph,double const c_obj[8]=NULL) const;
  virtual void getVelocity(double const pos[4], double vel[4]) ;



};

#endif
