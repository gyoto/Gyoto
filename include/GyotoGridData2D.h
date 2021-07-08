/**
 * \file GyotoGridData2D.h
 * \brief Base class for reading 2D gridded data
 * 
 *
 */

/*
  Copyright (c) 2019-2021 Frederic Vincent, Thibaut Paumard, Nicolas Aimar
  
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

#ifndef __GyotoGridData2D_H_
#define __GyotoGridData2D_H_ 

#include <iostream>
#include <fstream>
#include <iomanip>
#ifdef GYOTO_USE_CFITSIO
#include <fitsio.h>
#endif

namespace Gyoto {
  class GridData2D;
}

#ifdef GYOTO_USE_XERCES
#include <GyotoRegister.h>
#endif

/**
 * \class Gyoto::GridData2D
 * \brief Class for reading data stored in a 2D grid
 *        
 *        
 */
class Gyoto::GridData2D
{
 private:
  double phimin_; ///< Minimum phi in grid
  double phimax_; ///< Maximum phi in grid
  double dphi_; ///< &delta;&phi; between two grid columns
  size_t nphi_; ///< Grid size in the &phi; direction
  double dr_; ///< Radius step
  size_t nr_; ///< Grid size in the r direction
  double rmin_; ///< Minimum r in grid
  double rmax_; ///< Maximum r in grid
  double dt_; ///< Time step, if not constant would be ignored
  size_t nt_; ///< Grid size in the t direction
  double tmin_; ///< Minimum t in grid
  double tmax_; ///< Maximum t in grid
  //NB: phimin, phimax are always assumed to be 0, 2pi

 public:
  GridData2D(); ///< Constructor
  GridData2D(const GridData2D&); ///< Copy constructor
  virtual GridData2D* clone() const ;
  virtual ~GridData2D() ;        ///< Destructor

  // Accessors
  void rmin(double rmn);
  double rmin() const;
  void rmax(double rmx);
  double rmax() const;
  void nr(size_t nn);
  size_t nr() const;
  void dr(double dd);
  double dr() const;
  void phimin(double phimn);
  double phimin() const;
  void phimax(double phimx);
  double phimax() const;
  void dphi(double dd);
  double dphi() const;
  void tmin(double tmn);
  double tmin() const;
  void tmax(double tmx);
  double tmax() const;
  void nt(size_t nn);
  size_t nt() const;
  void nphi(size_t nn);
  size_t nphi() const;

#ifdef GYOTO_USE_CFITSIO


  virtual std::vector<size_t> fitsReadHDU(fitsfile* fptr,
					  std::string extname,
					  double *& dest,
					  size_t length = 0);

  /**
   * \brief Creates a FITS file with dummy primary HDU
   *
   *
   * Opens a new fits file referred to by a fitsfile pointer
   * and fills the primary HDU by a single pixel equal to 0.
   * Returns the fitsfile pointer to the new FITS file.
   *
   * \param filename Name of fits file to be created
   *
   */
  fitsfile* fitsCreate(std::string filename);

  /**
   * \brief Closes a fits file referred to by a fitsfile pointer
   *
   *
   * \param fptr fitsfile pointer to FITS file to be closed
   *
   */
  void fitsClose(fitsfile* fptr);

  /**
   * \brief Writes specific HDU in FITS files
   *
   *
   * \param fptr fitsfile pointer to FITS file
   * \param extname Name of extension to be written
   * \param src Data to be written in extension
   * \param length Data has shape {nr_,nphi_,nt_} if length is 0 (default; used for storing
   *        scalar data),
   *        or {nr_,nphi_,nt_,length} if length is not 0 (used for storing vector data,
   *        e.g.: length=2 - (v_r, v_phi) - for velocity)
   *
   */
  void fitsWriteHDU(fitsfile* fptr,
		    std::string extname,
  		    double* src,
		    size_t length = 0) ;


#endif

  void getIndices(size_t i[3], double const tt, double const phi, double const rr, double* const time_array=NULL) const ;
  double interpolate(double tt, double phi, double rr,
		     double* const array, double* const time_array=NULL) const ;
    


};

#endif
