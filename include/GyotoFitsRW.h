/**
 * \file GyotoFitsRW.h
 * \brief Class to read/write jnu and anu in FITS File
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

#ifndef __GyotoFitsRW_H_
#define __GyotoFitsRW_H_ 

#include <iostream>
#include <fstream>
#include <iomanip>
#ifdef GYOTO_USE_CFITSIO
#include <fitsio.h>
#endif

namespace Gyoto {
  class FitsRW;
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
class Gyoto::FitsRW
{
  private:

  size_t nnu_; //number of frequency lines in FITS File
  double numin_;
  double numax_;
  size_t nt_; //number of time steps in FITS File
  double dt_;
  double tmin_;
  double tmax_;

 public:
  FitsRW(); ///< Constructor
  FitsRW(const FitsRW&); ///< Copy constructor
  virtual FitsRW* clone() const ;
  virtual ~FitsRW() ;        ///< Destructor

  // Accessors
  void numin(double numn);
  double numin() const;
  void numax(double numx);
  double numax() const;
  void nnu(size_t nn);
  size_t nnu() const;
  void tmin(double tmn);
  double tmin() const;
  void tmax(double tmx);
  double tmax() const;
  void nt(size_t nn);
  size_t nt() const;
  void dt(double dd);
  double dt() const;

#ifdef GYOTO_USE_CFITSIO

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
   * 
   * \Brief Data has shape {nnu_,nt_}
   *
   */
  void fitsWriteHDU(fitsfile* fptr,
        std::string extname,
          double* src);

  void fitsWriteParams(fitsfile* fptr, double n_e, double theta, double kappa, double BB, double t_inj);

  virtual std::vector<size_t> fitsReadHDU(fitsfile* fptr,
            std::string extname,
            double *& dest);

  #endif

  void getIndices(size_t i[2], double const nu, double const tt, double* const freq_array) const ;
  double interpolate(double nu, double tt, double* const array, double* const freq_array) const ;

};


#endif