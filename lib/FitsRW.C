/*
    Copyright 2019-2021 Frederic Vincent, Thibaut Paumard, Nicolas Aimar

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

//Gyoto headers
#include "GyotoUtils.h"
#include "GyotoFitsRW.h"

#ifdef GYOTO_USE_CFITSIO
#define throwCfitsioError(status) \
    { fits_get_errstatus(status, ermsg); GYOTO_ERROR(ermsg); }
#endif

//Std headers
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <dirent.h>
#include <ctime>

using namespace std;
using namespace Gyoto;

FitsRW::FitsRW() {}

FitsRW::FitsRW(const FitsRW&o) {}

FitsRW* FitsRW::clone() const{
  return new FitsRW(*this);
}

FitsRW::~FitsRW() {}


#ifdef GYOTO_USE_CFITSIO
fitsfile* FitsRW::fitsCreate(string const filename) const{
  GYOTO_MSG << "FitsRW creating FITS file " << filename << endl;
  
  char*     fitsname   = const_cast<char*>(filename.c_str());
  int       status    = 0;
  char      ermsg[31] = ""; // ermsg is used in throwCfitsioError()
  long      fpixel[2]  = {1,1};
  fitsfile* fptr; // pointer to FITS file
  char *    CNULL     = NULL;
  
  if (remove(fitsname) == 0) {
    cout << "Existing file removed" << endl;
  } else {
    cout << "No existing file to remove" << endl;
  }
  
  fits_create_file(&fptr, fitsname, &status);
  long      naxes [2] = {1, 1};
  // Create the primary extension containing a zero pixel
  fits_create_img(fptr, DOUBLE_IMG, 1, naxes, &status);
  if (status) throwCfitsioError(status);
  short pixel_value = 0;
  // Write the pixel value to the image
  fits_write_pix(fptr, TSHORT, fpixel, 1, &pixel_value, &status);
  if (status) throwCfitsioError(status);

  return fptr;
}

fitsfile* FitsRW::fitsOpen(string const filename) const{
  GYOTO_DEBUG << "FitsRW opening FITS file " << filename << endl;

  fitsfile* fptr; // pointer to FITS file
  int       status    = 0;
  char      ermsg[31] = ""; // ermsg is used in throwCfitsioError()
  char*     pixfile   = const_cast<char*>(filename.c_str());

  fits_open_file(&fptr, pixfile, 0, &status);
  if (status) throwCfitsioError(status);

  return fptr;
}

void FitsRW::fitsClose(fitsfile* fptr) const{
  //GYOTO_MSG << "FitsRW Closing FITS file" << endl;

  int       status    = 0;
  char      ermsg[31] = ""; // ermsg is used in throwCfitsioError()
  int res=0;
  if (fits_close_file(fptr, &status)) throwCfitsioError(status) ;
  fptr = NULL;
}

void FitsRW::fitsWriteKey(fitsfile* const fptr, std::string const key, double value, std::string const hdu) const{
  int       status    = 0;
  char *    CNULL     = NULL;
  char      ermsg[31] = ""; // ermsg is used in throwCfitsioError()

  if (hdu=="PRIMARY"){
    fits_movabs_hdu(fptr, 1, NULL, &status);
  }else{
    fits_movnam_hdu(fptr, ANY_HDU, const_cast<char*>(hdu.c_str()), 0, &status);
  }
  fits_update_key(fptr, TDOUBLE, const_cast<char*>(key.c_str()), &value, CNULL, &status); // Updates if key already exists, appends if not
  if (status) throwCfitsioError(status) ;
}

void FitsRW::fitsWriteHDUData(fitsfile* const fptr, string const extname, double* const src, long const nelements) const{
  GYOTO_MSG << "FitsRW writing HDU " << extname << endl;
  if (!src) GYOTO_ERROR("FitsRW::fitsWrite: nothing to save!");
  int       status    = 0;
  int const ndim      = 1;
  long      naxes[ndim] = {nelements};
  long      fpixel[1] = {1};
  char      ermsg[31] = ""; // ermsg is used in throwCfitsioError()
  char *    CNULL     = NULL;

  ////// SAVE SRC IN APPROPRIATE HDU ///////
  fits_create_img(fptr, DOUBLE_IMG, ndim, naxes, &status);
  if (status) throwCfitsioError(status) ;
  fits_write_key(fptr, TSTRING, const_cast<char*>("EXTNAME"), const_cast<char*>(extname.c_str()), CNULL, &status);
  if (status) throwCfitsioError(status) ;
  fits_write_pix(fptr, TDOUBLE, fpixel, nelements, src, &status);
  if (status) throwCfitsioError(status) ;

}

double* FitsRW::fitsReadHDUData(fitsfile* const fptr, string const extname) const{
  int       status    = 0;
  int       anynul    = 0;
  int       ndim      = 1;
  long      naxes[ndim];
  long      fpixel[ndim]; fpixel[0] = {1};
  long      inc[ndim];    inc[0] = {1};
  char      ermsg[31] = ""; // ermsg is used in throwCfitsioError()
  
  ////// READ REQUIRED EXTENSION ///////
  GYOTO_DEBUG << "FitsRW::fitsReadHDUData(): search " << extname << " HDU" << endl;
  fits_movnam_hdu(fptr, ANY_HDU, const_cast<char*>(extname.c_str()), 0, &status);
  if (status) throwCfitsioError(status) ;

  GYOTO_DEBUG << "FitsRW::fitsReadHDUData(): get image size" << endl;
  fits_get_img_dim(fptr, &ndim, &status);
  if (ndim!=1)
    GYOTO_ERROR("FitsRW::fitsReadHDUData(): array dimension have to be 1.");
  fits_get_img_size(fptr, ndim, naxes, &status);
  if (status) throwCfitsioError(status) ;

  GYOTO_DEBUG << "FitsRW::fitsReadHDUData(): allocation array" << endl;
  double* dest = new double[naxes[0]];
  fits_read_subset(fptr, TDOUBLE, fpixel, naxes, inc, 0, dest, &anynul, &status);
  if (status) {
    GYOTO_DEBUG << " error reading array, freeing memory" << endl;
    delete [] dest; dest=NULL;
    throwCfitsioError(status) ;
  }
  GYOTO_DEBUG << " done." << endl;

  return dest;
}

double FitsRW::fitsReadKey(fitsfile* const fptr, string const key, std::string const hdu) const{
  int       status    = 0;
  double    tmpd;
  char      ermsg[31] = ""; // ermsg is used in throwCfitsioError()


  fits_movnam_hdu(fptr, ANY_HDU, const_cast<char*>(hdu.c_str()), 0, &status);
  if (status) throwCfitsioError(status) ;

  GYOTO_DEBUG << "FitsRW::fitsReadKey(): reading key " << key << endl;
  fits_read_key(fptr, TDOUBLE, const_cast<char*>(key.c_str()), &tmpd, NULL, &status);
  if (status) throwCfitsioError(status) ;

  return tmpd;

}

double FitsRW::fitsReadKey(fitsfile* const fptr, string const key, int const hdu_num) const{
  int       status    = 0;
  double    tmpd;
  char      ermsg[31] = ""; // ermsg is used in throwCfitsioError()
  int*      tmp;

  fits_movabs_hdu(fptr, hdu_num, tmp, &status);
  if (status) throwCfitsioError(status) ;

  GYOTO_DEBUG << "FitsRW::fitsReadKey(): reading key " << key << endl;
  fits_read_key(fptr, TDOUBLE, const_cast<char*>(key.c_str()), &tmpd, NULL, &status);
  if (status) throwCfitsioError(status) ;

  return tmpd;

}
#endif
