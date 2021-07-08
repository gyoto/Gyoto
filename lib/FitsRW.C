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

FitsRW::FitsRW() :
  nnu_(0), dt_(0.), nt_(0),
  numin_(0.), numax_(DBL_MAX),
  tmin_(-DBL_MAX), tmax_(DBL_MAX)
{
  GYOTO_DEBUG << endl;
}

FitsRW::FitsRW(const FitsRW&o):
  nnu_(o.nnu_), dt_(o.dt_), nt_(o.nt_),
  numin_(o.numin_), numax_(o.numax_),
  tmin_(o.tmin_), tmax_(o.tmax_)
{
  GYOTO_DEBUG << endl;
}

FitsRW* FitsRW::clone() const{
  GYOTO_DEBUG << endl;
  return new FitsRW(*this);
}

FitsRW::~FitsRW() 
{
  GYOTO_DEBUG<< endl;
}

void FitsRW::numin(double numn) {
  numin_ = numn;
}
double FitsRW::numin() const {return numin_;}

void FitsRW::numax(double numx) {
  numax_ = numx;
}
double FitsRW::numax() const {return numax_;}

void FitsRW::nnu(size_t nn) { nnu_ = nn;}
size_t FitsRW::nnu() const {return nnu_;}

void FitsRW::tmin(double tmn) {
  tmin_ = tmn;
  if (nt_>1) dt_ = (tmax_-tmin_) / double(nt_-1);
}
double FitsRW::tmin() const {return tmin_;}

void FitsRW::tmax(double tmx) {
  tmax_ = tmx;
  if (nt_>1) dt_ = (tmax_-tmin_) / double(nt_-1);
}
double FitsRW::tmax() const {return tmax_;}

void FitsRW::nt(size_t nn) { nt_ = nn;}
size_t FitsRW::nt() const {return nt_;}

void FitsRW::dt(double dd) { 
  dt_ = dd;
}
double FitsRW::dt() const {return dt_;}


#ifdef GYOTO_USE_CFITSIO
fitsfile* FitsRW::fitsCreate(string filename){
  //GYOTO_MSG << "FitsRW creating FITS file " << filename << endl; --> THIS leads to bug; CHECK; it seems that only this function, whcih doesn't returns a void, is affected; the same happens if I use cout; however cerr works:
  GYOTO_MSG << "FitsRW creating FITS file " << filename << endl;
  
  char*     pixfile   = const_cast<char*>(filename.c_str());
  int       status    = 0;
  char      ermsg[31] = ""; // ermsg is used in throwCfitsioError()
  long      fpixel[]  = {1};
  fitsfile* fptr; // pointer to FITS file
  char *    CNULL     = NULL;
  
  fits_create_file(&fptr, pixfile, &status);
  
  long      naxes [] = {1};
  // Create the primary extension containing a zero pixel
  fits_create_img(fptr, DOUBLE_IMG, 1, naxes, &status);
  double src[1] = {0.};
  fits_write_pix(fptr, TDOUBLE, fpixel, 1, src, &status);
  if (status) throwCfitsioError(status);



  // Create HDU containing the boundary of the arrays
  fits_create_img(fptr, DOUBLE_IMG, 1, naxes, &status);
  std::stringstream ss;
  ss << "GYOTO FitsRW " << "KEYS";
  fits_write_key(fptr, TSTRING, const_cast<char*>("EXTNAME"), const_cast<char*>(ss.str().c_str()), CNULL, &status);
  if (status) throwCfitsioError(status) ;

  ////// WRITE KEYWORDS (tmin_,tmax_,numin_,numax_) ///////
  if (tmin_>-DBL_MAX) {
    fits_write_key(fptr, TDOUBLE, const_cast<char*>("GYOTO FitsRW tmin"), &tmin_, CNULL, &status);
    if (status) throwCfitsioError(status) ;
  }else
    GYOTO_ERROR("In FitsRW::fitsCreate Please set tmin before calling this function");

  if (tmax_<DBL_MAX) {
    fits_write_key(fptr, TDOUBLE, const_cast<char*>("GYOTO FitsRW tmax"), &tmax_, CNULL, &status);
    if (status) throwCfitsioError(status) ;
  }else
    GYOTO_ERROR("In FitsRW::fitsCreate Please set tmax before calling this function");

  if (numin_>0.) {
    fits_write_key(fptr, TDOUBLE, const_cast<char*>("GYOTO FitsRW numin"), &numin_, CNULL, &status);
    if (status) throwCfitsioError(status) ;
  }else
    GYOTO_ERROR("In FitsRW::fitsCreate Please set numin before calling this function");

  if (numax_<DBL_MAX) {
    fits_write_key(fptr, TDOUBLE, const_cast<char*>("GYOTO FitsRW numax"), &numax_, CNULL, &status); 
    if (status) throwCfitsioError(status) ;
  }else
    GYOTO_ERROR("In FitsRW::fitsCreate Please set numax before calling this function");

  return fptr;
}

void FitsRW::fitsClose(fitsfile* fptr) {
  GYOTO_MSG << "FitsRW Closing FITS file" << endl;

  int       status    = 0;
  char      ermsg[31] = ""; // ermsg is used in throwCfitsioError()
  int res=0;
  if (fits_close_file(fptr, &status)) throwCfitsioError(status) ;
  fptr = NULL;
}

void FitsRW::fitsWriteParams(fitsfile* fptr, double n_e, double theta, double kappa, double BB, double t_inj){
  int       status    = 0;
  char *    CNULL     = NULL;
  char      ermsg[31] = ""; // ermsg is used in throwCfitsioError()
  string    ss        = "GYOTO FitsRW KEYS";


  fits_movnam_hdu(fptr, ANY_HDU, const_cast<char*>(ss.c_str()), 0, &status);

  if (n_e>0.){
    fits_write_key(fptr, TDOUBLE, const_cast<char*>("GYOTO FitsRW ne"), &n_e, CNULL, &status);
    if (status) throwCfitsioError(status) ;
  }else
    GYOTO_ERROR("In fitsWriteParams : n_e<0");

  if (theta>0.){
    fits_write_key(fptr, TDOUBLE, const_cast<char*>("GYOTO FitsRW theta"), &theta, CNULL, &status);
    if (status) throwCfitsioError(status) ;
  }else
    GYOTO_ERROR("In fitsWriteParams : theta<0");

  if (kappa>0.){
    fits_write_key(fptr, TDOUBLE, const_cast<char*>("GYOTO FitsRW kappa"), &kappa, CNULL, &status);
    if (status) throwCfitsioError(status) ;
  }else
    GYOTO_ERROR("In fitsWriteParams : kappa<0");

  if (BB>0.){
    fits_write_key(fptr, TDOUBLE, const_cast<char*>("GYOTO FitsRW BB"), &BB, CNULL, &status);
    if (status) throwCfitsioError(status) ;
  }else
    GYOTO_ERROR("In fitsWriteParams : BB<0");

  if (t_inj>0.){
    fits_write_key(fptr, TDOUBLE, const_cast<char*>("GYOTO FitsRW t_inj"), &t_inj, CNULL, &status);
    if (status) throwCfitsioError(status) ;
  }else
    GYOTO_ERROR("In fitsWriteParams : t_inj<0");
}

void FitsRW::fitsWriteHDU(fitsfile* fptr,
            string const extname,
            double* const src) {
  GYOTO_MSG << "FitsRW writing HDU " << extname << endl;
  if (!src) GYOTO_ERROR("FitsRW::fitsWrite: nothing to save!");
  int       status    = 0;
  long      fpixel [] = {1,1};
  long      ndim      = 2;
  long      naxes []  = {long(nnu_),long(nt_)};
  char      ermsg[31] = ""; // ermsg is used in throwCfitsioError()
  char *    CNULL     = NULL;

  if (nt_==0 || nnu_==0)
    GYOTO_ERROR("Axes lengths should be defined");

  //cout << "tmin, tmax, numin, numax= " << tmin_ << "," << tmax_ << "," << numin_ << "," << numax_ << endl;

  ////// SAVE SRC IN APPROPRIATE HDU ///////
  if (extname!="FREQUENCY"){
    fits_create_img(fptr, DOUBLE_IMG, ndim, naxes, &status);
    std::stringstream ss;
    ss << "GYOTO FitsRW " << extname;
    fits_write_key(fptr, TSTRING, const_cast<char*>("EXTNAME"),
       const_cast<char*>(ss.str().c_str()),
       CNULL, &status);
    fits_write_pix(fptr, TDOUBLE, fpixel, nnu_*nt_, src, &status);
    if (status) throwCfitsioError(status) ;
  }else{
    long naxe [] = {long(nnu_)}, fpix []={1};
    fits_create_img(fptr, DOUBLE_IMG, 1, naxe, &status);
    std::stringstream ss;
    ss << "GYOTO FitsRW " << extname;
    fits_write_key(fptr, TSTRING, const_cast<char*>("EXTNAME"),
     const_cast<char*>(ss.str().c_str()),
     CNULL, &status);
    fits_write_pix(fptr,TDOUBLE, fpix, nnu_, src, &status);
    if (status) throwCfitsioError(status) ;
  }

}

vector<size_t> FitsRW::fitsReadHDU(fitsfile* fptr,
               string extname,
               double *& dest) {
  //GYOTO_MSG << "FitsRW reading FITS extension " << extname << endl;

  int       status    = 0;
  int       anynul    = 0;
  double    tmpd;
  long      ndim = 2;
  long      naxes []  = {1,1};
  long      fpixel[]  = {1,1};
  long      inc   []  = {1,1};
  char      ermsg[31] = ""; // ermsg is used in throwCfitsioError()
  
  ////// READ REQUIRED EXTENSION ///////
  if (extname!="GYOTO FitsRW FREQUENCY"){
    GYOTO_DEBUG << "FitsRW::fitsRead(): search " << extname << " HDU" << endl;
    if (fits_movnam_hdu(fptr, ANY_HDU,
              const_cast<char*>(extname.c_str()),
              0, &status))
      throwCfitsioError(status) ;
    GYOTO_DEBUG << "FitsRW::fitsRead(): get image size" << endl;
    if (fits_get_img_size(fptr, ndim, naxes, &status)) throwCfitsioError(status) ;
    // update nt_, dt_
    nt_ = naxes[1];
    if (nt_>1) dt_ = (tmax_-tmin_)/double((nt_-1));

    // update nnu_
    nnu_ = naxes[0];

    if (dest) { delete [] dest; dest = NULL; }
    dest = new double[nnu_*nt_];
    for (int ii=0;ii<nnu_*nt_;ii++)
      dest[ii]=0.;
    if (fits_read_subset(fptr, TDOUBLE, fpixel, naxes, inc,
               0,dest,&anynul,&status)) {
      GYOTO_DEBUG << " error, trying to free pointer" << endl;
      delete [] dest; dest=NULL;
      throwCfitsioError(status) ;
    }
    GYOTO_DEBUG << " done." << endl;

    vector<size_t> dims(ndim, nnu_);
    dims[1] = nt_;

    return dims;
  }
  else{ // Read FREQUENCY
    long naxe[] = {1}, inc[]={1};
    GYOTO_DEBUG << "FitsRW::fitsRead(): search " << extname << " HDU" << endl;
    if (fits_movnam_hdu(fptr, ANY_HDU,
              const_cast<char*>(extname.c_str()),
              0, &status))
      throwCfitsioError(status) ;
    GYOTO_DEBUG << "FitsRW::fitsRead(): get image size" << endl;
    if (fits_get_img_size(fptr, 1, naxe, &status)) throwCfitsioError(status) ;
    nt_ = naxe[0];

    if (dest) { delete [] dest; dest = NULL; }
    dest = new double[nt_];
    for (int ii=0;ii<nnu_;ii++)
      dest[ii]=0.;
    if (fits_read_subset(fptr, TDOUBLE, fpixel, naxe, inc,
               0,dest,&anynul,&status)) {
      GYOTO_DEBUG << " error, trying to free pointer" << endl;
      delete [] dest; dest=NULL;
      throwCfitsioError(status) ;
    }
    GYOTO_DEBUG << " done." << endl;

    vector<size_t> dims(1, nnu_);
    return dims;
  }
}
#endif

void FitsRW::getIndices(size_t i[2], double const nu, double const tt, double* const freq_array) const {

  //cout << "tmin, tmax, tt, numin, numax, nu= " << tmin_ << "," << tmax_ << "," << tt << "," << numin_ << "," << numax_ << "," << nu << endl;
  if (tmin_>-DBL_MAX && tmax_<DBL_MAX && nt_>0. && dt_>0.) { // >1 time must be properly defined
    if (tt<tmin_) i[1]=0;         // assuming stationnarity before tmin_
    else if (tt>tmax_) i[1]=nt_-1;//                       and after tmax_
    else
      i[1] = size_t(floor((tt-tmin_)/dt_)); // index of closest grid point smaller than tt
  }else
    GYOTO_ERROR("In FitsRW::getIndices: time undefined!");

  if (nnu_>0.) {
    if (numin_ > 0. && numax_<DBL_MAX){
      if (numax_>nu && numin_<nu){
          i[0]=0;
        while (freq_array[i[0]+1]<nu){
          i[0]+=1;
        }
      }else
        GYOTO_ERROR("In FitsRW::getIndices: frequency out of boundary!");
    }else
      GYOTO_ERROR("In FitsRW::getIndices: frequencies badly defined!");
  } else {
    GYOTO_ERROR("In FitsRW::getIndices: frequencies undefined!");
  }

}

double FitsRW::interpolate(double nu, double tt, double* const array, double* const freq_array) const{
  if (!freq_array)
    GYOTO_ERROR("In FitsRW::interpolate freq_array not defined");

  size_t ind[2]; // {i_t, i_nu}
  getIndices(ind, nu , tt, freq_array);

  //cout << ind[0] << " " << ind[1] << endl;
  //cout << " TEST First slot= " << array[0] << " " << array[1] << " " << array[2] << endl;

  double array_interpo=0.;

  if (tt>tmax_)
    GYOTO_ERROR("In FitsRW::interpolate t>tmax");

  // Bilinear in nu,tt
  size_t inul=ind[0], inuu=ind[0]+1, itl=ind[1], itu=ind[1]+1;
  double array_ll = array[itl*nnu_+inul], // array_{nu,tt}, l=lower index, u=upper index
    array_lu = array[itl*nnu_+inuu],
    array_ul = array[itu*nnu_+inul],
    array_uu = array[itu*nnu_+inuu];
  double tl = tmin_+dt_*itl, tu = tmin_+dt_*itu,
    nul=freq_array[inul], nuu=freq_array[inuu];
  double rationu = (nu-nul)/(nuu-nul),
    ratiot = (tt-tl)/(tu-tl);
  array_interpo = array_ll+(array_ul-array_ll)*ratiot
    +(array_lu-array_ll)*rationu
    +(array_uu-array_lu-array_ul+array_ll)*rationu*ratiot;
  /*cout << "bilin interpo stuff: " << endl;
  cout << "T: " << tl << " " << tt << " " << tu << endl;
  cout << "NU: " << nul << " " << nu << " " << nuu << endl;
  cout << "ARRAY at 4 corners (tt,nu)=(ll,lu,ul,uu) + interpo: " << array_ll << " " << array_lu << " " << array_ul << " " << array_uu << " " << array_interpo << endl;*/
  return array_interpo;
}