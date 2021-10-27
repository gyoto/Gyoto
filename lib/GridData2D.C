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
#include "GyotoGridData2D.h"

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

GridData2D::GridData2D() :
  dphi_(0.), nphi_(0), dr_(0.), nr_(0.), dt_(0.), nt_(0),
  rmin_(0.), rmax_(DBL_MAX),
  phimin_(0.), phimax_(2.*M_PI),
  tmin_(-DBL_MAX), tmax_(DBL_MAX)
{
  GYOTO_DEBUG << endl;
}

GridData2D::GridData2D(const GridData2D&o):
  dphi_(o.dphi_), nphi_(o.nphi_), dr_(o.dr_), nr_(o.nr_), dt_(o.dt_), nt_(o.nt_),
  rmin_(o.rmin_), rmax_(o.rmax_),
  phimin_(o.phimin_), phimax_(o.phimax_),
  tmin_(o.tmin_), tmax_(o.tmax_)
{
  GYOTO_DEBUG << endl;
}

GridData2D* GridData2D::clone() const{
  GYOTO_DEBUG << endl;
  return new GridData2D(*this);
}

GridData2D::~GridData2D() 
{
  GYOTO_DEBUG<< endl;
}

void GridData2D::rmin(double rmn) {
  rmin_ = rmn;
  if (nr_>1) dr_ = (rmax_-rmin_) / double(nr_-1);
}
double GridData2D::rmin() const {return rmin_;}

void GridData2D::rmax(double rmx) {
  rmax_ = rmx;
  if (nr_>1) dr_ = (rmax_-rmin_) / double(nr_-1);
}
double GridData2D::rmax() const {return rmax_;}

void GridData2D::nr(size_t nn) { nr_ = nn;}
size_t GridData2D::nr() const {return nr_;}

void GridData2D::dr(double dd) { dr_ = dd;}
double GridData2D::dr() const {return dr_;}

void GridData2D::phimin(double phimn) {
  if (phimn<0. or phimn>2.*M_PI)
    throwError("In GridData2D::phimin: bad phimin");
  phimin_ = phimn;
  if (nphi_>1) dphi_ = (phimax_-phimin_) / double(nphi_-1);
}
double GridData2D::phimin() const {return phimin_;}

void GridData2D::phimax(double phimx) {
  if (phimx<0. or phimx>2.*M_PI)
    throwError("In GridData2D::phimax: bad phimax");
  phimax_ = phimx;
  if (nphi_>1) dphi_ = (phimax_-phimin_) / double(nphi_-1);
}
double GridData2D::phimax() const {return phimax_;}

void GridData2D::dphi(double dd) { dphi_ = dd;}
double GridData2D::dphi() const {return dphi_;}

void GridData2D::tmin(double tmn) {
  tmin_ = tmn;
  if (nt_>1) dt_ = (tmax_-tmin_) / double(nt_-1);
}
double GridData2D::tmin() const {return tmin_;}

void GridData2D::tmax(double tmx) {
  tmax_ = tmx;
  if (nt_>1) dt_ = (tmax_-tmin_) / double(nt_-1);
}
double GridData2D::tmax() const {return tmax_;}

void GridData2D::nt(size_t nn) { nt_ = nn;}
size_t GridData2D::nt() const {return nt_;}

void GridData2D::nphi(size_t nn) { nphi_ = nn;}
size_t GridData2D::nphi() const {return nphi_;}


#ifdef GYOTO_USE_CFITSIO
vector<size_t> GridData2D::fitsReadHDU(fitsfile* fptr,
				       string extname,
				       double *& dest,
				       size_t length) {
  GYOTO_MSG << "GridData2D reading FITS extension " << extname << endl;

  int       status    = 0;
  int       anynul    = 0;
  double    tmpd;
  long      ndim = length?4:3;
  long      naxes []  = {1,1,1,1};
  long      fpixel[]  = {1,1,1,1};
  long      inc   []  = {1,1,1,1};
  char      ermsg[31] = ""; // ermsg is used in throwCfitsioError()

  ////// READ REQUIRED EXTENSION ///////
  if (extname!="GYOTO GridData2D TIMEARRAY"){
    GYOTO_DEBUG << "GridData2D::fitsRead(): search " << extname << " HDU" << endl;
    if (fits_movnam_hdu(fptr, ANY_HDU,
    		      const_cast<char*>(extname.c_str()),
    		      0, &status))
      throwCfitsioError(status) ;
    GYOTO_DEBUG << "GridData2D::fitsRead(): get image size" << endl;
    if (fits_get_img_size(fptr, ndim, naxes, &status)) throwCfitsioError(status) ;
    // update nt_, dt_
    nt_ = naxes[2];
    if (nt_>1) dt_ = (tmax_-tmin_)/double((nt_-1));
    
    // update nphi_, dphi_
    nphi_ = naxes[1];
    if (nphi_>1) dphi_ = (phimax_ - phimin_)/double((nphi_-1));

    // update nr_, dr_
    nr_ = naxes[0];
    if (nr_>1) dr_ = (rmax_-rmin_) / double(nr_-1);

    if (dest) { delete [] dest; dest = NULL; }
    dest = new double[nt_ * nphi_ * nr_ * (length?length:1)];
    for (int ii=0;ii<nt_ * nphi_ * nr_ * (length?length:1);ii++)
      dest[ii]=0.;
    if (fits_read_subset(fptr, TDOUBLE, fpixel, naxes, inc,
    		       0,dest,&anynul,&status)) {
      GYOTO_DEBUG << " error, trying to free pointer" << endl;
      delete [] dest; dest=NULL;
      throwCfitsioError(status) ;
    }
    GYOTO_DEBUG << " done." << endl;

    vector<size_t> dims(ndim, nr_);
    dims[1] = nphi_;
    dims[2] = nt_;
    if (length) dims[3]=length;

    return dims;
  }
  else{ // Read TIME_ARRAY
    long naxe[] = {1}, inc[]={1};
    GYOTO_DEBUG << "GridData2D::fitsRead(): search " << extname << " HDU" << endl;
    if (fits_movnam_hdu(fptr, ANY_HDU,
              const_cast<char*>(extname.c_str()),
              0, &status))
      throwCfitsioError(status) ;
    GYOTO_DEBUG << "GridData2D::fitsRead(): get image size" << endl;
    if (fits_get_img_size(fptr, 1, naxe, &status)) throwCfitsioError(status) ;
    

    if (dest) { delete [] dest; dest = NULL; }
    dest = new double[nt_];
    for (int ii=0;ii<nt_;ii++)
      dest[ii]=0.;
    if (fits_read_subset(fptr, TDOUBLE, fpixel, naxe, inc,
               0,dest,&anynul,&status)) {
      GYOTO_DEBUG << " error, trying to free pointer" << endl;
      delete [] dest; dest=NULL;
      throwCfitsioError(status) ;
    }
    GYOTO_DEBUG << " done." << endl;

    vector<size_t> dims(1, nt_);
    return dims;
  }
}

fitsfile* GridData2D::fitsCreate(string filename){
  //				 size_t const naxes[3]) {
  //GYOTO_MSG << "GridData2D creating FITS file " << filename << endl; --> THIS leads to bug; CHECK; it seems that only this function, whcih doesn't returns a void, is affected; the same happens if I use cout; however cerr works:
  GYOTO_MSG << "GridData2D creating FITS file " << filename << endl;
  
  char*     pixfile   = const_cast<char*>(filename.c_str());
  int       status    = 0;
  char      ermsg[31] = ""; // ermsg is used in throwCfitsioError()
  long      fpixel[]  = {1};
  fitsfile* fptr; // pointer to FITS file
  
  fits_create_file(&fptr, pixfile, &status);
  
  long      naxes [] = {1};
  // Create the primary extension containing a zero pixel
  fits_create_img(fptr, DOUBLE_IMG, 1, naxes, &status);
  double src[1] = {0.};
  fits_write_pix(fptr, TDOUBLE, fpixel, 1, src, &status);
  if (status) throwCfitsioError(status) ;

  return fptr;
}

void GridData2D::fitsWriteHDU(fitsfile* fptr,
			      string const extname,
			      double* const src,
			      size_t length
			      ) {
  GYOTO_MSG << "GridData2D writing HDU " << extname << endl;
  if (!src) GYOTO_ERROR("GridData2D::fitsWrite: nothing to save!");
  int       status    = 0;
  long      fpixel [] = {1,1,1,1};
  long      ndim      = length?4:3;
  long      naxes []  = {long(nr_),long(nphi_),long(nt_),long(length)};
  char      ermsg[31] = ""; // ermsg is used in throwCfitsioError()
  char *    CNULL     = NULL;

  if (nt_==0 || nr_==0 || nphi_==0)
    GYOTO_ERROR("Axes lengths should be defined");


  
  ////// WRITE KEYWORDS (tmin_,tmax_,rmin_,rmax_) ///////

  if (extname=="VELOCITY") {
    // Writing primary HDU keywords (only once, when writing velocity)
    // Should ideally be placed in fitsCreate but leads to unclear code freezing
    // when calling fits_close_file
    if (tmin_>-DBL_MAX) {
      fits_write_key(fptr, TDOUBLE,
		     const_cast<char*>("GYOTO GridData2D tmin"),
		     &tmin_, CNULL, &status);
      if (status) throwCfitsioError(status) ;
    }

    if (tmax_<DBL_MAX) {
      fits_write_key(fptr, TDOUBLE,
		     const_cast<char*>("GYOTO GridData2D tmax"),
		     &tmax_, CNULL, &status);
      if (status) throwCfitsioError(status) ;
    }

    if (rmin_>0.) {
      fits_write_key(fptr, TDOUBLE,
		     const_cast<char*>("GYOTO GridData2D rmin"),
		     &rmin_, CNULL, &status);
      if (status) throwCfitsioError(status) ;
    }

    if (rmax_<DBL_MAX) {
      fits_write_key(fptr, TDOUBLE,
		     const_cast<char*>("GYOTO GridData2D rmax"),
		     &rmax_, CNULL, &status); 
      if (status) throwCfitsioError(status) ;
    }

    if (phimin_>0.) {
      fits_write_key(fptr, TDOUBLE,
		     const_cast<char*>("GYOTO GridData2D phimin"),
		     &phimin_, CNULL, &status);
      if (status) throwCfitsioError(status) ;
    }

    if (phimax_<2.*M_PI) {
      fits_write_key(fptr, TDOUBLE,
		     const_cast<char*>("GYOTO GridData2D phimax"),
		     &phimax_, CNULL, &status); 
      if (status) throwCfitsioError(status) ;
    }
  }

  ////// SAVE SRC IN APPROPRIATE HDU ///////
  if (extname!="TIMEARRAY"){
    fits_create_img(fptr, DOUBLE_IMG, ndim, naxes, &status);
    std::stringstream ss;
    ss << "GYOTO GridData2D " << extname;
    fits_write_key(fptr, TSTRING, const_cast<char*>("EXTNAME"),
       const_cast<char*>(ss.str().c_str()),
       CNULL, &status);
    fits_write_pix(fptr, TDOUBLE, fpixel, nt_*nphi_*nr_*(length?length:1), src, &status);
    if (status) throwCfitsioError(status) ;
  }else{
    long naxe [] = {long(nt_)}, fpix []={1};
    fits_create_img(fptr, DOUBLE_IMG, 1, naxe, &status);
    std::stringstream ss;
    ss << "GYOTO GridData2D " << extname;
  fits_write_key(fptr, TSTRING, const_cast<char*>("EXTNAME"),
     const_cast<char*>(ss.str().c_str()),
     CNULL, &status);
  fits_write_pix(fptr,TDOUBLE, fpix, nt_, src, &status);
  if (status) throwCfitsioError(status) ;
  }

}

void GridData2D::fitsClose(fitsfile* fptr) {
  GYOTO_MSG << "GridData2D Closing FITS file" << endl;

  int       status    = 0;
  char      ermsg[31] = ""; // ermsg is used in throwCfitsioError()
  int res=0;
  if (fits_close_file(fptr, &status)) throwCfitsioError(status) ;
  fptr = NULL;
}
#endif

void GridData2D::getIndices(size_t i[3], double const tt, double const phi, double const rr, double* const time_array) const {
  // rr is the radius projected in the equaotiral plane
  //cout << "in getInd R: " << rmin_ << " " << rmax_ << " " << nr_ << " " << dr_ << endl;

  if (rmin_>0. && rmax_<DBL_MAX && nr_>0. && dr_>0.) { // >1 radii must be properly defined
    i[2] = size_t(floor((rr-rmin_)/dr_)); // index of closest grid point smaller than rr
                                          // (and same for phi,t)
  } else {
    GYOTO_ERROR("In GridData2D::getIndices: radius undefined!");
  }
  // Then: rmin_ + dr_*i[2] < r < rmin_ + dr_*(i[2]+1)

  //cout << "in getInd PHI: " << phimin_ << " " << phimax_ << " " << nphi_ << " " << dphi_ << endl;
  if (nphi_>0.) {  // necessary condition
    if  (dphi_>0.) {// then >1 phi values provided
      if (phi<phimin_)
	      i[1]=nphi_-1; // then phimin_+dphi_*(nphi_-1) < phi < phimin_ (modulo 2pi)
      else if (phi>phimax_)
	      i[1]=nphi_-1; // idem
      else
	      i[1] = size_t(floor((phi-phimin_)/dphi_)); // then phimin_+dphi_*i[1] < phi < phimin_ +dphi_*(i[1]+1)
    }
    else          // then 1 phi value only, axisym
      i[1] = 0;
  } else {
    GYOTO_ERROR("In GridData2D::getIndices: phi undefined!");
  }

  //cout << "in getInd T: " << tmin_ << " " << tmax_ << " " << nt_ << " " << dt_ << endl;
  if (nt_>0.) {                                // necessary condition
    if (tmin_ > -DBL_MAX && tmax_<DBL_MAX) // >1 time values provided
      if (tt<tmin_) i[0]=0;        // assuming stationnarity before tmin_
      else if (tt>tmax_) i[0]=nt_-1;//                       and after tmax_
      else {
        if (!time_array){
          i[0] = size_t(floor((tt-tmin_)/dt_));
        }else{
          size_t i_t=0;
          while (tt>time_array[i_t] && tt>time_array[i_t+1] && i_t<nt_-1) { // search of i_t
            i_t+=1;
          }
          i[0]=i_t;
        }
      }
    else if (dt_==0. && tmin_==tmax_)          // only 1 time value, stationnary disk
      i[0]=0;
    else
      GYOTO_ERROR("In GridData2D::getIndices: time badly defined!");
  } else {
    GYOTO_ERROR("In GridData2D::getIndices: time undefined!");
  }

}

double GridData2D::interpolate(double tt, double phi, double rcyl,
			       double* const array, double* const time_array) const{
  size_t ind[3]; // {i_t, i_phi, i_r}
  getIndices(ind, tt, phi, rcyl, time_array);

  //cout << ind[0] << " " << ind[1] << " " << ind[2] << endl;
  //cout << " TEST First slot= " << array[0] << " " << array[1] << " " << array[2] << endl;

  double array_interpo=0.;

  if (nphi_==1) // axisym
    throwError("TBD axisym");

  // From here on, >1 phi values

  size_t iphil=-1, iphiu=-1;
  double phil=-1., phiu=-1.;
  // Special treatment for phi to deal with phi<phimin or >phimax (both valid)
  if (phi<phimin_){
    iphil=nphi_-1;
    iphiu=0;
    phil=phimax_-2*M_PI; // this a slightly <0 value
    phiu=phimin_;
    // Then: phil < phi < phiu as it should
  }else if(phi>phimax_){
    iphil=nphi_-1;
    iphiu=0;
    phil=phimax_; // this a slightly <0 value
    phiu=phimin_+2.*M_PI;
    // Then: phil < phi < phiu as it should
  }else{
    iphil=ind[1];
    iphiu=ind[1]+1;
    phil = phimin_+dphi_*iphil;
    phiu = phimin_+dphi_*iphiu;
  }
  
  if (nt_==1 || tt<=tmin_ || tt>=tmax_){
    // Bilinear in r,phi
    size_t irl=ind[2], iru=ind[2]+1, it=ind[0];
    double array_ll = array[it*nr_*nphi_+iphil*nr_+irl], // array_{phi,r}, l=lower index, u=upper index
      array_lu = array[it*nr_*nphi_+iphil*nr_+iru],
      array_ul = array[it*nr_*nphi_+iphiu*nr_+irl],
      array_uu = array[it*nr_*nphi_+iphiu*nr_+iru];
    double rl = rmin_+dr_*irl,
      ru = rmin_+dr_*iru;
    double ratiophi = (phi-phil)/(phiu-phil),
      ratior = (rcyl-rl)/(ru-rl);
    if (ratiophi<0. || ratior<0.)
      GYOTO_ERROR("One ratio is negative (over t, phi or r) !");

    array_interpo = array_ll+(array_ul-array_ll)*ratiophi
      +(array_lu-array_ll)*ratior
      +(array_uu-array_lu-array_ul+array_ll)*ratiophi*ratior;
    /*cout << "bilin interpo stuff: " << endl;
    cout << "T: " << tmin_ << " " << tt << " " << tmax_ << " " << it << endl;
    cout << "PHI: " << phil << " " << phi << " " << phiu << endl;
    cout << "R: " << rl << " " << rcyl << " " << ru << endl;
    cout << "ARRAY at 4 corners (phi,r)=(ll,lu,ul,uu) + interpo: " << array_ll << " " << array_lu << " " << array_ul << " " << array_uu << " " << array_interpo << endl;*/
  }else{
    // Trilinear in r,phi,t
    size_t irl=ind[2], iru=ind[2]+1, itl=ind[0], itu=ind[0]+1;
    double array_lll = array[itl*nr_*nphi_+iphil*nr_+irl], // array_{t,phi,r}
      array_llu = array[itl*nr_*nphi_+iphil*nr_+iru],
      array_lul = array[itl*nr_*nphi_+iphiu*nr_+irl],
      array_luu = array[itl*nr_*nphi_+iphiu*nr_+iru],
      array_ull = array[itu*nr_*nphi_+iphil*nr_+irl],
      array_ulu = array[itu*nr_*nphi_+iphil*nr_+iru],
      array_uul = array[itu*nr_*nphi_+iphiu*nr_+irl],
      array_uuu = array[itu*nr_*nphi_+iphiu*nr_+iru];

    double rl = rmin_+dr_*irl,
      ru = rmin_+dr_*iru,
      tl = 0, tu = 0;
    if (time_array){
      tl = time_array[itl],
      tu = time_array[itu];
    }else{
      tl = tmin_+dt_*itl;
      tu = rmin_+dt_*itu;
    }

    double ratiot = (tt-tl)/(tu-tl),
      ratiophi = (phi-phil)/(phiu-phil),
      ratior = (rcyl-rl)/(ru-rl);
    if (ratiot<0. || ratiophi<0. || ratior<0.)
      GYOTO_ERROR("One ratio is negative (over t, phi or r) !");

    array_interpo = array_lll
      + (array_ull-array_lll)*ratiot
      + (array_lul-array_lll)*ratiophi
      + (array_llu-array_lll)*ratior
      +(array_uul-array_lul-array_ull+array_lll)*ratiot*ratiophi
      +(array_luu-array_lul-array_llu+array_lll)*ratiophi*ratior
      +(array_ulu-array_llu-array_ull+array_lll)*ratiot*ratior
      +(array_uuu-array_luu-array_ulu-array_uul
	+array_ull+array_llu+array_lul-array_lll)*ratiot*ratiophi*ratior;
    /*cout << "trilin interpo stuff: " << endl;
    cout << "T: " << tl << " " << tt << " " << tu << endl;
    cout << "PHI: " << phil << " " << phi << " " << phiu << endl;
    cout << "R: " << rl << " " << rcyl << " " << ru << endl;
    cout << "ratiot,phi,r= " << ratiot << " " << ratiophi << " " << ratior << endl;
    cout << "ARRAY at 8 corners (t,phi,r)=(lll,llu,lul,luu,ull,ulu,uul,uuu) + interpo: " << array_lll << " " << array_llu << " " << array_lul << " " << array_luu << " " << array_ull << " " << array_ulu << " " << array_uul << " " << array_uuu << " " << array_interpo << endl;*/
  }
  
  return array_interpo;
}


