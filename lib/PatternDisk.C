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
#define throwCfitsioError(status) \
    { fits_get_errstatus(status, ermsg); throwError(ermsg); }

#include "GyotoPhoton.h"
#include "GyotoPatternDisk.h"
#include "GyotoUtils.h"
#include "GyotoFactoryMessenger.h"
#include "GyotoKerrBL.h"
#include "GyotoKerrKS.h"


#include <fitsio.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <fstream>
#include <cstring>
#include <cmath>
#include <limits>

using namespace std;
using namespace Gyoto;
using namespace Gyoto::Astrobj;

PatternDisk::PatternDisk() :
  ThinDisk("PatternDisk"), filename_(""),
  emission_(NULL), opacity_(NULL), velocity_(NULL), radius_(NULL),
  Omega_(0.), t0_(0.),
  dnu_(1.), nu0_(0), nnu_(0),
  dphi_(0.), phimin_(0.), 
  nphi_(0), phimax_(2*M_PI), repeat_phi_(1),
  dr_(0.), nr_(0)
{
  GYOTO_DEBUG << "PatternDisk Construction" << endl;
}

PatternDisk::PatternDisk(const PatternDisk& o) :
  ThinDisk(o), filename_(o.filename_),
  emission_(NULL), opacity_(NULL), velocity_(NULL), radius_(NULL),
  Omega_(o.Omega_), t0_(o.t0_),
  dnu_(o.dnu_), nu0_(o.nu0_), nnu_(o.nnu_),
  dphi_(o.dphi_), phimin_(o.phimin_),
  nphi_(o.nphi_), phimax_(o.phimax_), repeat_phi_(o.repeat_phi_),
  dr_(o.dr_), nr_(o.nr_)
{
  GYOTO_DEBUG << "PatternDisk Copy" << endl;
  size_t ncells = 0;
  if (o.emission_) {
    emission_ = new double[ncells = nnu_ * nphi_ * nr_];
    memcpy(emission_, o.emission_, ncells * sizeof(double));
  }
  if (o.opacity_) {
    opacity_ = new double[ncells = nnu_ * nphi_ * nr_];
    memcpy(opacity_, o.opacity_, ncells * sizeof(double));
  }
  if (o.velocity_) {
    velocity_ = new double[ncells = 2 * nphi_ * nr_];
    memcpy(velocity_, o.velocity_, ncells * sizeof(double));
  }
  if (o.radius_) {
    radius_ = new double[ncells = 2 * nphi_ * nr_];
    memcpy(radius_, o.radius_, ncells * sizeof(double));
  }
}
PatternDisk* PatternDisk::clone() const
{ return new PatternDisk(*this); }

PatternDisk::~PatternDisk() {
  GYOTO_DEBUG << "PatternDisk Destruction" << endl;
  if (emission_) delete [] emission_;
  if (opacity_) delete [] opacity_;
  if (velocity_) delete [] velocity_;
  if (radius_) delete [] radius_;
}

void PatternDisk::setEmission(double * pattern) {
  emission_ = pattern;
}

void PatternDisk::setVelocity(double * pattern) {
  velocity_ = pattern;
}

void PatternDisk::setRadius(double * pattern) {
  radius_ = pattern;
}

void PatternDisk::copyIntensity(double const *const pattern, size_t const naxes[3]) {
  GYOTO_DEBUG << endl;
  if (emission_) {
    GYOTO_DEBUG << "delete [] emission_;" << endl;
    delete [] emission_; emission_ = NULL;
  }
  if (pattern) {
    size_t nel;
    if (nnu_ != naxes[0]) {
      if (opacity_)  { delete [] opacity_; opacity_  = NULL; }
    }
    if (nphi_ != naxes[1]) {
      GYOTO_DEBUG <<"nphi_ changed, freeing velocity_" << endl;
      if (opacity_)  { delete [] opacity_; opacity_  = NULL; }
      if (velocity_) { delete [] velocity_; velocity_= NULL; }
    }
    if (nr_ != naxes[2]) {
      GYOTO_DEBUG <<"nr_ changed, freeing velocity_ and radius_" << endl;
      if (opacity_)  { delete [] opacity_;  opacity_ = NULL; }
      if (velocity_) { delete [] velocity_; velocity_= NULL; }
      if (radius_)   { delete [] radius_;   radius_  = NULL; }
    }
    if (!(nel=(nnu_ = naxes[0]) * (nphi_=naxes[1]) * (nr_=naxes[2])))
      throwError( "dimensions can't be null");
    if (nr_==1 || nphi_==1)
      throwError("In PatternDisk::copyIntensity: dimensions should be >1");
    dr_ = (rout_ - rin_) / double(nr_-1);
    if (repeat_phi_==0.)
      throwError("In PatternDisk::copyIntensity: repeat_phi is 0!");
    dphi_ = (phimax_-phimin_)/double((nphi_-1)*repeat_phi_);
    GYOTO_DEBUG << "allocate emission_;" << endl;
    emission_ = new double[nel];
    GYOTO_DEBUG << "pattern >> emission_" << endl;
    memcpy(emission_, pattern, nel*sizeof(double));
  }
}

double const * PatternDisk::getIntensity() const { return emission_; }
void PatternDisk::getIntensityNaxes( size_t naxes[3] ) const
{ naxes[0] = nnu_; naxes[1] = nphi_; naxes[2] = nr_; }

void PatternDisk::copyOpacity(double const *const opacity, size_t const naxes[3]) {
  GYOTO_DEBUG << endl;
  if (opacity_) {
    GYOTO_DEBUG << "delete [] opacity_;" << endl;
    delete [] opacity_; opacity_ = NULL;
    flag_radtransf_=0;
  }
  if (opacity) {
    if (nnu_ != naxes[0] || nphi_ != naxes[1] || nr_ != naxes[2])
      throwError("Please set intensity before opacity. "
		 "The two arrays must have the same dimensions.");
    GYOTO_DEBUG << "allocate opacity_;" << endl;
    opacity_ = new double[nnu_ * nphi_ * nr_];
    GYOTO_DEBUG << "opacity >> opacity_" << endl;
    memcpy(opacity_, opacity, nnu_ * nphi_ * nr_ * sizeof(double));
    flag_radtransf_=1;
  }
}

double const * PatternDisk::getOpacity() const { return opacity_; }

void PatternDisk::copyVelocity(double const *const velocity, size_t const naxes[2]) {
  GYOTO_DEBUG << endl;
  if (velocity_) {
    GYOTO_DEBUG << "delete [] velocity_;\n";
    delete [] velocity_; velocity_ = NULL;
  }
  if (velocity) {
    if (!emission_) throwError("Please use copyIntensity() before copyVelocity()");
    if (nphi_ != naxes[0] || nr_ != naxes[1])
      throwError("emission_ and velocity_ have inconsistent dimensions");
    GYOTO_DEBUG << "allocate velocity_;" << endl;
    velocity_ = new double[2*nphi_*nr_];
    GYOTO_DEBUG << "velocity >> velocity_" << endl;
    memcpy(velocity_, velocity, 2*nphi_*nr_*sizeof(double));
  }
}
double const * PatternDisk::getVelocity() const { return velocity_; }

void PatternDisk::copyGridRadius(double const *const radius, size_t nr) {
  GYOTO_DEBUG << endl;
  if (radius_) {
    GYOTO_DEBUG << "delete [] radius_;" << endl;
    delete [] radius_; radius_ = NULL;
  }
  if (radius) {
    if (!emission_) throwError("Please use copyIntensity() before copyGridRadius()");
    if (nr_ != nr)
      throwError("emission_ and radius_ have inconsistent dimensions");
    GYOTO_DEBUG << "allocate velocity_;" << endl;
    radius_ = new double[nr_];
    GYOTO_DEBUG << "velocity >> velocity_" << endl;
    memcpy(radius_, radius, nr_*sizeof(double));
    rin_=radius_[0];
    rout_=radius_[nr_-1];
  }
}
double const * PatternDisk::getGridRadius() const { return radius_; }

void PatternDisk::repeatPhi(size_t n) {
  repeat_phi_ = n;
  if ((nphi_-1)*repeat_phi_>0) 
    dphi_=(phimax_-phimin_)/double((nphi_-1)*repeat_phi_);
}
size_t PatternDisk::repeatPhi() const { return repeat_phi_; }

void PatternDisk::nu0(double freq) { nu0_ = freq; }
double PatternDisk::nu0() const { return nu0_; }

void PatternDisk::dnu(double dfreq) { dnu_ = dfreq; }
double PatternDisk::dnu() const { return dnu_; }

void PatternDisk::phimin(double phimn) {
  phimin_ = phimn;
  if (nphi_>1) dphi_ = (phimax_-phimin_) / double(nphi_-1);
}
double PatternDisk::phimin() const {return phimin_;}

void PatternDisk::phimax(double phimx) {
  phimax_ = phimx;
  if (nphi_>1) dphi_ = (phimax_-phimin_) / double(nphi_-1);
}
double PatternDisk::phimax() const {return phimax_;}

void PatternDisk::fitsRead(string filename) {
  GYOTO_MSG << "PatternDisk reading FITS file: " << filename << endl;

  filename_ = filename;
  int rin_set=0, rout_set=0;
  char*     pixfile   = const_cast<char*>(filename_.c_str());
  fitsfile* fptr      = NULL;
  int       status    = 0;
  int       anynul    = 0;
  long      tmpl;
  double    tmpd;
  long      naxes []  = {1, 1, 1};
  long      fpixel[]  = {1,1,1};
  long      inc   []  = {1,1,1};
  char      ermsg[31] = ""; // ermsg is used in throwCfitsioError()

  GYOTO_DEBUG << "PatternDisk::readFile(): opening file" << endl;
  if (fits_open_file(&fptr, pixfile, 0, &status)) throwCfitsioError(status) ;

  ////// READ FITS KEYWORDS COMMON TO ALL TABLES ///////
  //get Omega and t0;
  GYOTO_DEBUG << "PatternDisk::readFile(): read Omega_" << endl;
  fits_read_key(fptr, TDOUBLE, "GYOTO PatternDisk Omega", &tmpd, NULL, &status);
  if (status) {
    if (status == KEY_NO_EXIST) status = 0; // not fatal
    else throwCfitsioError(status) ;
  } else Omega_ = tmpd; // Omega_ found
  GYOTO_DEBUG << "PatternDisk::readFile(): read t0_" << endl;
  fits_read_key(fptr, TDOUBLE, "GYOTO PatternDisk t0", &tmpd, NULL, &status);
  if (status) {
    if (status == KEY_NO_EXIST) status = 0; // not fatal
    else throwCfitsioError(status) ;
  } else t0_ = tmpd; // T0_ found

  GYOTO_DEBUG << "PatternDisk::readFile(): read RepeatPhi_" << endl;
  fits_read_key(fptr, TLONG, "GYOTO PatternDisk RepeatPhi", &tmpl,
		NULL, &status);
  if (status) {
    if (status == KEY_NO_EXIST) status = 0; // not fatal
    else throwCfitsioError(status) ;
  } else repeat_phi_ = size_t(tmpl); // RepeatPhi found

  GYOTO_DEBUG << "PatternDisk::readFile(): read InnerRadius_" << endl;
  fits_read_key(fptr, TDOUBLE, "GYOTO ThinDisk InnerRadius", &tmpd,
		NULL, &status);
  if (status) {
    if (status == KEY_NO_EXIST) status = 0; // not fatal
    else throwCfitsioError(status) ;
  } else {
    rin_ = tmpd; // InnerRadius found
    rin_set=1;
  }
  GYOTO_DEBUG << "PatternDisk::readFile(): read OuterRadius_" << endl;
  fits_read_key(fptr, TDOUBLE, "GYOTO ThinDisk OuterRadius", &tmpd,
		NULL, &status);
  if (status) {
    if (status == KEY_NO_EXIST) status = 0; // not fatal
    else throwCfitsioError(status) ;
  } else {
    rout_ = tmpd; // OuterRadius found
    rout_set=1;
  }

  GYOTO_DEBUG << "PatternDisk::fitsRead(): read Phimin_" << endl;
  fits_read_key(fptr, TDOUBLE, "GYOTO PatternDisk Phimin", &tmpd,
		NULL, &status);
  if (status) {
    if (status == KEY_NO_EXIST) {
      status = 0; // not fatal
      phimin_=0.; //phimin defaults to 0
    }
    else throwCfitsioError(status) ;
  } else {
    phimin_ = tmpd; // Phimin found
  }
  GYOTO_DEBUG << "PatternDisk::fitsRead(): read Phimax_" << endl;
  fits_read_key(fptr, TDOUBLE, "GYOTO PatternDisk Phimax", &tmpd,
		NULL, &status);
  if (status) {
    if (status == KEY_NO_EXIST) {
      status = 0; // not fatal
      phimax_=2.*M_PI; //phimin defaults to 0
    }
    else throwCfitsioError(status) ;
  } else {
    phimax_ = tmpd; // Phimax found
  }

  ////// FIND MANDATORY EMISSION HDU, READ KWDS & DATA ///////
  GYOTO_DEBUG << "PatternDisk::readFile(): search emission HDU" << endl;
  if (fits_movnam_hdu(fptr, ANY_HDU,
		      const_cast<char*>("GYOTO PatternDisk emission"),
		      0, &status))
    throwCfitsioError(status) ;
  GYOTO_DEBUG << "PatternDisk::readFile(): get image size" << endl;
  if (fits_get_img_size(fptr, 3, naxes, &status)) throwCfitsioError(status) ;

  //update nu0_, nnu_, dnu_;
  nnu_ = naxes[0]; 
  double CRPIX1;
  GYOTO_DEBUG << "PatternDisk::readFile(): read CRPIX1, CRVAL1, CDELT1"
		    << endl;
  fits_read_key(fptr, TDOUBLE, "CRVAL1", &nu0_, NULL, &status);
  fits_read_key(fptr, TDOUBLE, "CDELT1", &dnu_, NULL, &status);
  fits_read_key(fptr, TDOUBLE, "CRPIX1", &CRPIX1, NULL, &status);
  if (status) throwCfitsioError(status) ;
  if (CRPIX1 != 1) nu0_ -= dnu_*(CRPIX1 - 1.);

  // update repeat_phi_, nphi_, dphi_
  nphi_ = naxes[1];
  if (nphi_==1 || repeat_phi_==0)
    throwError("In PatternDisk::fitsRead: nphi is 1 or repeat_phi is 0!");
  dphi_ = (phimax_-phimin_)/double((nphi_-1)*repeat_phi_);

  // update rin_, rout_, nr_, dr_
  nr_ = naxes[2];

  if (emission_) { delete [] emission_; emission_ = NULL; }
  emission_ = new double[nnu_ * nphi_ * nr_];
  if (debug())
    cerr << "PatternDisk::readFile(): read emission: "
	 << "nnu_=" << nnu_ << ", nphi_="<<nphi_ << ", nr_="<<nr_ << "...";
  if (fits_read_subset(fptr, TDOUBLE, fpixel, naxes, inc,
		       0, emission_,&anynul,&status)) {
    GYOTO_DEBUG << " error, trying to free pointer" << endl;
    delete [] emission_; emission_=NULL;
    throwCfitsioError(status) ;
  }
  GYOTO_DEBUG << " done." << endl;

  ////// FIND OPTIONAL OPACITY HDU ///////

  fits_movnam_hdu(fptr, ANY_HDU,
		  const_cast<char*>("GYOTO PatternDisk opacity"),
		  0, &status);
  if (status) {
    if (status == BAD_HDU_NUM) {
      // FITS file does not contain opacity information
      status = 0;
      if (opacity_) { delete [] opacity_; opacity_ = NULL; }
    } else throwCfitsioError(status) ;
  } else {
    if (fits_get_img_size(fptr, 3, naxes, &status)) throwCfitsioError(status) ;
    if (   size_t(naxes[0]) != nnu_
	|| size_t(naxes[1]) != nphi_
	|| size_t(naxes[2]) != nr_)
      throwError("PatternDisk::readFile(): opacity array not conformable");
    if (opacity_) { delete [] opacity_; opacity_ = NULL; }
    opacity_ = new double[nnu_ * nphi_ * nr_];
    if (fits_read_subset(fptr, TDOUBLE, fpixel, naxes, inc, 
			 0, opacity_,&anynul,&status)) {
      delete [] opacity_; opacity_=NULL;
      throwCfitsioError(status) ;
    }
  }

  ////// FIND OPTIONAL VELOCITY HDU ///////

  fits_movnam_hdu(fptr, ANY_HDU,
		  const_cast<char*>("GYOTO PatternDisk velocity"),
		  0, &status);
  if (status) {
    if (status == BAD_HDU_NUM) {
      // FITS file does not contain velocity information
      status = 0;
      if (velocity_) { delete [] velocity_; velocity_ = NULL; }
    } else throwCfitsioError(status) ;
  } else {
    if (fits_get_img_size(fptr, 3, naxes, &status)) throwCfitsioError(status) ;
    if (   size_t(naxes[0]) != size_t(2)
	|| size_t(naxes[1]) != nphi_
	|| size_t(naxes[2]) != nr_)
      throwError("PatternDisk::readFile(): velocity array not conformable");
    if (velocity_) { delete [] velocity_; velocity_ = NULL; }
    velocity_ = new double[2 * nphi_ * nr_];
    if (fits_read_subset(fptr, TDOUBLE, fpixel, naxes, inc, 
			 0, velocity_,&anynul,&status)) {
      delete [] velocity_; velocity_=NULL;
      throwCfitsioError(status) ;
    }
  }

  ////// FIND OPTIONAL RADIUS HDU ///////
  
  fits_movnam_hdu(fptr, ANY_HDU,
		  const_cast<char*>("GYOTO PatternDisk radius"),
		  0, &status);
   if (status) {
    if (status == BAD_HDU_NUM) {
      // FITS file does not contain explicit radius information
      status = 0;
      if (radius_) { delete [] radius_; radius_ = NULL; }
    } else throwCfitsioError(status) ;
  } else {
    if (fits_get_img_size(fptr, 1, naxes, &status)) throwCfitsioError(status) ;
    if (size_t(naxes[0]) != nr_)
      throwError("PatternDisk::readFile(): radius array not conformable");
    if (radius_) { delete [] radius_; radius_ = NULL; }
    radius_ = new double[nr_];
    if (fits_read_subset(fptr, TDOUBLE, fpixel, naxes, inc, 
			 0, radius_,&anynul,&status)) {
      delete [] radius_; radius_=NULL;
      throwCfitsioError(status) ;
    }
    if (!rin_set) rin_=radius_[0];
    if (!rout_set) rout_=radius_[nr_-1];
  }

   if (nr_ == 1.)
     throwError("In PatternDisk::fitsRead: nr_ should not be 0 here!");
   dr_ = (rout_-rin_) / double(nr_-1);

  if (fits_close_file(fptr, &status)) throwCfitsioError(status) ;
  fptr = NULL;
}

void PatternDisk::fitsWrite(string filename) {
  if (!emission_) throwError("PatternDisk::fitsWrite(filename): nothing to save!");
  filename_ = filename;
  char*     pixfile   = const_cast<char*>(filename_.c_str());
  fitsfile* fptr      = NULL;
  int       status    = 0;
  long      naxes []  = {long(nnu_), long(nphi_), long(nr_)};
  long      fpixel[]  = {1,1,1};
  char * CNULL=NULL;

  char      ermsg[31] = ""; // ermsg is used in throwCfitsioError()

  ////// CREATE FILE
  GYOTO_DEBUG << "creating file \"" << pixfile << "\"... ";
  fits_create_file(&fptr, pixfile, &status);
  if (debug()) cerr << "done." << endl;
  fits_create_img(fptr, DOUBLE_IMG, 3, naxes, &status);
  if (status) throwCfitsioError(status) ;

  ////// WRITE FITS KEYWORDS COMMON TO ALL TABLES ///////
  //set Omega and t0;
  if (Omega_!=0)
    fits_write_key(fptr, TDOUBLE,
		   const_cast<char*>("GYOTO PatternDisk Omega"),
		   &Omega_, CNULL, &status);
  if (t0_!=0)
    fits_write_key(fptr, TDOUBLE,
		   const_cast<char*>("GYOTO PatternDisk t0"),
		   &t0_, CNULL, &status);
  if (repeat_phi_!=1)
    fits_write_key(fptr, TLONG,
		   const_cast<char*>("GYOTO PatternDisk RepeatPhi"),
		   &repeat_phi_, CNULL, &status);

  if ((radius_ && rin_ != radius_[0]) || (!radius_ && rin_ != 0.))
    fits_write_key(fptr, TDOUBLE,
		   const_cast<char*>("GYOTO ThinDisk InnerRadius"),
		   &rin_, CNULL, &status);
  if ((radius_ && rout_ != radius_[nr_-1]) || (!radius_ && rout_ != DBL_MAX))
    fits_write_key(fptr, TDOUBLE,
		   const_cast<char*>("GYOTO ThinDisk OuterRadius"),
		   &rout_, CNULL, &status);

  if (phimin_ > -DBL_MAX)
    fits_write_key(fptr, TDOUBLE,
		   const_cast<char*>("GYOTO PatternDisk Phimin"),
		   &phimin_, CNULL, &status);

  if (phimax_ < DBL_MAX)
    fits_write_key(fptr, TDOUBLE,
		   const_cast<char*>("GYOTO PatternDisk Phimax"),
		   &phimax_, CNULL, &status);

  ////// SAVE EMISSION IN PRIMARY HDU ///////
  GYOTO_DEBUG << "saving emission_\n";
  fits_write_key(fptr, TSTRING,
		 const_cast<char*>("EXTNAME"),
		 const_cast<char*>("GYOTO PatternDisk emission"),
		 CNULL, &status);
  fits_write_key(fptr, TDOUBLE,
		 const_cast<char*>("CRVAL1"),
		 &nu0_, CNULL, &status);
  fits_write_key(fptr, TDOUBLE,
		 const_cast<char*>("CDELT1"),
		 &dnu_, CNULL, &status);
  double CRPIX1 = 1.;
  fits_write_key(fptr, TDOUBLE,
		 const_cast<char*>("CRPIX1"),
		 &CRPIX1, CNULL, &status);
  fits_write_pix(fptr, TDOUBLE, fpixel, nnu_*nphi_*nr_, emission_, &status);
  if (status) throwCfitsioError(status) ;

  ////// SAVE OPTIONAL OPACITY HDU ///////
  if (opacity_) {
    GYOTO_DEBUG << "saving opacity_\n";
    fits_create_img(fptr, DOUBLE_IMG, 3, naxes, &status);
    fits_write_key(fptr, TSTRING, const_cast<char*>("EXTNAME"),
		   const_cast<char*>("GYOTO PatternDisk opacity"),
		   CNULL, &status);
    fits_write_pix(fptr, TDOUBLE, fpixel, nnu_*nphi_*nr_, opacity_, &status);
    if (status) throwCfitsioError(status) ;
  }

  ////// SAVE OPTIONAL VELOCITY HDU ///////
  if (velocity_) {
    GYOTO_DEBUG << "saving velocity_\n";
    naxes[0]=2;
    fits_create_img(fptr, DOUBLE_IMG, 3, naxes, &status);
    fits_write_key(fptr, TSTRING, const_cast<char*>("EXTNAME"),
		   const_cast<char*>("GYOTO PatternDisk velocity"),
		   CNULL, &status);
    fits_write_pix(fptr, TDOUBLE, fpixel, 2*nphi_*nr_, velocity_, &status);
    if (status) throwCfitsioError(status) ;
  }

  ////// SAVE OPTIONAL RADIUS HDU ///////
  if (radius_) {
    GYOTO_DEBUG << "saving velocity_\n";
    fits_create_img(fptr, DOUBLE_IMG, 1, naxes+2, &status);
    fits_write_key(fptr, TSTRING, const_cast<char*>("EXTNAME"),
		   const_cast<char*>("GYOTO PatternDisk radius"),
		   CNULL, &status);
    fits_write_pix(fptr, TDOUBLE, fpixel, nr_, radius_, &status);
    if (status) throwCfitsioError(status) ;
  }

  ////// CLOSING FILE ///////
  GYOTO_DEBUG << "close FITS file\n";
  if (fits_close_file(fptr, &status)) throwCfitsioError(status) ;
  fptr = NULL;
}

void PatternDisk::getIndices(size_t i[3], double const co[4], double nu) const {
  GYOTO_DEBUG << "dnu_="<<dnu_<<", dphi_="<<dphi_<<", dr_="<<dr_<<endl;
  if (nu <= nu0_) i[0] = 0;
  else {
    i[0] = size_t(floor((nu-nu0_)/dnu_+0.5));
    if (i[0] >= nnu_) i[0] = nnu_-1;
  }
  double r = projectedRadius(co);
  double phi = sphericalPhi(co);
  double t = co[0];

  phi -= Omega_*(t-t0_);

  while (phi<0) phi += 2.*M_PI;
  if (dphi_==0.)
    throwError("In PatternDisk::getIndices: dphi_ should not be 0 here!");
  if (phi<phimin_) //Possible: any phi value is in the grid anyway
    i[1]=0;
  else if (phi>phimax_)
    i[1]=nphi_-1;
  else
    i[1] = size_t(floor((phi-phimin_)/dphi_)+1) % nphi_;

  /*
    With this definition:
    phimin_+(i[1]-1)*dphi_ <= phi < phimin_+i[1]*dphi_ 
    provided phi is not bigger than phimax_ nor smaller
    than phimin_
   */
  if (radius_) {
    GYOTO_DEBUG <<"radius_ != NULL" << endl;
    // if the radius_ vector is set, find closest value
    if (r >= radius_[nr_-1]) i[2] = nr_-1;
    else {
      for(i[2]=0; r > radius_[i[2]]; ++i[2]){}
      /*
	With this definition:
	radius_[i[2]-1] <= r < radius_[i[2]]
	provided r<rmax (i[2] is always at least 1)
       */
      /*if (i[2]>0 && r-radius_[i[2]-1] < radius_[i[2]]) {
	--i[2];
	}*/
    }
  } else {
    GYOTO_DEBUG <<"radius_ == NULL, dr_==" << dr_ << endl;
    // radius_ is not set: assume linear repartition
    if (dr_==0.)
      throwError("In PatternDisk::getIndices: dr_ should not be 0 here!");
    i[2] = size_t(floor((r-rin_)/dr_+0.5));
    if (i[2] >= nr_) i[2] = nr_ - 1;
  }

}

void PatternDisk::getVelocity(double const pos[4], double vel[4]) {
  if (velocity_) {
    if (dir_ != 1)
      throwError("PatternDisk::getVelocity(): "
		 "dir_ should be 1 if velocity_ is provided");
    size_t i[3]; // {i_nu, i_phi, i_r}
    getIndices(i, pos);

    double rr = projectedRadius(pos);
    double phi = sphericalPhi(pos);

    double phiprime=0., rprime=0.;
    if (i[1]==0 || i[1]==nphi_-1 || i[2]==nr_-1){
      // Extreme cases no interpolation
      phiprime=velocity_[i[2]*(nphi_*2)+i[1]*2+0];
      rprime=velocity_[i[2]*(nphi_*2)+i[1]*2+1];
    }else{
      // Bilinear interpolation
      double phip0=velocity_[(i[2]-1)*(nphi_*2)+(i[1]-1)*2+0];
      double phip1=velocity_[(i[2]-1)*(nphi_*2)+i[1]*2+0];
      double phip2=velocity_[i[2]*(nphi_*2)+i[1]*2+0];
      double phip3=velocity_[i[2]*(nphi_*2)+(i[1]-1)*2+0];
      
      double rp0=velocity_[(i[2]-1)*(nphi_*2)+(i[1]-1)*2+1];
      double rp1=velocity_[(i[2]-1)*(nphi_*2)+i[1]*2+1];
      double rp2=velocity_[i[2]*(nphi_*2)+i[1]*2+1];
      double rp3=velocity_[i[2]*(nphi_*2)+(i[1]-1)*2+1];

      double rinf=radius_[i[2]-1], rsup=radius_[i[2]],
	phiinf=phimin_+(i[1]-1)*dphi_, phisup=phiinf+dphi_;

      if (phi<phiinf || phi>phisup || rr<rinf || rr>rsup){
	throwError("In PatternDisk::getVelocity: "
		   "bad interpolation");
      }

      double cr = (rr-rinf)/(rsup-rinf),
	cp = (phi-phiinf)/(phisup-phiinf);

      rprime=(1-cr)*(1-cp)*rp0
	+cr*(1-cp)*rp1
	+cr*cp*rp2
	+(1-cr)*cp*rp3;

      phiprime=(1-cr)*(1-cp)*phip0
	+cr*(1-cp)*phip1
	+cr*cp*phip2
	+(1-cr)*cp*phip3;
    }
	
    switch (gg_->coordKind()) {
    case GYOTO_COORDKIND_SPHERICAL:
      {
	double pos2[4] = {pos[0], pos[1], pos[2], pos[3]};
	pos2[1] = radius_ ? radius_[i[2]] : rin_+double(i[2])*dr_;
	vel[1] = rprime;
	vel[2] = 0.;
	vel[3] = phiprime;
	vel[0] = gg_->SysPrimeToTdot(pos2, vel+1);
	vel[1] *= vel[0];
	vel[3] *= vel[0];
      }
      break;
    case GYOTO_COORDKIND_CARTESIAN:
      throwError("PatternDisk::getVelocity(): metric must be in "
		 "spherical coordinaites if velocity field is provided");
      break;
    default:
      throwError("PatternDisk::getVelocity(): unknown COORDKIND");
    } 
  }else ThinDisk::getVelocity(pos, vel);
}

double PatternDisk::emission(double nu, double dsem,
				    double *,
				    double co[8]) const{
  //See Page & Thorne 74 Eqs. 11b, 14, 15. This is F(r).
  GYOTO_DEBUG << endl;
  size_t i[3]; // {i_nu, i_phi, i_r}
  getIndices(i, co, nu);

    double rr = projectedRadius(co);
    double phi = sphericalPhi(co);

    double Iem=0.;
    if (i[1]==0 || i[1]==nphi_-1 || i[2]==nr_-1){
      // Extreme cases no interpolation
      Iem=emission_[i[2]*(nphi_*nnu_)+i[1]*nnu_+i[0]];
    }else{
      // Bilinear interpolation
      double Iem0=emission_[(i[2]-1)*(nphi_*nnu_)+(i[1]-1)*nnu_+i[0]];
      double Iem1=emission_[(i[2]-1)*(nphi_*nnu_)+i[1]*nnu_+i[0]];
      double Iem2=emission_[i[2]*(nphi_*nnu_)+i[1]*nnu_+i[0]];
      double Iem3=emission_[i[2]*(nphi_*nnu_)+(i[1]-1)*nnu_+i[0]];

      double rinf=radius_[i[2]-1], rsup=radius_[i[2]],
	phiinf=phimin_+(i[1]-1)*dphi_, phisup=phiinf+dphi_;
      
      if (phi<phiinf || phi>phisup || rr<rinf || rr>rsup){
	throwError("In PatternDisk::emission: "
		   "bad interpolation");
      }

      double cr = (rr-rinf)/(rsup-rinf),
	cp = (phi-phiinf)/(phisup-phiinf);

      Iem = (1-cr)*(1-cp)*Iem0
	+cr*(1-cp)*Iem1
	+cr*cp*Iem2
	+(1-cr)*cp*Iem3;
    }

  if (!flag_radtransf_) return Iem;
  double thickness;
  // NB: thickness is not interpolated so far
  if (opacity_ && (thickness=opacity_[i[2]*(nphi_*nnu_)+i[1]*nnu_+i[0]]*dsem))
    return Iem * (1. - exp (-thickness)) ;
  return 0.;
}

double PatternDisk::transmission(double nu, double dsem, double*co) const {
  GYOTO_DEBUG << endl;
  if (!flag_radtransf_) return 0.;
  if (!opacity_) return 1.;
  size_t i[3]; // {i_nu, i_phi, i_r}
  getIndices(i, co, nu);
  // NB: opacity is not interpolated so far
  double opacity = opacity_[i[2]*(nphi_*nnu_)+i[1]*nnu_+i[0]];
  GYOTO_DEBUG << "nu="<<nu <<", dsem="<<dsem << ", opacity="<<opacity <<endl;
  if (!opacity) return 1.;
  return exp(-opacity*dsem);
}

void PatternDisk::setInnerRadius(double rin) {
  ThinDisk::setInnerRadius(rin);
  if (nr_>1 && !radius_) dr_ = (rout_-rin_) / double(nr_-1);
}

void PatternDisk::setOuterRadius(double rout) {
  ThinDisk::setOuterRadius(rout);
  if (nr_>1 && !radius_) dr_ = (rout_-rin_) / double(nr_-1);
}

void PatternDisk::setPatternVelocity(double omega) { Omega_ = omega; }
double PatternDisk::getPatternVelocity() { return Omega_; }

int PatternDisk::setParameter(std::string name,
			      std::string content,
			      std::string unit) {
  if      (name == "File")          fitsRead( content );
  else if (name=="PatternVelocity") setPatternVelocity(atof(content.c_str()));
  else return ThinDisk::setParameter(name, content, unit);
  return 0;
}

#ifdef GYOTO_USE_XERCES
void PatternDisk::fillElement(FactoryMessenger *fmp) const {
  fmp->setParameter("File", (filename_.compare(0,1,"!") ?
			     filename_ :
			     filename_.substr(1)));
  if (Omega_) fmp->setParameter("PatternVelocity", Omega_);
  ThinDisk::fillElement(fmp);
}

void PatternDisk::setParameters(FactoryMessenger* fmp) {
  string name, content, unit;
  metric(fmp->metric());
  while (fmp->getNextParameter(&name, &content, &unit)) {
    if  (name == "File") setParameter(name, fmp -> fullPath(content), unit);
    else setParameter(name, content, unit);
  }
}
#endif
