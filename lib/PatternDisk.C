/*
    Copyright 2011-2014, 2017-2019 Thibaut Paumard, Frederic Vincent

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
#include "GyotoPhoton.h"
#include "GyotoPatternDisk.h"
#include "GyotoProperty.h"
#include "GyotoUtils.h"
#include "GyotoFactoryMessenger.h"
#include "GyotoKerrBL.h"
#include "GyotoKerrKS.h"

#ifdef GYOTO_USE_CFITSIO
#include <fitsio.h>
#define throwCfitsioError(status) \
    { fits_get_errstatus(status, ermsg); GYOTO_ERROR(ermsg); }
#endif
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

/// Properties

GYOTO_PROPERTY_START(PatternDisk)
GYOTO_PROPERTY_FILENAME(PatternDisk, File, file)
GYOTO_PROPERTY_DOUBLE(PatternDisk, PatternVelocity, patternVelocity)
GYOTO_PROPERTY_END(PatternDisk, ThinDisk::properties)

void PatternDisk::fillProperty(Gyoto::FactoryMessenger *fmp,
			       Property const &p) const {
  if (p.name == "File")
    fmp->setParameter("File", (filename_.compare(0,1,"!") ?
			       filename_ :
			       filename_.substr(1)) );
  else ThinDisk::fillProperty(fmp, p);
}

///

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

void PatternDisk::radius(double * pattern) {
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
      GYOTO_ERROR( "dimensions can't be null");
    if (nr_==1)
      GYOTO_ERROR("In PatternDisk::copyIntensity: "
		 "radial dimension should be >1");
    dr_ = (rout_ - rin_) / double(nr_-1);
    if (repeat_phi_==0.)
      GYOTO_ERROR("In PatternDisk::copyIntensity: repeat_phi is 0!");
    if (nphi_>1)
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

void PatternDisk::copyOpacity(double const *const opac, size_t const naxes[3]) {
  GYOTO_DEBUG << endl;
  if (opacity_) {
    GYOTO_DEBUG << "delete [] opacity_;" << endl;
    delete [] opacity_; opacity_ = NULL;
    flag_radtransf_=0;
  }
  if (opac) {
    if (nnu_ != naxes[0] || nphi_ != naxes[1] || nr_ != naxes[2])
      GYOTO_ERROR("Please set intensity before opacity. "
		 "The two arrays must have the same dimensions.");
    GYOTO_DEBUG << "allocate opacity_;" << endl;
    opacity_ = new double[nnu_ * nphi_ * nr_];
    GYOTO_DEBUG << "opacity >> opacity_" << endl;
    memcpy(opacity_, opac, nnu_ * nphi_ * nr_ * sizeof(double));
    flag_radtransf_=1;
  }
}

double const * PatternDisk::opacity() const { return opacity_; }

void PatternDisk::copyVelocity(double const *const velocity, size_t const naxes[2]) {
  GYOTO_DEBUG << endl;

  if (velocity_) {
    GYOTO_DEBUG << "delete [] velocity_;\n";
    delete [] velocity_; velocity_ = NULL;
  }
  if (velocity) {
    if (!emission_) GYOTO_ERROR("Please use copyIntensity() before copyVelocity()");
    if (nphi_ != naxes[0] || nr_ != naxes[1])
      GYOTO_ERROR("emission_ and velocity_ have inconsistent dimensions");
    GYOTO_DEBUG << "allocate velocity_;" << endl;
    velocity_ = new double[2*nphi_*nr_];
    GYOTO_DEBUG << "velocity >> velocity_" << endl;
    memcpy(velocity_, velocity, 2*nphi_*nr_*sizeof(double));
  }
}
double const * PatternDisk::getVelocity() const { return velocity_; }

void PatternDisk::copyGridRadius(double const *const rad, size_t nr) {
  GYOTO_DEBUG << endl;
  if (radius_) {
    GYOTO_DEBUG << "delete [] radius_;" << endl;
    delete [] radius_; radius_ = NULL;
  }
  if (rad) {
    if (!emission_) GYOTO_ERROR("Please use copyIntensity() before copyGridRadius()");
    if (nr_ != nr)
      GYOTO_ERROR("emission_ and radius_ have inconsistent dimensions");
    GYOTO_DEBUG << "allocate velocity_;" << endl;
    radius_ = new double[nr_];
    GYOTO_DEBUG << "velocity >> velocity_" << endl;
    memcpy(radius_, rad, nr_*sizeof(double));
    rin_=radius_[0];
    rout_=radius_[nr_-1];
    dr_ = (rout_ - rin_) / double(nr_-1);
  }
}
double const * PatternDisk::getGridRadius() const { return radius_; }

void PatternDisk::repeatPhi(size_t n) {
  repeat_phi_ = n;
  if ((nphi_-1)*repeat_phi_>0) 
    dphi_=(phimax_-phimin_)/double((nphi_-1)*repeat_phi_);
  GYOTO_WARNING << "PatternDisk: not tested for repeat_phi_>1; "
    "check your results" << endl;
}
size_t PatternDisk::repeatPhi() const { return repeat_phi_; }

void PatternDisk::nu0(double freq) { nu0_ = freq; }
double PatternDisk::nu0() const { return nu0_; }

void PatternDisk::dnu(double dfreq) { dnu_ = dfreq; }
double PatternDisk::dnu() const { return dnu_; }

void PatternDisk::phimin(double phimn) {
  phimin_ = phimn;
  if (nphi_>1) dphi_ = (phimax_-phimin_) / double((nphi_-1)*repeat_phi_);
}
double PatternDisk::phimin() const {return phimin_;}

void PatternDisk::phimax(double phimx) {
  phimax_ = phimx;
  if (nphi_>1) dphi_ = (phimax_-phimin_) / double((nphi_-1)*repeat_phi_);
}
double PatternDisk::phimax() const {return phimax_;}

#ifdef GYOTO_USE_CFITSIO
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
  if (repeat_phi_==0)
    GYOTO_ERROR("In PatternDisk::fitsRead: repeat_phi is 0!");
  if (nphi_>1)
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
      GYOTO_ERROR("PatternDisk::readFile(): opacity array not conformable");
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
      GYOTO_ERROR("PatternDisk::readFile(): velocity array not conformable");
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
      GYOTO_ERROR("PatternDisk::readFile(): radius array not conformable");
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
     GYOTO_ERROR("In PatternDisk::fitsRead: nr_ should not be 0 here!");
   dr_ = (rout_-rin_) / double(nr_-1);

  if (fits_close_file(fptr, &status)) throwCfitsioError(status) ;
  fptr = NULL;
}

void PatternDisk::fitsWrite(string filename) {
  if (!emission_) GYOTO_ERROR("PatternDisk::fitsWrite(filename): nothing to save!");
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
#endif

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

  if (repeat_phi_>1) phi=fmod((phi-phimin_),(phimax_-phimin_)/double(repeat_phi_))+phimin_;
  // wrap phi if the disk is periodic with period repeat_phi_

  if (nphi_>1){
    if (phi<phimin_) //Possible: any phi value is in the grid anyway
      i[1]=0;
    else if (phi>phimax_)
      i[1]=nphi_;
    else {
      i[1] = size_t(floor((phi-phimin_)/dphi_)+1); // % (nphi_-1);
      if (i[1]==0 || i[1]==nphi_){
	cerr << "iphi stuff= " << phi << " " << dphi_ << " " << nphi_ << " " << floor((phi-phimin_)/dphi_) << " " << i[1] << endl;
	GYOTO_ERROR("In PatternDisk:getIndices: bad i[1]");
      }				
    }
  }else{
    i[1] = 0;
  }

  /*
    With this definition:
    phimin_+(i[1]-1)*dphi_ <= phi < phimin_+i[1]*dphi_ 
    provided phi is not <phimin_ or >phimax_
  */
  //cerr << "phi stuff: " << phimin_ << " " << phimax_ << " " << dphi_ << " " << nphi_ << endl;
  //cerr << "in indice --PHI: " << phimin_+(i[1]-1)*dphi_ << " " << phi << " " << phimin_+i[1]*dphi_  << endl;

  if (radius_) {
    GYOTO_DEBUG <<"radius_ != NULL" << endl;
    // if the radius_ vector is set, find closest value
    if (r >= radius_[nr_-1]) i[2] = nr_;
    else {
      for(i[2]=0; r > radius_[i[2]]; ++i[2]){}
      //cerr << "in indice --RAD: " << radius_[i[2]-1] << " " << r << " " << radius_[i[2]] << endl;
      /*
	With this definition:
	radius_[i[2]-1] <= r < radius_[i[2]]
	Special case: i[2]=0 if r < radius_[0]
       */
    }
  } else {
    GYOTO_DEBUG <<"radius_ == NULL, dr_==" << dr_ << endl;
    // radius_ is not set: assume linear repartition
    if (dr_==0.)
      GYOTO_ERROR("In PatternDisk::getIndices: dr_ should not be 0 here!");
    i[2] = size_t(floor((r-rin_)/dr_)+1);
    if (i[2] >= nr_) i[2] = nr_ - 1;
    //cerr << "in indice no radius --RAD: " << rin_ + (i[2]-1)*dr_ << " " << r << " " <<  rin_ + (i[2])*dr_ << endl;
  }

}

void PatternDisk::getVelocity(double const pos[4], double vel[4]) {
  //cerr << "In pattern get vel" << endl;
  if (velocity_) {
    if (dir_ != 1)
      GYOTO_ERROR("PatternDisk::getVelocity(): "
		 "dir_ should be 1 if velocity_ is provided");
    size_t i[3]; // {i_nu, i_phi, i_r}
    //cout << "in pattern getvel go to getIndices" << endl;
    getIndices(i, pos);
    //cout << "in pattern getvel after getIndices" << endl;
    
    double rr = projectedRadius(pos);
    double phi = sphericalPhi(pos);

    if (repeat_phi_>1) phi=fmod(phi-phimin_,(phimax_-phimin_)/double(repeat_phi_))+phimin_;

    if (rr < rin_ || rr > rout_){
      // Outside radial grid, emission is zero: put arbitrary velocity
      vel[0]=1;vel[1]=0;vel[2]=0;vel[3]=0;
      return;
    }

    //cout << "in velo r phi= " << rr << " " << phi << endl;
    //cout << "and indices= " << i[0] << " " << i[1] << " " << i[2] << endl;
    
    double phiprime=0., rprime=0.;
    if (nphi_==1){
      // If axisym: 1D interpo in radius only
      double rprimelow=velocity_[i[2]-1],
	rprimehigh=velocity_[i[2]],
	phiprimelow=velocity_[nr_+i[2]-1],
	phiprimehigh=velocity_[nr_+i[2]];

      double radlow, radhigh;
      if (radius_){
	radlow = radius_[i[2]-1];
	radhigh = radius_[i[2]];
      }else{
	radlow = rin_ + double(i[2]-1)*dr_;
	radhigh = rin_ + double(i[2])*dr_;
      }

      if (rr<radlow || rr>radhigh){
	//cout << "r= " << i[2] << " " <<  radlow << " " << rr << " " << radhigh << endl;
	GYOTO_ERROR("In PatternDisk::getVelocity: "
		   "bad radial interpolation");
      }

      rprime = rprimelow + (rr-radlow)/(radhigh-radlow)*(rprimehigh-rprimelow);
      phiprime = phiprimelow
	+ (rr-radlow)/(radhigh-radlow)*(phiprimehigh-phiprimelow);
    }else{
      // Bilinear interpolation in r,phi
      // Notation: X_{phi,r}
      int iphil, iphiu;
      double philow, phihigh;
      if ((i[1]==0 || i[1]==nphi_) && repeat_phi_==1){
	// then phi is below phimin or above phimax,
	// i.e. phimax < phi[2pi] < phimin+2pi
	// [I don't want to code the repeat_phi_>1 case...]
	iphil = nphi_-1;
	philow = phimin_+double(iphil)*dphi_;
	iphiu = 0;
	phihigh = phimin_ + 2.*M_PI;
	if (phi<phimin_)
	  phi+=2*M_PI; // to get phimax < phi < phimin+2pi
      }else{ // standard case
        iphil = i[1]-1;
        philow = phimin_+double(iphil)*dphi_;
        iphiu = i[1];
        phihigh = phimin_+double(iphiu)*dphi_;
      }

      double radlow, radhigh;
      if (radius_){
	radlow = radius_[i[2]-1];
	radhigh = radius_[i[2]];
      }else{
	radlow = rin_ + double(i[2]-1)*dr_;
	radhigh = rin_ + double(i[2])*dr_;
      }
	    
      double phip00=velocity_[nr_*nphi_+iphil*nr_+(i[2]-1)];
      double phip10=velocity_[nr_*nphi_+iphiu*nr_+(i[2]-1)];
      double phip11=velocity_[nr_*nphi_+iphiu*nr_+i[2]];
      double phip01=velocity_[nr_*nphi_+iphil*nr_+i[2]];
      
      double rp00=velocity_[iphil*nr_+(i[2]-1)];
      double rp10=velocity_[iphiu*nr_+(i[2]-1)];
      double rp11=velocity_[iphiu*nr_+i[2]];
      double rp01=velocity_[iphil*nr_+i[2]];

      //cout << "velo bilin interpo, phip, rp: " << phip00 << " " << rp00 << " " << phip01 << " " << rp01 << " " << phip10 << " " << rp10 << " " << phip11 << " " << rp11 << endl;


      //cout << "rin sup, phi inf sup, r phi: " << radlow << " " << radhigh << " " << philow << " " << phihigh << " " << rr << " " << phi << endl;

      if (phi<philow || phi>phihigh || rr<radlow || rr>radhigh){
	//cout << "r, phis= " << i[2] << " " <<  radlow << " " << rr << " " << radhigh << " " << philow << " " << phi << " " << phihigh << endl;
	GYOTO_ERROR("In PatternDisk::getVelocity: "
		   "bad interpolation");
      }

      double cr = (rr-radlow)/(radhigh-radlow),
	cp = (phi-philow)/(phihigh-philow);

      rprime=rp00 + cp*(rp10-rp00) + cr*(rp01-rp00)
	+ cr*cp*(rp11-rp01+rp00-rp10);

      phiprime=phip00 + cp*(phip10-phip00) + cr*(phip01-phip00)
	+ cr*cp*(phip11-phip01+phip00-phip10);

      //cout << "interpol velo phi, rp= " << phiprime << " " << rprime << endl;
    }
	
    switch (gg_->coordKind()) {
    case GYOTO_COORDKIND_SPHERICAL:
      {
	vel[1] = rprime;
	vel[2] = 0.;
	vel[3] = phiprime;
	//cout << "pos and vel= " << pos[1] << " " << pos[2] << " " << pos[3] << " " << vel[1] << " " << vel[2] << " " << vel[3] << endl;
	vel[0] = gg_->SysPrimeToTdot(pos, vel+1);
	vel[1] *= vel[0];
	vel[3] *= vel[0];
      }
      break;
    case GYOTO_COORDKIND_CARTESIAN:
      GYOTO_ERROR("PatternDisk::getVelocity(): metric must be in "
		 "spherical coordinaites if velocity field is provided");
      break;
    default:
      GYOTO_ERROR("PatternDisk::getVelocity(): unknown COORDKIND");
    } 
  }else ThinDisk::getVelocity(pos, vel);
}

double PatternDisk::emission(double nu, double dsem,
				    state_t const &,
				    double const co[8]) const{
  //See Page & Thorne 74 Eqs. 11b, 14, 15. This is F(r).
  GYOTO_DEBUG << endl;
  size_t i[3]; // {i_nu, i_phi, i_r}
  getIndices(i, co, nu);

  double rr = projectedRadius(co);
  double phi = sphericalPhi(co);

  if (repeat_phi_>1) phi=fmod(phi-phimin_,(phimax_-phimin_)/double(repeat_phi_))+phimin_;

  //cout << "in emission r, phi= " << rr << " " << phi << endl;
  //cout << "and indices= " << i[0] << " " << i[1] << " " << i[2] << endl;
  //cout << "dphi= " << dphi_ << endl;
  // No emission outside radius limits
  //cerr << "r checks: " << rin_ << " " << rout_ << endl;
  if (rr < rin_ || rr > rout_) return 0.;
  double Iem=0.;

  if (nnu_>1) GYOTO_ERROR("In PatternDisk: multifrequency case not implemented");
  
  if (nphi_==1){
    // If axisym: 1D interpo in radius only
    double Iemlow = emission_[i[2]-1],
      Iemhigh = emission_[i[2]];

    double radlow, radhigh;
    if (radius_){
      radlow = radius_[i[2]-1];
      radhigh = radius_[i[2]];
    }else{
      radlow = rin_ + double(i[2]-1)*dr_;
      radhigh = rin_ + double(i[2])*dr_;
    }

    if (rr<radlow || rr>radhigh){
      //cout << "r= " << i[2] << " " <<  radlow << " " << rr << " " << radhigh << endl;
      GYOTO_ERROR("In PatternDisk::emission: "
		 "bad radial interpolation");
    }
    
    Iem = Iemlow + (rr-radlow)/(radhigh-radlow)*(Iemhigh-Iemlow);
    
    //cout << "In emission only rad: " << radlow << " " << radhigh << " " << Iemlow << " " << Iemhigh << " " << Iem << endl;
  }else{
    // Bilinear interpolation
    // Notation I_{phi,r}
    //cout << "entering bilin em in pattern" << endl;
    int iphil, iphiu;
    double philow, phihigh;
    if ((i[1]==0 || i[1]==nphi_) && repeat_phi_==1){
      // then phi is below phimin or above phimax,
      // i.e. phimax < phi[2pi] < phimin+2pi
      // [I don't want to code the repeat_phi_>1 case...]
      iphil = nphi_-1;
      philow = phimin_+double(iphil)*dphi_;
      iphiu = 0;
      phihigh = phimin_ + 2.*M_PI;
      if (phi<phimin_)
	phi+=2*M_PI; // to get phimax < phi < phimin+2pi
    }else{ // standard case
      iphil = i[1]-1;
      philow = phimin_+double(iphil)*dphi_;
      iphiu = i[1];
      phihigh = phimin_+double(iphiu)*dphi_;
    }
    //cout << "iphi: " << i[1] << " " << iphil << " " << iphiu << endl;

    double radlow, radhigh;
    if (radius_){
      radlow = radius_[i[2]-1];
      radhigh = radius_[i[2]];
    }else{
      radlow = rin_ + double(i[2]-1)*dr_;
      radhigh = rin_ + double(i[2])*dr_;
    }
    
    double I00=emission_[iphil*nr_+(i[2]-1)];
    double I10=emission_[iphiu*nr_+(i[2]-1)];
    double I11=emission_[iphiu*nr_+i[2]];
    double I01=emission_[iphil*nr_+i[2]];
    //cout << "In emission Iem grid= " << I00 << " " << I01 << " " << I10 << " " << I11 << endl;
    
    //cout << " In emission rin sup, phi inf sup, r phi: " << radlow << " " << radhigh << " " << philow << " " << phihigh << " " << rr << " " << phi << endl;
    
    if (phi<philow || phi>phihigh || rr<radlow || rr>radhigh){
      cout << "phi: " << philow << " " << phi << " " << phihigh << endl;
      cout << "r: " << radlow << " " << rr << " " << radhigh << endl;
      GYOTO_ERROR("In PatternDisk::emission: "
		 "bad interpolation");
    }
    
    double cr = (rr-radlow)/(radhigh-radlow),
      cp = (phi-philow)/(phihigh-philow);
    
    Iem = I00 + cp*(I10-I00) + cr*(I01-I00) + cr*cp*(I11-I01+I00-I10);
    //cout << "In emission I interpo= " << Iem << endl;
  }
  
  if (!flag_radtransf_) return Iem;
  double thickness;
  // NB: thickness is not interpolated so far
  if (opacity_ && (thickness=opacity_[i[1]*nr_+i[2]]*dsem))
    return Iem * (1. - exp (-thickness)) ;
  return 0.;
}

double PatternDisk::transmission(double nu, double dsem, state_t const &, double const co[8]) const {
  GYOTO_DEBUG << endl;
  if (!flag_radtransf_) return 0.;
  if (!opacity_) return 1.;
  size_t i[3]; // {i_nu, i_phi, i_r}
  getIndices(i, co, nu);
  // NB: opacity is not interpolated so far
  double opac = opacity_[i[1]*nr_+i[2]];
  GYOTO_DEBUG << "nu="<<nu <<", dsem="<<dsem << ", opacity="<<opac <<endl;
  if (!opac) return 1.;
  return exp(-opac*dsem);
}

void PatternDisk::innerRadius(double rin) {
  ThinDisk::innerRadius(rin);
  if (nr_>1 && !radius_) dr_ = (rout_-rin_) / double(nr_-1);
}

void PatternDisk::outerRadius(double rout) {
  ThinDisk::outerRadius(rout);
  if (nr_>1 && !radius_) dr_ = (rout_-rin_) / double(nr_-1);
}

void PatternDisk::patternVelocity(double omega) { Omega_ = omega; }
double PatternDisk::patternVelocity() const { return Omega_; }

void PatternDisk::file(std::string const &f) {
# ifdef GYOTO_USE_CFITSIO
  fitsRead(f);
# else
  GYOTO_ERROR("This Gyoto has no FITS i/o");
# endif
}

std::string PatternDisk::file() const {
  return filename_;
}
