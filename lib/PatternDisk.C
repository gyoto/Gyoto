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
  emission_(NULL), velocity_(NULL), radius_(NULL),
  Omega_(0.), t0_(0.),
  dnu_(0.), nu0_(0), nnu_(0),
  dphi_(0.), nphi_(0), repeat_phi_(1),
  dr_(0.), nr_(0)
{
  if (debug()) cerr << "DEBUG: PatternDisk Construction" << endl;
}

PatternDisk::PatternDisk(const PatternDisk& o) :
  ThinDisk(o), filename_(o.filename_),
  emission_(NULL), velocity_(NULL), radius_(NULL),
  Omega_(o.Omega_), t0_(o.t0_),
  dnu_(o.dnu_), nu0_(o.nu0_), nnu_(o.nnu_),
  dphi_(o.dphi_), nphi_(o.nphi_), repeat_phi_(o.repeat_phi_),
  dr_(o.dr_), nr_(o.nr_)
{
  if (debug()) cerr << "DEBUG: PatternDisk Copy" << endl;
  size_t ncells = 0;
  if (o.emission_) {
    emission_ = new double[ncells = nnu_ * nphi_ * nr_];
    memcpy(emission_, o.emission_, ncells * sizeof(double));
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
  if (debug()) cerr << "DEBUG: PatternDisk Destruction" << endl;
  if (emission_) delete [] emission_;
  if (velocity_) delete [] velocity_;
  if (radius_) delete [] radius_;
}

void PatternDisk::readFile(string filename) {
  if (verbose() >= GYOTO_QUIET_VERBOSITY)
    cout << "PatternDisk reading FITS file: " << filename << endl;

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
#define throwCfitsioError(status) \
    { fits_get_errstatus(status, ermsg); throwError(ermsg); }

  if (debug()) cerr << "PatternDisk::readFile(): opening file" << endl;
  if (fits_open_file(&fptr, pixfile, 0, &status)) throwCfitsioError(status) ;

  ////// READ FITS KEYWORDS COMMON TO ALL TABLES ///////
  //get Omega and t0;
  if (debug()) cerr << "PatternDisk::readFile(): read Omega_" << endl;
  fits_read_key(fptr, TDOUBLE, "GYOTO PatternDisk Omega", &tmpd, NULL, &status);
  if (status) {
    if (status == KEY_NO_EXIST) status = 0; // not fatal
    else throwCfitsioError(status) ;
  } else Omega_ = tmpd; // Omega_ found
  if (debug()) cerr << "PatternDisk::readFile(): read t0_" << endl;
  fits_read_key(fptr, TDOUBLE, "GYOTO PatternDisk t0", &tmpd, NULL, &status);
  if (status) {
    if (status == KEY_NO_EXIST) status = 0; // not fatal
    else throwCfitsioError(status) ;
  } else t0_ = tmpd; // T0_ found

  if (debug()) cerr << "PatternDisk::readFile(): read RepeatPhi_" << endl;
  fits_read_key(fptr, TLONG, "GYOTO PatternDisk RepeatPhi", &tmpl,
		NULL, &status);
  if (status) {
    if (status == KEY_NO_EXIST) status = 0; // not fatal
    else throwCfitsioError(status) ;
  } else repeat_phi_ = size_t(tmpl); // RepeatPhi found

  if (debug()) cerr << "PatternDisk::readFile(): read InnerRadius_" << endl;
  fits_read_key(fptr, TDOUBLE, "GYOTO ThinDisk InnerRadius", &tmpd,
		NULL, &status);
  if (status) {
    if (status == KEY_NO_EXIST) status = 0; // not fatal
    else throwCfitsioError(status) ;
  } else {
    rin_ = tmpd; // InnerRadius found
    rin_set=1;
  }
  if (debug()) cerr << "PatternDisk::readFile(): read OuterRadius_" << endl;
  fits_read_key(fptr, TDOUBLE, "GYOTO ThinDisk OuterRadius", &tmpd,
		NULL, &status);
  if (status) {
    if (status == KEY_NO_EXIST) status = 0; // not fatal
    else throwCfitsioError(status) ;
  } else {
    rout_ = tmpd; // OuterRadius found
    rout_set=1;
  }

  ////// FIND MANDATORY EMISSION HDU, READ KWDS & DATA ///////
  if (debug()) cerr << "PatternDisk::readFile(): search emission HDU" << endl;
  if (fits_movnam_hdu(fptr, ANY_HDU,
		      const_cast<char*>("GYOTO PatternDisk emission"),
		      0, &status))
    throwCfitsioError(status) ;
  if (debug()) cerr << "PatternDisk::readFile(): get image size" << endl;
  if (fits_get_img_size(fptr, 3, naxes, &status)) throwCfitsioError(status) ;

  //update nu0_, nnu_, dnu_;
  nnu_ = naxes[0]; 
  double CRPIX1;
  if (debug()) cerr << "PatternDisk::readFile(): read CRPIX1, CRVAL1, CDELT1"
		    << endl;
  fits_read_key(fptr, TDOUBLE, "CRVAL1", &nu0_, NULL, &status);
  fits_read_key(fptr, TDOUBLE, "CDELT1", &dnu_, NULL, &status);
  fits_read_key(fptr, TDOUBLE, "CRPIX1", &CRPIX1, NULL, &status);
  if (status) throwCfitsioError(status) ;
  if (CRPIX1 != 1) nu0_ -= dnu_*(CRPIX1 - 1.);

  // update repeat_phi_, nphi_, dphi_
  nphi_ = naxes[1];
  dphi_ = 2.*M_PI/double(nphi_*repeat_phi_);

  // update rin_, rout_, nr_, dr_
  nr_ = naxes[2];
  dr_ = (rout_-rin_) / nr_;

  if (emission_) delete [] emission_; emission_=NULL;
  emission_ = new double[nnu_ * nphi_ * nr_];
  if (debug())
    cerr << "PatternDisk::readFile(): read emission: "
	 << "nnu_=" << nnu_ << ", nphi_="<<nphi_ << ", nr_="<<nr_ << "...";
  if (fits_read_subset(fptr, TDOUBLE, fpixel, naxes, inc,
		       0, emission_,&anynul,&status)) {
    if (debug()) cerr << " error, trying to free pointer" << endl;
    delete [] emission_; emission_=NULL;
    throwCfitsioError(status) ;
  }
  if (debug()) cerr << " done." << endl;

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
    dr_=0.;
    if (!rin_set) rin_=radius_[0];
    if (!rout_set) rout_=radius_[nr_-1];
  }
 

  if (fits_close_file(fptr, &status)) throwCfitsioError(status) ;
  fptr = NULL;
}

void PatternDisk::getIndices(size_t i[3], double const co[4], double nu) const {
  if (nu <= nu0_) i[0] = 0;
  else {
    i[0] = size_t((nu-nu0_)/dnu_);
    if (i[0] >= nnu_) i[0] = nnu_-1;
  }

  double r = projectedRadius(co);
  double phi = sphericalPhi(co);
  double t = co[0];

  i[1] = size_t((phi-Omega_*(t-t0_))/dphi_) % nphi_;

  if (radius_) {
    // if the radius_ vector is set, find closest value
    if (r >= radius_[nr_-1]) i[2] = nr_-1;
    else {
      for(i[2]=0; r > radius_[i[2]]; ++i[2]){}
      if (i[2]>0 && r-radius_[i[2]-1] < radius_[i[2]]) --i[2];
    }
  } else {
    // radius_ is not set: assume linear repartition
    i[2] = size_t((r-rin_)/dr_);
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
    double phiprime=velocity_[i[2]*(nphi_*2)+i[1]*2+0];
    double rprime=velocity_[i[2]*(nphi_*2)+i[1]*2+1];
    switch (gg_->getCoordKind()) {
    case GYOTO_COORDKIND_SPHERICAL:
      {
	double pos2[4] = {pos[0], pos[1], pos[2], pos[3]};
	pos2[1] = radius_ ? radius_[i[2]] : rin_+i[2]*dr_;
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
  }else gg_->circularVelocity(pos, vel);
}

double PatternDisk::emission(double nu, double dsem,
				    double *,
				    double co[8]) const{
  //See Page & Thorne 74 Eqs. 11b, 14, 15. This is F(r).
  size_t i[3]; // {i_nu, i_phi, i_r}
  getIndices(i, co, nu);
  double Iem = emission_[i[2]*(nphi_*nnu_)+i[1]*nnu_+i[0]];

  if (flag_radtransf_) Iem *= dsem;

  return Iem;
  
}

void PatternDisk::setInnerRadius(double rin) {
  ThinDisk::setInnerRadius(rin);
  if (nr_ && !radius_) dr_ = (rout_-rin_) / nr_;
}

void PatternDisk::setOuterRadius(double rout) {
  ThinDisk::setOuterRadius(rout);
  if (nr_ && !radius_) dr_ = (rout_-rin_) / nr_;
}

void PatternDisk::setPatternVelocity(double omega) { Omega_ = omega; }

#ifdef GYOTO_USE_XERCES
void PatternDisk::fillElement(FactoryMessenger *fmp) const {
  fmp->setParameter("File", filename_);
  if (Omega_) fmp->setParameter("PatternVelocity", Omega_);
  ThinDisk::fillElement(fmp);
}

SmartPointer<Astrobj::Generic> PatternDisk::Subcontractor(FactoryMessenger* fmp) {
  string name, content;
  SmartPointer<PatternDisk> ao = new PatternDisk();
  ao -> setMetric(fmp->getMetric());

  while (fmp->getNextParameter(&name, &content)) {
    if      (name == "File") ao -> readFile( fmp -> fullPath(content) );
    else if (name=="PatternVelocity")
      ao -> setPatternVelocity(atof(content.c_str()));
    else ao -> setGenericParameter(name, content);
  }

  return ao;
}

void PatternDisk::Init() {
  Astrobj::Register("PatternDisk", &Subcontractor);
}
#endif
