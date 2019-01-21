/*
  Copyright 2017, 2018 Frederic Vincent, Thibaut Paumard

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
#include "GyotoPhoton.h"
#include "GyotoXillverReflection.h"
#include "GyotoFactoryMessenger.h"
#include "GyotoKerrBL.h"
#include "GyotoKerrKS.h"

//Std headers
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <cmath>
#include <limits>
#include <cstring>
#include <sstream>

#ifdef GYOTO_USE_CFITSIO
#include <fitsio.h>
#define throwCfitsioError(status) \
    { fits_get_errstatus(status, ermsg); GYOTO_ERROR(ermsg); }
#endif

using namespace std;
using namespace Gyoto;
using namespace Gyoto::Astrobj;

/// Properties
#include "GyotoProperty.h"
GYOTO_PROPERTY_START(XillverReflection,
		     "Xillver reflection accretion disk.")
GYOTO_PROPERTY_FILENAME(XillverReflection, FileIllumination, fileillumination)
GYOTO_PROPERTY_FILENAME(XillverReflection, FileReflection, filereflection)
GYOTO_PROPERTY_DOUBLE(XillverReflection, LampRadius, lampradius)
GYOTO_PROPERTY_DOUBLE(XillverReflection, TimeLampPhiZero, timelampphizero)
GYOTO_PROPERTY_BOOL(XillverReflection,
		    AverageOverAngle, DontAverageOverAngle,
		    averageOverAngle)
GYOTO_PROPERTY_END(XillverReflection, ThinDisk::properties)

void XillverReflection::fillProperty(Gyoto::FactoryMessenger *fmp,
				     Property const &p) const {
  if (p.name == "FileIllumination")
    fmp->setParameter("FileIllumination", (filenameIllum_.compare(0,1,"!") ?
					   filenameIllum_ :
					   filenameIllum_.substr(1)) );
  else if (p.name == "FileReflection")
    fmp->setParameter("FileReflection", (filenameRefl_.compare(0,1,"!") ?
					 filenameRefl_ :
					 filenameRefl_.substr(1)) );
  else ThinDisk::fillProperty(fmp, p);
}

XillverReflection::XillverReflection() :
  ThinDisk("XillverReflection"), filenameIllum_(""), filenameRefl_(""),
  lampradius_(0), timelampphizero_(0.),
  aa_(0.),
  illumination_(NULL), reflection_(NULL),
  radius_(NULL), phi_(NULL),
  logxi_(NULL), incl_(NULL), freq_(NULL),
  nnu_(0), ni_(0), nxi_(0), nr_(0), nphi_(0),
  average_over_angle_(0) {
  
  GYOTO_DEBUG << endl;

}

XillverReflection::XillverReflection(const XillverReflection& o) :
  ThinDisk(o), filenameIllum_(o.filenameIllum_), filenameRefl_(o.filenameRefl_),
  lampradius_(o.lampradius_), timelampphizero_(o.timelampphizero_),
  aa_(o.aa_),
  illumination_(NULL), reflection_(NULL),
  radius_(NULL), phi_(NULL),
  logxi_(NULL), incl_(NULL), freq_(NULL),
  nnu_(o.nnu_), ni_(o.ni_), nxi_(o.nxi_),
  nr_(o.nr_), nphi_(o.nphi_), 
  average_over_angle_(o.average_over_angle_)
{
  GYOTO_DEBUG << endl;
  size_t ncells = 0;
  if (o.illumination_) {
    illumination_ = new double[ncells = nr_ * nphi_];
    memcpy(illumination_, o.illumination_, ncells * sizeof(double));
  }
  if (o.reflection_) {
    reflection_ = new double[ncells = nnu_ * ni_ * nxi_];
    memcpy(reflection_, o.reflection_, ncells * sizeof(double));
  }
  if (o.freq_) {
    freq_ = new double[ncells = nnu_];
    memcpy(freq_, o.freq_, ncells * sizeof(double));
  }
  if (o.incl_) {
    incl_ = new double[ncells = ni_];
    memcpy(incl_, o.incl_, ncells * sizeof(double));
  }
  if (o.logxi_) {
    logxi_ = new double[ncells = nxi_];
    memcpy(logxi_, o.logxi_, ncells * sizeof(double));
  }
  if (o.radius_) {
    radius_ = new double[ncells = nr_];
    memcpy(radius_, o.radius_, ncells * sizeof(double));
  }
  if (o.phi_) {
    phi_ = new double[ncells = nphi_];
    memcpy(phi_, o.phi_, ncells * sizeof(double));
  }

}

XillverReflection * XillverReflection::clone() const {
  return new XillverReflection(*this); }

XillverReflection::~XillverReflection() {
  GYOTO_DEBUG << endl;
  if (illumination_) delete [] illumination_;
  if (reflection_) delete [] reflection_;
  if (freq_) delete [] freq_;
  if (incl_) delete [] incl_;
  if (logxi_) delete [] logxi_;
  if (radius_) delete [] radius_;
  if (phi_) delete [] phi_;
}

void XillverReflection::metric(SmartPointer<Metric::Generic> gg) {
  if (gg_) gg_->unhook(this);
  string kin = gg->kind();
  if (kin != "KerrBL" && kin != "KerrKS")
    GYOTO_ERROR
      ("Xillver::metric(): metric must be KerrBL or KerrKS");
  Generic::metric(gg);
  updateSpin();
  gg->hook(this);
}

void XillverReflection::updateSpin() {
  if (!gg_) return;
  switch (gg_->coordKind()) {
  case GYOTO_COORDKIND_SPHERICAL:
    aa_ = static_cast<SmartPointer<Metric::KerrBL> >(gg_) -> spin();
    break;
  case GYOTO_COORDKIND_CARTESIAN:
    aa_ = static_cast<SmartPointer<Metric::KerrKS> >(gg_) -> spin();
    break;
  default:
    GYOTO_ERROR("Xillver::updateSpin(): unknown COORDKIND");
  }
}

void XillverReflection::tell(Hook::Teller* msg) {
  if (msg==gg_) updateSpin();
}

// Next 2 function are probably useless
void XillverReflection::setReflection(double * pattern) {
  reflection_ = pattern;
}
void XillverReflection::setIllumination(double * pattern) {
  illumination_ = pattern;
}

void XillverReflection::copyReflection(double const *const pattern,
				       size_t const naxes[3]) {
  GYOTO_DEBUG << endl;
  if (reflection_) {
    GYOTO_DEBUG << "delete [] reflection_;" << endl;
    delete [] reflection_; reflection_ = NULL;
  }
  if (pattern) {
    size_t nel;
    if (nnu_ != naxes[0]) {
      GYOTO_DEBUG <<"nnu_ changed, freeing freq_" << endl;
      if (freq_)  { delete [] freq_; freq_  = NULL; }
    }
    if (ni_ != naxes[1]) {
      GYOTO_DEBUG <<"ni_ changed, freeing freq_ and incl_" << endl;
      if (freq_)  { delete [] freq_; freq_  = NULL; }
      if (incl_) { delete [] incl_; incl_= NULL; }
    }
    if (nxi_ != naxes[2]) {
      GYOTO_DEBUG <<"nxi_ changed, freeing freq_, incl_ and logxi_" << endl;
      if (freq_)  { delete [] freq_; freq_  = NULL; }
      if (incl_) { delete [] incl_; incl_= NULL; }
      if (logxi_)   { delete [] logxi_;   logxi_  = NULL; }
    }
    if (!(nel=(nnu_ = naxes[0]) * (ni_=naxes[1]) * (nxi_=naxes[2])))
      GYOTO_ERROR( "dimensions can't be null");
    GYOTO_DEBUG << "allocate reflection_;" << endl;
    reflection_ = new double[nel];
    GYOTO_DEBUG << "pattern >> reflection_" << endl;
    memcpy(reflection_, pattern, nel*sizeof(double));
  }
}
double const * XillverReflection::getReflection() const {
  return reflection_; }
void XillverReflection::getReflectionNaxes( size_t naxes[3] ) const
{ naxes[0] = nnu_; naxes[1] = ni_; naxes[2] = nxi_; }

void XillverReflection::copyIllumination(double const *const pattern,
					 size_t const naxes[2]) {
  GYOTO_DEBUG << endl;
  if (illumination_) {
    GYOTO_DEBUG << "delete [] illumination_;" << endl;
    delete [] illumination_; illumination_ = NULL;
  }
  if (pattern) {
    size_t nel;
    if (nr_ != naxes[0]) {
      GYOTO_DEBUG <<"nr_ changed, freeing radius_" << endl;
      if (radius_)  { delete [] radius_; radius_  = NULL; }
    }
    if (nphi_ != naxes[1]) {
      GYOTO_DEBUG <<"nphi_ changed, freeing radius_ and phi_" << endl;
      if (radius_)  { delete [] radius_; radius_  = NULL; }
      if (phi_) { delete [] phi_; phi_= NULL; }
    }

    if (!(nel=(nr_ = naxes[0]) * (nphi_=naxes[1])))
      GYOTO_ERROR( "dimensions can't be null");
    GYOTO_DEBUG << "allocate illumination_;" << endl;
    illumination_ = new double[nel];
    GYOTO_DEBUG << "pattern >> illumination_" << endl;
    memcpy(illumination_, pattern, nel*sizeof(double));
  }
}
double const * XillverReflection::getIllumination() const {
  return illumination_; }
void XillverReflection::getIlluminationNaxes( size_t naxes[2] ) const
{ naxes[0] = nr_; naxes[1] = nphi_; }

void XillverReflection::copyGridReflLogxi(double const *const lxi,
					  size_t nxi) {
  GYOTO_DEBUG << endl;
  if (logxi_) {
    GYOTO_DEBUG << "delete [] logxi_;" << endl;
    delete [] logxi_; logxi_ = NULL;
  }
  if (lxi) {
    if (!reflection_) 
      GYOTO_ERROR("Please use copyReflection() before copyGridReflLogxi()");
    if (nxi_ != nxi)
      GYOTO_ERROR("reflection_ and logxi_ have inconsistent dimensions");
    GYOTO_DEBUG << "allocate logxi_;" << endl;
    logxi_ = new double[nxi_];
    GYOTO_DEBUG << "logxi >> logxi_" << endl;
    memcpy(logxi_, lxi, nxi_*sizeof(double));
  }
}
double const * XillverReflection::getGridReflLogxi() const {
  return logxi_; }

void XillverReflection::copyGridReflIncl(double const *const incl, size_t ni) {
  GYOTO_DEBUG << endl;
  if (incl_) {
    GYOTO_DEBUG << "delete [] incl_;" << endl;
    delete [] incl_; incl_ = NULL;
  }
  if (incl) {
    if (!reflection_) 
      GYOTO_ERROR("Please use copyReflection() before copyGridReflIncl()");
    if (ni_ != ni)
      GYOTO_ERROR("reflection_ and incl_ have inconsistent dimensions");
    GYOTO_DEBUG << "allocate incl_;" << endl;
    incl_ = new double[ni_];
    GYOTO_DEBUG << "incl >> incl_" << endl;
    memcpy(incl_, incl, ni_*sizeof(double));
  }
}
double const * XillverReflection::getGridReflIncl() const { return incl_; }

void XillverReflection::copyGridReflFreq(double const *const freq,
					 size_t nnu) {
  GYOTO_DEBUG << endl;
  if (freq_) {
    GYOTO_DEBUG << "delete [] freq_;" << endl;
    delete [] freq_; freq_ = NULL;
  }
  if (freq) {
    if (!reflection_) 
      GYOTO_ERROR("Please use copyReflection() before copyGridReflFreq()");
    if (nnu_ != nnu)
      GYOTO_ERROR("reflection_ and freq_ have inconsistent dimensions");
    GYOTO_DEBUG << "allocate freq_;" << endl;
    freq_ = new double[nnu_];
    GYOTO_DEBUG << "freq >> freq_" << endl;
    memcpy(freq_, freq, nnu_*sizeof(double));
  }
}
double const * XillverReflection::getGridReflFreq() const { return freq_; }

void XillverReflection::copyGridIllumRadius(double const *const radius,
					    size_t nr) {
  GYOTO_DEBUG << endl;
  if (radius_) {
    GYOTO_DEBUG << "delete [] radius_;" << endl;
    delete [] radius_; radius_ = NULL;
  }
  if (radius) {
    if (!illumination_) 
      GYOTO_ERROR("Please use copyIllumination() before copyGridIllumRadius()");
    if (nr_ != nr)
      GYOTO_ERROR("illumination_ and radius_ have inconsistent dimensions");
    GYOTO_DEBUG << "allocate radius_;" << endl;
    radius_ = new double[nr_];
    GYOTO_DEBUG << "radius >> radius_" << endl;
    memcpy(radius_, radius, nr_*sizeof(double));
  }
}
double const * XillverReflection::getGridIllumRadius() const { return radius_; }

void XillverReflection::copyGridIllumPhi(double const *const phi,
					 size_t nphi) {
  GYOTO_DEBUG << endl;
  if (phi_) {
    GYOTO_DEBUG << "delete [] phi_;" << endl;
    delete [] phi_; phi_ = NULL;
  }
  if (phi) {
    if (!illumination_) 
      GYOTO_ERROR("Please use copyIllumination() before copyGridIllumPhi()");
    if (nphi_ != nphi)
      GYOTO_ERROR("illumination_ and phi_ have inconsistent dimensions");
    GYOTO_DEBUG << "allocate phi_;" << endl;
    phi_ = new double[nphi_];
    GYOTO_DEBUG << "phi >> phi_" << endl;
    memcpy(phi_, phi, nphi_*sizeof(double));
  }
}
double const * XillverReflection::getGridIllumPhi() const { return phi_; }

void XillverReflection::averageOverAngle(bool t) {
  average_over_angle_=t;}
bool XillverReflection::averageOverAngle()const {
  return average_over_angle_;}

void XillverReflection::fileillumination(std::string const &f) {
# ifdef GYOTO_USE_CFITSIO
  fitsReadIllum(f);
# else
  GYOTO_ERROR("This Gyoto has no FITS i/o");
# endif
}

double XillverReflection::timelampphizero() const{
  return timelampphizero_;
}

void XillverReflection::timelampphizero(double tt){
  if (lampradius_==0.)
    {
      GYOTO_ERROR("In Xillver::timelempphizero: "
		 "update lampradius before timelampphizero.");
    }
  double lampperiod = 2*M_PI*(pow(lampradius_,1.5)+aa_);
  timelampphizero_=fmod(tt,lampperiod);
}

double XillverReflection::lampradius() const{
  return lampradius_;
}

void XillverReflection::lampradius(double rr){
  lampradius_=rr;
}

void XillverReflection::filereflection(std::string const &f) {
# ifdef GYOTO_USE_CFITSIO
  fitsReadRefl(f);
# else
  GYOTO_ERROR("This Gyoto has no FITS i/o");
# endif
}

std::string XillverReflection::fileillumination() const {
  return filenameIllum_;
}

std::string XillverReflection::filereflection() const {
  return filenameRefl_;
}

#ifdef GYOTO_USE_CFITSIO
void XillverReflection::fitsReadIllum(string filenameIllum) {
  GYOTO_MSG << "XillverReflection reading FITS files: " <<
    filenameIllum << endl;

  filenameIllum_ = filenameIllum;
  char*     pixfileI   = const_cast<char*>(filenameIllum_.c_str());
  fitsfile* fptrI      = NULL;
  int       statusI    = 0;
  int       anynulI    = 0;
  long      naxesI []  = {1,1};
  long      fpixelI[]  = {1,1};
  long      incI   []  = {1,1};
  char      ermsg[31] = ""; // ermsg is used in throwCfitsioError()

  GYOTO_DEBUG << "XillverReflection::readFile(): opening illum file" << endl;
  if (fits_open_file(&fptrI, pixfileI, 0, &statusI)) throwCfitsioError(statusI) ;
  ////// FIND MANDATORY ILLUMINATION HDU, READ KWDS & DATA ///////
  GYOTO_DEBUG << "XillverReflection::readFile(): "
    "search illumination HDU" << endl;
  if (fits_movnam_hdu(fptrI, ANY_HDU,
		      const_cast<char*>("GYOTO XillverReflection illumination"),
		      0, &statusI))
    throwCfitsioError(statusI) ;
  GYOTO_DEBUG << "XillverReflection::readFile(): get image size" << endl;
  if (fits_get_img_size(fptrI, 2, naxesI, &statusI)) throwCfitsioError(statusI) ;

  //update nr_, nphi_
  nr_ = naxesI[0]; 
  nphi_  = naxesI[1];
  
  if (illumination_) { delete [] illumination_; illumination_ = NULL; }
  illumination_ = new double[nr_ * nphi_];
  if (debug())
    cerr << "XillverReflection::readFile(): read illumination: "
	 << "nr_=" << nr_ << ", nphi_="<<nphi_ << "...";
  if (fits_read_subset(fptrI, TDOUBLE, fpixelI, naxesI, incI,
		       0, illumination_,&anynulI,&statusI)) {
    GYOTO_DEBUG << " error, trying to free pointer" << endl;
    delete [] illumination_; illumination_=NULL;
    throwCfitsioError(statusI) ;
  }
  GYOTO_DEBUG << " done." << endl;
  
  // double minemission=DBL_MAX, maxemission=DBL_MIN;
  // for (int myi=0;myi<nnu_ * ni_ * nsg_-1;myi++){
  //   if (emission_[myi]<minemission) minemission=emission_[myi];
  //   if (emission_[myi]>maxemission) maxemission=emission_[myi];
  // }
  //cout << "In XillverRefl::fitsRead: Min and max emission= " <<
  //  minemission << " " << maxemission << endl;
  
  ////// FIND MANDATORY ILLUM::RADIUS HDU ///////
  
  if (fits_movnam_hdu(fptrI, ANY_HDU,
		      const_cast<char*>("GYOTO XillverReflection radius"),
		      0, &statusI))
    throwCfitsioError(statusI) ;
  if (fits_get_img_size(fptrI, 1, naxesI, &statusI)) throwCfitsioError(statusI) ;
  if (size_t(naxesI[0]) != nr_)
    GYOTO_ERROR("XillverReflection::readFile(): radius array not conformable");
  if (radius_) { delete [] radius_; radius_ = NULL; }
  radius_ = new double[nr_];
  if (fits_read_subset(fptrI, TDOUBLE, fpixelI, naxesI, incI, 
		       0, radius_,&anynulI,&statusI)) {
    delete [] radius_; radius_=NULL;
    throwCfitsioError(statusI) ;
  }
  
  ////// FIND MANDATORY ILLUM::PHI HDU ///////
  
  if (fits_movnam_hdu(fptrI, ANY_HDU,
		      const_cast<char*>("GYOTO XillverReflection phi"),
		      0, &statusI))
    throwCfitsioError(statusI) ;
  if (fits_get_img_size(fptrI, 1, naxesI, &statusI)) throwCfitsioError(statusI) ;
  if (size_t(naxesI[0]) != nphi_)
    GYOTO_ERROR("XillverReflection::readFile(): phi array not conformable");
  if (phi_) { delete [] phi_; phi_ = NULL; }
  phi_ = new double[nphi_];
  if (fits_read_subset(fptrI, TDOUBLE, fpixelI, naxesI, incI, 
		       0, phi_,&anynulI,&statusI)) {
    delete [] phi_; phi_=NULL;
    throwCfitsioError(statusI) ;
  }
  
  ////// CLOSING FITS /////////
  
  if (fits_close_file(fptrI, &statusI)) throwCfitsioError(statusI) ;
  fptrI = NULL;
}

void XillverReflection::fitsReadRefl(string filenameRefl) {
  GYOTO_MSG << "XillverReflection reading FITS files: " <<
    filenameRefl << endl;

  filenameRefl_ = filenameRefl;
  char*     pixfileR   = const_cast<char*>(filenameRefl_.c_str());
  fitsfile* fptrR      = NULL;
  int       statusR    = 0;
  int       anynulR    = 0;
  long      naxesR []  = {1,1,1};
  long      fpixelR[]  = {1,1,1};
  long      incR   []  = {1,1,1};
  char      ermsg[31] = ""; // ermsg is used in throwCfitsioError()

  GYOTO_DEBUG << "XillverReflection::readFile(): opening refl file" << endl;
  if (fits_open_file(&fptrR, pixfileR, 0, &statusR)) throwCfitsioError(statusR) ;

  ////// FIND MANDATORY REFLECTION HDU, READ KWDS & DATA ///////
  GYOTO_DEBUG << "XillverReflection::readFile(): "
    "search reflection HDU" << endl;
  if (fits_movnam_hdu(fptrR, ANY_HDU,
		      const_cast<char*>("GYOTO XillverReflection reflection"),
		      0, &statusR))
    throwCfitsioError(statusR) ;
  GYOTO_DEBUG << "XillverReflection::readFile(): get image size" << endl;
  if (fits_get_img_size(fptrR, 3, naxesR, &statusR)) throwCfitsioError(statusR) ;

  //update nnu_, ni_, nxi_
  nnu_ = naxesR[0]; 
  ni_  = naxesR[1];
  nxi_  = naxesR[2];

  if (reflection_) { delete [] reflection_; reflection_ = NULL; }
  reflection_ = new double[nnu_ * ni_ * nxi_];
  if (debug())
    cerr << "XillverReflection::readFile(): read reflection: "
	 << "nnu_=" << nnu_ << ", ni_="<<ni_ << ", nxi_="<<nxi_ << "...";
  if (fits_read_subset(fptrR, TDOUBLE, fpixelR, naxesR, incR,
		       0, reflection_,&anynulR,&statusR)) {
    GYOTO_DEBUG << " error, trying to free pointer" << endl;
    delete [] reflection_; reflection_=NULL;
    throwCfitsioError(statusR) ;
  }
  GYOTO_DEBUG << " done." << endl;

  // double minemission=DBL_MAX, maxemission=DBL_MIN;
  // for (int myi=0;myi<nnu_ * ni_ * nsg_-1;myi++){
  //   if (emission_[myi]<minemission) minemission=emission_[myi];
  //   if (emission_[myi]>maxemission) maxemission=emission_[myi];
  // }
  //cout << "In NSModelAtm::fitsRead: Min and max emission= " <<
  //  minemission << " " << maxemission << endl;
  
  
  ////// FIND MANDATORY REFL::FREQ HDU ///////
  
  if (fits_movnam_hdu(fptrR, ANY_HDU,
		      const_cast<char*>("GYOTO XillverReflection freq"),
		      0, &statusR))
    throwCfitsioError(statusR) ;
  if (fits_get_img_size(fptrR, 1, naxesR, &statusR)) throwCfitsioError(statusR) ;
  if (size_t(naxesR[0]) != nnu_)
    GYOTO_ERROR("XillverReflection::readFile(): freq array not conformable");
  if (freq_) { delete [] freq_; freq_ = NULL; }
  freq_ = new double[nnu_];
  if (fits_read_subset(fptrR, TDOUBLE, fpixelR, naxesR, incR, 
		       0, freq_,&anynulR,&statusR)) {
    delete [] freq_; freq_=NULL;
    throwCfitsioError(statusR) ;
  }
  
  ////// FIND MANDATORY REFL::INCL HDU ///////
  
  if (fits_movnam_hdu(fptrR, ANY_HDU,
		      const_cast<char*>("GYOTO XillverReflection incl"),
		      0, &statusR))
    throwCfitsioError(statusR) ;
  if (fits_get_img_size(fptrR, 1, naxesR, &statusR)) throwCfitsioError(statusR) ;
  if (size_t(naxesR[0]) != ni_)
    GYOTO_ERROR("XillverReflection::readFile(): incl array not conformable");
  if (incl_) { delete [] incl_; incl_ = NULL; }
  incl_ = new double[ni_];
  if (fits_read_subset(fptrR, TDOUBLE, fpixelR, naxesR, incR, 
		       0, incl_,&anynulR,&statusR)) {
    delete [] incl_; incl_=NULL;
    throwCfitsioError(statusR) ;
  }
  
  ////// FIND MANDATORY REFL::LOGXI HDU ///////
  
  if (fits_movnam_hdu(fptrR, ANY_HDU,
		      const_cast<char*>("GYOTO XillverReflection logxi"),
		      0, &statusR))
    throwCfitsioError(statusR) ;
  if (fits_get_img_size(fptrR, 1, naxesR, &statusR)) throwCfitsioError(statusR) ;
  if (size_t(naxesR[0]) != nxi_)
    GYOTO_ERROR("XillverReflection::readFile(): logxi array not conformable");
  if (logxi_) { delete [] logxi_; logxi_ = NULL; }
  logxi_ = new double[nxi_];
  if (fits_read_subset(fptrR, TDOUBLE, fpixelR, naxesR, incR, 
		       0, logxi_,&anynulR,&statusR)) {
    delete [] logxi_; logxi_=NULL;
    throwCfitsioError(statusR) ;
  }
  
  ////// CLOSING FITS /////////
  
  if (fits_close_file(fptrR, &statusR)) throwCfitsioError(statusR) ;
  fptrR = NULL;
}

void XillverReflection::fitsWriteIllum(string filenameIllumination
				       ) {
  GYOTO_DEBUG_EXPR(illumination_);
  if (!illumination_) GYOTO_ERROR("XillverReflection::fitsWrite(filename): no illumination to save!");
  
  filenameIllum_ = filenameIllumination;
  char*     pixfileI   = const_cast<char*>(filenameIllum_.c_str());
  fitsfile* fptrI      = NULL;
  int       statusI    = 0;
  long      naxesI []  = {long(nr_), long(nphi_)};
  long      fpixelI[]  = {1,1};
  char *    CNULLI     = NULL;
  char      ermsg[31] = ""; // ermsg is used in throwCfitsioError()

  ////// CREATE FILE
  GYOTO_DEBUG << "creating illum file \"" << pixfileI << "\"... ";
  fits_create_file(&fptrI, pixfileI, &statusI);
  if (debug()) cerr << "done." << endl;
  fits_create_img(fptrI, DOUBLE_IMG, 2, naxesI, &statusI);
  if (statusI) throwCfitsioError(statusI) ;

  ////// SAVE ILLUMINATION IN PRIMARY HDU ///////
  GYOTO_DEBUG << "saving illumination_\n";
  fits_write_key(fptrI, TSTRING,
		 const_cast<char*>("EXTNAME"),
		 const_cast<char*>("GYOTO XillverReflection illumination"),
		 CNULLI, &statusI);
  fits_write_pix(fptrI, TDOUBLE, fpixelI, nr_*nphi_, illumination_,
		 &statusI);
  if (statusI) throwCfitsioError(statusI) ;

  ////// SAVE ILLUM::RADIUS HDU ///////
  if (!radius_) GYOTO_ERROR("XillverReflection::fitsWrite(filename): "
			   "no radius to save!");
  GYOTO_DEBUG << "saving radius_\n";
  fits_create_img(fptrI, DOUBLE_IMG, 1, naxesI, &statusI);
  fits_write_key(fptrI, TSTRING, const_cast<char*>("EXTNAME"),
		 const_cast<char*>("GYOTO XillverReflection radius"),
		 CNULLI, &statusI);
  fits_write_pix(fptrI, TDOUBLE, fpixelI, nr_, radius_, &statusI);
  if (statusI) throwCfitsioError(statusI) ;
  
  ////// SAVE ILLUM::PHI HDU ///////
  if (!phi_) GYOTO_ERROR("XillverReflection::fitsWrite(filename): "
			"no phi to save!");
  GYOTO_DEBUG << "saving phi_\n";
  fits_create_img(fptrI, DOUBLE_IMG, 1, naxesI+1, &statusI);
  fits_write_key(fptrI, TSTRING, const_cast<char*>("EXTNAME"),
		 const_cast<char*>("GYOTO XillverReflection phi"),
		 CNULLI, &statusI);
  fits_write_pix(fptrI, TDOUBLE, fpixelI, nphi_, phi_, &statusI);
  if (statusI) throwCfitsioError(statusI) ;
  
  ////// CLOSING FILE ///////
  GYOTO_DEBUG << "close FITS file\n";
  if (fits_close_file(fptrI, &statusI)) throwCfitsioError(statusI) ;
  fptrI = NULL;
}

void XillverReflection::fitsWriteRefl(
				      string filenameReflection) {
  GYOTO_DEBUG_EXPR(reflection_);
  
  if (!reflection_) GYOTO_ERROR("XillverReflection::fitsWrite(filename): no reflection to save!");
  
  filenameRefl_ = filenameReflection;
  char*     pixfileR   = const_cast<char*>(filenameRefl_.c_str());
  fitsfile* fptrR      = NULL;
  int       statusR    = 0;
  long      naxesR []  = {long(nnu_), long(ni_), long(nxi_)};
  long      fpixelR[]  = {1,1,1};
  char *    CNULLR     = NULL;
  char      ermsg[31] = ""; // ermsg is used in throwCfitsioError()

  ////// CREATE FILE
  GYOTO_DEBUG << "creating refl file \"" << pixfileR << "\"... ";
  fits_create_file(&fptrR, pixfileR, &statusR);
  if (debug()) cerr << "done." << endl;
  fits_create_img(fptrR, DOUBLE_IMG, 3, naxesR, &statusR);
  if (statusR) throwCfitsioError(statusR) ;

  ////// SAVE REFLECTION IN PRIMARY HDU ///////
  GYOTO_DEBUG << "saving reflection_\n";
  fits_write_key(fptrR, TSTRING,
		 const_cast<char*>("EXTNAME"),
		 const_cast<char*>("GYOTO XillverReflection reflection"),
		 CNULLR, &statusR);
  fits_write_pix(fptrR, TDOUBLE, fpixelR, nnu_*ni_*nxi_, reflection_,
		 &statusR);
  if (statusR) throwCfitsioError(statusR) ;

  ////// SAVE REFL::FREQ HDU ///////
  if (!freq_) GYOTO_ERROR("XillverReflection::fitsWrite(filename): "
			 "no freq to save!");
  GYOTO_DEBUG << "saving freq_\n";
  fits_create_img(fptrR, DOUBLE_IMG, 1, naxesR, &statusR);
  fits_write_key(fptrR, TSTRING, const_cast<char*>("EXTNAME"),
		 const_cast<char*>("GYOTO XillverReflection freq"),
		 CNULLR, &statusR);
  fits_write_pix(fptrR, TDOUBLE, fpixelR, nnu_, freq_, &statusR);
  if (statusR) throwCfitsioError(statusR) ;
  
  ////// SAVE REFL::INCL HDU ///////
  if (!incl_) GYOTO_ERROR("XillverReflection::fitsWrite(filename): "
			 "no incl to save!");
  GYOTO_DEBUG << "saving incl_\n";
  fits_create_img(fptrR, DOUBLE_IMG, 1, naxesR+1, &statusR);
  fits_write_key(fptrR, TSTRING, const_cast<char*>("EXTNAME"),
		 const_cast<char*>("GYOTO XillverReflection incl"),
		 CNULLR, &statusR);
  fits_write_pix(fptrR, TDOUBLE, fpixelR, ni_, incl_, &statusR);
  if (statusR) throwCfitsioError(statusR) ;
  
  ////// SAVE REFL::LOGXI HDU ///////
  if (!logxi_) GYOTO_ERROR("XillverReflection::fitsWrite(filename): "
			  "no logxi to save!");
  GYOTO_DEBUG << "saving logxi_\n";
  fits_create_img(fptrR, DOUBLE_IMG, 1, naxesR+2, &statusR);
  fits_write_key(fptrR, TSTRING, const_cast<char*>("EXTNAME"),
		 const_cast<char*>("GYOTO XillverReflection logxi"),
		 CNULLR, &statusR);
  fits_write_pix(fptrR, TDOUBLE, fpixelR, nxi_, logxi_, &statusR);
  if (statusR) throwCfitsioError(statusR) ;
  
  ////// CLOSING FILE ///////
  GYOTO_DEBUG << "close FITS file\n";
  if (fits_close_file(fptrR, &statusR)) throwCfitsioError(statusR) ;
  fptrR = NULL;
}
#endif

void XillverReflection::getIndicesRefl(size_t i[3], double const co[4], 
				       double logxi, double incl, double freq) const {
  if (logxi_) {
    // Here we know that logxi is within the bounds of logxi_
    for(i[2]=0; logxi > logxi_[i[2]]; ++i[2]){}
    /*
      With this definition:
      logxi_[i[2]-1] <= r < logxi_[i[2]]
      and i[2]>0
    */
  } else {
    GYOTO_ERROR("In XillverReflection::getIndicesRefl: logxi undefined!");
  }

  if (incl_) {
    if (incl >= incl_[ni_-1]) i[1] = ni_-1;
    else{
      for(i[1]=0; incl > incl_[i[1]]; ++i[1]){}
      /*
	With this definition:
	incl_[i[1]-1] <= r < incl_[i[1]]
	and i[1]=0 means incl<incl_[0] (allowed)
      */
    }
  } else {
    GYOTO_ERROR("In XillverReflection::getIndicesRefl: incl undefined!");
  }

  if (freq_) {
    // Here we know that freq is within the bounds of freq_
    for(i[0]=0; freq > freq_[i[0]]; ++i[0]){}
    /*
      With this definition:
      freq_[i[0]-1] <= r < freq_[i[0]]
      and i[0]>0
    */
  } else {
    GYOTO_ERROR("In XillverReflection::getIndicesRefl: freq undefined!");
  }
}

void XillverReflection::getIndicesIllum(size_t i[2], double const co[4])
  const {
  double rr = projectedRadius(co),
    phi = co[3];
  if (phi<=0. or phi>2*M_PI) GYOTO_ERROR("In Xillver::getIndicesIllum: "
				       "phi value not in 0,2pi");
  if (phi<phi_[0]) phi+=2*M_PI; // this is to conveniently deal
        // with phi values between phimax and phimin, modulo 2pi
  
  if (radius_) {
    // Here we know that rr is within the bounds of radius_
    for(i[0]=0; rr > radius_[i[0]]; ++i[0]){}
    /*
      With this definition:
      radius_[i[0]-1] <= r < radius_[i[0]]	
      The case i[0]=0 is impossible here.
    */
  } else {
    GYOTO_ERROR("In XillverReflection::getIndicesIllum: radius undefined!");
  }

  if (phi_) {
    if (phi >= phi_[nphi_-1]) {
      i[1] = nphi_-1;
    } else {
      for(i[1]=0; phi > phi_[i[1]]; ++i[1]){
      }
      /*
	phi_[i[1]-1] <= phi < phi_[i[1]]
	The redefinition of phi at the begining
	makes it impossible to have i[1]=0.
      */
    }
  } else {
    GYOTO_ERROR("In XillverReflection::getIndicesIllum: phi undefined!");
  }
  
}

double XillverReflection::emission(double nu, double,
				   state_t const &cp,
				   double const co[8]) const{
  // ************* ILLUMINATION *************

  double lampperiod = 2*M_PI*(pow(lampradius_,1.5)+aa_),
    timehit = cp[0];
  while (timehit<0.){
    // we only care about the hit time modulo the period time
    // and want a positive hit time
    timehit+=lampperiod;
  }
  double timerescale = fmod(timehit,lampperiod); // hit time modulo lamp period
  //cout << "times= " << lampperiod << " " << timehit << " " << timerescale << endl; 
  // No illumination and no reflection outside
  // the radius range
  double rr=co[1], phi=co[3];
  if (rr<=radius_[0] || rr>=radius_[nr_-1]) return 0.;
  double philamp = (timerescale-timelampphizero_)/lampperiod*2*M_PI; // lamp phi position at hit
  double dphi = phi - philamp; // phi difference between lamp and disk hit pos
  // this is the quantity with which to interpolate the illum data computed
  // only for philamp=0
  // At this stage we have -2pi < dphi < 4pi (checked)
  
  if (dphi<-2*M_PI or dphi>4*M_PI) {
    cout << "tresc, tlamp, period, phi, philamp, dphi= "
	 << timerescale << " " << timelampphizero_ << " " << lampperiod
	 << " " << phi << " " << philamp << " " << dphi << endl;
    GYOTO_ERROR("In Xillver::emission: "
	       "bad dphi");
  }
  while (dphi<0.){
    dphi+=2.*M_PI;
  }
  while (dphi>2*M_PI){
    dphi-=2.*M_PI;
  }
  // Here 0<dphi<2pi, ready for interpolation
  //cout << "r, phi, philamp, dphi= " << rr << " " << phi << " " << philamp << " " << dphi << endl;
  if (dphi<0. or dphi>2*M_PI) {
    GYOTO_ERROR("In Xillver::emission: bad dphi after correction");
  }
  // Illumination indices of the current closest grid point
  size_t indIllum[2]; // {i_r, i_phi}
  double col[4]={co[0],co[1],co[2],dphi};
  getIndicesIllum(indIllum, col);
  
  // **** Bilinear interpo for illumination
  size_t iru=indIllum[0];
  if (iru==0) GYOTO_ERROR("In Xillver::emission: bad radius index");
  // indeed, iru=0 means rr<radius_[0] which is impossible here
  size_t irl=iru-1;
  size_t iphiu=indIllum[1];
  if (iphiu==0) GYOTO_ERROR("In Xillver==emission: bad phi index");
  size_t iphil=iphiu-1;
  if (dphi>phi_[nphi_-1]){
    iphil = nphi_-1;
    iphiu = 0;
    // indeed, in this case, dphi is above phimax,
    // and below phimin+2pi. So interpolate between
    // the last and first values of phi_.
  }
  //cout << "r indices= " << irl << " " << iru << " " << radius_[irl] << " " << rr << " " << radius_[iru] << endl;
  //cout << "phi indices= " << iphil << " " << iphiu << " " << phi_[iphil] << " " << dphi << " " << phi_[iphiu] << endl;
  double F00 = illumination_[nphi_*irl+iphil], // F_{r,phi}
    F01 = illumination_[nphi_*irl+iphiu],
    F10 = illumination_[nphi_*iru+iphil],
    F11 = illumination_[nphi_*iru+iphiu];
  double ratior = (rr-radius_[irl])/(radius_[iru]-radius_[irl]),
    ratiophi = (dphi-phi_[iphil])/(phi_[iphiu]-phi_[iphil]);
  if (dphi>phi_[nphi_-1]) // add 2pi to phi_[iphiu] if dphi not within bounds
    ratiophi = (dphi-phi_[iphil])/(2*M_PI+phi_[iphiu]-phi_[iphil]);
  double fluxillum = F00
    +(F10-F00)*ratior
    +(F01-F00)*ratiophi
    +(F11-F01-F10+F00)*ratior*ratiophi;
  //cout << "Illum bilin= " << F00 << " " << F01 << " " << F10 << " " << F11 << " " << fluxillum << endl;

  // ************* REFLECTION ************* 

  // **** Compute ionization param
  double ne = 1e15; // electron number density in cm^-3 assumed constant
  double logxi = log10(4.*M_PI*fluxillum/ne); // log of ionization param
  // logxi should be within the bounds of logxi_
  //cout << "logxi= " << logxi << endl;
  if (logxi<logxi_[0] or logxi>logxi_[nxi_-1]) {
    cout << "Illum and logxi= " << fluxillum << " " << logxi << endl;
    GYOTO_ERROR("In Xillver::emission: logxi out of bounds");
  }

  // **** Emission angle
  double normal[4]={0.,0.,-1.,0.}; // parallel to -d_theta (upwards)
  double normal_norm=gg_->ScalarProd(&cp[0],normal,normal);
  if (normal_norm<=0.) GYOTO_ERROR("In XillverReflection::emission"
				  " normal should be spacelike");
  normal_norm=sqrt(normal_norm);
  double np = 1./normal_norm*gg_->ScalarProd(&cp[0],normal,&cp[4]),
    up = gg_->ScalarProd(&cp[0],co+4,&cp[4]);
  // cos between unit normal n and tangent to photon p
  // is equal -n.p/u.p (u being the emitter's 4-vel);
  // fabs because assuming plane symmetry
  double cosi = fabs(-np/up),
    incl = acos(cosi)*180./M_PI;
  //double tolcos = 0.005;
  //if (cosi>1.){
  //  if (fabs(cosi-1)>tolcos) GYOTO_ERROR("In XillverReflection: bad cos!");
  //  cosi=1.;
  //}

  // Frequency should not be outisde the range
  if (nu<=freq_[0] || nu>=freq_[nnu_-1]) {
    cout << "nu= " << nu << endl;
    GYOTO_ERROR("In Xillver::emission: freq outside range");
  }
					   
  // Reflection indices of the current closest grid point
  size_t ind[3]; // {i_nu, i_incl, i_xi}
  getIndicesRefl(ind, co, logxi, incl, nu);
  size_t inuu = ind[0];
  if (inuu==0) GYOTO_ERROR("In Xillver::emission: bad nu index");
  size_t inul = inuu-1;
  size_t iiu = ind[1];
  size_t ixiu = ind[2];
  //cout << "illum, logxi: " << ixiu << " " << logxi << " " << logxi_[ixiu] << endl;
  if (ixiu==0) GYOTO_ERROR("In Xillver::emission: bad logxi index");
  size_t ixil = ixiu-1;
  /*cout << "logxi: " << ixiu << " " << logxi_[ixil] << " " << logxi << " " << logxi_[ixiu] << endl;
  cout << "nu: " << inuu << " " << freq_[inul] << " " << nu << " " << freq_[inuu] << endl;
  if (iiu>0){
    cout << "incl: " << iiu << " " << incl_[iiu-1] << " " << incl << " " << incl_[iiu] << endl;
    }*/

  double reflectedintens=0.;
  if (incl<incl_[0] or incl>incl_[ni_-1]){
    // Bilinear interpo if incl is not within bounds
    size_t ii=iiu; // unique index for incl, no interpo
    // bilin interpo in nu, logxi:
    double I00 = reflection_[inul*ni_*nxi_+ii*nxi_+ixil], // I_{nu,xi}
      I01 = reflection_[inul*ni_*nxi_+ii*nxi_+ixiu],
      I10 = reflection_[inuu*ni_*nxi_+ii*nxi_+ixil],
      I11 = reflection_[inuu*ni_*nxi_+ii*nxi_+ixiu];
    double rationu = (nu-freq_[inul])/(freq_[inuu]-freq_[inul]),
      ratioxi = (logxi-logxi_[ixil])/(logxi_[ixiu]-logxi_[ixil]);
    reflectedintens = I00
      +(I10-I00)*rationu
      +(I01-I00)*ratioxi
      +(I11-I01-I10+I00)*rationu*ratioxi;
    //cout << "Refl bilin: " << I00 << " " << I01 << " " << I10 << " " << I11 << " " << reflectedintens << endl;
  }else{
    // Trilinear interpo general case
    if (iiu==0) GYOTO_ERROR("In Xillver::emission: bad incl index");
    size_t iil = iiu-1;

    double I000 = reflection_[inul*ni_*nxi_+iil*nxi_+ixil], // I_{nu,incl,xi}
      I100      = reflection_[inuu*ni_*nxi_+iil*nxi_+ixil],
      I110      = reflection_[inuu*ni_*nxi_+iiu*nxi_+ixil], 
      I010      = reflection_[inul*ni_*nxi_+iiu*nxi_+ixil],
      I001      = reflection_[inul*ni_*nxi_+iil*nxi_+ixiu], 
      I101      = reflection_[inuu*ni_*nxi_+iil*nxi_+ixiu],
      I111      = reflection_[inuu*ni_*nxi_+iiu*nxi_+ixiu],
      I011      = reflection_[inul*ni_*nxi_+iiu*nxi_+ixiu];
    double rationu = (nu-freq_[inul])/(freq_[inuu]-freq_[inul]),
      ratioi = (incl-incl_[iil])/(incl_[iiu]-incl_[iil]),
      ratioxi = (logxi-logxi_[ixil])/(logxi_[ixiu]-logxi_[ixil]);
    reflectedintens = I000
	+ (I100-I000)*rationu
	+ (I010-I000)*ratioi
	+ (I001-I000)*ratioxi
	+ (I110-I010-I100+I000)*rationu*ratioi
	+ (I011-I010-I001+I000)*ratioi*ratioxi
	+ (I101-I001-I100+I000)*rationu*ratioxi
	+ (I111-I011-I101-I110+I100+I001+I010-I000)*rationu*ratioi*ratioxi;
    //cout << "Refl trilin: " << I000 << " " << I100 << " " << I110 << " " << I010 << " " << I001 << " " << I101 << " " << I111 << " " << I011 << " " << reflectedintens << endl;
  }
  
  return reflectedintens/1e3; // 1e3 factor translates from cgs to SI,
                              // gyoto speaks in SI
								 
}
