/*
    Copyright 2014, 2018 Frederic Vincent, Thibaut Paumard

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
#include "GyotoDirectionalDisk.h"
#include "GyotoUtils.h"
#include "GyotoFactoryMessenger.h"
#include "GyotoProperty.h"
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

//// Properties:

GYOTO_PROPERTY_START(DirectionalDisk)
GYOTO_PROPERTY_FILENAME(DirectionalDisk, File, file)
GYOTO_PROPERTY_DOUBLE(DirectionalDisk, LampAltitude, lampaltitude)
GYOTO_PROPERTY_VECTOR_DOUBLE(DirectionalDisk, LampCutOffsIneV, lampcutoffsinev)
GYOTO_PROPERTY_BOOL(DirectionalDisk,
		    AverageOverAngle, DontAverageOverAngle,
		    averageOverAngle)
GYOTO_PROPERTY_END(DirectionalDisk, ThinDisk::properties)

void DirectionalDisk::fillProperty(Gyoto::FactoryMessenger *fmp,
			       Property const &p) const {
  if (p.name == "File")
    fmp->setParameter("File", (filename_.compare(0,1,"!") ?
			       filename_ :
			       filename_.substr(1)) );
  else ThinDisk::fillProperty(fmp, p);
}

////

DirectionalDisk::DirectionalDisk() :
  ThinDisk("DirectionalDisk"), filename_(""),
  emission_(NULL), radius_(NULL), cosi_(NULL), freq_(NULL),
  nnu_(0), ni_(0), nr_(0),
  lampaltitude_(10.),
  minfreq_computed_(DBL_MIN), maxfreq_computed_(DBL_MAX),
  minfreq_lampframe_(1.), maxfreq_lampframe_(1e30),
  average_over_angle_(0)
{
  GYOTO_DEBUG << "DirectionalDisk Construction" << endl;
}

DirectionalDisk::DirectionalDisk(const DirectionalDisk& o) :
  ThinDisk(o), filename_(o.filename_),
  emission_(NULL), radius_(NULL), cosi_(NULL), freq_(NULL),
  nnu_(o.nnu_), ni_(o.ni_), nr_(o.nr_),
  lampaltitude_(o.lampaltitude_),
  minfreq_computed_(o.minfreq_computed_), 
  maxfreq_computed_(o.maxfreq_computed_),
  minfreq_lampframe_(o.minfreq_lampframe_), 
  maxfreq_lampframe_(o.maxfreq_lampframe_), 
  average_over_angle_(o.average_over_angle_)
{
  GYOTO_DEBUG << "DirectionalDisk Copy" << endl;
  size_t ncells = 0;
  if (o.emission_) {
    emission_ = new double[ncells = nnu_ * ni_ * nr_];
    memcpy(emission_, o.emission_, ncells * sizeof(double));
  }
  if (o.freq_) {
    freq_ = new double[ncells = nnu_];
    memcpy(freq_, o.freq_, ncells * sizeof(double));
  }
  if (o.cosi_) {
    cosi_ = new double[ncells = ni_];
    memcpy(cosi_, o.cosi_, ncells * sizeof(double));
  }
  if (o.radius_) {
    radius_ = new double[ncells = nr_];
    memcpy(radius_, o.radius_, ncells * sizeof(double));
  }
}
DirectionalDisk* DirectionalDisk::clone() const
{ return new DirectionalDisk(*this); }

DirectionalDisk::~DirectionalDisk() {
  GYOTO_DEBUG << "DirectionalDisk Destruction" << endl;
  if (emission_) delete [] emission_;
  if (radius_) delete [] radius_;
  if (cosi_) delete [] cosi_;
  if (freq_) delete [] freq_;

}

void DirectionalDisk::setEmission(double * pattern) {
  emission_ = pattern;
}

void DirectionalDisk::radius(double * pattern) {
  radius_ = pattern;
}

void DirectionalDisk::copyIntensity(double const *const pattern, size_t const naxes[3]) {
  GYOTO_DEBUG << endl;
  if (emission_) {
    GYOTO_DEBUG << "delete [] emission_;" << endl;
    delete [] emission_; emission_ = NULL;
  }
  if (pattern) {
    size_t nel;
    if (nnu_ != naxes[0]) {
      GYOTO_DEBUG <<"nnu_ changed, freeing freq_" << endl;
      if (freq_)  { delete [] freq_; freq_  = NULL; }
    }
    if (ni_ != naxes[1]) {
      GYOTO_DEBUG <<"ni_ changed, freeing freq_ and cosi_" << endl;
      if (freq_)  { delete [] freq_; freq_  = NULL; }
      if (cosi_) { delete [] cosi_; cosi_= NULL; }
    }
    if (nr_ != naxes[2]) {
      GYOTO_DEBUG <<"nr_ changed, freeing freq_, cosi_ and radius_" << endl;
      if (freq_)  { delete [] freq_; freq_  = NULL; }
      if (cosi_) { delete [] cosi_; cosi_= NULL; }
      if (radius_)   { delete [] radius_;   radius_  = NULL; }
    }
    if (!(nel=(nnu_ = naxes[0]) * (ni_=naxes[1]) * (nr_=naxes[2])))
      GYOTO_ERROR( "dimensions can't be null");
    GYOTO_DEBUG << "allocate emission_;" << endl;
    emission_ = new double[nel];
    GYOTO_DEBUG << "pattern >> emission_" << endl;
    memcpy(emission_, pattern, nel*sizeof(double));
  }
}

double const * DirectionalDisk::getIntensity() const { return emission_; }
void DirectionalDisk::getIntensityNaxes( size_t naxes[3] ) const
{ naxes[0] = nnu_; naxes[1] = ni_; naxes[2] = nr_; }

void DirectionalDisk::copyGridRadius(double const *const rad, size_t nr) {
  GYOTO_DEBUG << endl;
  if (radius_) {
    GYOTO_DEBUG << "delete [] radius_;" << endl;
    delete [] radius_; radius_ = NULL;
  }
  if (rad) {
    if (!emission_) 
      GYOTO_ERROR("Please use copyIntensity() before copyGridRadius()");
    if (nr_ != nr)
      GYOTO_ERROR("emission_ and radius_ have inconsistent dimensions");
    GYOTO_DEBUG << "allocate radius_;" << endl;
    radius_ = new double[nr_];
    GYOTO_DEBUG << "radius >> radius_" << endl;
    memcpy(radius_, rad, nr_*sizeof(double));
  }
}
double const * DirectionalDisk::getGridRadius() const { return radius_; }

void DirectionalDisk::copyGridCosi(double const *const cosi, size_t ni) {
  GYOTO_DEBUG << endl;
  if (cosi_) {
    GYOTO_DEBUG << "delete [] cosi_;" << endl;
    delete [] cosi_; cosi_ = NULL;
  }
  if (cosi) {
    if (!emission_) 
      GYOTO_ERROR("Please use copyIntensity() before copyGridCosi()");
    if (ni_ != ni)
      GYOTO_ERROR("emission_ and cosi_ have inconsistent dimensions");
    GYOTO_DEBUG << "allocate cosi_;" << endl;
    cosi_ = new double[ni_];
    GYOTO_DEBUG << "cosi >> cosi_" << endl;
    memcpy(cosi_, cosi, ni_*sizeof(double));
  }
}
double const * DirectionalDisk::getGridCosi() const { return cosi_; }

void DirectionalDisk::copyGridFreq(double const *const freq, size_t nnu) {
  GYOTO_DEBUG << endl;
  if (freq_) {
    GYOTO_DEBUG << "delete [] freq_;" << endl;
    delete [] freq_; freq_ = NULL;
  }
  if (freq) {
    if (!emission_) 
      GYOTO_ERROR("Please use copyIntensity() before copyGridFreq()");
    if (nnu_ != nnu)
      GYOTO_ERROR("emission_ and freq_ have inconsistent dimensions");
    GYOTO_DEBUG << "allocate freq_;" << endl;
    freq_ = new double[nnu_];
    GYOTO_DEBUG << "freq >> freq_" << endl;
    memcpy(freq_, freq, nnu_*sizeof(double));
  }
}
double const * DirectionalDisk::getGridFreq() const { return freq_; }

void DirectionalDisk::averageOverAngle(bool t) {average_over_angle_=t;}
bool DirectionalDisk::averageOverAngle()const {return average_over_angle_;}

void DirectionalDisk::file(std::string const &f) {
# ifdef GYOTO_USE_CFITSIO
  fitsRead(f);
# else
  GYOTO_ERROR("This Gyoto has no FITS i/o");
# endif
}

std::string DirectionalDisk::file() const {
  return filename_;
}

void DirectionalDisk::lampaltitude(double zz) {
  lampaltitude_ = zz;
}

double DirectionalDisk::lampaltitude() const {
  return lampaltitude_;
}

void DirectionalDisk::lampcutoffsinev(std::vector<double> const &v) {
  if (v.size() != 2)
    GYOTO_ERROR("In DirectionalDisk: Only 2 arguments to define lamp energy range");
  minfreq_lampframe_ = v[0]*GYOTO_eV2Hz;
  maxfreq_lampframe_ = v[1]*GYOTO_eV2Hz;
}

std::vector<double> DirectionalDisk::lampcutoffsinev() const {
  std::vector<double> v (2, 0.);
  v[0]=minfreq_lampframe_; v[1]=maxfreq_lampframe_;
  return v;
}

#ifdef GYOTO_USE_CFITSIO
void DirectionalDisk::fitsRead(string filename) {
  GYOTO_MSG << "DirectionalDisk reading FITS file: " << filename << endl;

  filename_ = filename;
  char*     pixfile   = const_cast<char*>(filename_.c_str());
  fitsfile* fptr      = NULL;
  int       status    = 0;
  int       anynul    = 0;
  long      naxes []  = {1, 1, 1};
  long      fpixel[]  = {1,1,1};
  long      inc   []  = {1,1,1};
  char      ermsg[31] = ""; // ermsg is used in throwCfitsioError()

  GYOTO_DEBUG << "DirectionalDisk::readFile(): opening file" << endl;
  if (fits_open_file(&fptr, pixfile, 0, &status)) throwCfitsioError(status) ;

  ////// FIND MANDATORY EMISSION HDU, READ KWDS & DATA ///////
  GYOTO_DEBUG << "DirectionalDisk::readFile(): search emission HDU" << endl;
  if (fits_movnam_hdu(fptr, ANY_HDU,
		      const_cast<char*>("GYOTO DirectionalDisk emission"),
		      0, &status))
    throwCfitsioError(status) ;
  GYOTO_DEBUG << "DirectionalDisk::readFile(): get image size" << endl;
  if (fits_get_img_size(fptr, 3, naxes, &status)) throwCfitsioError(status) ;

  //update nnu_, ni_, nr_
  nnu_ = naxes[0]; 
  ni_  = naxes[1];
  nr_  = naxes[2];

  if (emission_) { delete [] emission_; emission_ = NULL; }
  emission_ = new double[nnu_ * ni_ * nr_];
  if (debug())
    cerr << "DirectionalDisk::readFile(): read emission: "
	 << "nnu_=" << nnu_ << ", ni_="<<ni_ << ", nr_="<<nr_ << "...";
  if (fits_read_subset(fptr, TDOUBLE, fpixel, naxes, inc,
		       0, emission_,&anynul,&status)) {
    GYOTO_DEBUG << " error, trying to free pointer" << endl;
    delete [] emission_; emission_=NULL;
    throwCfitsioError(status) ;
  }
  GYOTO_DEBUG << " done." << endl;

  ////// FIND MANDATORY FREQ HDU ///////
  
   if (fits_movnam_hdu(fptr, ANY_HDU,
		       const_cast<char*>("GYOTO DirectionalDisk freq"),
		       0, &status))
     throwCfitsioError(status) ;
   if (fits_get_img_size(fptr, 1, naxes, &status)) throwCfitsioError(status) ;
   if (size_t(naxes[0]) != nnu_)
     GYOTO_ERROR("DirectionalDisk::readFile(): freq array not conformable");
   if (freq_) { delete [] freq_; freq_ = NULL; }
   freq_ = new double[nnu_];
   if (fits_read_subset(fptr, TDOUBLE, fpixel, naxes, inc, 
			0, freq_,&anynul,&status)) {
     delete [] freq_; freq_=NULL;
     throwCfitsioError(status) ;
   }

   // Computing min and max of freq_

   minfreq_computed_ = DBL_MAX;
   maxfreq_computed_ = DBL_MIN;
   for (int ii=0;ii<nnu_;ii++){
     if (freq_[ii]<minfreq_computed_) minfreq_computed_=freq_[ii];
     if (freq_[ii]>maxfreq_computed_) maxfreq_computed_=freq_[ii];
   }
   
   GYOTO_DEBUG << "Min, max freq= " << minfreq_computed_ << " " << maxfreq_computed_ << endl;

  ////// FIND MANDATORY COSI HDU ///////
  
   if (fits_movnam_hdu(fptr, ANY_HDU,
		       const_cast<char*>("GYOTO DirectionalDisk cosi"),
		       0, &status))
     throwCfitsioError(status) ;
   if (fits_get_img_size(fptr, 1, naxes, &status)) throwCfitsioError(status) ;
   if (size_t(naxes[0]) != ni_)
     GYOTO_ERROR("DirectionalDisk::readFile(): cosi array not conformable");
   if (cosi_) { delete [] cosi_; cosi_ = NULL; }
   cosi_ = new double[ni_];
   if (fits_read_subset(fptr, TDOUBLE, fpixel, naxes, inc, 
			0, cosi_,&anynul,&status)) {
     delete [] cosi_; cosi_=NULL;
     throwCfitsioError(status) ;
   }

  ////// FIND MANDATORY RADIUS HDU ///////
  
   if (fits_movnam_hdu(fptr, ANY_HDU,
		       const_cast<char*>("GYOTO DirectionalDisk radius"),
		       0, &status))
     throwCfitsioError(status) ;
   if (fits_get_img_size(fptr, 1, naxes, &status)) throwCfitsioError(status) ;
   if (size_t(naxes[0]) != nr_)
     GYOTO_ERROR("DirectionalDisk::readFile(): radius array not conformable");
   if (radius_) { delete [] radius_; radius_ = NULL; }
   radius_ = new double[nr_];
   if (fits_read_subset(fptr, TDOUBLE, fpixel, naxes, inc, 
			0, radius_,&anynul,&status)) {
     delete [] radius_; radius_=NULL;
     throwCfitsioError(status) ;
   }

   ////// CLOSING FITS /////////

  if (fits_close_file(fptr, &status)) throwCfitsioError(status) ;
  fptr = NULL;
}

void DirectionalDisk::fitsWrite(string filename) {
  if (!emission_) GYOTO_ERROR("DirectionalDisk::fitsWrite(filename): nothing to save!");
  filename_ = filename;
  char*     pixfile   = const_cast<char*>(filename_.c_str());
  fitsfile* fptr      = NULL;
  int       status    = 0;
  long      naxes []  = {long(nnu_), long(ni_), long(nr_)};
  long      fpixel[]  = {1,1,1};
  char * CNULL=NULL;

  char      ermsg[31] = ""; // ermsg is used in throwCfitsioError()

  ////// CREATE FILE
  GYOTO_DEBUG << "creating file \"" << pixfile << "\"... ";
  fits_create_file(&fptr, pixfile, &status);
  if (debug()) cerr << "done." << endl;
  fits_create_img(fptr, DOUBLE_IMG, 3, naxes, &status);
  if (status) throwCfitsioError(status) ;

  ////// SAVE EMISSION IN PRIMARY HDU ///////
  GYOTO_DEBUG << "saving emission_\n";
  fits_write_key(fptr, TSTRING,
		 const_cast<char*>("EXTNAME"),
		 const_cast<char*>("GYOTO DirectionalDisk emission"),
		 CNULL, &status);
  fits_write_pix(fptr, TDOUBLE, fpixel, nnu_*ni_*nr_, emission_, &status);
  if (status) throwCfitsioError(status) ;

  ////// SAVE FREQ HDU ///////
  if (!freq_) GYOTO_ERROR("DirectionalDisk::fitsWrite(filename): no freq to save!");
  GYOTO_DEBUG << "saving freq_\n";
  fits_create_img(fptr, DOUBLE_IMG, 1, naxes, &status);
  fits_write_key(fptr, TSTRING, const_cast<char*>("EXTNAME"),
		 const_cast<char*>("GYOTO DirectionalDisk freq"),
		 CNULL, &status);
  fits_write_pix(fptr, TDOUBLE, fpixel, nnu_, freq_, &status);
  if (status) throwCfitsioError(status) ;
  
  ////// SAVE COSI HDU ///////
  if (!cosi_) GYOTO_ERROR("DirectionalDisk::fitsWrite(filename): no cosi to save!");
  GYOTO_DEBUG << "saving cosi_\n";
  fits_create_img(fptr, DOUBLE_IMG, 1, naxes+1, &status);
  fits_write_key(fptr, TSTRING, const_cast<char*>("EXTNAME"),
		 const_cast<char*>("GYOTO DirectionalDisk cosi"),
		 CNULL, &status);
  fits_write_pix(fptr, TDOUBLE, fpixel, ni_, cosi_, &status);
  if (status) throwCfitsioError(status) ;
  
  ////// SAVE RADIUS HDU ///////
  if (!radius_) GYOTO_ERROR("DirectionalDisk::fitsWrite(filename): no radius to save!");
    GYOTO_DEBUG << "saving radius_\n";
    fits_create_img(fptr, DOUBLE_IMG, 1, naxes+2, &status);
    fits_write_key(fptr, TSTRING, const_cast<char*>("EXTNAME"),
		   const_cast<char*>("GYOTO DirectionalDisk radius"),
		   CNULL, &status);
    fits_write_pix(fptr, TDOUBLE, fpixel, nr_, radius_, &status);
    if (status) throwCfitsioError(status) ;

  ////// CLOSING FILE ///////
  GYOTO_DEBUG << "close FITS file\n";
  if (fits_close_file(fptr, &status)) throwCfitsioError(status) ;
  fptr = NULL;
}
#endif

void DirectionalDisk::getIndices(size_t i[3], double const co[4], 
				 double cosi, double nu) const {
  double rr = projectedRadius(co);
  if (radius_) {
    if (rr >= radius_[nr_-1]) i[2] = nr_-1; // emission will be 0
    else {
      for(i[2]=0; rr > radius_[i[2]]; ++i[2]){}
      /*
	With this definition:
	radius_[i[2]-1] <= r < radius_[i[2]]
	
	The case i[2]=0 (if r<radius_[0]) is dealt
	with later on, it returns 0
      */
    }
  } else {
    GYOTO_ERROR("In DirectionalDisk::getIndices: radius undefined!");
  }

  if (cosi_) {
    if (cosi >= cosi_[ni_-1]) i[1] = ni_-1;
    else {
      for(i[1]=0; cosi > cosi_[i[1]]; ++i[1]){}
      /*
	cosi_[i[1]-1] <= cosi < cosi_[i[1]]
      */
    }
  } else {
    GYOTO_ERROR("In DirectionalDisk::getIndices: cosi undefined!");
  }

  if (freq_) {
    if (nu <= freq_[nnu_-1]) i[0] = nnu_-1;
    else { 
      for(i[0]=nnu_-1; nu > freq_[i[0]]; --i[0]){}
      /*
	Caution: freq is ordered decreasingly!
	freq_[i[0]+1] <= nu < freq_[i[0]]
      */
    }
  } else {
    GYOTO_ERROR("In DirectionalDisk::getIndices: freq undefined!");
  }

}

double DirectionalDisk::emission(double nu, double,
				 state_t const &cp,
				 double const co[8]) const{
  GYOTO_DEBUG << endl;
  // Checking whether the current freq is outside of 
  // the redshifted illumination range
  double aa = static_cast<SmartPointer<Metric::KerrBL> >(gg_) -> spin();
  double zz = lampaltitude_;
  double rr = co[1];
  double gg_lampdisk = (pow(rr,1.5)+aa)/sqrt(rr*rr*rr+2*aa*pow(rr,1.5)-3*rr*rr)
    *sqrt((zz*zz+aa*aa-2.*zz)/(zz*zz+aa*aa)); // this is the redshift factor
  // linking the lamp-frame and the disk-frame, it has nothing to do with the
  // redshift factor between the disk and the far-away-observer frames.
  double minfreq_diskframe = minfreq_lampframe_*gg_lampdisk,
    maxfreq_diskframe = maxfreq_lampframe_*gg_lampdisk;
  //cout << "Limits computed= " << minfreq_computed_ << " " << maxfreq_computed_<< endl;
  //cout << "Limits lamp= " << minfreq_lampframe_ << " " << maxfreq_lampframe_<< endl;
  //cout << "Limits disk= " << minfreq_diskframe << " " << maxfreq_diskframe << endl;
  //cout << "Local nu= " << nu << endl;
  if (minfreq_diskframe < minfreq_computed_ || maxfreq_diskframe > maxfreq_computed_){
    GYOTO_ERROR("In DirectionalDisk::emission(): "
	       "bad freq value ; update LampCutOffsIneV in XML");
  }

  // Cut-offs in disk frame: if the local freq is not inside the
  // redshifted illumination band, no signal
  if (nu < minfreq_diskframe || nu > maxfreq_diskframe) return 0.;

  // Compute angle between photon direction and normal
  double normal[4]={0.,0.,-1.,0.}; // parallel to -d_theta (upwards)
  double normal_norm=gg_->ScalarProd(&cp[0],normal,normal);
  if (normal_norm<=0.) GYOTO_ERROR("In DirectionalDisk::emission"
				  " normal should be spacelike");
  normal_norm=sqrt(normal_norm);
  double np = 1./normal_norm*gg_->ScalarProd(&cp[0],normal,&cp[4]),
    up = gg_->ScalarProd(&cp[0],co+4,&cp[4]);
  double cosi = fabs(-np/up);
  double tolcos = 0.005;
  GYOTO_DEBUG_EXPR(cosi);
  GYOTO_DEBUG_ARRAY(co, 8);
  GYOTO_DEBUG_ARRAY(cp, 8);
  if (cosi>1.){
    if (fabs(cosi-1)>tolcos) GYOTO_ERROR("In DirectionalDisk: bad cos!");
    cosi=1.;
  }
  //cout << "cosi= " << cosi << endl;
  // Don't put a "return cosi" here, see later

  // cos between unit normal n and tangent to photon p
  // is equal -n.p/u.p (u being the emitter's 4-vel);
  // fabs because assuming plane symmetry

  // Indices of the current closest grid point
  size_t ind[3]; // {i_nu, i_cosi, i_r}
  getIndices(ind, co, cosi, nu);

  //cout << "r, i2, nr= " << co[1] << " " << ind[2] << " " << nr_ <<endl;

  //if (ind[2]==nr_) return 0.; // 0 emission outside simulation scope

  // Specific intensity emitted at the current location
  // No emission outside radius and frequency data range
  if (rr<=radius_[0] || rr>=radius_[nr_-1]) return 0.;
  if (nu<=freq_[nnu_-1] || nu>=freq_[0]) return 0.;
  // So here, ind[2] should be >0 and ind[0]<nnu_-1
  if (ind[2]==0 || ind[0]==nnu_-1){
    GYOTO_ERROR("In DirectionalDisk::emission "
	       "bad {nu,r} indices");
  }

  //return acos(cosi)*180./M_PI; // TEST!!! Don't forget to impose redshift to 1

  //cout << "nu(eV), r(rS), cosi= " << nu/GYOTO_eV2Hz << " " << rr/2. << " " << cosi << endl;
  double Iem=0.;
  size_t i0l=ind[0]+1, i0u=ind[0], 
    i2l=ind[2]-1, i2u=ind[2]; // Correct: i0 is freq, ordered decreasingly,
                              // i2 is radius ordered increasingly

  //  cout << "ind_cosi=, ni= " << ind[1] << " " << ni_ << endl;
  //cout << "min max r= " << radius_[0] << " " << radius_[nr_-1] << endl;

  if (!average_over_angle_){
    if (cosi <= cosi_[0] || cosi >= cosi_[ni_-1]){
      // If cosi is out of the cosi_ range, bilinear interpol in nu,r
      size_t i1=ind[1];
      //cout << "cos value unique= " << cosi_[i1] << endl;
      double I00 = emission_[i2l*(ni_*nnu_)+i1*nnu_+i0l], // I_{nu,r}
	I01 = emission_[i2u*(ni_*nnu_)+i1*nnu_+i0l],
	I10 = emission_[i2l*(ni_*nnu_)+i1*nnu_+i0u],
	I11 = emission_[i2u*(ni_*nnu_)+i1*nnu_+i0u];
      //cout << "bilin dir: " << I00 << " " << I01 << " " << I10 << " " << I11 << endl;
      double rationu = (nu-freq_[i0l])/(freq_[i0u]-freq_[i0l]),
	ratior = (rr-radius_[i2l])/(radius_[i2u]-radius_[i2l]);
      Iem = I00+(I10-I00)*rationu
	+(I01-I00)*ratior
	+(I11-I01-I10+I00)*rationu*ratior;
      //cout << "I interp= " << Iem << endl;
    }else{
      // Trilinear interpol
      if (ind[1]==0){
	GYOTO_ERROR("In DirectionalDisk::emission "
		   "bad cosi indice");
      }
      size_t i1l=ind[1]-1, i1u=ind[1];
      double I000 = emission_[i2l*(ni_*nnu_)+i1l*nnu_+i0l], // I_{nu,cosi,r}
	I100 = emission_[i2l*(ni_*nnu_)+i1l*nnu_+i0u],
	I110 = emission_[i2l*(ni_*nnu_)+i1u*nnu_+i0u], 
	I010 = emission_[i2l*(ni_*nnu_)+i1u*nnu_+i0l],
	I001 = emission_[i2u*(ni_*nnu_)+i1l*nnu_+i0l], 
	I101 = emission_[i2u*(ni_*nnu_)+i1l*nnu_+i0u],
	I111 = emission_[i2u*(ni_*nnu_)+i1u*nnu_+i0u],
	I011 = emission_[i2u*(ni_*nnu_)+i1u*nnu_+i0l];
      //cout << "trilin dir: " << I000 << " " << I100 << " " << I110 << " " << I010 << " " << I001 << " " << I101 << " " << I111 << " " << I011 << endl;
      double rationu = (nu-freq_[i0l])/(freq_[i0u]-freq_[i0l]),
	ratioi = (cosi-cosi_[i1l])/(cosi_[i1u]-cosi_[i1l]),
	ratior = (rr-radius_[i2l])/(radius_[i2u]-radius_[i2l]);
      Iem = I000
	+ (I100-I000)*rationu
	+ (I010-I000)*ratioi
	+ (I001-I000)*ratior
	+ (I110-I010-I100+I000)*rationu*ratioi
	+ (I011-I010-I001+I000)*ratioi*ratior
	+ (I101-I001-I100+I000)*rationu*ratior
	+ (I111-I011-I101-I110+I100+I001+I010-I000)*rationu*ratioi*ratior;
      //cout << "I interp= " << Iem << endl;
    }
  }else{
    // Average over cosi values
    // with bilinear interpol in nu,r
    double I00=0., I01=0., I10=0., I11=0.;
    double I00min=DBL_MAX, I00max=DBL_MIN, I01min=DBL_MAX, I01max=DBL_MIN, I10min=DBL_MAX, I10max=DBL_MIN, I11min=DBL_MAX, I11max=DBL_MIN;
    /* Using trapezoidal rule, I_integ = \int I(mu)*dmu, mu=cos(i)
       NB: in Garcia+14, they compute a flux because they don't raytrace,
       so they use F = 1/4pi * \int I(i) cos(i) di = 1/2 * \int I(mu) mu dmu,
       here we are not interested in the same quantity */
    double dcostot = 0.; // will contain \int d\mu (~1 but not exactly)
    for (size_t ii=0; ii<ni_-1; ++ii){
      double dcos = cosi_[ii+1]-cosi_[ii];
      I00 += 0.5*dcos*
	(emission_[i2l*(ni_*nnu_)+(ii+1)*nnu_+i0l]
	 +emission_[i2l*(ni_*nnu_)+ii*nnu_+i0l]);
      I01 += 0.5*dcos*
	(emission_[i2u*(ni_*nnu_)+(ii+1)*nnu_+i0l]
	 +emission_[i2u*(ni_*nnu_)+ii*nnu_+i0l]);
      I10 += 0.5*dcos*
	(emission_[i2l*(ni_*nnu_)+(ii+1)*nnu_+i0u]
	 +emission_[i2l*(ni_*nnu_)+ii*nnu_+i0u]);
      I11 += 0.5*dcos*
	(emission_[i2u*(ni_*nnu_)+(ii+1)*nnu_+i0u]
	 +emission_[i2u*(ni_*nnu_)+ii*nnu_+i0u]);
      dcostot+=dcos;

      /*
	// CHECK MINIMUM
	if (emission_[i2l*(ni_*nnu_)+(ii+1)*nnu_+i0l]<I00min) I00min=emission_[i2l*(ni_*nnu_)+(ii+1)*nnu_+i0l];
      if (emission_[i2l*(ni_*nnu_)+(ii)*nnu_+i0l]<I00min) I00min=emission_[i2l*(ni_*nnu_)+(ii)*nnu_+i0l];
      if (emission_[i2u*(ni_*nnu_)+(ii+1)*nnu_+i0l]<I01min) I01min=emission_[i2u*(ni_*nnu_)+(ii+1)*nnu_+i0l];
      if (emission_[i2u*(ni_*nnu_)+(ii)*nnu_+i0l]<I01min) I01min=emission_[i2u*(ni_*nnu_)+(ii)*nnu_+i0l];
      if (emission_[i2l*(ni_*nnu_)+(ii+1)*nnu_+i0u]<I10min) I10min=emission_[i2l*(ni_*nnu_)+(ii+1)*nnu_+i0u];
      if (emission_[i2l*(ni_*nnu_)+(ii)*nnu_+i0u]<I10min) I10min=emission_[i2l*(ni_*nnu_)+(ii)*nnu_+i0u];
      if (emission_[i2u*(ni_*nnu_)+(ii+1)*nnu_+i0u]<I11min) I11min=emission_[i2u*(ni_*nnu_)+(ii+1)*nnu_+i0u];
      if (emission_[i2u*(ni_*nnu_)+(ii)*nnu_+i0u]<I11min) I11min=emission_[i2u*(ni_*nnu_)+(ii)*nnu_+i0u];
      */

      /*
      // CHECK MAXIMUM
      if (emission_[i2l*(ni_*nnu_)+(ii+1)*nnu_+i0l]>I00max) I00max=emission_[i2l*(ni_*nnu_)+(ii+1)*nnu_+i0l];
      if (emission_[i2l*(ni_*nnu_)+(ii)*nnu_+i0l]>I00max) I00max=emission_[i2l*(ni_*nnu_)+(ii)*nnu_+i0l];
      if (emission_[i2u*(ni_*nnu_)+(ii+1)*nnu_+i0l]>I01max) I01max=emission_[i2u*(ni_*nnu_)+(ii+1)*nnu_+i0l];
      if (emission_[i2u*(ni_*nnu_)+(ii)*nnu_+i0l]>I01max) I01max=emission_[i2u*(ni_*nnu_)+(ii)*nnu_+i0l];
      if (emission_[i2l*(ni_*nnu_)+(ii+1)*nnu_+i0u]>I10max) I10max=emission_[i2l*(ni_*nnu_)+(ii+1)*nnu_+i0u];
      if (emission_[i2l*(ni_*nnu_)+(ii)*nnu_+i0u]>I10max) I10max=emission_[i2l*(ni_*nnu_)+(ii)*nnu_+i0u];
      if (emission_[i2u*(ni_*nnu_)+(ii+1)*nnu_+i0u]>I11max) I11max=emission_[i2u*(ni_*nnu_)+(ii+1)*nnu_+i0u];
      if (emission_[i2u*(ni_*nnu_)+(ii)*nnu_+i0u]>I11max) I11max=emission_[i2u*(ni_*nnu_)+(ii)*nnu_+i0u];
      */
      
      
      //cout << "Raw data 1 for I00= " << radius_[i2l]/2. << " " << radius_[i2u]/2. << " " << freq_[i0l]/GYOTO_eV2Hz << " " << freq_[i0u]/GYOTO_eV2Hz << " " << cosi_[ii] << " " << cosi_[ii+1] << endl;
      //cout << "Raw data 2 for I00= " << emission_[i2l*(ni_*nnu_)+(ii+1)*nnu_+i0l] << " " << emission_[i2l*(ni_*nnu_)+ii*nnu_+i0l] << endl;

      //cout << "IO in avg i= " << ii << " and I0= " << I00 << endl;

      
    } 

    // Normalizing (int d co(i) is very close to 1 but not exactly 1)
    I00/=dcostot;
    I01/=dcostot;
    I10/=dcostot;
    I11/=dcostot;

    //cout << "bilin avg: " << I00 << " " << I01 << " " << I10 << " " << I11 << endl;

    //if (I00<I00min || I01<I01min || I10<I10min || I11<I11min) GYOTO_ERROR("test");
    //if (I00>I00max || I01>I01max || I10>I10max || I11>I11max) GYOTO_ERROR("test");
    double rationu = (nu-freq_[i0l])/(freq_[i0u]-freq_[i0l]),
      ratior = (rr-radius_[i2l])/(radius_[i2u]-radius_[i2l]);
    Iem = I00+(I10-I00)*rationu
      +(I01-I00)*ratior
      +(I11-I01-I10+I00)*rationu*ratior;
    //cout << "I interp= " << Iem << endl;
  }
  //cout << "return= " << Iem << endl;
  return Iem;
}

void DirectionalDisk::metric(SmartPointer<Metric::Generic> gg) {
  //Metric must be KerrBL (see emission function)
  string kin = gg->kind();
  if (kin != "KerrBL")
    GYOTO_ERROR("DirectionalDisk::metric(): metric must be KerrBL");
  ThinDisk::metric(gg);
}
