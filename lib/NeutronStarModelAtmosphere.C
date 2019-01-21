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

// Lorene headers
#include "metric.h"
#include "nbr_spx.h"
#include "utilitaires.h"
#include "graphique.h"

//Gyoto headers
#include "GyotoUtils.h"
#include "GyotoPhoton.h"
#include "GyotoNeutronStarModelAtmosphere.h"
#include "GyotoFactoryMessenger.h"


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

#define LORENE_UNIT_ACCEL GYOTO_C*GYOTO_C/1e4

using namespace std;
using namespace Gyoto;
using namespace Gyoto::Astrobj;
using namespace Lorene;

/// Properties
#include "GyotoProperty.h"
GYOTO_PROPERTY_START(NeutronStarModelAtmosphere,
		     "Neutron star emitting at its surface.")
GYOTO_PROPERTY_FILENAME(NeutronStarModelAtmosphere, File, file)
GYOTO_PROPERTY_BOOL(NeutronStarModelAtmosphere,
		    AverageOverAngle, DontAverageOverAngle,
		    averageOverAngle)
GYOTO_PROPERTY_END(NeutronStarModelAtmosphere, NeutronStar::properties)

void NeutronStarModelAtmosphere::fillProperty(Gyoto::FactoryMessenger *fmp,
			       Property const &p) const {
  if (p.name == "File")
    fmp->setParameter("File", (filename_.compare(0,1,"!") ?
			       filename_ :
			       filename_.substr(1)) );
  else NeutronStar::fillProperty(fmp, p);
}

NeutronStarModelAtmosphere::NeutronStarModelAtmosphere() :
NeutronStar("NeutronStarModelAtmosphere"),
  emission_(NULL), surfgrav_(NULL), cosi_(NULL), freq_(NULL),
  nnu_(0), ni_(0), nsg_(0), average_over_angle_(0) {
  
  GYOTO_DEBUG << endl;

}

NeutronStarModelAtmosphere::NeutronStarModelAtmosphere(const NeutronStarModelAtmosphere& o) :
  NeutronStar(o),
  emission_(NULL), surfgrav_(NULL), cosi_(NULL), freq_(NULL),
  nnu_(o.nnu_), ni_(o.ni_), nsg_(o.nsg_),
  average_over_angle_(o.average_over_angle_)
{
  GYOTO_DEBUG << endl;
  size_t ncells = 0;
  if (o.emission_) {
    emission_ = new double[ncells = nnu_ * ni_ * nsg_];
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
  if (o.surfgrav_) {
    surfgrav_ = new double[ncells = nsg_];
    memcpy(surfgrav_, o.surfgrav_, ncells * sizeof(double));
  }

}
NeutronStarModelAtmosphere * NeutronStarModelAtmosphere::clone() const {
  return new NeutronStarModelAtmosphere(*this); }

NeutronStarModelAtmosphere::~NeutronStarModelAtmosphere() {
  GYOTO_DEBUG << endl;
  if (emission_) delete [] emission_;
  if (surfgrav_) delete [] surfgrav_;
  if (cosi_) delete [] cosi_;
  if (freq_) delete [] freq_;
}

void NeutronStarModelAtmosphere::setEmission(double * pattern) {
  emission_ = pattern;
}

void NeutronStarModelAtmosphere::surfgrav(double * pattern) {
  surfgrav_ = pattern;
}

void NeutronStarModelAtmosphere::copyIntensity(double const *const pattern, size_t const naxes[3]) {
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
    if (nsg_ != naxes[2]) {
      GYOTO_DEBUG <<"nsg_ changed, freeing freq_, cosi_ and surfgrav_" << endl;
      if (freq_)  { delete [] freq_; freq_  = NULL; }
      if (cosi_) { delete [] cosi_; cosi_= NULL; }
      if (surfgrav_)   { delete [] surfgrav_;   surfgrav_  = NULL; }
    }
    if (!(nel=(nnu_ = naxes[0]) * (ni_=naxes[1]) * (nsg_=naxes[2])))
      GYOTO_ERROR( "dimensions can't be null");
    GYOTO_DEBUG << "allocate emission_;" << endl;
    emission_ = new double[nel];
    GYOTO_DEBUG << "pattern >> emission_" << endl;
    memcpy(emission_, pattern, nel*sizeof(double));
  }
}

double const * NeutronStarModelAtmosphere::getIntensity() const {
  return emission_; }
void NeutronStarModelAtmosphere::getIntensityNaxes( size_t naxes[3] ) const
{ naxes[0] = nnu_; naxes[1] = ni_; naxes[2] = nsg_; }

void NeutronStarModelAtmosphere::copyGridSurfgrav(double const *const sg,
						size_t nsg) {
  GYOTO_DEBUG << endl;
  if (surfgrav_) {
    GYOTO_DEBUG << "delete [] surfgrav_;" << endl;
    delete [] surfgrav_; surfgrav_ = NULL;
  }
  if (sg) {
    if (!emission_) 
      GYOTO_ERROR("Please use copyIntensity() before copyGridSurfgrav()");
    if (nsg_ != nsg)
      GYOTO_ERROR("emission_ and surfgrav_ have inconsistent dimensions");
    GYOTO_DEBUG << "allocate surfgrav_;" << endl;
    surfgrav_ = new double[nsg_];
    GYOTO_DEBUG << "surfgrav >> surfgrav_" << endl;
    memcpy(surfgrav_, sg, nsg_*sizeof(double));
  }
}
double const * NeutronStarModelAtmosphere::getGridSurfgrav() const {
  return surfgrav_; }

void NeutronStarModelAtmosphere::copyGridCosi(double const *const cosi, size_t ni) {
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
double const * NeutronStarModelAtmosphere::getGridCosi() const { return cosi_; }

void NeutronStarModelAtmosphere::copyGridFreq(double const *const freq,
					      size_t nnu) {
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
double const * NeutronStarModelAtmosphere::getGridFreq() const { return freq_; }

void NeutronStarModelAtmosphere::averageOverAngle(bool t) {
  average_over_angle_=t;}
bool NeutronStarModelAtmosphere::averageOverAngle()const {
  return average_over_angle_;}

void NeutronStarModelAtmosphere::file(std::string const &f) {
# ifdef GYOTO_USE_CFITSIO
  fitsRead(f);
# else
  GYOTO_ERROR("This Gyoto has no FITS i/o");
# endif
}

std::string NeutronStarModelAtmosphere::file() const {
  return filename_;
}

#ifdef GYOTO_USE_CFITSIO
void NeutronStarModelAtmosphere::fitsRead(string filename) {
  GYOTO_MSG << "NeutronStarModelAtmosphere reading FITS file: " <<
    filename << endl;

  filename_ = filename;
  char*     pixfile   = const_cast<char*>(filename_.c_str());
  fitsfile* fptr      = NULL;
  int       status    = 0;
  int       anynul    = 0;
  long      naxes []  = {1, 1, 1};
  long      fpixel[]  = {1,1,1};
  long      inc   []  = {1,1,1};
  char      ermsg[31] = ""; // ermsg is used in throwCfitsioError()

  GYOTO_DEBUG << "NeutronStarModelAtmosphere::readFile(): opening file" << endl;
  if (fits_open_file(&fptr, pixfile, 0, &status)) throwCfitsioError(status) ;

  ////// FIND MANDATORY EMISSION HDU, READ KWDS & DATA ///////
  GYOTO_DEBUG << "NeutronStarModelAtmosphere::readFile(): search emission HDU" << endl;
  if (fits_movnam_hdu(fptr, ANY_HDU,
		      const_cast<char*>("GYOTO NeutronStarModelAtmosphere emission"),
		      0, &status))
    throwCfitsioError(status) ;
  GYOTO_DEBUG << "NeutronStarModelAtmosphere::readFile(): get image size" << endl;
  if (fits_get_img_size(fptr, 3, naxes, &status)) throwCfitsioError(status) ;

  //update nnu_, ni_, nsg_
  nnu_ = naxes[0]; 
  ni_  = naxes[1];
  nsg_  = naxes[2];

  if (emission_) { delete [] emission_; emission_ = NULL; }
  emission_ = new double[nnu_ * ni_ * nsg_];
  if (debug())
    cerr << "NeutronStarModelAtmosphere::readFile(): read emission: "
	 << "nnu_=" << nnu_ << ", ni_="<<ni_ << ", nsg_="<<nsg_ << "...";
  if (fits_read_subset(fptr, TDOUBLE, fpixel, naxes, inc,
		       0, emission_,&anynul,&status)) {
    GYOTO_DEBUG << " error, trying to free pointer" << endl;
    delete [] emission_; emission_=NULL;
    throwCfitsioError(status) ;
  }
  GYOTO_DEBUG << " done." << endl;

  double minemission=DBL_MAX, maxemission=DBL_MIN;
  for (int myi=0;myi<nnu_ * ni_ * nsg_-1;myi++){
    if (emission_[myi]<minemission) minemission=emission_[myi];
    if (emission_[myi]>maxemission) maxemission=emission_[myi];
  }
  //cout << "In NSModelAtm::fitsRead: Min and max emission= " <<
  //  minemission << " " << maxemission << endl;

  ////// FIND MANDATORY FREQ HDU ///////
  
   if (fits_movnam_hdu(fptr, ANY_HDU,
		       const_cast<char*>("GYOTO NeutronStarModelAtmosphere freq"),
		       0, &status))
     throwCfitsioError(status) ;
   if (fits_get_img_size(fptr, 1, naxes, &status)) throwCfitsioError(status) ;
   if (size_t(naxes[0]) != nnu_)
     GYOTO_ERROR("NeutronStarModelAtmosphere::readFile(): freq array not conformable");
   if (freq_) { delete [] freq_; freq_ = NULL; }
   freq_ = new double[nnu_];
   if (fits_read_subset(fptr, TDOUBLE, fpixel, naxes, inc, 
			0, freq_,&anynul,&status)) {
     delete [] freq_; freq_=NULL;
     throwCfitsioError(status) ;
   }

  ////// FIND MANDATORY COSI HDU ///////
  
   if (fits_movnam_hdu(fptr, ANY_HDU,
		       const_cast<char*>("GYOTO NeutronStarModelAtmosphere cosi"),
		       0, &status))
     throwCfitsioError(status) ;
   if (fits_get_img_size(fptr, 1, naxes, &status)) throwCfitsioError(status) ;
   if (size_t(naxes[0]) != ni_)
     GYOTO_ERROR("NeutronStarModelAtmosphere::readFile(): cosi array not conformable");
   if (cosi_) { delete [] cosi_; cosi_ = NULL; }
   cosi_ = new double[ni_];
   if (fits_read_subset(fptr, TDOUBLE, fpixel, naxes, inc, 
			0, cosi_,&anynul,&status)) {
     delete [] cosi_; cosi_=NULL;
     throwCfitsioError(status) ;
   }

  ////// FIND MANDATORY SURFGRAV HDU ///////
  
   if (fits_movnam_hdu(fptr, ANY_HDU,
		       const_cast<char*>("GYOTO NeutronStarModelAtmosphere surfgrav"),
		       0, &status))
     throwCfitsioError(status) ;
   if (fits_get_img_size(fptr, 1, naxes, &status)) throwCfitsioError(status) ;
   if (size_t(naxes[0]) != nsg_)
     GYOTO_ERROR("NeutronStarModelAtmosphere::readFile(): surfgrav array not conformable");
   if (surfgrav_) { delete [] surfgrav_; surfgrav_ = NULL; }
   surfgrav_ = new double[nsg_];
   if (fits_read_subset(fptr, TDOUBLE, fpixel, naxes, inc, 
			0, surfgrav_,&anynul,&status)) {
     delete [] surfgrav_; surfgrav_=NULL;
     throwCfitsioError(status) ;
   }

   ////// CLOSING FITS /////////

  if (fits_close_file(fptr, &status)) throwCfitsioError(status) ;
  fptr = NULL;
}

void NeutronStarModelAtmosphere::fitsWrite(string filename) {
  GYOTO_DEBUG_EXPR(emission_);
  if (!emission_) GYOTO_ERROR("NeutronStarModelAtmosphere::fitsWrite(filename): nothing to save!");
  filename_ = filename;
  char*     pixfile   = const_cast<char*>(filename_.c_str());
  fitsfile* fptr      = NULL;
  int       status    = 0;
  long      naxes []  = {long(nnu_), long(ni_), long(nsg_)};
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
		 const_cast<char*>("GYOTO NeutronStarModelAtmosphere emission"),
		 CNULL, &status);
  fits_write_pix(fptr, TDOUBLE, fpixel, nnu_*ni_*nsg_, emission_, &status);
  if (status) throwCfitsioError(status) ;

  ////// SAVE FREQ HDU ///////
  if (!freq_) GYOTO_ERROR("NeutronStarModelAtmosphere::fitsWrite(filename): no freq to save!");
  GYOTO_DEBUG << "saving freq_\n";
  fits_create_img(fptr, DOUBLE_IMG, 1, naxes, &status);
  fits_write_key(fptr, TSTRING, const_cast<char*>("EXTNAME"),
		 const_cast<char*>("GYOTO NeutronStarModelAtmosphere freq"),
		 CNULL, &status);
  fits_write_pix(fptr, TDOUBLE, fpixel, nnu_, freq_, &status);
  if (status) throwCfitsioError(status) ;
  
  ////// SAVE COSI HDU ///////
  if (!cosi_) GYOTO_ERROR("NeutronStarModelAtmosphere::fitsWrite(filename): no cosi to save!");
  GYOTO_DEBUG << "saving cosi_\n";
  fits_create_img(fptr, DOUBLE_IMG, 1, naxes+1, &status);
  fits_write_key(fptr, TSTRING, const_cast<char*>("EXTNAME"),
		 const_cast<char*>("GYOTO NeutronStarModelAtmosphere cosi"),
		 CNULL, &status);
  fits_write_pix(fptr, TDOUBLE, fpixel, ni_, cosi_, &status);
  if (status) throwCfitsioError(status) ;
  
  ////// SAVE SURFGRAV HDU ///////
  if (!surfgrav_) GYOTO_ERROR("NeutronStarModelAtmosphere::fitsWrite(filename): no surfgrav to save!");
    GYOTO_DEBUG << "saving surfgrav_\n";
    fits_create_img(fptr, DOUBLE_IMG, 1, naxes+2, &status);
    fits_write_key(fptr, TSTRING, const_cast<char*>("EXTNAME"),
		   const_cast<char*>("GYOTO NeutronStarModelAtmosphere surfgrav"),
		   CNULL, &status);
    fits_write_pix(fptr, TDOUBLE, fpixel, nsg_, surfgrav_, &status);
    if (status) throwCfitsioError(status) ;

  ////// CLOSING FILE ///////
  GYOTO_DEBUG << "close FITS file\n";
  if (fits_close_file(fptr, &status)) throwCfitsioError(status) ;
  fptr = NULL;
}
#endif

void NeutronStarModelAtmosphere::getIndices(size_t i[3], double const co[4], 
				 double cosi, double nu) const {
  const Vector& a_i = *(gg_->getAccel_tab()[0]);
  double rr=co[1], th=co[2], phi=co[3];
  if (rr==0.) GYOTO_ERROR("In NeutronStarModelAtm.C::getIndices r is 0!");
  double rsinth = rr*sin(th);
  if (rsinth==0.) GYOTO_ERROR("In NeutronStarModelAtm.C::getIndices on z axis!");
  double rm1 = 1./rr, rm2 = rm1*rm1, sm1 = 1./sin(th),
    sm2 = sm1*sm1;
  double a_r = a_i(1).val_point(rr,th,phi),
    a_t = rr*a_i(2).val_point(rr,th,phi),
    a_p = rr*sin(th)*a_i(3).val_point(rr,th,phi);
  if (a_p!=0.) {GYOTO_ERROR("In NeutronStarModelAtm::getIndices: "
			   "For axisym spacetime phi-compo should be zero");}
  const Sym_tensor& g_up_ij = *(gg_->getGamcon_tab()[0]);
  double grr=g_up_ij(1,1).val_point(rr,th,phi), 
    gtt=rm2*g_up_ij(2,2).val_point(rr,th,phi);
  double ar = a_r*grr, at = a_t*gtt; //contravariant 3-accel

  double accelvecNorm2 = grr*a_r*a_r + gtt*a_t*a_t; // squared norm of accel vector
  if (accelvecNorm2<=0.) GYOTO_ERROR("In NeutronStarModelAtmosphere::getIndices"
				    " accel vector should be spacelike");
  double accelvecNorm = sqrt(accelvecNorm2);

  double sgloc = accelvecNorm*LORENE_UNIT_ACCEL*100.; // LORENE speaks in SI, the 100 translates to cgs

  //cout << "Accel vec compo, gs= " << a_r << " " << a_t << " " << a_p << " " << sgloc << endl;
  if (surfgrav_) { 
    if (nsg_==1){ // Only one value of surfgrav, put some value, won't be used
      i[2]=1; // don't put 0, see later why: it would return an error
    }else{
      if (sgloc >= surfgrav_[nsg_-1]) i[2] = nsg_-1; // emission will be 0
      else {
	for(i[2]=0; sgloc > surfgrav_[i[2]]; ++i[2]){}
	//cout << "In indices sg: " << i[2] << " " << surfgrav_[i[2]-1] << " " << sgloc << " " << surfgrav_[i[2]] << endl;
	/*
	  With this definition:
	  surfgrav_[i[2]-1] <= sgloc < surfgrav_[i[2]]
	  
	  The case i[2]=0 (if r<surfgrav_[0]) is dealt
	  with later on, it returns 0
	*/
      }
    }
  } else {
    GYOTO_ERROR("In NeutronStarModelAtmosphere::getIndices: surfgrav undefined!");
  }

  if (cosi_) {
    if (cosi >= cosi_[ni_-1]) i[1] = ni_-1;
    else {
      for(i[1]=0; cosi > cosi_[i[1]]; ++i[1]){}
      //cout << "In indices cos: " << i[1] << " " << cosi_[i[1]-1] << " " << cosi << " " << cosi_[i[1]] << endl;
      /*
	cosi_[i[1]-1] <= cosi < cosi_[i[1]]
      */
    }
  } else {
    GYOTO_ERROR("In NeutronStarModelAtmosphere::getIndices: cosi undefined!");
  }

  if (freq_) {
    if (nu <= freq_[nnu_-1]) i[0] = nnu_-1;
    else {
      for(i[0]=nnu_-1; nu > freq_[i[0]]; --i[0]){}
      //cout << "In indices nu: " << i[0] << " " << freq_[i[0]+1]/GYOTO_eV2Hz << " " << nu/GYOTO_eV2Hz << " " << freq_[i[0]]/GYOTO_eV2Hz << endl;
      /*
	Caution: freq is ordered decreasingly!
	freq_[i[0]+1] <= nu < freq_[i[0]]
      */
    }
  } else {
    GYOTO_ERROR("In NeutronStarModelAtmosphere::getIndices: freq undefined!");
  }

}

double NeutronStarModelAtmosphere::emission(double nu, double,
					    state_t const &cp,
					    double const co[8]) const{
  /*
    Important remarks on the precision: the variable GYOTO_T_TOL
    defined in GyotoDefs.h is important as it tunes the precision
    with which Gyoto will fine the star's surface. GYOTO_T_TOL=1e-4
    e.g. leads to an error on r_emission of approx 1e-4 as well.
    This for instance leads to small changes of Iobs when varying
    the inclination for a non-rotating star (although the emission
    should be indep of i). There is another limitation specifically
    for photons that hit the star tangentially. These guys should
    be integrated with care to find a precise (r_emission,theta_emission),
    and thus get a precise photon tangent vector at emission, that in
    turn manages the precision of cosi. To ensure this, a good
    precaution is to decrease DeltaMaxOverR in the xml, eg to 0.1
    to get something very precise. [Tested on August 2017 for a 30*30
    map: actually some few pixels even need 0.01 to get the same
    value for i=1° and i=90°! This would be really crazy long for
    a full map...]

    There is a second important limitation: the precision of Lorene.
    The metrics given by Michal in May 2017 e.g. lead to a value of
    r_star constant with theta to within approx 1e-8. 

    So to have the most precise calcualtion (to machine prec), take
    GYOTO_T_TOL=machine prec, make sure that Lorene works at
    machine prec as well, and use a small DeltaMAxOverR. Of course
    this is only needed to get crazy high precision (definitely not
    a problem to fit observations e.g.)
   */

  GYOTO_DEBUG << endl;
  //cout << "In emission NSatm, intens test= " << emission_[0] << endl;
  const Vector& a_i = *(gg_->getAccel_tab()[0]);
  double rr=co[1], th=co[2], phi=co[3];
  //cout << "r,th,phi in emiss= " << setprecision(10) << rr << " " << th << " " << phi << endl;

  // First, check that we are on the star surface (not obvious, photon
  // could be a bit inside, see StandardAstrobj.C). If not, return 0.
  // This is important coz if not present, sgloc can be computed inside
  // the star and be out of the range computed in the grid, leading to error.
  Valeur* ns_surf = gg_->getNssurf_tab()[0];
  ns_surf->std_base_scal();
  double rstar = ns_surf->val_point(0,0.,th,phi);
  //cout << "rstar= " << rstar << endl;
  double rtol = 1e-4; // should be such that GYOTO_T_TOL ensures a
                      // a convergence to rstar better than rtol
  if (fabs(rstar-rr)>rtol) return 0.;
  
  if (rr==0.) GYOTO_ERROR("In NeutronStarModelAtm.C::emission r is 0!");
  double rsinth = rr*sin(th);
  if (rsinth==0.) GYOTO_ERROR("In NeutronStarModelAtm.C::emission on z axis!");
  double rm1 = 1./rr, rm2 = rm1*rm1, sm1 = 1./sin(th),
    sm2 = sm1*sm1;

  // Finding acceleration 4-vector
  /*
    Let us call a^\alpha the acceleration 4-vector (normal to the surface).
    Because we assume circular fluid motion, the covariant
    acceleration reads
    a_\alpha = - \nabla (ln u_t) = - \partial_\alpha (ln u_t)
    Thus, because of stationarity + axisymetry: a_t = a_phi = 0
    Thus a^t = g^{tt} a_t + g^{tphi} a_phi = 0 and the same for a^phi.
    Thus only a^r and a^theta remain.
   */
  
  double a_r = a_i(1).val_point(rr,th,phi),
    a_t = rr*a_i(2).val_point(rr,th,phi),
    a_p = rr*sin(th)*a_i(3).val_point(rr,th,phi);
  if (a_p!=0.) {GYOTO_ERROR("In NeutronStarModelAtm::emission: "
			   "For axisym spacetime phi-compo should be zero");}
  const Sym_tensor& g_up_ij = *(gg_->getGamcon_tab()[0]);
  double grr=g_up_ij(1,1).val_point(rr,th,phi), 
    gtt=rm2*g_up_ij(2,2).val_point(rr,th,phi);
    //gpp=rm2*sm2*g_up_ij(3,3).val_point(rr,th,phi); // here gpp is gamma^{phi,phi} ; it is useless as a_p is zero
  double ar = a_r*grr, at = a_t*gtt; //contravariant 3-accel compo

  double accelvec[4]={0.,ar,at,0.}; // acceleration 4-vector, normal to surf

  //cout << "accel= " << a_r  << " " << a_t << " " << a_p << endl;

  double accelvecNorm2 = grr*a_r*a_r + gtt*a_t*a_t; // squared norm of accel vector
  if (accelvecNorm2<=0.) GYOTO_ERROR("In NeutronStarModelAtmosphere::emission"
				    " accel vector should be spacelike");
  double accelvecNorm = sqrt(accelvecNorm2);
  // Surface gravity is that quantity, scaled to cgs units:
  double sgloc = accelvecNorm*LORENE_UNIT_ACCEL*100.; // LORENE speaks in SI, the 100 translates to cgs

  //cout << "r,rstar,sgloc=" << rr << " " << rstar << " " << sgloc << endl;
  
  //cout << "r, sg, nsg= " << rr << " " << sgloc << " " << nsg_ << endl; 
  
  //cout << "accel vector= " << ar << " " << at << endl;
  //cout << "photon vector= " << cp[4] << " " << cp[5] << " " << cp[6] << " " << cp[7] << endl;
    
  // Compute angle between photon direction and normal
  double np = 1./accelvecNorm*gg_->ScalarProd(&cp[0],accelvec,&cp[4]),
    up = gg_->ScalarProd(&cp[0],co+4,&cp[4]);
  //cout << "accel and velo= " << accelvec[0] << " " << accelvec[1] << " " << accelvec[2] << " " << accelvec[3] << " " << cp[4] << " " << cp[5] << " " << cp[6] << " " << cp[7] << " " << endl;
  //cout << "scalar prods= " << np << " " << up << endl;
  double p_r=cp[5]/grr, p_t=cp[6]/gtt;
  double myscalprod = 1./accelvecNorm*(grr*a_r*p_r + gtt*a_t*p_t);
  //cout << "myscalprod= " << myscalprod << endl;
  //cout << "gmunu= " << grr << " " << gtt << endl;
  //cout << "prods of ap= " << a_r*p_r << " " << a_t*p_t << endl;
  //cout << "parts of scal prod= " << grr*a_r*p_r << " " << gtt*a_t*p_t << endl;
  //cout << "ar, pr, arpr= " << a_r << " " << p_r << " " << a_r*p_r << endl;

  // cos between unit normal n and tangent to photon p
  // is equal -n.p/u.p (u being the emitter's 4-vel);
  double cosi = -np/up;
  //cout << "cosi= " << cosi << endl;
  double tolcos = 0.005;
  if (cosi>1.){
    if (fabs(cosi-1)>tolcos){
      cout << "Bad cosi= " << cosi << endl;
      GYOTO_ERROR("In NeutronStarModelAtmosphere: bad cos!");
    }
    cosi=1.;
  }
  if (cosi<0.){
    // cosi should be >0, the photon cannot come from
    // inside the star!
    if (fabs(cosi)>tolcos){
      cout << "Bad cosi= " << cosi << endl;
      //GYOTO_ERROR("In NeutronStarModelAtmosphere: bad cos!");
    }
    cosi=0.;
  }
  //cout << "cosi= " << cosi << endl;
  // Don't put a "return cosi" here, see later

  // Indices of the current closest grid point
  size_t ind[3]; // {i_nu, i_cosi, i_surfgrav}
  getIndices(ind, co, cosi, nu);

  //cout << "sg, isg, nsg= " << sgloc << " " << ind[2] << " " << nsg_ <<endl;

  //if (ind[2]==nsg_) return 0.; // 0 emission outside simulation scope

  // Error if current surfgrav is not in provided range
  if (nsg_>1 && (sgloc<=surfgrav_[0] || sgloc>=surfgrav_[nsg_-1])){
    cout << "With surf grav= " << sgloc << endl;
    GYOTO_ERROR("In NeutronStarModelAtmosphere: bad value of surface gravity");
  }
  // No emission outside freq range
  if (nu<=freq_[nnu_-1] || nu>=freq_[0]) {
    //cout << "bad freq nu= " << nu << " " << freq_[0] << " " << freq_[nnu_-1] << endl;
    return 0.;
  }
  
  // So here, ind[2] should be >0 and ind[0]<nnu_-1
  if (ind[2]==0 || ind[0]==nnu_-1){
    GYOTO_ERROR("In NeutronStarModelAtmosphere::emission "
	       "bad {nu,r} indices");
  }

  //return acos(cosi)*180./M_PI; // TEST!!! Don't forget to impose redshift to 1

  //cout << setprecision(10) << "nu(eV), surfgrav, cosi= " << nu/GYOTO_eV2Hz << " " << sgloc << " " << cosi << endl;
  //cout << "indices= " << ind[0] << " " << ind[1] << " " << ind[2] << endl;
  double Iem=0.;
  size_t inul=ind[0]+1, inuu=ind[0], 
    isgl=ind[2]-1, isgu=ind[2]; // Correct: inu is freq, ordered decreasingly,
                              // isg is surfgrav ordered increasingly

  if (nsg_==1){
    // Only one value of surfgrav, i.e. non-rotating star
    // put surfgrav indices to zero, no interpolation in this direction
    isgl=0;
    isgu=0;
  }

  //  cout << "ind_cosi=, ni= " << ind[1] << " " << ni_ << endl;
  //cout << "min max sg= " << surfgrav_[0] << " " << surfgrav_[nsg_-1] << endl;

  /* 
     How emission_ is organized:

     [
     (nu=0,cos=0,sg=0),(nu=0,cos=0,sg=1),...,(nu=0,cos=0,sg=nsg-1),
     (nu=0,cos=1,sg=0),(nu=0,cos=1,sg=1),...,(nu=0,cos=1,sg=nsg-1),
     ...
     (nu=0,cos=ni-1,sg=0),(nu=0,cos=ni-1,sg=1),...,(nu=0,cos=ni-1,sg=nsg-1),
     (nu=1,cos=0,sg=0),(nu=1,cos=0,sg=1),...,(nu=1,cos=0,sg=nsg-1),
     ...
     ]

   */
  if (!average_over_angle_){
    if (cosi <= cosi_[0] || cosi >= cosi_[ni_-1]){
      // If cosi is out of the cosi_ range, bilinear interpol in nu,sg
      size_t icos=ind[1];
      //cout << "Bilin cos value unique= " << cosi_[icos] << endl;
      double I00 = emission_[inul*ni_*nsg_+icos*nsg_+isgl], // I_{nu,sg}
	I01 = emission_[inul*ni_*nsg_+icos*nsg_+isgu],
	I10 = emission_[inuu*ni_*nsg_+icos*nsg_+isgl],
	I11 = emission_[inuu*ni_*nsg_+icos*nsg_+isgu];
      //cout << "bilin dir: " << I00 << " " << I01 << " " << I10 << " " << I11 << endl;
      double rationu = (nu-freq_[inul])/(freq_[inuu]-freq_[inul]),
	ratiosg = (sgloc-surfgrav_[isgl])/(surfgrav_[isgu]-surfgrav_[isgl]);
      if (nsg_==1) ratiosg=0.; // no interpolation in sg
      Iem = I00+(I10-I00)*rationu
	+(I01-I00)*ratiosg
	+(I11-I01-I10+I00)*rationu*ratiosg;
      //cout << "I interp= " << Iem << endl;
    }else{
      // Trilinear interpol
      if (ind[1]==0){
	GYOTO_ERROR("In NeutronStarModelAtmosphere::emission "
		   "bad cosi indice");
      }
      size_t icosl=ind[1]-1, icosu=ind[1];
      //cout << "Trilin interpol indices= " << inul << " " << inuu << " " << icosl << " " << icosu << " " << isgl << " " << isgu << endl;
      double I000 = emission_[inul*ni_*nsg_+icosl*nsg_+isgl], // I_{nu,cosi,sg}
	I100 = emission_[inuu*ni_*nsg_+icosl*nsg_+isgl],
	I110 = emission_[inuu*ni_*nsg_+icosu*nsg_+isgl], 
	I010 = emission_[inul*ni_*nsg_+icosu*nsg_+isgl],
	I001 = emission_[inul*ni_*nsg_+icosl*nsg_+isgu], 
	I101 = emission_[inuu*ni_*nsg_+icosl*nsg_+isgu],
	I111 = emission_[inuu*ni_*nsg_+icosu*nsg_+isgu],
	I011 = emission_[inul*ni_*nsg_+icosu*nsg_+isgu];
      //cout << setprecision(10) << "trilin dir: " << I000 << " " << I100 << " " << I110 << " " << I010 << " " << I001 << " " << I101 << " " << I111 << " " << I011 << endl;
      double rationu = (nu-freq_[inul])/(freq_[inuu]-freq_[inul]),
	ratioi = (cosi-cosi_[icosl])/(cosi_[icosu]-cosi_[icosl]),
	ratiosg = (sgloc-surfgrav_[isgl])/(surfgrav_[isgu]-surfgrav_[isgl]);
      if (nsg_==1) ratiosg=0.; // no interpolation in sg
      Iem = I000
	+ (I100-I000)*rationu
	+ (I010-I000)*ratioi
	+ (I001-I000)*ratiosg
	+ (I110-I010-I100+I000)*rationu*ratioi
	+ (I011-I010-I001+I000)*ratioi*ratiosg
	+ (I101-I001-I100+I000)*rationu*ratiosg
	+ (I111-I011-I101-I110+I100+I001+I010-I000)*rationu*ratioi*ratiosg;
      //cout << "I interp= " << Iem << endl;
    }
  }else{
    // Average over cosi values
    // with bilinear interpol in nu,sg
    double I00=0., I01=0., I10=0., I11=0.;
    /* Using trapezoidal rule, I_integ = \int I(mu)*dmu, mu=cos(i)
       NB: in Garcia+14, they compute a flux because they don't raytrace,
       so they use F = 1/4pi * \int I(i) cos(i) di = 1/2 * \int I(mu) mu dmu,
       here we are not interested in the same quantity */
    double dcostot = 0.; // will contain \int d\mu (~1 but not exactly)
    for (size_t ii=0; ii<ni_-1; ++ii){
      double dcos = cosi_[ii+1]-cosi_[ii];
      I00 += 0.5*dcos*
	(emission_[inul*ni_*nsg_+(ii+1)*nsg_+isgl]
	 +emission_[inul*ni_*nsg_+ii*nsg_+isgl]);
      I01 += 0.5*dcos*
	(emission_[inul*ni_*nsg_+(ii+1)*nsg_+isgu]
	 +emission_[inul*ni_*nsg_+ii*nsg_+isgu]);
      I10 += 0.5*dcos*
	(emission_[inuu*ni_*nsg_+(ii+1)*nsg_+isgl]
	 +emission_[inuu*ni_*nsg_+ii*nsg_+isgl]);
      I11 += 0.5*dcos*
	(emission_[inuu*ni_*nsg_+(ii+1)*nsg_+isgu]
	 +emission_[inuu*ni_*nsg_+ii*nsg_+isgu]);
      dcostot+=dcos;
      
    } 

    // Normalizing (int d cos(i) is very close to 1 but not exactly 1)
    I00/=dcostot;
    I01/=dcostot;
    I10/=dcostot;
    I11/=dcostot;

    //cout << "\int dcos, and I bilin avg: " << dcostot << " " << I00 << " " << I01 << " " << I10 << " " << I11 << endl;

    double rationu = (nu-freq_[inul])/(freq_[inuu]-freq_[inul]),
      ratiosg = (sgloc-surfgrav_[isgl])/(surfgrav_[isgu]-surfgrav_[isgl]);
    if (nsg_==1) ratiosg=0.; // no interpolation in sg
    Iem = I00+(I10-I00)*rationu
      +(I01-I00)*ratiosg
      +(I11-I01-I10+I00)*rationu*ratiosg;
    //cout << "I interp= " << Iem << endl;
  }
  //cout << "return= " << Iem << endl;
  return Iem/1e3; // 1e3 factor translates from cgs to SI, gyoto speaks in SI

}
