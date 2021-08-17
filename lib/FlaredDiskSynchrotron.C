/*
    Copyright 2019-2021 Frederic Vincent & Thibaut Paumard & Nicolas Aimar

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
#include "GyotoFlaredDiskSynchrotron.h"
#include "GyotoProperty.h"
#include "GyotoUtils.h"
#include "GyotoFactoryMessenger.h"

#ifdef GYOTO_USE_CFITSIO
#define throwCfitsioError(status) \
    { fits_get_errstatus(status, ermsg); GYOTO_ERROR(ermsg); }
#endif

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <cmath>
#include <limits>
#include <cstring>

using namespace std;
using namespace Gyoto;
using namespace Gyoto::Astrobj;

GYOTO_PROPERTY_START(FlaredDiskSynchrotron)
GYOTO_PROPERTY_FILENAME(FlaredDiskSynchrotron, File, file,
			"File name of FITS file containing data")
GYOTO_PROPERTY_DOUBLE(FlaredDiskSynchrotron, TimeTranslation_inMunit,
		      timeTranslation_inMunit,
		      "Shift simulation times by this amount, in GM/c3 unit")
GYOTO_PROPERTY_DOUBLE(FlaredDiskSynchrotron, HoverR, hoverR,
		      "Aspect ratio H/r of flared disk")
GYOTO_PROPERTY_DOUBLE_UNIT(FlaredDiskSynchrotron, NumberDensityMax, numberDensityMax,
		      "Maximum value of nb density in SI")
GYOTO_PROPERTY_DOUBLE(FlaredDiskSynchrotron, TemperatureMax, temperatureMax,
		      "Maximum value of temperature in K")
GYOTO_PROPERTY_DOUBLE(FlaredDiskSynchrotron, BetaAtMax, betaAtMax, 
          "Value of Beta at Maximum nb density")
GYOTO_PROPERTY_DOUBLE(FlaredDiskSynchrotron, MagnetizationParameter,
		      magnetizationParameter,
		      "Standard magnetization parameter (B^2/4pi) / (rho*c^2) "
		      "where rho is mass density")
GYOTO_PROPERTY_DOUBLE(FlaredDiskSynchrotron, KappaIndex, kappaIndex)
GYOTO_PROPERTY_DOUBLE(FlaredDiskSynchrotron, PolytropicIndex, polytropicIndex)
GYOTO_PROPERTY_END(FlaredDiskSynchrotron, Standard::properties)

FlaredDiskSynchrotron::FlaredDiskSynchrotron() :
Standard("FlaredDiskSynchrotron"), GridData2D(),
  filename_(""), hoverR_(0.), time_array_(NULL),
  density_(NULL), velocity_(NULL), Bvector_(NULL),
  numberDensityMax_cgs_(0.), temperatureMax_(0.),
  BMax_cgs_(0.), beta_(1.), magnetizationParameter_(1.), 
  deltat_(0.), flag_(false), gamm1_(5./3.)
{
  GYOTO_DEBUG << endl;
  spectrumKappaSynch_ = new Spectrum::KappaDistributionSynchrotron();
}

FlaredDiskSynchrotron::FlaredDiskSynchrotron(const FlaredDiskSynchrotron& o) :
  Standard(o), GridData2D(o),
  filename_(o.filename_), hoverR_(o.hoverR_), time_array_(NULL),
  density_(NULL), velocity_(NULL), Bvector_(NULL),
  numberDensityMax_cgs_(o.numberDensityMax_cgs_), temperatureMax_(o.temperatureMax_), beta_(o.beta_),
  magnetizationParameter_(o.magnetizationParameter_), deltat_(o.deltat_), flag_(o.flag_),
  gamm1_(o.gamm1_), BMax_cgs_(o.BMax_cgs_)
{
  GYOTO_DEBUG << endl;
  size_t ncells = 0;
  size_t nt=GridData2D::nt(), nphi=GridData2D::nphi(), nr=GridData2D::nr();
  if (o.density_) {
    density_ = new double[ncells = nt * nphi * nr];
    memcpy(density_, o.density_, ncells * sizeof(double));
  }
  if (o.velocity_) {
    velocity_ = new double[ncells = 2 * nt * nphi * nr];
    memcpy(velocity_, o.velocity_, ncells * sizeof(double));
  }
  if (o.Bvector_) {
    Bvector_= new double[ncells = 4* nt * nphi * nr];
    memcpy(Bvector_,o.Bvector_, ncells * sizeof(double));
  }
  if (o.time_array_) {
    time_array_ = new double[nt];
    memcpy(time_array_,o.time_array_, nt*sizeof(double));
  }
  if (o.spectrumKappaSynch_()) spectrumKappaSynch_=o.spectrumKappaSynch_->clone();
}
FlaredDiskSynchrotron* FlaredDiskSynchrotron::clone() const
{ return new FlaredDiskSynchrotron(*this); }

FlaredDiskSynchrotron::~FlaredDiskSynchrotron() {
  GYOTO_DEBUG << endl;
  if (density_) delete [] density_;
  if (velocity_) delete [] velocity_;
  if (time_array_) delete [] time_array_;
  if (Bvector_) delete [] Bvector_;
}

void FlaredDiskSynchrotron::file(std::string const &f) {
# ifdef GYOTO_USE_CFITSIO
  fitsRead(f);
# else
  GYOTO_ERROR("This Gyoto has no FITS i/o");
# endif
}

std::string FlaredDiskSynchrotron::file() const {
  return filename_;
}

void FlaredDiskSynchrotron::hoverR(double const hor) {
  double hmin=1e-4;
  if (hor < hmin){
    cerr << " " << endl;
    cerr << "***!!WARNING!!*** In FlaredDiskSynchrotron::hoverR: "
      "H/R very small, you might not resolve your disk; "
      "increase H/R or decrease GYOTO_T_TOL." << endl;
    cerr << " " << endl;
  }
  hoverR_ = hor;
}

double FlaredDiskSynchrotron::hoverR() const {
  return hoverR_;
}

void FlaredDiskSynchrotron::timeTranslation_inMunit(double const dt) {
  double tmin=GridData2D::tmin(), tmax=GridData2D::tmax();
  GridData2D::tmin(tmin-deltat_+dt);
  GridData2D::tmax(tmax-deltat_+dt);
  deltat_=dt;
  if (GridData2D::nt()==0)
    GYOTO_ERROR("In FlaredDiskSynchrotron::timeTranslation nt not yet defined");
  int nt = GridData2D::nt();
  if (!time_array_)
    GYOTO_ERROR("In FlaredDiskSynchrotron::timeTranslation time_array_ not defined. Please use FlaredDiskSynchrotron::file(string) before this function");
  for (int ii=0;ii<nt;ii++){
    time_array_[ii]+=dt;
  }
  if (GridData2D::tmin()>0)
    cout << "\nWARNING : tmin is positive, in most cases the stationnary boundary condition will be applied. You should decrease more timeTranslation_inMunit until at least " << -tmin << "\n" << endl;
}

double FlaredDiskSynchrotron::timeTranslation_inMunit() const {
  return deltat_;
}

double FlaredDiskSynchrotron::numberDensityMax() const {
  // Converts internal cgs dens to SI
  double dens=numberDensityMax_cgs_;
# ifdef HAVE_UDUNITS
  dens = Units::Converter("cm-3", "m-3")(dens);
# else
  GYOTO_WARNING << "Units ignored, please recompile Gyoto with --with-udunits"
		<< endl ;
# endif
  return dens; }
double FlaredDiskSynchrotron::numberDensityMax(string const &unit) const
{
  double dens = numberDensityMax();
  if (unit != "") {
# ifdef HAVE_UDUNITS
    dens = Units::Converter("m-3", unit)(dens);
# else
    GYOTO_WARNING << "Units ignored, please recompile Gyoto with --with-udunits"
		  << endl ;
# endif
  }
  return dens;
}
void FlaredDiskSynchrotron::numberDensityMax(double dens) {
# ifdef HAVE_UDUNITS
  dens = Units::Converter("m-3", "cm-3")(dens);
# else
  GYOTO_WARNING << "Units ignored, please recompile Gyoto with --with-udunits"
		<< endl ;
# endif
  numberDensityMax_cgs_=dens;
  BMax_cgs_=sqrt(8.*M_PI*numberDensityMax_cgs_*GYOTO_BOLTZMANN_CGS*temperatureMax_/beta_);
}
void FlaredDiskSynchrotron::numberDensityMax(double dens, string const &unit) {
  if (unit != "") {
# ifdef HAVE_UDUNITS
    dens = Units::Converter(unit, "m-3")(dens);
# else
    GYOTO_WARNING << "Units ignored, please recompile Gyoto with --with-udunits"
		  << endl ;
# endif
  }
  numberDensityMax(dens);
}

void FlaredDiskSynchrotron::temperatureMax(double tt) {
  temperatureMax_=tt;
  BMax_cgs_=sqrt(8.*M_PI*numberDensityMax_cgs_*GYOTO_BOLTZMANN_CGS*temperatureMax_/beta_);
}

double FlaredDiskSynchrotron::temperatureMax() const{return temperatureMax_;}

void FlaredDiskSynchrotron::polytropicIndex(double gamma) {gamm1_=gamma-1;}

double FlaredDiskSynchrotron::polytropicIndex() const {return gamm1_+1;}

void FlaredDiskSynchrotron::betaAtMax(double beta){
  if (beta<=0.)
    GYOTO_ERROR("In betaAtMax: beta must be >0!");
  beta_=beta;
  BMax_cgs_=sqrt(8.*M_PI*numberDensityMax_cgs_*GYOTO_BOLTZMANN_CGS*temperatureMax_/beta_);
}

double FlaredDiskSynchrotron::betaAtMax() const {return beta_;}

void FlaredDiskSynchrotron::magnetizationParameter(double rr) {
  magnetizationParameter_=rr;}

double FlaredDiskSynchrotron::magnetizationParameter()const{
  return magnetizationParameter_;}

void FlaredDiskSynchrotron::kappaIndex(double index) {
  spectrumKappaSynch_->kappaindex(index);
}

double FlaredDiskSynchrotron::kappaIndex()const{
  return spectrumKappaSynch_->kappaindex();
}

void FlaredDiskSynchrotron::copyDensity(double const *const density,
					size_t const naxes[3]) {
  GYOTO_DEBUG << endl;
  if (density_) {
    GYOTO_DEBUG << "delete [] density_;" << endl;
    delete [] density_; density_ = NULL;
  }
  size_t nt=GridData2D::nt(), nphi=GridData2D::nphi(), nr=GridData2D::nr();
  if (density) {
    if (nt != naxes[2] || nphi != naxes[1] || nr != naxes[0]) {
      GYOTO_DEBUG <<"grid dims changed, freeing velocity_" << endl;
      if (velocity_) { delete [] velocity_; velocity_= NULL; }
    }
    
    size_t nel;
    GridData2D::nt(naxes[2]);
    GridData2D::nphi(naxes[1]);
    GridData2D::nr(naxes[0]);
    //cout << naxes[0] << "," << naxes[1] << "," << naxes[2] << endl;
    if (!(nel=naxes[0] * naxes[1] * naxes[2]))
      GYOTO_ERROR( "dimensions can't be null");

    // NB: not updating dr_ contrary to PD
    GYOTO_DEBUG << "allocate density_;" << endl;
    density_ = new double[nel];
    GYOTO_DEBUG << "density >> density_" << endl;
    memcpy(density_, density, nel*sizeof(double));
  }

  //cout << "density stored= " << endl;
  //for (int ii=0;ii<30;ii++) cerr << density_[ii] << " " ;
  //cout << endl;

}

double const * FlaredDiskSynchrotron::getDensity() const { return density_; }


void FlaredDiskSynchrotron::copyVelocity(double const *const velocity,
					 size_t const naxes[3]) {
  GYOTO_DEBUG << endl;
  
  if (velocity_) {
    GYOTO_DEBUG << "delete [] velocity_;\n";
    delete [] velocity_; velocity_ = NULL;
  }
  size_t nt=GridData2D::nt(), nphi=GridData2D::nphi(), nr=GridData2D::nr();
  if (velocity) {
    if (!density_) GYOTO_ERROR("Please use copyDensity() before copyVelocity()");
    if (nt != naxes[2] || nphi != naxes[1] || nr != naxes[0])
      GYOTO_ERROR("density_ and velocity_ have inconsistent dimensions");
    GYOTO_DEBUG << "allocate velocity_;" << endl;
    size_t nel = 2*nt*nphi*nr;
    velocity_ = new double[nel];
    GYOTO_DEBUG << "velocity >> velocity_" << endl;
    memcpy(velocity_, velocity, nel*sizeof(double));
  }
  
  //cout << "velo stored= " << endl;
  //for (int ii=0;ii<60;ii++) cerr << velocity_[ii] << " " ;
  //cout << endl;
}

double const * FlaredDiskSynchrotron::getVelocity() const { return velocity_; }


void FlaredDiskSynchrotron::copyBvector(double const *const Bvector,
           size_t const naxes[3]) {
  GYOTO_DEBUG << endl;
  
  if (Bvector_) {
    GYOTO_DEBUG << "delete [] Bvector_;\n";
    delete [] Bvector_; Bvector_ = NULL;
  }
  size_t nt=GridData2D::nt(), nphi=GridData2D::nphi(), nr=GridData2D::nr();
  if (Bvector) {
    if (!density_) GYOTO_ERROR("Please use copyDensity() before copyBvector()");
    if (nt != naxes[2] || nphi != naxes[1] || nr != naxes[0])
      GYOTO_ERROR("density_ and Bvector_ have inconsistent dimensions");
    GYOTO_DEBUG << "allocate Bvector_;" << endl;
    size_t nel = 4*nt*nphi*nr;
    Bvector_ = new double[nel];
    GYOTO_DEBUG << "Bvector >> Bvector_" << endl;
    memcpy(Bvector_, Bvector, nel*sizeof(double));
    flag_=true;
  }
  
  //cout << "velo stored= " << endl;
  //for (int ii=0;ii<60;ii++) cerr << velocity_[ii] << " " ;
  //cout << endl;
}

double const * FlaredDiskSynchrotron::getBvector() const { return Bvector_; }

void FlaredDiskSynchrotron::copyTimeArray(double const *const time_array, size_t const ntimes) {
  GYOTO_DEBUG << endl;
  
  if (time_array_) {
    GYOTO_DEBUG << "delete [] time_array_;\n";
    delete [] time_array_; time_array_ = NULL;
  }
  size_t nt=GridData2D::nt();
  if (time_array) {
    if (nt != ntimes)
      GYOTO_ERROR("the given ntimes and nt from FITS file are inconsistent");
    GYOTO_DEBUG << "allocate time_array_;" << endl;
    time_array_ = new double[ntimes];
    GYOTO_DEBUG << "time_array >> time_array_" << endl;
    memcpy(time_array_, time_array, ntimes*sizeof(double));
  }
}

double const * FlaredDiskSynchrotron::getTimeArray() const { return time_array_; }

#ifdef GYOTO_USE_CFITSIO
vector<size_t> FlaredDiskSynchrotron::fitsRead(string filename) {
  // Remove first char if it is "!"
  if (filename.substr(0,1)=="!")
    filename.erase(0,1);

  GYOTO_MSG << "FlaredDiskSynchrotron reading FITS file: " << filename << endl;

  filename_ = filename;
  char*     pixfile   = const_cast<char*>(filename_.c_str());
  fitsfile* fptr      = NULL;
  int       status    = 0;
  double    tmpd;
  char      ermsg[31] = ""; // ermsg is used in throwCfitsioError()

  GYOTO_DEBUG << "FlaredDiskSynchrotron::fitsRead: opening file" << endl;
  if (fits_open_file(&fptr, pixfile, 0, &status)) throwCfitsioError(status) ;

  ////// READ FITS KEYWORDS COMMON TO ALL TABLES ///////
  // These are: tmin, tmax, rmin, rmax, phimin, phimax

  GYOTO_DEBUG << "FlaredDiskSynchrotron::fitsRead(): read tmin_" << endl;
  fits_read_key(fptr, TDOUBLE, "GYOTO GridData2D tmin", &tmpd,
  		NULL, &status);
  if (status) {
    if (status == KEY_NO_EXIST) status = 0; // not fatal
    else throwCfitsioError(status) ;
  } else GridData2D::tmin(tmpd+deltat_); // tmin_ found

  GYOTO_DEBUG << "FlaredDiskSynchrotron::fitsRead(): read tmax_" << endl;
  fits_read_key(fptr, TDOUBLE, "GYOTO GridData2D tmax", &tmpd,
  		NULL, &status);
  if (status) {
    if (status == KEY_NO_EXIST) status = 0; // not fatal
    else throwCfitsioError(status) ;
  } else GridData2D::tmax(tmpd+deltat_); // tmax_ found

  GYOTO_DEBUG << "FlaredDiskSynchrotron::fitsRead(): read rmin_" << endl;
  fits_read_key(fptr, TDOUBLE, "GYOTO GridData2D rmin", &tmpd,
  		NULL, &status);
  if (status) {
    if (status == KEY_NO_EXIST) status = 0; // not fatal
    else throwCfitsioError(status) ;
  } else GridData2D::rmin(tmpd); // rmin_ found
  
  GYOTO_DEBUG << "FlaredDiskSynchrotron::fitsRead(): read rmax_" << endl;
  fits_read_key(fptr, TDOUBLE, "GYOTO GridData2D rmax", &tmpd,
  		NULL, &status);
  if (status) {
    if (status == KEY_NO_EXIST) status = 0; // not fatal
    else throwCfitsioError(status) ;
  } else GridData2D::rmax(tmpd); // rmax_ found

  GYOTO_DEBUG << "FlaredDiskSynchrotron::fitsRead(): read phimin_" << endl;
  fits_read_key(fptr, TDOUBLE, "GYOTO GridData2D phimin", &tmpd,
  		NULL, &status);
  if (status) {
    if (status == KEY_NO_EXIST) status = 0; // not fatal
    else throwCfitsioError(status) ;
  } else GridData2D::phimin(tmpd); // phimin_ found
  
  GYOTO_DEBUG << "FlaredDiskSynchrotron::fitsRead(): read phimax_" << endl;
  fits_read_key(fptr, TDOUBLE, "GYOTO GridData2D phimax", &tmpd,
  		NULL, &status);
  if (status) {
    if (status == KEY_NO_EXIST) status = 0; // not fatal
    else throwCfitsioError(status) ;
  } else GridData2D::phimax(tmpd); // phimax_ found

  // READ EXTENSIONS

  // Density
  vector<size_t> naxes_dens = GridData2D::fitsReadHDU(fptr,"GYOTO GridData2D DENSITY",
						      density_);

  //cout << "density read= " << endl;
  //for (int ii=0;ii<30;ii++) cerr << density_[ii] << " " ;
  //cout << endl;

  // Velocity
  size_t length=2; // velocity is a 2-vector
  vector<size_t> naxes_velo = GridData2D::fitsReadHDU(fptr,"GYOTO GridData2D VELOCITY",
						      velocity_,
						      length);

  if (naxes_dens[0]!=naxes_velo[0] ||
      naxes_dens[1]!=naxes_velo[1] ||
      naxes_dens[2]!=naxes_velo[2])
    throwError("In FlaredDiskSynchro: density and velocity, dimensions "
         "do not agree");

  // 4-vector B
  string name="GYOTO GridData2D BVECTOR";
  fits_movnam_hdu(fptr, ANY_HDU, const_cast<char*>(name.c_str()), 0, &status);
  if (status==0){ // read only if the HDU exist
    flag_=true;
    size_t lengthB=4; // Bvector is a 4-vector
    vector<size_t> naxes_Bvec = GridData2D::fitsReadHDU(fptr,"GYOTO GridData2D BVECTOR",
                    Bvector_,
                    lengthB);

    if (naxes_dens[0]!=naxes_velo[0] || naxes_dens[0]!=naxes_Bvec[0] ||
        naxes_dens[1]!=naxes_velo[1] || naxes_dens[1]!=naxes_Bvec[1] ||
        naxes_dens[2]!=naxes_velo[2] || naxes_dens[2]!=naxes_Bvec[2])
      throwError("In FlaredDiskSynchro: density and B4vector dimensions "
  	       "do not agree");

  /*cout << "B4vector read= " << endl;
  for (int ii=0;ii<4*naxes_Bvec[1]*naxes_Bvec[2];ii++) cerr << Bvector_[ii] << " " ;
  cout << endl;*/
  }       

  // Time array
  vector<size_t> naxes_time = GridData2D::fitsReadHDU(fptr,"GYOTO GridData2D TIMEARRAY",
                  time_array_);
  
  if (naxes_time[0]!=naxes_dens[2])
    GYOTO_ERROR("In FlaredDiskSynchro: ntimes differ from density array and time_array");
  /*cout << "velo read= " << endl;
  for (int ii=0;ii<60;ii++) cerr << velocity_[ii] << " " ;
  cout << endl;*/

  GridData2D::nr(naxes_dens[0]);
  GridData2D::nphi(naxes_dens[1]);
  GridData2D::nt(naxes_dens[2]);
  //cout << "axes dens: " << naxes_dens[0] << " " << naxes_dens[1] << " " << naxes_dens[2] << endl;
  //cout << "axes velo: " << naxes_velo[0] << " " << naxes_velo[1] << " " << naxes_velo[2] << " " << naxes_velo[3] << endl;

  return naxes_velo;

}
#endif

double FlaredDiskSynchrotron::operator()(double const coord[4]) {
  // zpos: modulus of altitude above equatorial plane
  // rproj: radius projected in the equatorial plane
  double zpos=0., rproj=0.;
  
  switch (gg_ -> coordKind()) {
  case GYOTO_COORDKIND_SPHERICAL:
    rproj  = coord[1]*sin(coord[2]);
    zpos  = fabs(coord[1]*cos(coord[2]));
    break;
  case GYOTO_COORDKIND_CARTESIAN:
    zpos  = fabs(coord[3]);
    rproj  = sqrt(coord[1]*coord[1]+coord[2]*coord[2]);
    break;
  default:
    GYOTO_ERROR("FlaredDiskSynchrotron::operator(): unknown COORDKIND");
  }
  
  //cout << coord[0] << ", " << coord[1] << ", "  << coord[2] << ", " << coord[3] << endl;

  if (rproj < GridData2D::rmin() or rproj > GridData2D::rmax())
    return 1.; // outside disk
  
  double zdisk = rproj*hoverR_; // altitude of flared disk at this rproj
  return zpos - zdisk; // >0 outside, <0 inside flared disk
}

void FlaredDiskSynchrotron::radiativeQ(double Inu[], // output
		     double Taunu[], // output
		     double const nu_ems[], size_t nbnu, // input
		     double dsem,
		     state_t const &coord_ph,
		     double const coord_obj[8]) const {

  double rcyl=0.; // cylindrical radius
  double zz=0.; // height, z coord
  switch (gg_->coordKind()) {
  case GYOTO_COORDKIND_SPHERICAL:
    rcyl = coord_ph[1]*sin(coord_ph[2]);
    zz   = coord_ph[1]*cos(coord_ph[2]);
    break;
  case GYOTO_COORDKIND_CARTESIAN:
    rcyl = pow(coord_ph[1]*coord_ph[1]+coord_ph[2]*coord_ph[2], 0.5);
    zz   = coord_ph[3];
    break;
  default:
    GYOTO_ERROR("In FlaredDiskSynchrotron::radiativeQ: Unknown coordinate system kind");
  }

  double tt = coord_ph[0], phi = coord_ph[3];
  
  if (rcyl<GridData2D::rmin() || rcyl>GridData2D::rmax())
    throwError("In FlaredDiskSynchrotron::radiativeQ: r is not in grid!");
  if (phi<0. or phi>2.*M_PI)
    throwError("In FlaredDiskSynchrotron::radiativeQ: phi is not in 0,2pi!");
  // NB: phi is always in grid, and t might be outside, assuming stationnary
  // disk at t<tmin_ and t>tmax_

  //cout << "CALLING INTERPO FOR RHO" << endl;

  // Interpolating the density_ table, which contains
  // Sigma/r (surface density divided by radius), normalized.
  double Sigma_over_r_interpo=GridData2D::interpolate(tt,phi,rcyl,density_,time_array_);

  // Polytrop: p = kappa*rho^gamma
  double HH = hoverR_*rcyl; // Height of the flared disk at the local rcyl
  // 3D number density as a function of 2D surface density,
  // for a polytropic disk. See notes for details:
  double zfactor = 1.-zz*zz/(HH*HH),
    number_density = Sigma_over_r_interpo*numberDensityMax_cgs_
    *pow(zfactor,1./gamm1_);
  
  // 3D temperature now, see notes for details.
  // Perfect gas: p = (rho/mp)*k*T
  // Thus: T \propto rho^{gamma-1}
  double temperature = pow(Sigma_over_r_interpo,gamm1_)*temperatureMax_*zfactor;

  //cout << "stuff: " << Sigma_over_r_interpo << " " << hoverR_ << " " << rcyl << " " << zz << " " << HH << " " << zfactor << endl;

  double thetae = GYOTO_BOLTZMANN_CGS*temperature
    /(GYOTO_ELECTRON_MASS_CGS*GYOTO_C2_CGS);

  double hypergeom = Gyoto::hypergeom(kappaIndex(), thetae);

  double BB=0;
  double theta_mag=0.;

  if (flag_){ // magnetic field defined in FITS data
    size_t nr = GridData2D::nr(), nphi = GridData2D::nphi(),
    nt = GridData2D::nt(), nel = nt*nphi*nr-1;

    //cout << "CALLING INTERPO FOR B" << endl;
    double Bt=GridData2D::interpolate(tt, phi, rcyl, Bvector_, time_array_),
    Br=GridData2D::interpolate(tt, phi, rcyl, Bvector_+nel+1, time_array_),
    Btheta=GridData2D::interpolate(tt, phi, rcyl, Bvector_+2*(nel+1), time_array_),
    Bphi=GridData2D::interpolate(tt, phi, rcyl, Bvector_+3*(nel+1), time_array_);

    double b4vec[4]={Bt,Br,Btheta,Bphi}; // B 4-vector in BL frame
    double vel[4]; // 4-velocity of emitter
    const_cast<FlaredDiskSynchrotron*>(this)->getVelocity(coord_obj, vel);
    gg_->projectFourVect(&coord_ph[0],b4vec,vel); //Projection of the 4-vector B to 4-velocity to be in the rest frame of the emitter
    double photon_emframe[4]; // photon tgt vector projected in comoving frame
    for (int ii=0;ii<4;ii++){
      photon_emframe[ii]=coord_ph[ii+4]+vel[ii]*gg_->ScalarProd(&coord_ph[0],&coord_ph[4],vel);
    }
    double bnorm = gg_->norm(&coord_ph[0],b4vec);
    double lnorm = gg_->norm(&coord_ph[0],photon_emframe);
    double lscalb = gg_->ScalarProd(&coord_ph[0],photon_emframe,b4vec);
    theta_mag = acos(lscalb/(lnorm*bnorm));

    double sth = sin(theta_mag);//, cth = cos(theta_mag);
    if (sth==0.) GYOTO_ERROR("In FlaredDiskSynchrotron::radiativeQ: "
          "theta_mag is zero leads to undefined emission");
    if (sth!=sth) GYOTO_ERROR("In FlaredDiskSynchrotron::radiativeQ: "
          "theta_mag is nan");
    spectrumKappaSynch_->angle_averaged(0);

    //cout << "4 composant of B: " << b4vec[0] << "," << b4vec[1] << "," << b4vec[2] << "," << b4vec[3] << endl;
    BB = BMax_cgs_*sqrt(pow(b4vec[1],2.)+pow(b4vec[2],2.)+pow(b4vec[3],2.)); // norm of the 3-vector B
    
    
  }
  else{
    spectrumKappaSynch_->angle_averaged(1); // impose angle-averaging

    BB = sqrt(4.*M_PI*magnetizationParameter_
       *GYOTO_PROTON_MASS_CGS * GYOTO_C_CGS * GYOTO_C_CGS
       *number_density);
  }

  double nu0 = GYOTO_ELEMENTARY_CHARGE_CGS*BB
    /(2.*M_PI*GYOTO_ELECTRON_MASS_CGS*GYOTO_C_CGS); // cyclotron freq

  //cout << "jet stuff= " << coord_ph[1] << " " << coord_ph[2] << " " << zz << " " << rcyljetbase << " " << rcyl << " " << number_density << " " << thetae << " " << temperatureSlope_ << " " << nu0 << endl;
  //cout << "jet zz,rcyl,th,ph,ne,Te= " <<  zz << " " << rcyl << " " << coord_ph[2] << " " << coord_ph[3] << " " << number_density << " " << temperature << endl;
  // Use that line for Compton study:
  //cout <<  zz << " " << rcyl << " " << number_density << " " << temperature << endl;

  // KAPPA-DISTRIB SYNCHROTRON
  double jnu_synch_kappa[nbnu], anu_synch_kappa[nbnu];
  for (size_t ii=0; ii<nbnu; ++ii){
    // Initializing to <0 value to create errors if not updated
    jnu_synch_kappa[ii]=-1.;
    anu_synch_kappa[ii]=-1.;
  }

  spectrumKappaSynch_->numberdensityCGS(number_density);
  spectrumKappaSynch_->cyclotron_freq(nu0);
  spectrumKappaSynch_->thetae(thetae);
  spectrumKappaSynch_->angle_B_pem(theta_mag); // so we don't care about angle
  spectrumKappaSynch_->hypergeometric(hypergeom);
  //cout << "jet stuff for kappa: " << nu_ems[0] << " " << number_density << " " << nu0 << " " << thetae << " " << BB << " " << temperature << " " << hypergeom << endl;
  spectrumKappaSynch_->radiativeQ(jnu_synch_kappa,anu_synch_kappa,
				  nu_ems,nbnu);


  // RETURNING TOTAL INTENSITY AND TRANSMISSION
  for (size_t ii=0; ii<nbnu; ++ii){

    double jnu_tot = jnu_synch_kappa[ii],
      anu_tot = anu_synch_kappa[ii];
      //cout << jnu_tot << ", " << anu_tot << endl;

    //cout << "in jet stuff: " << number_density << " " << nu0 << " " << thetae << " " << hypergeom << " " << jnu_tot << " " << anu_tot << " " << dsem << endl;

    //cout << "at r,th= " << coord_ph[1] << " " << coord_ph[2] << endl;
    //cout << "jet jnu anu kappa= " << jnu_tot << " " << anu_tot << endl; //x" " << jnu_tot/anu_tot << " " << dsem << endl;

    // expm1 is a precise implementation of exp(x)-1
    double em1=std::expm1(-anu_tot * dsem * gg_->unitLength());
    Taunu[ii] = em1+1.;
    Inu[ii] = anu_tot == 0. ? jnu_tot * dsem * gg_->unitLength() :
      -jnu_tot / anu_tot * em1;

    if (Inu[ii]<0.)
      GYOTO_ERROR("In FlaredDiskSynchrotron::radiativeQ: Inu<0");
    if (Inu[ii]!=Inu[ii] or Taunu[ii]!=Taunu[ii]){
      GYOTO_ERROR("In FlaredDiskSynchrotron::radiativeQ: Inu or Taunu is nan");
    }
    if (Inu[ii]==Inu[ii]+1. or Taunu[ii]==Taunu[ii]+1.)
      GYOTO_ERROR("In FlaredDiskSynchrotron::radiativeQ: Inu or Taunu is infinite");

  }
}

void FlaredDiskSynchrotron::getVelocity(double const pos[4], double vel[4]){

  double rcyl=0.; // cylindrical radius
  double zz=0.; // height, z coord
  switch (gg_->coordKind()) {
  case GYOTO_COORDKIND_SPHERICAL:
    rcyl = pos[1]*sin(pos[2]);
    zz   = pos[1]*cos(pos[2]);
    break;
  case GYOTO_COORDKIND_CARTESIAN:
    rcyl = pow(pos[1]*pos[1]+pos[2]*pos[2], 0.5);
    zz   = pos[3];
    break;
  default:
    GYOTO_ERROR("In FlaredDiskSynchrotron::getVelocity: "
		" Unknown coordinate system kind");
  }

  double tt = pos[0], phi = pos[3];
  
  if (rcyl<GridData2D::rmin() || rcyl>GridData2D::rmax())
    throwError("In FlaredDiskSynchrotron::getVelocity: r is not in grid!");
  if (phi<0. or phi>2*M_PI)
    throwError("In FlaredDiskSynchrotron::getVelocity phi is not in 0;2pi!");
  // NB: phi is always in grid, and t might be outside, assuming stationnary
  // disk at t<tmin_ and t>tmax_

  size_t nr = GridData2D::nr(), nphi = GridData2D::nphi(),
    nt = GridData2D::nt(), nel = nt*nphi*nr-1;

  // first half od velocity_ contains all values of dr/dt
  // second hald contains all values of dphi/dt
  //cout << "In VELO R: " << velocity_[0] << " " << velocity_[1] << " " << velocity_[nr-1] << " " << velocity_[nr] << " " << velocity_[nr*nphi-1] << endl;
  //cout << "In VELO PHI: " << velocity_[0+nr*nphi] << " " << velocity_[1+nr*nphi] << " " << velocity_[nr-1+nr*nphi] << " " << velocity_[nr+nr*nphi] << " " << velocity_[nr*nphi-1+nr*nphi] << endl;
  //cout << "CALLING INTERPO FOR DR/DT" << endl;
  double drdt_interpo=GridData2D::interpolate(tt,phi,rcyl,velocity_,time_array_);
  //cout << "CALLING INTERPO FOR DPHI/DT" << endl;
  double dphidt_interpo=GridData2D::interpolate(tt,phi,rcyl,velocity_+nel+1,time_array_);

  switch (gg_->coordKind()) {
  case GYOTO_COORDKIND_SPHERICAL:
    {
      vel[1] = drdt_interpo;
      vel[2] = 0.;
      vel[3] = dphidt_interpo;
      //cout << "IN FLARED: pos and vel= " << pos[1] << " " << pos[2] << " " << pos[3] << " " << vel[1] << " " << vel[2] << " " << vel[3] << endl;
      vel[0] = gg_->SysPrimeToTdot(pos, vel+1);
      vel[1] *= vel[0];
      vel[3] *= vel[0];
    }
    break;
  case GYOTO_COORDKIND_CARTESIAN:
    GYOTO_ERROR("FlaredDiskSynchro::getVelocity(): metric must be in "
		"spherical coordinates if velocity field is provided");
    break;
  default:
    GYOTO_ERROR("FlaredDiskSynchro::getVelocity(): unknown COORDKIND");
  }
  
}

bool FlaredDiskSynchrotron::isThreadSafe() const {
  return Standard::isThreadSafe();
}

