/*
    Copyright 2019-2021 Frederic Vincent & Thibaut Paumard

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
#include "GyotoThinDiskGridIntensity.h"
#include "GyotoProperty.h"
#include "GyotoUtils.h"
#include "GyotoFactoryMessenger.h"
#include "GyotoKerrBL.h"

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

GYOTO_PROPERTY_START(ThinDiskGridIntensity)
GYOTO_PROPERTY_FILENAME(ThinDiskGridIntensity, File, file,
			"File name of FITS file containing data")
GYOTO_PROPERTY_DOUBLE(ThinDiskGridIntensity, TimeTranslation_inMunit,
		      timeTranslation_inMunit,
		      "Shift simulation times by this amount, in GM/c3 unit")
//GYOTO_PROPERTY_DOUBLE(ThinDiskGridIntensity, dt, dt)
GYOTO_PROPERTY_END(ThinDiskGridIntensity, ThinDisk::properties)

ThinDiskGridIntensity::ThinDiskGridIntensity() :
ThinDisk("ThinDiskGridIntensity"), GridData2D(),
  filename_(""), time_array_(NULL),
  intensity_(NULL), 
  deltat_(0.)
{
  GYOTO_DEBUG << endl;
}

ThinDiskGridIntensity::ThinDiskGridIntensity(const ThinDiskGridIntensity& o) :
  ThinDisk(o), GridData2D(o),
  filename_(o.filename_), time_array_(NULL),
  intensity_(NULL), 
  deltat_(o.deltat_)
{
  GYOTO_DEBUG << endl;
  size_t ncells = 0;
  size_t nt=GridData2D::nt(), nphi=GridData2D::nphi(), nr=GridData2D::nr();
  if (o.intensity_) {
    intensity_ = new double[ncells = nt * nphi * nr];
    memcpy(intensity_, o.intensity_, ncells * sizeof(double));
  }
  if (o.time_array_) {
    time_array_ = new double[nt];
    memcpy(time_array_,o.time_array_, nt*sizeof(double));
  }
}
ThinDiskGridIntensity* ThinDiskGridIntensity::clone() const
{ return new ThinDiskGridIntensity(*this); }

ThinDiskGridIntensity::~ThinDiskGridIntensity() {
  GYOTO_DEBUG << endl;
  if (intensity_) delete [] intensity_;
  if (time_array_) delete [] time_array_;
}

void ThinDiskGridIntensity::file(std::string const &f) {
# ifdef GYOTO_USE_CFITSIO
  fitsRead(f);
# else
  GYOTO_ERROR("This Gyoto has no FITS i/o");
# endif
}

std::string ThinDiskGridIntensity::file() const {
  return filename_;
}

void ThinDiskGridIntensity::timeTranslation_inMunit(double const dt) {
  double tmin=GridData2D::tmin(), tmax=GridData2D::tmax();
  GridData2D::tmin(tmin-deltat_+dt);
  GridData2D::tmax(tmax-deltat_+dt);
  deltat_=dt;
  if (GridData2D::nt()==0)
    GYOTO_ERROR("In ThinDiskGridIntensity::timeTranslation nt not yet defined");
  int nt = GridData2D::nt();
  if (!time_array_)
    GYOTO_ERROR("In ThinDiskGridIntensity::timeTranslation time_array_ not defined. Please use ThinDiskGridIntensity::file(string) before this function");
  for (int ii=0;ii<nt;ii++){
    time_array_[ii]+=dt;
  }
  if (GridData2D::tmin()>0)
    cout << "\nWARNING : tmin is positive, in most cases the stationnary boundary condition will be applied. You should decrease more timeTranslation_inMunit until at least " << -tmin << "\n" << endl;
}

double ThinDiskGridIntensity::timeTranslation_inMunit() const {
  return deltat_;
}

/*void ThinDiskGridIntensity::dt(double dd) {
  GridData2D::dt(dd);
}*/

void ThinDiskGridIntensity::copyIntensity(double const *const intensity,
					  size_t const naxes[3]) {
  GYOTO_DEBUG << endl;
  if (intensity_) {
    GYOTO_DEBUG << "delete [] intensity_;" << endl;
    delete [] intensity_; intensity_ = NULL;
  }
  size_t nt=GridData2D::nt(), nphi=GridData2D::nphi(), nr=GridData2D::nr();
  if (intensity) {
    size_t nel;
    GridData2D::nt(naxes[2]);
    GridData2D::nphi(naxes[1]);
    GridData2D::nr(naxes[0]);
    //cout << naxes[0] << "," << naxes[1] << "," << naxes[2] << endl;
    if (!(nel=naxes[0] * naxes[1] * naxes[2]))
      GYOTO_ERROR( "dimensions can't be null");

    // NB: not updating dr_ contrary to PD
    GYOTO_DEBUG << "allocate intensity_;" << endl;
    intensity_ = new double[nel];
    GYOTO_DEBUG << "intensity >> intensity_" << endl;
    memcpy(intensity_, intensity, nel*sizeof(double));
  }

  //cout << "intensity stored= " << endl;
  //for (int ii=0;ii<30;ii++) cerr << intensity_[ii] << " " ;
  //cout << endl;

}

double const * ThinDiskGridIntensity::getIntensity() const {
  return intensity_; }


void ThinDiskGridIntensity::copyTimeArray(double const *const time_array,
					  size_t const ntimes) {
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

double const * ThinDiskGridIntensity::getTimeArray() const {
  return time_array_; }

#ifdef GYOTO_USE_CFITSIO
vector<size_t> ThinDiskGridIntensity::fitsRead(string filename) {
  // Remove first char if it is "!"
  if (filename.substr(0,1)=="!")
    filename.erase(0,1);

  GYOTO_MSG << "ThinDiskGridIntensity reading FITS file: " << filename << endl;

  filename_ = filename;
  char*     pixfile   = const_cast<char*>(filename_.c_str());
  fitsfile* fptr      = NULL;
  int       status    = 0;
  double    tmpd;
  char      ermsg[31] = ""; // ermsg is used in throwCfitsioError()

  GYOTO_DEBUG << "ThinDiskGridIntensity::fitsRead: opening file" << endl;
  if (fits_open_file(&fptr, pixfile, 0, &status)) throwCfitsioError(status) ;

  ////// READ FITS KEYWORDS COMMON TO ALL TABLES ///////
  // These are: tmin, tmax, rmin, rmax, phimin, phimax

  GYOTO_DEBUG << "ThinDiskGridIntensity::fitsRead(): read tmin_" << endl;
  fits_read_key(fptr, TDOUBLE, "GYOTO GridData2D tmin", &tmpd,
  		NULL, &status);
  if (status) {
    if (status == KEY_NO_EXIST) status = 0; // not fatal
    else throwCfitsioError(status) ;
  } else GridData2D::tmin(tmpd+deltat_); // tmin_ found

  GYOTO_DEBUG << "ThinDiskGridIntensity::fitsRead(): read tmax_" << endl;
  fits_read_key(fptr, TDOUBLE, "GYOTO GridData2D tmax", &tmpd,
  		NULL, &status);
  if (status) {
    if (status == KEY_NO_EXIST) status = 0; // not fatal
    else throwCfitsioError(status) ;
  } else GridData2D::tmax(tmpd+deltat_); // tmax_ found

  GYOTO_DEBUG << "ThinDiskGridIntensity::fitsRead(): read rmin_" << endl;
  fits_read_key(fptr, TDOUBLE, "GYOTO GridData2D rmin", &tmpd,
  		NULL, &status);
  if (status) {
    if (status == KEY_NO_EXIST) status = 0; // not fatal
    else throwCfitsioError(status) ;
  } else GridData2D::rmin(tmpd); // rmin_ found
  
  GYOTO_DEBUG << "ThinDiskGridIntensity::fitsRead(): read rmax_" << endl;
  fits_read_key(fptr, TDOUBLE, "GYOTO GridData2D rmax", &tmpd,
  		NULL, &status);
  if (status) {
    if (status == KEY_NO_EXIST) status = 0; // not fatal
    else throwCfitsioError(status) ;
  } else GridData2D::rmax(tmpd); // rmax_ found

  GYOTO_DEBUG << "ThinDiskGridIntensity::fitsRead(): read phimin_" << endl;
  fits_read_key(fptr, TDOUBLE, "GYOTO GridData2D phimin", &tmpd,
  		NULL, &status);
  if (status) {
    if (status == KEY_NO_EXIST) status = 0; // not fatal
    else throwCfitsioError(status) ;
  } else GridData2D::phimin(tmpd); // phimin_ found
  
  GYOTO_DEBUG << "ThinDiskGridIntensity::fitsRead(): read phimax_" << endl;
  fits_read_key(fptr, TDOUBLE, "GYOTO GridData2D phimax", &tmpd,
  		NULL, &status);
  if (status) {
    if (status == KEY_NO_EXIST) status = 0; // not fatal
    else throwCfitsioError(status) ;
  } else GridData2D::phimax(tmpd); // phimax_ found

  // READ EXTENSIONS

  // Intensity
  vector<size_t> naxes_intens =
    GridData2D::fitsReadHDU(fptr,
			    "GYOTO GridData2D INTENSITY",
			    intensity_);

  //cout << "intensity read= " << endl;
  //for (int ii=0;ii<30;ii++) cerr << intensity_[ii] << " " ;
  //cout << endl;

  // Time array
  vector<size_t> naxes_time =
    GridData2D::fitsReadHDU(fptr,
			    "GYOTO GridData2D TIMEARRAY",
			    time_array_);
  
  if (naxes_time[0]!=naxes_intens[2])
    GYOTO_ERROR("In FlaredDiskSynchro: ntimes differ from intensity array and time_array");
  /*cout << "velo read= " << endl;
  for (int ii=0;ii<60;ii++) cerr << velocity_[ii] << " " ;
  cout << endl;*/

  GridData2D::nr(naxes_intens[0]);
  GridData2D::nphi(naxes_intens[1]);
  GridData2D::nt(naxes_intens[2]);
  //cout << "axes intens: " << naxes_intens[0] << " " << naxes_intens[1] << " " << naxes_intens[2] << endl;
  //cout << "axes velo: " << naxes_velo[0] << " " << naxes_velo[1] << " " << naxes_velo[2] << " " << naxes_velo[3] << endl;

  return naxes_intens;

}
#endif

double ThinDiskGridIntensity::emission(double nu, double,
			    state_t const &coord_ph,
			    double const coord_obj[8]) const{

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
    GYOTO_ERROR("In ThinDiskGridIntensity::radiativeQ: Unknown coordinate system kind");
  }

  double tt = coord_ph[0], phi = coord_ph[3];
  
  if (rcyl<GridData2D::rmin() || rcyl>GridData2D::rmax())
    return 0.;
  if (phi<0. or phi>2.*M_PI)
    throwError("In ThinDiskGridIntensity::radiativeQ: phi is not in 0,2pi!");
  // NB: phi is always in grid, and t might be outside, assuming stationnary
  // disk at t<tmin_ and t>tmax_

  //cout << "CALLING INTERPO FOR RHO" << endl;

  //cout << "tmin max t phi rcyl: " << GridData2D::tmin() << " " << GridData2D::tmax() << " " << tt << " " << phi << " " << rcyl << endl;

  // Interpolating the intensity_ table
  double intensity_interpo=GridData2D::interpolate(tt,phi,rcyl,intensity_,time_array_);

  //cout << "interpo intens= "<< intensity_interpo << endl;

  return intensity_interpo;
}

void ThinDiskGridIntensity::getVelocity(double const pos[4], double vel[4]){

  string kin = gg_->kind();
  if (kin != "KerrBL")
    GYOTO_ERROR("ThinDiskGridIntensity: KerrBL needed!");
  double SPIN = static_cast<SmartPointer<Metric::KerrBL> >(gg_) -> spin();
  
  double risco = gg_->getRms(); // prograde Kerr ISCO

  double rr = pos[1];
  if (rr > risco){
    // Keplerian velocity above ISCO
    gg_ -> circularVelocity(pos, vel, 1);
  }else{
    // See formulas in Gralla, Lupsasca & Marrone 2020, Eqs B8-B14
    // initally from Cunnigham 1975
    double lambda_ms = (risco*risco - 2.*SPIN*sqrt(risco) + SPIN*SPIN)/(pow(risco,1.5) - 2.*sqrt(risco) + SPIN),
      gamma_ms = sqrt(1.-2./(3.*risco)),
      delta = rr*rr - 2.*rr + SPIN*SPIN,
      hh = (2.*rr - SPIN*lambda_ms)/delta;

    vel[0] = gamma_ms*(1.+2./rr*(1.+hh)); // this is: -Ems*g^{tt} + Lms*g^{tp}
    vel[1] = -sqrt(2./(3.*risco))*pow(risco/rr-1.,1.5); // this is: -sqrt{(-1 - g_{tt}*u^t - g_{pp}*u^p - 2*g_{tp}*u^t*u^p)/grr}
    vel[2] = 0.;
    vel[3] = gamma_ms/(rr*rr)*(lambda_ms+SPIN*hh);

    //cout << "u2 = " << gg_->ScalarProd(pos,vel,vel) << endl;
  }

  
}

bool ThinDiskGridIntensity::isThreadSafe() const {
  return ThinDisk::isThreadSafe();
}

