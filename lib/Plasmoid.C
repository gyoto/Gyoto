/*
    Copyright 2019, 2020 Frederic Vincent, Thibaut Paumard, Nicolas Aimar

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

#include "GyotoUtils.h"
#include "GyotoPlasmoid.h"
#include "GyotoPhoton.h"
#include "GyotoWorldline.h"
#include "GyotoFactoryMessenger.h"

#include <iostream>
#include <cmath>
#include <string>
#include <cstdlib>
#include <float.h>
#include <sstream>
#include <string.h>

#ifdef GYOTO_USE_CFITSIO
#define throwCfitsioError(status) \
    { fits_get_errstatus(status, ermsg); GYOTO_ERROR(ermsg); }
#endif

using namespace std;
using namespace Gyoto;
using namespace Gyoto::Astrobj;

/// Properties
#include "GyotoProperty.h"
GYOTO_PROPERTY_START(Plasmoid, "Synchrotron-emitting orbiting plasmoid heated by magnetic reconnection")
GYOTO_PROPERTY_VECTOR_DOUBLE(Plasmoid, InitPosition, initPosition,
              "(t,r,theta,phi) initial position of plasmoid")
GYOTO_PROPERTY_VECTOR_DOUBLE(Plasmoid, InitVelocity, initVelocity,
              "(dr/dt,dtheta/dt,dphi/dt) initial 3-velocity "
              "of plasmoid")
GYOTO_PROPERTY_DOUBLE(Plasmoid, RadiusMax, radiusMax,
		          "Maximun radius of the Plasmoid")
GYOTO_PROPERTY_END(Plasmoid, UniformSphere::properties)

Plasmoid::Plasmoid() : 
  FitsRW(), 
  UniformSphere("Plasmoid"),
  flag_("None"),
  posSet_(false),
  posIni_(NULL),
  fourveldt_(NULL),
  radiusMax_(1.),
  varyRadius_("None"),
  filename_("None"),
  jnu_array_(NULL),
  anu_array_(NULL),
  freq_array_(NULL),
  t_inj_(1.)
{
  kind_="Plasmoid";
# ifdef GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << "done." << endl;
# endif

  posIni_= new double[4];
  fourveldt_= new double[4];
}

Plasmoid::Plasmoid(const Plasmoid& orig) :
  FitsRW(orig),
  UniformSphere(orig),
  flag_(orig.flag_),
  posSet_(orig.posSet_),
  posIni_(NULL),
  fourveldt_(NULL),
  radiusMax_(orig.radiusMax_),
  varyRadius_(orig.varyRadius_),
  filename_(orig.filename_),
  jnu_array_(NULL),
  anu_array_(NULL),
  freq_array_(NULL),
  t_inj_(orig.t_inj_)
{

  if(orig.posIni_){
	  posIni_= new double[4];
	  memcpy(posIni_,orig.posIni_, 4*sizeof(double));
  }

  if(orig.fourveldt_){
	  fourveldt_= new double[4];
	  memcpy(fourveldt_,orig.fourveldt_, 4*sizeof(double));
  }

  size_t ncells=0;
  size_t nnu=FitsRW::nnu(), nt=FitsRW::nt();
  ncells=nnu*nt;
  if (orig.jnu_array_){
    jnu_array_ = new double[ncells];
    memcpy(jnu_array_,orig.jnu_array_, ncells*sizeof(double));
  }
  if (orig.anu_array_){
    anu_array_ = new double[ncells];
    memcpy(anu_array_,orig.anu_array_, ncells*sizeof(double));
  }
  if (orig.freq_array_){
    freq_array_ = new double[nnu];
    memcpy(freq_array_,orig.freq_array_, nnu*sizeof(double));
  }
}


Plasmoid* Plasmoid::clone() const { return new Plasmoid(*this); }

Plasmoid::~Plasmoid() {
  if (debug()) cerr << "DEBUG: Plasmoid::~Plasmoid()\n";
  if (jnu_array_) delete [] jnu_array_;
  if (anu_array_) delete [] anu_array_;
  if (freq_array_) delete [] freq_array_;
}

string Plasmoid::className() const { return  string("Plasmoid"); }
string Plasmoid::className_l() const { return  string("Plasmoid"); }


void Plasmoid::radiativeQ(double Inu[], // output
              double Taunu[], // output
              double const nu_ems[], size_t nbnu, // input
              double dsem,
              state_t const &coord_ph,
              double const coord_obj[8]) const {

# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << endl;
# endif
  if (filename_=="None")
      GYOTO_ERROR("In Plamsoid RadiativeQ : filename_ not defined, please use file(string)");

  double tcur=coord_ph[0]*GYOTO_G_OVER_C_SQUARE*gg_->mass()/GYOTO_C/60.; // in min
  double t0 = posIni_[0]*GYOTO_G_OVER_C_SQUARE*gg_->mass()/GYOTO_C/60.;  // t0 in min

  // Defining jnus, anus
  double jnu[nbnu];
  double anu[nbnu];
  
  for (size_t ii=0; ii<nbnu; ++ii){
    // Initializing to <0 value to create errors if not updated
    // [ exp(-anu*ds) will explose ]
    jnu[ii]=-1.;
    anu[ii]=-1.;
  }

  
  // COMPUTE VALUES IN FUNCTION OF PHASE
  if (tcur<=t0){ // HEATING PHASE
    for (size_t ii=0; ii<nbnu; ++ii){
      jnu[ii]=0;
      anu[ii]=0;
    }
  }
  else{ // COOLING PHASE
    double tt=(tcur-t0)*60.; // in sec
    for (size_t ii=0; ii<nbnu; ++ii){
      jnu[ii]=FitsRW::interpolate(nu_ems[ii], tt, jnu_array_, freq_array_);
      anu[ii]=FitsRW::interpolate(nu_ems[ii], tt, anu_array_, freq_array_);
    }
  }


  // RETURNING TOTAL INTENSITY AND TRANSMISSION
  for (size_t ii=0; ii<nbnu; ++ii){
    double jnu_tot = jnu[ii], anu_tot = anu[ii];
    // expm1 is a precise implementation of exp(x)-1
    double em1=std::expm1(-anu_tot * dsem * gg_->unitLength());
    Taunu[ii] = em1+1.;
    Inu[ii] = anu_tot == 0. ? jnu_tot * dsem * gg_->unitLength() :
      -jnu_tot / anu_tot * em1;
    
    if (Inu[ii]<0.)
      GYOTO_ERROR("In Plasmoid::radiativeQ: Inu<0");
    if (Inu[ii]!=Inu[ii] or Taunu[ii]!=Taunu[ii])
      GYOTO_ERROR("In Plasmoid::radiativeQ: Inu or Taunu is nan");
    if (Inu[ii]==Inu[ii]+1. or Taunu[ii]==Taunu[ii]+1.)
      GYOTO_ERROR("In Plasmoid::radiativeQ: Inu or Taunu is infinite");
    
  }

}

void Plasmoid::motionType(std::string const type){
  if (type=="Helical" || type=="Equatorial")
  {
    flag_=type;
  }
  else
    GYOTO_ERROR("In Plasmoid::motonType: motion not recognized, please enter a valid motion type (Helical or Equatorial)");
}

SmartPointer<Metric::Generic> Plasmoid::metric() const { return gg_; }

void Plasmoid::metric(SmartPointer<Metric::Generic> gg) {
  UniformSphere::metric(gg);
  gg_=gg;
}

void Plasmoid::initPosition(std::vector<double> const &v) {
  posIni_[0] = v[0];
  posIni_[1] = v[1];
  posIni_[2] = v[2];
  posIni_[3] = v[3];
  posSet_=true;
}

std::vector<double> Plasmoid::initPosition() const {
  std::vector<double> v (4, 0.);
  v[0] = posIni_[0];
  v[1] = posIni_[1];
  v[2] = posIni_[2];
  v[3] = posIni_[3];
  return v;
}

void Plasmoid::initVelocity(std::vector<double> const &v) {
  if (!posSet_)
  	GYOTO_ERROR("In Plasmoid::initVelocity initial Position not defined");
  fourveldt_[1] = v[0];
  fourveldt_[2] = v[1];
  fourveldt_[3] = v[2];
  fourveldt_[0] = 1.;

  double sum = 0;
  double g[4][4];

  gg_->gmunu(g, posIni_);

  for (int i=0;i<4;++i) {
    for (int j=0;j<4;++j) {
      sum+=g[i][j]*fourveldt_[i]*fourveldt_[j];
    }
  }
  if (sum>=0)
 	GYOTO_ERROR("In Plasmoid::initVelocity Initial Velocity over C");

}

std::vector<double> Plasmoid::initVelocity() const {
  std::vector<double> v (3, 0.);
  v[0] = fourveldt_[1];
  v[1] = fourveldt_[2];
  v[2] = fourveldt_[3];
  return v;
}

void Plasmoid::initCoord(std::vector<double> const &v) {
  posIni_[0] = v[0];
  posIni_[1] = v[1];
  posIni_[2] = v[2];
  posIni_[3] = v[3];
  fourveldt_[0] = v[4];
  fourveldt_[1] = v[5];
  fourveldt_[2] = v[6];
  fourveldt_[3] = v[7];
}

std::vector<double> Plasmoid::initCoord() const {
  std::vector<double> v (8, 0.);
  v[0] = posIni_[0];
  v[1] = posIni_[1];
  v[2] = posIni_[2];
  v[3] = posIni_[3];
  v[4] = fourveldt_[0];
  v[5] = fourveldt_[1];
  v[6] = fourveldt_[2];
  v[7] = fourveldt_[3];
  return v;
}

void Plasmoid::radiusMax(double rr) {
	if (rr<0.2)
		GYOTO_ERROR("In Plasmoid::radiusMax radiusMax<0.2 (minimum value)");
	radiusMax_=rr;
}

double Plasmoid::radiusMax() const {
	return radiusMax_;
}

void Plasmoid::Radius(std::string vary) {
  if (vary=="Constant" || vary=="Varying") varyRadius_=vary;
  else
    GYOTO_ERROR("In Plasmoid::Radius operation on radius not recognized, please enter a valid operation (Constant or Varying)");
}

void Plasmoid::getCartesian(double const * const dates, size_t const n_dates,
          double * const x, double * const y, double * const z, 
          double * const xprime, double * const yprime, double * const zprime){
  // this yields the position of the center of the UnifSphere
  // at time t
  // fourveldt_ is the initial 3-velocity dxi/dt
  // vel is the 4-velocity dxnu/dtau

  if (n_dates!=1)
    GYOTO_ERROR("In Plasmoid::getCartesian n_dates!=1");

  if (flag_=="None")
      GYOTO_ERROR("In Plasmoid::getCartesian Motion not defined; motionType('Helical' or 'Equatorial'");

  double tt=dates[0];
  
  double r, theta, phi; // spherical coordinates
  double vel[4];
  
  if (flag_=="Helical") // Helical ejection
  {
    r = posIni_[1]+fourveldt_[1]*(tt-posIni_[0]);
    theta = posIni_[2];
    phi = posIni_[3] + posIni_[1]*posIni_[1]*fourveldt_[3]/fourveldt_[1]*(pow(posIni_[1],-1.)-pow(r,-1.)); // result of integrale of vphi over time
    //cout << "t, r, theta, phi = " << tt << ", " << r << ", " << theta << ", " << phi << endl;

  }
  else // Equatorial motion (Keplerian orbit)
  {
    if (posIni_[2]!=M_PI/2.)
      cout << "Warning input theta value incompatible with 'Equatorial' motion. Theta fixed to pi/2." << endl;
    getVelocity(posIni_, vel);

    r = posIni_[1];
    theta = M_PI/2.;
    phi = posIni_[3] + vel[3]/vel[0]*(tt-posIni_[0]);

  }
  // Convertion into cartesian coordinates
  x[0] = r*sin(theta)*cos(phi);
  y[0] = r*sin(theta)*sin(phi);
  z[0] = r*cos(theta);

  if (xprime!=NULL && yprime!=NULL && zprime!=NULL)
  {
    xprime[0] = r*sin(theta)*sin(phi)*vel[2];
    yprime[0] = -r*sin(theta)*cos(phi)*vel[2];
    zprime[0] = 0.;
  }
}

void Plasmoid::getVelocity(double const pos[4], double vel[4]){
  if (!gg_)
    GYOTO_ERROR("In Plasmoid::getVelocity Metric not set");
  if (flag_=="None")
    GYOTO_ERROR("In Plasmoid::getVelocity Motion not defined; motionType('Helical' or 'Equatorial'");
  
  if (flag_=="Helical") // Helical case
  {
  	vel[0] = 1.;
	vel[1] = fourveldt_[1];
	vel[2] = 0.;
	vel[3] = fourveldt_[3]*pow(posIni_[1]/pos[1],2.); // conservation of the Newtonian angular momentum [Ball et al. 2020]
	gg_->normalizeFourVel(pos, vel);

  }
  else // Equatorial case
  {
    gg_->circularVelocity(pos, vel);
  }
}


int Plasmoid::Impact(Photon* ph, size_t index, Properties *data){
	// Overload function of StandardAstrobj::Impact
	// This function update the radius of the plasmoid 
	// which increase linearly during the injection phase
	// before calling the StandardAstrobj function

	double radiusMin = 0.2;
  double t0 = posIni_[0]*GYOTO_G_OVER_C_SQUARE*gg_->mass()/GYOTO_C/60.;  // t0 in min

  size_t sz = ph -> parallelTransport()?16:8;
  state_t p1(sz);
  ph->getCoord(index, p1);
  double tcur = p1[0]*GYOTO_G_OVER_C_SQUARE*gg_->mass()/GYOTO_C/60.; //tcur in min

  
  if (varyRadius_== "Varying")
  {
    if (tcur<=t0) radius(radiusMin);
    else if (tcur<=t0+t_inj_) radius(radiusMin+(radiusMax_-radiusMin)*(tcur-t0)/t_inj_);
    else radius(radiusMax_);
  }
	else if (varyRadius_== "Constant") radius(radiusMax_);
  else{
    GYOTO_ERROR("In Plasmoid::Impact operation on radius not recognized. Use Radius('Constant' or 'Varying')");
  }

	return Standard::Impact(ph, index, data);
}




void Plasmoid::file(std::string const &f) {
  # ifdef GYOTO_USE_CFITSIO
    fitsRead(f);
  # else
    GYOTO_ERROR("This Gyoto has no FITS i/o");
  # endif
}

#ifdef GYOTO_USE_CFITSIO
vector<size_t> Plasmoid::fitsRead(string filename) {
  // Remove first char if it is "!"
  if (filename.substr(0,1)=="!")
    filename.erase(0,1);

  GYOTO_MSG << "Plasmoid reading FITS file: " << filename << endl;

  filename_ = filename;
  char*     pixfile   = const_cast<char*>(filename_.c_str());
  fitsfile* fptr      = NULL;
  int       status    = 0;
  double    tmpd;
  char      ermsg[31] = ""; // ermsg is used in throwCfitsioError()

  GYOTO_DEBUG << "Plasmoid::fitsRead: opening file" << endl;
  if (fits_open_file(&fptr, pixfile, 0, &status)) throwCfitsioError(status);

  ////// READ FITS KEYWORDS COMMON TO ALL TABLES ///////
  // These are: tmin, tmax, numin, numax
  
  string extname = "GYOTO FitsRW KEYS";
  fits_movnam_hdu(fptr, ANY_HDU, const_cast<char*>(extname.c_str()), 0, &status);

  GYOTO_DEBUG << "FitsRW::fitsRead(): read tmin_" << endl;
  fits_read_key(fptr, TDOUBLE, "GYOTO FitsRW tmin", &tmpd,
    NULL, &status);
  if (status) {
      if (status == KEY_NO_EXIST) status = 0; // not fatal
      else throwCfitsioError(status) ;
  } else FitsRW::tmin(tmpd); // tmin_ found

  GYOTO_DEBUG << "FitsRW::fitsRead(): read tmax_" << endl;
  fits_read_key(fptr, TDOUBLE, "GYOTO FitsRW tmax", &tmpd,
    NULL, &status);
  if (status) {
      if (status == KEY_NO_EXIST) status = 0; // not fatal
      else throwCfitsioError(status) ;
  } else FitsRW::tmax(tmpd); // tmax_ found

  GYOTO_DEBUG << "FitsRW::fitsRead(): read numin_" << endl;
  fits_read_key(fptr, TDOUBLE, "GYOTO FitsRW numin", &tmpd,
    NULL, &status);
  if (status) {
      if (status == KEY_NO_EXIST) status = 0; // not fatal
      else throwCfitsioError(status) ;
  } else FitsRW::numin(tmpd); // numin_ found

  GYOTO_DEBUG << "FitsRW::fitsRead(): read numax_" << endl;
  fits_read_key(fptr, TDOUBLE, "GYOTO FitsRW numax", &tmpd,
    NULL, &status);
  if (status) {
      if (status == KEY_NO_EXIST) status = 0; // not fatal
      else throwCfitsioError(status) ;
  } else FitsRW::numax(tmpd); // rmax_ found

  GYOTO_DEBUG << "FitsRW::fitsRead(): read t_inj" << endl;
  fits_read_key(fptr, TDOUBLE, "GYOTO FitsRW t_inj", &tmpd,
    NULL, &status);
  if (status){
    throwCfitsioError(status) ;
  }
  else t_inj_=tmpd; // t_inj_ found

  // READ EXTENSIONS
  vector<size_t> naxes_jnu = FitsRW::fitsReadHDU(fptr,"GYOTO FitsRW Jnu",
                  jnu_array_);

  vector<size_t> naxes_anu = FitsRW::fitsReadHDU(fptr,"GYOTO FitsRW Anu",
                  anu_array_);

  if (naxes_jnu[0]!=naxes_anu[0] || naxes_jnu[1]!=naxes_anu[1])
    throwError("In Plasmoid: jnu_array_ and anu_array_ dimensions "
         "do not agree");

  // freq array
  vector<size_t> naxes_freq = FitsRW::fitsReadHDU(fptr,"GYOTO FitsRW FREQUENCY",
                  freq_array_);
  
  if (naxes_freq[0]!=naxes_jnu[0])
    GYOTO_ERROR("In Plasmoid: nnu differ from jnu_array_ and freq_array");
  
  /*cout << "jnu read= " << endl;
  for (int ii=0;ii<60;ii++) cerr << jnu_array_[ii] << " " ;
  cout << endl;*/

  FitsRW::nnu(naxes_jnu[0]);
  FitsRW::nt(naxes_jnu[1]);

  return naxes_anu;

}
#endif