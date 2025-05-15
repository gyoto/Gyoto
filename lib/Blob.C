/*
    Copyright 2019, 2020 Frederic Vincent, Thibaut Paumard

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
#include "GyotoBlob.h"
#include "GyotoPhoton.h"
#include "GyotoFactoryMessenger.h"
#include "GyotoKerrBL.h"

#include <iostream>
#include <cmath>
#include <string>
#include <cstdlib>
#include <float.h>
#include <sstream>
#include <string.h>

using namespace std;
using namespace Gyoto;
using namespace Gyoto::Astrobj;
using namespace Eigen;

/// Properties
#include "GyotoProperty.h"
GYOTO_PROPERTY_START(Blob, "Synchrotron-emitting orbiting blob of plasma")
GYOTO_PROPERTY_BOOL(Blob,
		    TimeGaussianModulated, NoTimeGaussianModulated,
		    timeGaussianModulated)
GYOTO_PROPERTY_BOOL(Blob,
		    SpaceGaussianModulated, NoSpaceGaussianModulated,
		    spaceGaussianModulated)
GYOTO_PROPERTY_VECTOR_DOUBLE(Blob, Init4Coord, init4Coord,
			     "Initial 4-position vector of the Blob, eg (t,r,theta,phi)")
GYOTO_PROPERTY_VECTOR_DOUBLE(Blob, Init3Velo, init3Velo,
			     "Initial coordinate 3-velocity of the Blob, eg (dr/dt,dtheta/dt,dphi/dt).")
GYOTO_PROPERTY_STRING(Blob, BlobMotionType, blobMotionType,
		      "\"Equatorial\" (default), \"HelicalConical\", or \"HelicalCylindrical\".")
GYOTO_PROPERTY_DOUBLE_UNIT(Blob, NumberDensity, numberDensity,
			   "cgs number density, constant through blob")
GYOTO_PROPERTY_DOUBLE(Blob, Temperature, temperature,
		      "temperature, constant through blob")
GYOTO_PROPERTY_DOUBLE_UNIT(Blob, TimeRef, timeRef,
			   "time of max of Gaussian evolution "
			   "of blob density and temperature")
GYOTO_PROPERTY_DOUBLE_UNIT(Blob, TimeSigma, timeSigma,
			   "temporal sigma of Gaussian evolution "
			   "of blob density and temperature")
GYOTO_PROPERTY_DOUBLE(Blob, MagnetizationParameter,
		      magnetizationParameter,
		      "magnetization parameter")
GYOTO_PROPERTY_DOUBLE(Blob, KappaIndex, kappaIndex,
		      "PL index of kappa-synchrotron")
GYOTO_PROPERTY_STRING(Blob, ElectronDistribution, electronDistribution,
		      "\"Thermal\" (default), or \"Kappa\".")
GYOTO_PROPERTY_END(Blob, UniformSphere::properties)

#define USE_IPOLE_FORMALISM 0

Blob::Blob() :
UniformSphere("Blob"),
  time_gauss_modulated_(false),
  space_gauss_modulated_(false),
  init4Coord_(NULL),
  init3Velo_(NULL),
  blobMotionType_("Equatorial"),
  numberDensity_cgs_(1.),
  temperature_(1.),
  timeRef_M_(1.),
  timeSigma_M_(1.),
  magnetizationParameter_(1.),
  kappaIndex_(1.),
  spectrumKappaSynch_(NULL),
  spectrumPLSynch_(NULL),
  spectrumThermalSynch_(NULL),
  electronDistrib_("Thermal")
{
  kind_="Blob";
# ifdef GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << "done." << endl;
# endif
  spectrumKappaSynch_ = new Spectrum::KappaDistributionSynchrotron();
  spectrumPLSynch_ = new Spectrum::PowerLawSynchrotron();
  spectrumThermalSynch_ = new Spectrum::ThermalSynchrotron();

  init4Coord_ = new double[4];
  init3Velo_ = new double[3];
  for (int ii=0;ii<3;ii++){
    init4Coord_[ii]=0.;
    init3Velo_[ii]=0.;
  }
  init4Coord_[3]=0.;
}

Blob::Blob(const Blob& orig) :
  UniformSphere(orig),
  time_gauss_modulated_(orig.time_gauss_modulated_),
  space_gauss_modulated_(orig.space_gauss_modulated_),
  init4Coord_(NULL),
  init3Velo_(NULL),
  blobMotionType_(orig.blobMotionType_),
  numberDensity_cgs_(orig.numberDensity_cgs_),
  temperature_(orig.temperature_),
  timeRef_M_(orig.timeRef_M_),
  timeSigma_M_(orig.timeSigma_M_),
  kappaIndex_(orig.kappaIndex_),
  magnetizationParameter_(orig.magnetizationParameter_),
  spectrumKappaSynch_(NULL),
  spectrumPLSynch_(NULL),
  spectrumThermalSynch_(NULL),
  electronDistrib_(orig.electronDistrib_)
{
  if (orig.spectrumKappaSynch_()) spectrumKappaSynch_=orig.spectrumKappaSynch_->clone();
  if (orig.spectrumPLSynch_()) spectrumPLSynch_=orig.spectrumPLSynch_->clone();
  if (orig.spectrumThermalSynch_()) spectrumThermalSynch_=orig.spectrumThermalSynch_->clone();

  init4Coord_ = new double[4];
  init3Velo_ = new double[3];
  for (int ii=0;ii<3;ii++){
    init4Coord_[ii]=orig.init4Coord_[ii];
    init3Velo_[ii]=orig.init3Velo_[ii];
  }
  init4Coord_[3]=orig.init4Coord_[3];
}

Blob* Blob::clone() const { return new Blob(*this); }

Blob::~Blob() {
  if (debug()) cerr << "DEBUG: Blob::~Blob()\n";
}

bool Blob::timeGaussianModulated() const
{return time_gauss_modulated_;}
void Blob::timeGaussianModulated(bool timemod)
{
  time_gauss_modulated_=timemod;
}

bool Blob::spaceGaussianModulated() const
{return space_gauss_modulated_;}
void Blob::spaceGaussianModulated(bool spacemod)
{
  space_gauss_modulated_=spacemod;
}

void Blob::init4Coord(std::vector<double> const &v) {
  size_t n = v.size();
  if (n!=4)
    GYOTO_ERROR("Initial coordinate should have 4 entries.");
  for (size_t i=0; i<n; ++i) {
    init4Coord_[i]=v[i];
  }
  if (blobMotionType_=="HelicalCylindrical" or blobMotionType_=="Equatorial"){
    if (init4Coord_[2] != M_PI/2.)
      GYOTO_ERROR("Equatorial and cylindrical motions must start at equator");
  }
}
std::vector<double> Blob::init4Coord() const {
  std::vector<double> v(4, 0.);
  for (size_t i=0; i<4; ++i) v[i]=init4Coord_[i];
  return v;
}

void Blob::init3Velo(std::vector<double> const &v) {
  size_t n = v.size();
  if (n!=3)
    GYOTO_ERROR("Initial velocity should have 3 entries.");
  for (size_t i=0; i<n; ++i) {
    init3Velo_[i]=v[i];
  }
}
std::vector<double> Blob::init3Velo() const {
  std::vector<double> v(3, 0.);
  for (size_t i=0; i<3; ++i) v[i]=init3Velo_[i];
  return v;
}

void Blob::blobMotionType(const string &kind) {
  if(kind == "Equatorial")
    blobMotionType_ = "Equatorial";
  else if (kind == "HelicalConical")
    blobMotionType_ = "HelicalConical";
  else if (kind == "HelicalCylindrical")
    blobMotionType_ = "HelicalCylindrical";
  else
    throwError("unknown blob motion type!");
}
string Blob::blobMotionType() const {
  return blobMotionType_;
}

void Blob::electronDistribution(const string &kind) {
  if(kind == "Thermal")
    electronDistrib_ = "Thermal";
  else if(kind == "Kappa")
    electronDistrib_ = "Kappa";
  else if (kind == "PL")
    electronDistrib_ = "PL";
  else
    throwError("unknown electron distribution!");
}
string Blob::electronDistribution() const {
  return electronDistrib_;
}

string Blob::className() const { return  string("Blob"); }
string Blob::className_l() const { return  string("blob"); }

double Blob::numberDensity() const {
  // Converts internal cgs central enthalpy to SI
  double dens=numberDensity_cgs_;
# ifdef HAVE_UDUNITS
  dens = Units::Converter("cm-3", "m-3")(dens);
# else
  GYOTO_WARNING << "Units ignored, please recompile Gyoto with --with-udunits"
		<< endl ;
# endif
  return dens; }
double Blob::numberDensity(string const &unit) const
{
  double dens = numberDensity();
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
void Blob::numberDensity(double dens) {
# ifdef HAVE_UDUNITS
  dens = Units::Converter("m-3", "cm-3")(dens);
# else
  GYOTO_WARNING << "Units ignored, please recompile Gyoto with --with-udunits"
		<< endl ;
# endif
  numberDensity_cgs_=dens;
}
void Blob::numberDensity(double dens, string const &unit) {
  if (unit != "") {
# ifdef HAVE_UDUNITS
    dens = Units::Converter(unit, "m-3")(dens);
# else
    GYOTO_WARNING << "Units ignored, please recompile Gyoto with --with-udunits"
		  << endl ;
# endif
  }
  numberDensity(dens);
}

double Blob::temperature() const { return temperature_; }
void Blob::temperature(double tt) { temperature_ = tt; }

double Blob::timeRef() const {
  // Converts internal M-unit time to SI
  double tt=timeRef_M_;
# ifdef HAVE_UDUNITS
  if (gg_)
    tt = Units::ToSeconds(tt,"geometrical_time",gg_);
  else
    GYOTO_SEVERE << "Cannot convert to seconds as metric is not set!" << endl;
# else
  GYOTO_WARNING << "Units ignored, please recompile Gyoto with --with-udunits"
		<< endl ;
# endif
  return tt; }
double Blob::timeRef(string const &unit) const
{
  double tt = timeRef();
  if (unit != "") {
# ifdef HAVE_UDUNITS
    tt = Units::Converter("s", unit)(tt);
# else
    GYOTO_WARNING << "Units ignored, please recompile Gyoto with --with-udunits"
		  << endl ;
# endif
  }
  return tt;
}
void Blob::timeRef(double tt) {
# ifdef HAVE_UDUNITS
  tt = Units::ToGeometricalTime(tt, "s", gg_);
# else
  GYOTO_WARNING << "Units ignored, please recompile Gyoto with --with-udunits"
		<< endl ;
# endif
  timeRef_M_ = tt; }
void Blob::timeRef(double tt, string const &unit) {
  if (unit != "") {
# ifdef HAVE_UDUNITS
    if (gg_)
      tt = Units::ToSeconds(tt,unit,gg_);
  else
    GYOTO_SEVERE << "Cannot convert to seconds as metric is not set!" << endl;
# else
    GYOTO_WARNING << "Units ignored, please recompile Gyoto with --with-udunits"
		  << endl ;
# endif
  }
  timeRef(tt);
}

double Blob::timeSigma() const {
  // Converts internal M-unit time to SI
  double tt=timeSigma_M_;
# ifdef HAVE_UDUNITS
  if (gg_)
    tt = Units::ToSeconds(tt,"geometrical_time",gg_);
  else
    GYOTO_SEVERE << "Cannot convert to seconds as metric is not set!" << endl;
# else
  GYOTO_WARNING << "Units ignored, please recompile Gyoto with --with-udunits"
		<< endl ;
# endif
  return tt; }
double Blob::timeSigma(string const &unit) const
{
  double tt = timeSigma();
  if (unit != "") {
# ifdef HAVE_UDUNITS
    tt = Units::Converter("s", unit)(tt);
# else
    GYOTO_WARNING << "Units ignored, please recompile Gyoto with --with-udunits"
		  << endl ;
# endif
  }
  return tt;
}
void Blob::timeSigma(double tt) {
# ifdef HAVE_UDUNITS
  tt = Units::ToGeometricalTime(tt, "s", gg_);
# else
  GYOTO_WARNING << "Units ignored, please recompile Gyoto with --with-udunits"
		<< endl ;
# endif
  timeSigma_M_ = tt; }
void Blob::timeSigma(double tt, string const &unit) {
  if (unit != "") {
# ifdef HAVE_UDUNITS
    tt = Units::ToSeconds(tt,unit,gg_);
# else
    GYOTO_WARNING << "Units ignored, please recompile Gyoto with --with-udunits"
		  << endl ;
# endif
  }
  timeSigma(tt);
}

void Blob::magnetizationParameter(double rr) {
  magnetizationParameter_=rr;}
double Blob::magnetizationParameter()const{
  return magnetizationParameter_;}

double Blob::kappaIndex() const { return kappaIndex_; }
void Blob::kappaIndex(double ind) { kappaIndex_ = ind; }

//////////////////////////////////////////
// NON-POLARIZED RADIATIVEQ //////////////
//////////////////////////////////////////

void Blob::radiativeQ(double Inu[], // output
		      double Taunu[], // output
		      double const nu_ems[], size_t nbnu, // input
		      double dsem,
		      state_t const &coord_ph,
		      double const co[8]) const {

# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << endl;
# endif

  double rr, rcyl, theta, phi, xx, yy, zz=0.;
  switch (gg_->coordKind()) {
  case GYOTO_COORDKIND_SPHERICAL:
    rr = coord_ph[1];
    theta = coord_ph[2];
    phi = coord_ph[3];
    
    rcyl = rr*sin(theta);
    xx = rcyl*cos(phi);
    yy = rcyl*sin(phi);
    zz   = rr*cos(theta);
    break;
  case GYOTO_COORDKIND_CARTESIAN:
    xx = coord_ph[1];
    yy = coord_ph[2];
    zz   = coord_ph[3];
    break;
  default:
    GYOTO_ERROR("In Blob::radiativeQ(): Unknown coordinate system kind");
  }

  ////////////// COMPUTE MODULATIONS /////////////////////
  double time_modulation = 1., space_modulation = 1.;
  if (time_gauss_modulated_==true){ // TIME GAUSSIAN MODULATION
    double tcur=coord_ph[0];
    time_modulation = exp(-pow((tcur-timeRef_M_)/timeSigma_M_,2));
  }

  if (space_gauss_modulated_==true){ // SPACE GAUSSIAN MODULATION
    double coord_spot[4]={co[0]};
    const_cast<Blob*>(this)
      ->getCartesian(coord_spot, 1, coord_spot+1, coord_spot+2, coord_spot+3);
    //above: nasty trick to deal with constness of emission
    double xspot=coord_spot[1], yspot=coord_spot[2], zspot=coord_spot[3];
    
    double difx=(xx-xspot), dify=(yy-yspot), difz=(zz-zspot);
    double d2 = difx*difx+dify*dify+difz*difz; // square coord distance between photon and blob's center
    
    double blobsize=radius_/3.; // NB: radius_ contains the "blob extension" ie the maximum extension over which Photon::hit returns 1 and emission is computed. This is assumed to coincide with the 3sigma extension of the Gaussian from the blob's center. So the blob's radius stricto sensu is radius_/3.
    double ds2=blobsize*blobsize;
    
    space_modulation=exp(-d2/(2.*ds2));
  }
  double total_modulation = time_modulation*space_modulation;
  ////////////// END COMPUTE MODULATIONS /////////////////
  
  double temperature = total_modulation*temperature_,
    number_density = total_modulation*numberDensity_cgs_;
  //cout << "spot tcur, time_ref, time_sigma, time_modulation, number_density=" << tcur << " " << timeRef_M_ << " " << timeSigma_M_ << " " << time_modulation << " " << numberDensity_cgs_ << " " << temperature_ << " " << number_density << " " << temperature << " " << kappaIndex_ << " " << magnetizationParameter_ << endl;
  double thetae = GYOTO_BOLTZMANN_CGS*temperature
    /(GYOTO_ELECTRON_MASS_CGS*GYOTO_C2_CGS);
  
  double BB = sqrt(4.*M_PI*magnetizationParameter_
		   *GYOTO_PROTON_MASS_CGS * GYOTO_C_CGS * GYOTO_C_CGS
		   *number_density);
  double nu0 = GYOTO_ELEMENTARY_CHARGE_CGS*BB
    /(2.*M_PI*GYOTO_ELECTRON_MASS_CGS*GYOTO_C_CGS); // cyclotron freq
  
  // Defining jnus, anus
  double jnu[nbnu], anu[nbnu];
  
  for (size_t ii=0; ii<nbnu; ++ii){
    // Initializing to <0 value to create errors if not updated
    // [ exp(-anu*ds) will explode ]
    jnu[ii]=-1.;
    anu[ii]=-1.;
  }

    if (electronDistrib_=="Thermal"){
    // THERMAL SYNCHROTRON
    double besselK2 = bessk(2, 1./thetae);
    //cout << "In Blob: ne, temperature, BB, nu0, besselK2, theta_mag: " << number_density << " " << temperature << " " << BB << " " << nu0 << " " << besselK2 << " " << theta_mag << endl;
    spectrumThermalSynch_->temperature(temperature);
    spectrumThermalSynch_->numberdensityCGS(number_density);
    spectrumThermalSynch_->angle_averaged(1); 
    spectrumThermalSynch_->angle_B_pem(0.); // avg so we don't care
    //spectrumThermalSynch_->angle_B_pem(0.785); // TEST!!
    spectrumThermalSynch_->cyclotron_freq(nu0);
    spectrumThermalSynch_->besselK2(besselK2);
    //cout << "for anu jnu: " << coord_ph[1] << " " << zz << " " << temperature << " " << number_density << " " << nu0 << " " << thetae << " " << besselK2 << endl;
    spectrumThermalSynch_->radiativeQ(jnu, anu, nu_ems, nbnu);   
  }else if (electronDistrib_=="Kappa"){
    // KAPPA SYNCHRO
    //double hypergeom = Gyoto::hypergeom(kappaIndex_, 10.); // TEST
    double hypergeom = Gyoto::hypergeom(kappaIndex_, thetae);
    //cout << "In Blob: ne, temperature, BB, nu0, besselK2, theta_mag: " << number_density << " " << temperature << " " << BB << " " << nu0 << " " << hypergeom << " " << theta_mag << endl;
    spectrumKappaSynch_->kappaindex(kappaIndex_);
    spectrumKappaSynch_->numberdensityCGS(number_density);
    spectrumKappaSynch_->angle_averaged(1); // TEST!!
    //cout << "In Blob: sin theta (k,B) = " << sin(theta_mag) << endl;
    spectrumKappaSynch_->angle_B_pem(0.); // avg so we don't care
    //spectrumKappaSynch_->angle_B_pem(0.785); // TEST!!
    spectrumKappaSynch_->cyclotron_freq(nu0);
    spectrumKappaSynch_->thetae(thetae);
    spectrumKappaSynch_->hypergeometric(hypergeom);
    
    spectrumKappaSynch_->radiativeQ(jnu, anu, nu_ems, nbnu);   
  }else if (electronDistrib_ == "PL"){
    spectrumPLSynch_->numberdensityCGS(number_density);
    spectrumPLSynch_->angle_averaged(1);
    //cout << "In Blob: sin theta (k,B) = " << sin(theta_mag) << endl;
    spectrumPLSynch_->angle_B_pem(0.); // avg so we don't care
    //spectrumPLSynch_->angle_B_pem(0.785); // TEST!!
    spectrumPLSynch_->cyclotron_freq(nu0);
    spectrumPLSynch_->PLindex(kappaIndex_-1);
    
    spectrumPLSynch_->radiativeQ(jnu, anu, nu_ems, nbnu);   
  }else{
    GYOTO_ERROR("Unknown electron distribution");
  }
  

  // RETURNING TOTAL INTENSITY AND TRANSMISSION
  for (size_t ii=0; ii<nbnu; ++ii){
    double jnu_tot = jnu[ii],
      anu_tot = anu[ii];
    //cout << "At r,th= " << coord_ph[1] << " " << coord_ph[2] << endl;
    //cout << "in unif stuff: " << number_density << " " << nu0 << " " << thetae << " " << hypergeom << " " << jnu_tot << " " << anu_tot << " " << dsem << endl;

    // expm1 is a precise implementation of exp(x)-1
    double em1=std::expm1(-anu_tot * dsem * gg_->unitLength());
    Taunu[ii] = em1+1.;
    Inu[ii] = anu_tot == 0. ? jnu_tot * dsem * gg_->unitLength() :
      -jnu_tot / anu_tot * em1;
    
    if (Inu[ii]<0.)
      GYOTO_ERROR("In Blob::radiativeQ: Inu<0");
    if (Inu[ii]!=Inu[ii] or Taunu[ii]!=Taunu[ii])
      GYOTO_ERROR("In Blob::radiativeQ: Inu or Taunu is nan");
    if (Inu[ii]==Inu[ii]+1. or Taunu[ii]==Taunu[ii]+1.)
      GYOTO_ERROR("In Blob::radiativeQ: Inu or Taunu is infinite");
    
  }

}

//////////////////////////////////////////
// POLARIZED RADIATIVEQ //////////////////
//////////////////////////////////////////

void Blob::radiativeQ(double *Inu, double *Qnu, double *Unu,
           double *Vnu,
           Eigen::Matrix4d *Onu,
           double const *nuem , size_t nbnu,
           double dsem,
           state_t const &coord_ph,
           double const *co) const {

  // polarized radiativeQ

  //cout << "In radQ pos vel= " << co[0] << " " << co[1] << " " << co[2] << " " << co[3] << " -- " << co[5]/co[4] << " " << co[7]/co[4] << endl;
  //cout << "ds in Blob= " << dsem << endl;
  //cout << "in blob photon tgt in orthon frame= " << coord_ph[4] << " " << coord_ph[5] << " " << coord_ph[1]*coord_ph[6] << " " <<  coord_ph[1]*abs(sin(coord_ph[2]))*coord_ph[7] << endl;
  
  double rr, rcyl, theta, phi, xx, yy, zz=0.;
  switch (gg_->coordKind()) {
  case GYOTO_COORDKIND_SPHERICAL:
    rr = coord_ph[1];
    theta = coord_ph[2];
    phi = coord_ph[3];
    
    rcyl = rr*sin(theta);
    xx = rcyl*cos(phi);
    yy = rcyl*sin(phi);
    zz   = rr*cos(theta);
    break;
  case GYOTO_COORDKIND_CARTESIAN:
    rcyl = pow(coord_ph[1]*coord_ph[1]+coord_ph[2]*coord_ph[2], 0.5);
    rr = sqrt(coord_ph[1]*coord_ph[1]+coord_ph[2]*coord_ph[2]
        +coord_ph[3]*coord_ph[3]);
    theta   = acos(coord_ph[3]/rr);
    xx = coord_ph[1];
    yy = coord_ph[2];
    zz   = coord_ph[3];
    break;
  default:
    GYOTO_ERROR("In Blob::radiativeQ(): Unknown coordinate system kind");
  }

  //cout << "Blob rcyl= " << rcyl << endl;
  //cout << "in Blob z= " << blobMotionType_ << " " << zz << endl;

  double vel[4]; // 4-velocity of emitter (needed for mf definition below)
  for (int ii=0;ii<4;ii++){
    vel[ii]=co[ii+4];
  }
  //cout << "***blob velo in orthonorm frame= " << vel[0] << " " << vel[1] << " " << coord_ph[1]*vel[2] << " " << coord_ph[1]*abs(sin(coord_ph[2]))*vel[3] << endl;

  ////////////// COMPUTE MODULATIONS /////////////////////
  double time_modulation = 1., space_modulation = 1.;
  if (time_gauss_modulated_==true){ // TIME GAUSSIAN MODULATION
    double tcur=coord_ph[0];
    time_modulation = exp(-pow((tcur-timeRef_M_)/timeSigma_M_,2));
    //cout << "timecur ref mod= " << tcur << " " << timeRef_M_ << " " << timeSigma_M_ << " " << time_modulation << endl;
  }

  if (space_gauss_modulated_==true){ // SPACE GAUSSIAN MODULATION
    double coord_spot[4]={co[0]};
    const_cast<Blob*>(this)
      ->getCartesian(coord_spot, 1, coord_spot+1, coord_spot+2, coord_spot+3);
    //above: nasty trick to deal with constness of emission
    double xspot=coord_spot[1], yspot=coord_spot[2], zspot=coord_spot[3];

    //cout << "in blob center at t=" << co[0] << " is at xyz=" << xspot << " " << yspot << " " << zspot << endl;
    
    double difx=(xx-xspot), dify=(yy-yspot), difz=(zz-zspot);
    double d2 = difx*difx+dify*dify+difz*difz; // square coord distance between photon and blob's center
    //cout << "In radQ dist to blob center at tcur; rcur zcur= " << sqrt(d2) << " " << rr << " " << zz << endl;
    
    double blobsize=radius_/3.; // NB: radius_ contains the "blob extension" ie the maximum extension over which Photon::hit returns 1 and emission is computed. This is assumed to coincide with the 3sigma extension of the Gaussian from the blob's center. So the blob's radius stricto sensu is radius_/3.
    double ds2=blobsize*blobsize;
    
    space_modulation=exp(-d2/(2.*ds2));
  }
  double total_modulation = time_modulation*space_modulation;
  ////////////// END COMPUTE MODULATIONS /////////////////

  // Computing density, temperature, mf magnitude
  double number_density, temperature, BB, Theta;
  if (USE_IPOLE_FORMALISM==0){
    /*******************/
    /* GYOTO FORMALISM */
    /*******************/
    // Modulate density, temperature and mf with this Gaussian
    number_density=numberDensity_cgs_*total_modulation;
    temperature = temperature_*total_modulation;
    BB = sqrt(4.*M_PI*magnetizationParameter_
	      *GYOTO_PROTON_MASS_CGS * GYOTO_C_CGS * GYOTO_C_CGS
	      *number_density); // *total_modulation --> already in ne!; // equipartition mf
  

    
    Theta = GYOTO_BOLTZMANN_CGS*temperature
      /(GYOTO_ELECTRON_MASS_CGS*GYOTO_C2_CGS);
  }else{
    /*******************/
    /* IPOLE FORMALISM */
    /*******************/
    //number_density=6e6/3.*total_modulation;
    double n0=6e6, theta0=200, B0=100;
    number_density=n0*total_modulation;
    Theta = theta0*total_modulation; // bug here before, was theta!
    BB=B0*total_modulation;
    //BB=50*total_modulation;
    temperature=GYOTO_ELECTRON_MASS_CGS*GYOTO_C2_CGS*theta/GYOTO_BOLTZMANN_CGS;
  }

  double nu0 = GYOTO_ELEMENTARY_CHARGE_CGS*BB
    /(2.*M_PI*GYOTO_ELECTRON_MASS_CGS*GYOTO_C_CGS); // cyclotron freq
      
  /////////////////// CHOOSE BFIELD GEOMETRY ////////////////////
  // Note: Bfield is simply a normalized spacelike vector, we only
  // need its direction, the norm is computed independently.

  double B4vect[4]={0.,0.,0.,0.};

  if (USE_IPOLE_FORMALISM==0){
    /*********************/
    /* GYOTO B FORMALISM */
    /*********************/

    computeB4vect(B4vect, magneticConfig_, co, coord_ph);
    
  }else{
    /*********************/
    /* IPOLE B FORMALISM */
    /*********************/

    string kin = gg_->kind();
    if (kin != "KerrBL") GYOTO_ERROR("Blob should be in Kerr for ipole formalism!");
    double spin = static_cast<SmartPointer<Metric::KerrBL> >(gg_) -> spin();

    computeB4vect_ipole(B4vect, magneticConfig_, co, coord_ph, spin);

  }
  
  //cout << "***BL mf in orthonorm frame= "<< B4vect[0]  << " " << B4vect[1]  << " " << coord_ph[1]*B4vect[2]  << " " << coord_ph[1]*abs(sin(coord_ph[2]))*B4vect[3]  << endl;

  //cout << "***k photon in orthonorm frame= " << coord_ph[4] << " " << coord_ph[5] << " " << coord_ph[1]*coord_ph[6] << " " << coord_ph[1]*abs(sin(coord_ph[2]))*coord_ph[7] << endl;
  
  double norm=sqrt(gg_->ScalarProd(&coord_ph[0], B4vect, B4vect));
  //cout << "norm= " << norm << endl;
  
  if (USE_IPOLE_FORMALISM==0){
    // Gyoto B formalism gives a unit vector
    if (fabs(norm-1.)>GYOTO_DEFAULT_ABSTOL) GYOTO_ERROR("Bad mf normalization");
  }else{
    // Ipole B formalism does not
    //cout << "renormalize ipole" << endl;
    gg_->multiplyFourVect(B4vect,1./norm);
  }

  double Chi=getChi(B4vect, coord_ph, vel); // this is EVPA
  //cout << endl;
  //cout << "At r,phi,x,y,z= " << coord_ph[1] << " " << coord_ph[3] << " " << coord_ph[1]*sin(coord_ph[2])*cos(coord_ph[3]) << " " << coord_ph[1]*sin(coord_ph[2])*sin(coord_ph[3]) << " " << coord_ph[1]*cos(coord_ph[2]) << " ; Chi=" << Chi << ", tan 2chi= " << tan(2.*Chi) << endl;
  //cout << endl;

  double theta_mag = get_theta_mag(B4vect, coord_ph, vel);
  //cout << "in Blob thetaB=" << theta_mag << endl;
  //cout << "at tcur= " << coord_ph[0] << "M = " << coord_ph[0]*21.18/60. << "min" << endl;
  //cout << endl;
  
  Eigen::Matrix4d Omat;
  Omat << 1, 0, 0, 0,
          0, 1, 0, 0,
          0, 0, 1, 0,
          0, 0, 0, 1;

 // Defining emission, absoprtion and rotation coefficients for the transmission matrix
  double jInu[nbnu], jQnu[nbnu], jUnu[nbnu], jVnu[nbnu];
  double aInu[nbnu], aQnu[nbnu], aUnu[nbnu], aVnu[nbnu];
  double rotQnu[nbnu], rotUnu[nbnu], rotVnu[nbnu];
  
  for (size_t ii=0; ii<nbnu; ++ii){
    // Initialze them to -1 to create error if not updated
    jInu[ii]=-1.;
    jQnu[ii]=-1.;
    jUnu[ii]=-1.;
    jVnu[ii]=-1.;
    aInu[ii]=-1.;
    aQnu[ii]=-1.;
    aUnu[ii]=-1.;
    aVnu[ii]=-1.;
    rotQnu[ii]=-1.;
    rotUnu[ii]=-1.;
    rotVnu[ii]=-1.;
  }

  if (electronDistrib_=="Thermal"){
    // THERMAL SYNCHROTRON
    double besselK2 = bessk(2, 1./Theta);
    //cout << "In Blob: ne, temperature, BB, nu0, besselK2, theta_mag: " << number_density << " " << temperature << " " << BB << " " << nu0 << " " << besselK2 << " " << theta_mag << endl;
    spectrumThermalSynch_->temperature(temperature);
    spectrumThermalSynch_->numberdensityCGS(number_density);
    spectrumThermalSynch_->angle_averaged(0); //  no angle avg of course
    spectrumThermalSynch_->angle_B_pem(theta_mag);
    //spectrumThermalSynch_->angle_B_pem(0.785); // TEST!!
    spectrumThermalSynch_->cyclotron_freq(nu0);
    spectrumThermalSynch_->besselK2(besselK2);
    //cout << "for anu jnu: " << coord_ph[1] << " " << zz << " " << temperature << " " << number_density << " " << nu0 << " " << thetae << " " << besselK2 << endl;
    spectrumThermalSynch_->radiativeQ(jInu, jQnu, jUnu, jVnu,
				      aInu, aQnu, aUnu, aVnu,
				      rotQnu, rotUnu, rotVnu, nuem, nbnu);
  }else if (electronDistrib_=="Kappa"){
    // KAPPA SYNCHRO
    //double hypergeom = Gyoto::hypergeom(kappaIndex_, 10.); // TEST
    double hypergeom = Gyoto::hypergeom(kappaIndex_, Theta);
    //cout << "In Blob: ne, temperature, BB, nu0, besselK2, theta_mag: " << number_density << " " << temperature << " " << BB << " " << nu0 << " " << hypergeom << " " << theta_mag << endl;
    spectrumKappaSynch_->kappaindex(kappaIndex_);
    spectrumKappaSynch_->numberdensityCGS(number_density);
    spectrumKappaSynch_->angle_averaged(0); // TEST!!
    //cout << "In Blob: sin theta (k,B) = " << sin(theta_mag) << endl;
    spectrumKappaSynch_->angle_B_pem(theta_mag);
    //spectrumKappaSynch_->angle_B_pem(0.785); // TEST!!
    spectrumKappaSynch_->cyclotron_freq(nu0);
    spectrumKappaSynch_->thetae(Theta);
    spectrumKappaSynch_->hypergeometric(hypergeom);
    
    spectrumKappaSynch_->radiativeQ(jInu, jQnu, jUnu, jVnu,
				    aInu, aQnu, aUnu, aVnu,
				    rotQnu, rotUnu, rotVnu, nuem, nbnu);
  }else if (electronDistrib_ == "PL"){
    spectrumPLSynch_->numberdensityCGS(number_density);
    spectrumPLSynch_->angle_averaged(0);
    //cout << "In Blob: sin theta (k,B) = " << sin(theta_mag) << endl;
    spectrumPLSynch_->angle_B_pem(theta_mag);
    //spectrumPLSynch_->angle_B_pem(0.785); // TEST!!
    spectrumPLSynch_->cyclotron_freq(nu0);
    spectrumPLSynch_->PLindex(kappaIndex_-1);
    
    spectrumPLSynch_->radiativeQ(jInu, jQnu, jUnu, jVnu,
				 aInu, aQnu, aUnu, aVnu,
				 rotQnu, rotUnu, rotVnu, nuem, nbnu);   
  }else{
    GYOTO_ERROR("Unknown electron distribution");
  }

  // RETURNING TOTAL INTENSITY AND TRANSMISSION
  for (size_t ii=0; ii<nbnu; ++ii) {
    
    // if (jInu[ii] > 2.5e-17){
    //   cout << "In Blob: ne, temperature, BB, theta_mag, Gauss modul: " << number_density << " " << temperature << " " << BB << " " << theta_mag << " " << total_modulation << endl;
    //   cout << "In Blob nuem jnu= " << 1.36e14/nuem[ii] << " " << jInu[ii]  << endl;
    // }
    //cout << "In Blob nuem jnu= " << 1.36e14/nuem[ii] << " " << jInu[ii]  << endl;
    //cout << "In Blob: jInu, jQnu, jUnu, jVnu: " << jInu[ii] << ", " << jQnu[ii] << ", " << jUnu[ii] << ", " << jVnu[ii] << endl;
    //cout << "In Blob: aInu, aQnu, aUnu, aVnu: " << aInu[ii] << ", " << aQnu[ii] << ", " << aUnu[ii] << ", " << aVnu[ii] << endl;
    //cout << "In Blob: rQnu, rUnu, rVnu: " << rotQnu[ii] << ", " << rotUnu[ii] << ", " << rotVnu[ii] << endl;
    Eigen::Vector4d Jstokes=rotateJs(jInu[ii], jQnu[ii], jUnu[ii], jVnu[ii], Chi)*dsem*gg_->unitLength();
    //cout << Jstokes << endl;
    Omat = Omatrix(aInu[ii], aQnu[ii], aUnu[ii], aVnu[ii], rotQnu[ii], rotUnu[ii], rotVnu[ii], Chi, dsem);
    //cout << Omat << endl;
    // Computing the increment of the Stokes parameters. Equivalent to dInu=exp(-anu*dsem)*jnu*dsem in the non-polarised case.
    Eigen::Vector4d Stokes=Omat*Jstokes;
    //cout << Stokes << endl;
    Inu[ii] = Stokes(0);
    Qnu[ii] = Stokes(1);
    Unu[ii] = Stokes(2);
    //cout << "Q and U and U/Q= " << Qnu[ii] << " " << Unu[ii] << " " << Unu[ii]/Qnu[ii] << endl;
    //cout << endl;
    
    Vnu[ii] = Stokes(3);
    Onu[ii] = Omat;

    //cout << "In Blob: r,th,ph, Inu, Qnu, Unu, Vnu, dsem, LP: " << rr << " " << theta << " " << phi << " " << Inu[ii] << ", " << Qnu[ii] << ", " << Unu[ii] << ", " << Vnu[ii] << ", " << dsem << ", " << pow(Qnu[ii]*Qnu[ii]+Unu[ii]*Unu[ii],0.5)/Inu[ii] << endl;

    if (Inu[ii]<0.)
      GYOTO_ERROR("In Blob::radiativeQ(): Inu<0");
    if (Inu[ii]!=Inu[ii] or Onu[ii](0,0)!=Onu[ii](0,0))
      GYOTO_ERROR("In Blob::radiativeQ(): Inu or Taunu is nan");
    if (Inu[ii]==Inu[ii]+1. or Onu[ii](0,0)==Onu[ii](0,0)+1.)
      GYOTO_ERROR("In Blob::radiativeQ(): Inu or Taunu is infinite");
  }
}

void Blob::getCartesian(double const * const dates, size_t const n_dates,
			double * const x, double * const y, double * const z, 
			double * const xprime, double * const yprime, double * const zprime){

  // this yields the position of the center of the UnifSphere
  // at time t
  // fourveldt_ is the initial 3-velocity dxi/dt
  // vel is the 4-velocity dxnu/dtau

  if (n_dates!=1)
    GYOTO_ERROR("In Blob::getCartesian n_dates!=1");

  double tt=dates[0];
  
  double r, theta, phi; // spherical coordinates
  double vel[4];
  
  if (blobMotionType_=="HelicalConical") // Helical conical ejection
    {
      r = init4Coord_[1]+init3Velo_[0]*(tt-init4Coord_[0]); // careful: init4Coord_[1] is r, but init3Velo_[0] is dr/dt
      theta = init4Coord_[2];
      phi = init4Coord_[3] + init4Coord_[1]*init4Coord_[1]*init3Velo_[2]/init3Velo_[0]*(pow(init4Coord_[1],-1.)-pow(r,-1.)); // result of integral of vphi over time
      //cout << "t0, r0, vr= " << init4Coord_[0] << " " << init4Coord_[1] << " " <<  init3Velo_[0] << endl;
      //cout << "In getCart t, r, theta, phi = " << tt << ", " << r << ", " << theta << ", " << phi << endl;
      
    }
  else if (blobMotionType_=="HelicalCylindrical")
    {
      double tmp = init3Velo_[1] * (tt-init4Coord_[0]); // vth0/r0*t = dthetadt_0 * t (vth0 = r0 * dthetadt_0)
      r = init4Coord_[1]/2. * (exp(tmp) + exp(-tmp));
      theta = 2.*atan(exp(tmp));
      phi = init4Coord_[3] + init3Velo_[2]*(tt-init4Coord_[0]);
    }
  else if (blobMotionType_=="Equatorial") // Equatorial orbit
    {
      if (init4Coord_[2]!=M_PI/2.)
	cout << "Warning input theta value incompatible with 'Equatorial' motion. Theta fixed to pi/2." << endl;
      
      r = init4Coord_[1];
      theta = M_PI/2.;
      phi = init4Coord_[3] + init3Velo_[2]*(tt-init4Coord_[0]);
      
    }
  else{
    GYOTO_ERROR("Unrecognized type of motion.");
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

void Blob::getVelocity(double const pos[4], double vel[4]){

  if (!gg_)
    GYOTO_ERROR("In Blob::getVelocity Metric not set");

  vel[0] = 1.; // put u^t to 1 as initilization, it will be computed in normalizeFourVel below

  if (blobMotionType_=="Equatorial"){
    vel[1] = 0.;
    vel[2] = 0.;
    vel[3] = init3Velo_[2]; // dphi/dt here, will be multiplied by dt/dtau in normalizeFourVel below
  }else if (blobMotionType_=="HelicalConical"){
    vel[1] = init3Velo_[0];
    vel[2] = 0.;
    double r0 = init4Coord_[1], rr = pos[1];
    //cout << "in blob velo t, r0, r= " << pos[0] << " " << r0 << " " << rr << endl;
    vel[3] = init3Velo_[2] * r0*r0 / (rr*rr); // assuming conservation of Newtonian ang mom
    //vel[3] = init3Velo_[2]; // TEST!!!
  }else if (blobMotionType_=="HelicalCylindrical"){
    // See FV 2024 notes for details of eqs
    double tmp = init3Velo_[1] * (pos[0] - init4Coord_[0]); // vth0/r0*t = dthetadt_0 * t
    vel[1] = init4Coord_[1]*init3Velo_[1]/2. * (exp(tmp) - exp(-tmp)); // dr/dt
    // Careful: init4Coord_[1] is r ; init3Velo_[1] is dtheta/dt
    vel[2] = init4Coord_[1]/pos[1]*init3Velo_[1]; // vth = r*dth/dt = cst -> dth/dt = r_0/r * dth/dt_0
    vel[3] = init3Velo_[2]; // vph = r*sinth*Omega = cst and Omega = cst
    //cout << "In Blob pos, velos, Velo square= " << pos[0] << " " << pos[1] << " " << pos[2] << " " << pos[3] << " " << pos[1]*cos(pos[2]) << " " << vel[1] << " " << pos[1]*vel[2] << " " <<  pos[1]*sin(pos[2])*vel[3] << " " << init4Coord_[1]*sin(init4Coord_[2])*init3Velo_[2] << " " << vel[1]*vel[1] + pos[1]*pos[1]*vel[2]*vel[2] + pos[1]*pos[1]*sin(pos[2])*sin(pos[2])*vel[3]*vel[3] << endl;
  }

  //cout << "At t, r, theta= " << pos[0] << " " << pos[1] << " " << pos[2] << ", blob 3vel= " << vel[1] << " " << vel[3] << endl;
  
  gg_->normalizeFourVel(pos, vel);
  //cout << "in blob 4u= " << vel[0] << " " << vel[1] << " " << vel[2] << " " << vel[3] << endl;

}
