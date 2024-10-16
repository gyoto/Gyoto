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
GYOTO_PROPERTY_END(Blob, Star::properties)

#define USE_IPOLE_FORMALISM 0

Blob::Blob() :
  Star(),
  time_gauss_modulated_(false),
  space_gauss_modulated_(false),
  numberDensity_cgs_(1.),
  temperature_(1.),
  timeRef_M_(1.),
  timeSigma_M_(1.),
  magnetizationParameter_(1.),
  kappaIndex_(1.),
  magneticConfig_("None"),
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
}

Blob::Blob(const Blob& orig) :
  Star(orig),
  time_gauss_modulated_(orig.time_gauss_modulated_),
  space_gauss_modulated_(orig.space_gauss_modulated_),
  numberDensity_cgs_(orig.numberDensity_cgs_),
  temperature_(orig.temperature_),
  timeRef_M_(orig.timeRef_M_),
  timeSigma_M_(orig.timeSigma_M_),
  kappaIndex_(orig.kappaIndex_),
  magnetizationParameter_(orig.magnetizationParameter_),
  magneticConfig_(orig.magneticConfig_),
  spectrumKappaSynch_(NULL),
  spectrumPLSynch_(NULL),
  spectrumThermalSynch_(NULL),
  electronDistrib_(orig.electronDistrib_)
{
  if (orig.spectrumKappaSynch_()) spectrumKappaSynch_=orig.spectrumKappaSynch_->clone();
  if (orig.spectrumPLSynch_()) spectrumPLSynch_=orig.spectrumPLSynch_->clone();
  if (orig.spectrumThermalSynch_()) spectrumThermalSynch_=orig.spectrumThermalSynch_->clone();
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

void Blob::magneticConfiguration(string config){
  magneticConfig_=config;
}

string Blob::magneticConfiguration() const{
  return magneticConfig_;
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

  string kin = gg_->kind();
  if (kin != "KerrBL" and kin != "Minkowski") GYOTO_ERROR("Blob should be in Kerr or Minko!");
  double spin = 0.;
  if (kin == "KerrBL"){
    // Check that Kerr spin is 0 (mf formulas below are so far in Sch only)
    spin = static_cast<SmartPointer<Metric::KerrBL> >(gg_) -> spin();
    if (spin!=0.) GYOTO_ERROR("Blob should be in Schwarzschild!");
  }

  // polarized radiativeQ
  
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
    
    double difx=(xx-xspot), dify=(yy-yspot), difz=(zz-zspot);
    double d2 = difx*difx+dify*dify+difz*difz; // square coord distance between photon and blob's center
    
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
    
    double gtt = gg_->gmunu(&coord_ph[0],0,0),
      grr = gg_->gmunu(&coord_ph[0],1,1),
      gthth = gg_->gmunu(&coord_ph[0],2,2),
      gtp = gg_->gmunu(&coord_ph[0],0,3),
      gpp = gg_->gmunu(&coord_ph[0],3,3);
    
    // SCHWARZSCHILD EXPRESSION SO FAR!!!
    if (magneticConfig_=="Vertical"){
      double Br = cos(coord_ph[2])/sqrt(grr),
	Bth = -sin(coord_ph[2])/sqrt(gthth); // along +ez
      
      B4vect[1]=Br;
      B4vect[2]=Bth;
    }else if(magneticConfig_=="Radial"){
      double Br = 1./sqrt(grr); // along +er
      
      B4vect[1]=Br;
    }else if(magneticConfig_=="Toroidal"){
      if (vel[0]==0.) GYOTO_ERROR("Undefined 4-velocity for toroidal mf");
      double omega=vel[3]/vel[0], omega2 = omega*omega;
      double Bt2 = gpp/gtt*omega2/(gtt+gpp*omega2),
	Bp2 = gtt/gpp*1./(gtt+gpp*omega2);
      //cout << "Btor stuff: " << omega2 << " " << Bt2 << " " << Bt2/Bp2 << endl;
      if (Bt2<0. or Bp2<0.) GYOTO_ERROR("Bad configuration for toroidal mf");
      B4vect[0]=sqrt(Bt2);
      B4vect[3]=sqrt(Bp2);
    }else{
      GYOTO_ERROR("Not implemented Bfield orientation");
    }
  }else{
    /*********************/
    /* IPOLE B FORMALISM */
    /*********************/
    
    double B_1=0.,B_2=0.,B_3=0;
    //double spin = 0.; // ONLY VALID FOR SCH
    double gtt = gg_->gmunu(&coord_ph[0],0,0),
      grr = gg_->gmunu(&coord_ph[0],1,1),
      gthth = gg_->gmunu(&coord_ph[0],2,2),
      gpp = gg_->gmunu(&coord_ph[0],3,3);
    double dx1=0.025,
      dx2=0.025;
    
    if (magneticConfig_=="None")
      GYOTO_ERROR("Specify the magnetic field configuration");
    if (magneticConfig_=="Vertical"){
      double g_det = sqrt(M_PI*M_PI*pow(rr,6)*pow(sin(theta),2));
      
      double F11 = exp(log(rr)-dx1)*sin(theta-dx2*M_PI),
	F12 = exp(log(rr)-dx1)*sin(theta+dx2*M_PI),
	F21 = exp(log(rr)+dx1)*sin(theta-dx2*M_PI),
	F22 = exp(log(rr)+dx1)*sin(theta+dx2*M_PI);
      B_1 = -(F11-F12+F21-F22)/(2.*dx2*g_det);
      B_2 =  (F11+F12-F21-F22)/(2.*dx1*g_det);
      B_3 = 0.;
    }
    else if (magneticConfig_=="Radial"){
      double g_det = sqrt(M_PI*M_PI*pow(rr,6)*pow(sin(theta),2));
      double F11 = 1.-cos(theta-dx2*M_PI),
	F12 = 1.-cos(theta+dx2*M_PI),
	F21 = 1.-cos(theta-dx2*M_PI),
	F22 = 1.-cos(theta+dx2*M_PI);
      B_1 = -(F11-F12+F21-F22)/(2.*dx2*g_det),
	B_2 =  (F11+F12-F21-F22)/(2.*dx1*g_det),
	B_3 = 0.;
    }
    else if (magneticConfig_=="Toroidal"){
      /*double gtt=gg_->gmunu(&coord_ph[0],0,0),
	gpp=gg_->gmunu(&coord_ph[0],3,3),
	Bp=1.,
	Bt=(-gpp*Bp*vel[3])/(gtt*vel[0]);
	B4vect[0]=Bt;
	B4vect[3]=Bp;*/
      B_1 = 0.;
      B_2 = 0.;
      B_3 = 1.;
    }
    else
      GYOTO_ERROR("Unknown magnetic field configuration");

    // compute contravariant velocity in KS' from BL
    double dtKS_drBL   = 2. * rr / (rr*rr - 2.*rr + spin*spin);
    double dphiKS_drBL = spin / (rr*rr - 2.*rr + spin*spin);
    double Ucon_KSm[4]={0.,0.,0.,0.};
    Ucon_KSm[0]=vel[0]+vel[1]*dtKS_drBL;
    Ucon_KSm[1]=vel[1]/rr;
    Ucon_KSm[2]=vel[2]/M_PI;
    Ucon_KSm[3]=vel[3]+vel[1]*dphiKS_drBL;
    
    // Compute KS' metric
    double gcov_ksm[4][4];
    double sin2=sin(theta)*sin(theta), rho2=rr*rr+spin*spin*cos(theta)*cos(theta);
    double gcov_ks[4][4];
    for(int mu=0;mu<4;mu++)
      for(int nu=0;nu<4;nu++)
	gcov_ks[mu][nu]=0.;
    
    gcov_ks[0][0] = -1. + 2. * rr / rho2 ;
    gcov_ks[0][1] = 2. * rr / rho2 ;
    gcov_ks[0][3] = -2. * spin * rr * sin(theta)*sin(theta) / rho2;
    gcov_ks[1][0] = gcov_ks[0][1];
    gcov_ks[1][1] = 1. + 2. * rr / rho2 ;
    gcov_ks[1][3] = -spin * sin(theta)*sin(theta) * (1. + 2. * rr / rho2);
    gcov_ks[2][2] = rho2 ;
    gcov_ks[3][0] = gcov_ks[0][3];
    gcov_ks[3][1] = gcov_ks[1][3];
    gcov_ks[3][3] = sin(theta)*sin(theta) * (rho2 + spin * spin * sin(theta)*sin(theta) * (1. + 2. * rr / rho2));
    
    // convert from ks metric to a modified one using Jacobian
    double dxdX[4][4];
    double hslope=0.;
    for(int mu=0;mu<4;mu++)
      for(int nu=0;nu<4;nu++)
	dxdX[mu][nu]=0.;
    
    dxdX[0][0] = 1.;
    dxdX[1][1] = rr;
    dxdX[2][2] = M_PI  + hslope*2*M_PI*cos(2*theta); 
    dxdX[3][3] = 1.;
    
    for(int mu=0;mu<4;mu++){
      for(int nu=0;nu<4;nu++){
	gcov_ksm[mu][nu] = 0;
	for (int lam = 0; lam < 4; ++lam) {
	  for (int kap = 0; kap < 4; ++kap) {
	    gcov_ksm[mu][nu] += gcov_ks[lam][kap] * dxdX[lam][mu] * dxdX[kap][nu];
	  }
	}
      }
    }
    
    // Compute covariant velocity in KS'
    double Ucov_KSm[4]={0.,0.,0.,0.};
    for(int mu=0;mu<4;mu++){
      for(int nu=0;nu<4;nu++){
	Ucov_KSm[mu] += gcov_ksm[mu][nu]*Ucon_KSm[nu];
      }
    }
    
    // Copute Magnetic field in KS'
    //cout << "r sth, velBL, ucov KS'= " << co[1] << " " << sin(co[2]) << " " << vel[0] << " " << vel[3] << " " << Ucov_KSm[1] << " " << Ucov_KSm[2] << " " << Ucov_KSm[3] << endl;
    //throwError("test disk");
    double B0=B_1*Ucov_KSm[1]+B_2*Ucov_KSm[2]+B_3*Ucov_KSm[3],
      B1=(B_1+B0*Ucon_KSm[1])/Ucon_KSm[0],
      B2=(B_2+B0*Ucon_KSm[2])/Ucon_KSm[0],
      B3=(B_3+B0*Ucon_KSm[3])/Ucon_KSm[0];
    
    // Conversion Magnetic field from KS' -> BL
    double Delta = pow(rr,2)-2.*rr+pow(spin,2.);
    B4vect[0]=B0-B1*2.*pow(rr,2)/Delta;
    B4vect[1]=B1*rr;
    B4vect[2]=B2*M_PI;
    B4vect[3]=B3-B1*spin*rr/Delta;
    // end of ipole Bfield formalism    
  }

  /////////////////// END CHOOSE BFIELD GEOMETRY ////////////////////
    
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

  // Computing the angle theta_mag between the magnetic field vector and photon tgt vector in the rest frame of the emitter
  gg_->projectFourVect(&coord_ph[0],B4vect,vel); //Projection of the 4-vector B to 4-velocity to be in the rest frame of the emitter
  double photon_emframe[4]; // photon tgt vector projected in comoving frame
  for (int ii=0;ii<4;ii++){
    photon_emframe[ii]=coord_ph[ii+4];
  }

  // Angle between mf and K in emitter's frame. B is already in this frame,
  // it is by construction normal to u. We still need to project k
  // normal to u, this is K = k + (k.u) u. Then:
  // cos(thetaB) = (K / |K|) . (B / |B|)
  gg_->projectFourVect(&coord_ph[0],photon_emframe,vel);
  //cout << "***K photon proj in orthonorm frame= " << photon_emframe[0] << " " << photon_emframe[1] << " " << coord_ph[1]*photon_emframe[2] << " " << coord_ph[1]*abs(sin(coord_ph[2]))*photon_emframe[3] << endl;
  //cout << endl;
  double bnorm = gg_->norm(&coord_ph[0],B4vect);
  double lnorm = gg_->norm(&coord_ph[0],photon_emframe);
  double lscalb = gg_->ScalarProd(&coord_ph[0],photon_emframe,B4vect);
  double theta_mag = acos(lscalb/(lnorm*bnorm)); // acos is in 0,pi, which is appropriate for theta_mag.
  //cout << "thetaB=" << theta_mag << endl;
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

