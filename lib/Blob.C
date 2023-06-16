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
GYOTO_PROPERTY_END(Blob, Star::properties)

Blob::Blob() :
  Star(),
  numberDensity_cgs_(1.),
  temperature_(1.),
  timeRef_M_(1.),
  timeSigma_M_(1.),
  magnetizationParameter_(1.),
  kappaIndex_(1.),
  magneticConfig_("None"),
  spectrumKappaSynch_(NULL),
  spectrumThermalSynch_(NULL)
{
  kind_="Blob";
# ifdef GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << "done." << endl;
# endif
  spectrumKappaSynch_ = new Spectrum::KappaDistributionSynchrotron();
  spectrumThermalSynch_ = new Spectrum::ThermalSynchrotron();
}

Blob::Blob(const Blob& orig) :
  Star(orig),
  numberDensity_cgs_(orig.numberDensity_cgs_),
  temperature_(orig.temperature_),
  timeRef_M_(orig.timeRef_M_),
  timeSigma_M_(orig.timeSigma_M_),
  kappaIndex_(orig.kappaIndex_),
  magnetizationParameter_(orig.magnetizationParameter_),
  magneticConfig_(orig.magneticConfig_),
  spectrumKappaSynch_(NULL),
  spectrumThermalSynch_(NULL)
{
  if (orig.spectrumKappaSynch_()) spectrumKappaSynch_=orig.spectrumKappaSynch_->clone();
  if (orig.spectrumThermalSynch_()) spectrumThermalSynch_=orig.spectrumThermalSynch_->clone();
}

Blob* Blob::clone() const { return new Blob(*this); }

Blob::~Blob() {
  if (debug()) cerr << "DEBUG: Blob::~Blob()\n";
}

string Blob::className() const { return  string("Blob"); }
string Blob::className_l() const { return  string("blob"); }

void Blob::radiativeQ(double Inu[], // output
		      double Taunu[], // output
		      double const nu_ems[], size_t nbnu, // input
		      double dsem,
		      state_t const &coord_ph,
		      double const coord_obj[8]) const {

# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << endl;
# endif
  double tcur=coord_ph[0]; //*GMoc3/60.; // in min
  double modulation = 1; //exp(-pow((tcur-timeRef_M_)/timeSigma_M_,2));
  double temperature = modulation*temperature_,
    number_density = modulation*numberDensity_cgs_;
  //cout << "spot tcur, time_ref, time_sigma, modulation, number_density=" << tcur << " " << timeRef_M_ << " " << timeSigma_M_ << " " << modulation << " " << numberDensity_cgs_ << " " << temperature_ << " " << number_density << " " << temperature << " " << kappaIndex_ << " " << magnetizationParameter_ << endl;
  double thetae = GYOTO_BOLTZMANN_CGS*temperature
    /(GYOTO_ELECTRON_MASS_CGS*GYOTO_C2_CGS);
  
  double hypergeom = Gyoto::hypergeom(kappaIndex_, thetae);
  
  double BB = sqrt(4.*M_PI*magnetizationParameter_
		   *GYOTO_PROTON_MASS_CGS * GYOTO_C_CGS * GYOTO_C_CGS
		   *number_density);
  double nu0 = GYOTO_ELEMENTARY_CHARGE_CGS*BB
    /(2.*M_PI*GYOTO_ELECTRON_MASS_CGS*GYOTO_C_CGS); // cyclotron freq
  
  // Defining jnus, anus
  double jnu_synch_kappa[nbnu], anu_synch_kappa[nbnu];
  
  for (size_t ii=0; ii<nbnu; ++ii){
    // Initializing to <0 value to create errors if not updated
    // [ exp(-anu*ds) will explose ]
    jnu_synch_kappa[ii]=-1.;
    anu_synch_kappa[ii]=-1.;
  }
  
  // THERMAL SYNCHRO
  spectrumKappaSynch_->kappaindex(kappaIndex_);
  spectrumKappaSynch_->numberdensityCGS(number_density);
  spectrumKappaSynch_->angle_averaged(1);
  spectrumKappaSynch_->angle_B_pem(0.); // avg so we don't care
  spectrumKappaSynch_->cyclotron_freq(nu0);
  spectrumKappaSynch_->thetae(thetae);
  spectrumKappaSynch_->hypergeometric(hypergeom);

  spectrumKappaSynch_->radiativeQ(jnu_synch_kappa,anu_synch_kappa,
				  nu_ems,nbnu);

  // RETURNING TOTAL INTENSITY AND TRANSMISSION
  for (size_t ii=0; ii<nbnu; ++ii){
    double jnu_tot = 1., // TEST jnu_synch_kappa[ii],
      anu_tot = 0.; // TEST anu_synch_kappa[ii];
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

void Blob::radiativeQ(double *Inu, double *Qnu, double *Unu,
           double *Vnu,
           Eigen::Matrix4d *Onu,
           double const *nuem , size_t nbnu,
           double dsem,
           state_t const &coord_ph,
           double const *co) const {
  
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

  double vel[4]; // 4-velocity of emitter
  for (int ii=0;ii<4;ii++){
    vel[ii]=co[ii+4];
  }
  //cout << "blob velo= " << vel[0] << " " << vel[1] << " " << vel[2] << " " << vel[3] << endl;

  double coord_spot[4]={co[0]};
  const_cast<Blob*>(this)
    ->getCartesian(coord_spot, 1, coord_spot+1, coord_spot+2, coord_spot+3);
  //above: nasty trick to deal with constness of emission
  double xspot=coord_spot[1], yspot=coord_spot[2], zspot=coord_spot[3];

  double difx=(xx-xspot), dify=(yy-yspot), difz=(zz-zspot);
  double d2 = difx*difx+dify*dify+difz*difz; // square coord distance between photon and blob's center
  
  double blobsize=radius_/3.; // NB: radius_ contains the "blob extension" ie the maximum extension over which Photon::hit returns 1 and emission is computed. This is assumed to coincide with the 3sigma extension of the Gaussian from the blob's center. So the blob's radius stricto sensu is radius_/3.
  double ds2=blobsize*blobsize;

  double expo_fact=exp(-d2/(2.*ds2)); // Gaussian modulation with sigma corresponding to radius_/3

  // Modulate density, temperature and mf with this Gaussian
  double number_density=numberDensity_cgs_*expo_fact;
  double temperature = temperature_*expo_fact;
  double BB = sqrt(4.*M_PI*magnetizationParameter_
		   *GYOTO_PROTON_MASS_CGS * GYOTO_C_CGS * GYOTO_C_CGS
		   *number_density)*expo_fact; // equipartition mf
  
  double nu0 = GYOTO_ELEMENTARY_CHARGE_CGS*BB
    /(2.*M_PI*GYOTO_ELECTRON_MASS_CGS*GYOTO_C_CGS); // cyclotron freq

  double Theta = GYOTO_BOLTZMANN_CGS*temperature
    /(GYOTO_ELECTRON_MASS_CGS*GYOTO_C2_CGS);

  double hypergeom = Gyoto::hypergeom(kappaIndex_, Theta);
  
  // CHOOSE BFIELD GEOMETRY
  // Note: Bfield is simply a normalized spacelike vector, we only
  // need its direction, the norm is computed independently.

  double B4vect[4]={0.,0.,0.,0.};
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
    if (Bt2<0. or Bp2<0.) GYOTO_ERROR("Bad configuration for toroidal mf");
    B4vect[0]=sqrt(Bt2);
    B4vect[3]=sqrt(Bp2);
  }else{
    GYOTO_ERROR("Not implemented Bfield orientation");
  }
    
  //cout << "***BL mf in orthonorm frame= "<< B4vect[0]  << " " << B4vect[1]  << " " << coord_ph[1]*B4vect[2]  << " " << coord_ph[1]*abs(sin(coord_ph[2]))*B4vect[3]  << endl;

  //cout << "***k photon in orthonorm frame= " << coord_ph[4] << " " << coord_ph[5] << " " << coord_ph[1]*coord_ph[6] << " " << coord_ph[1]*abs(sin(coord_ph[2]))*coord_ph[7] << endl;
  
  double norm=sqrt(gg_->ScalarProd(&coord_ph[0], B4vect, B4vect));
  //cout << "norm= " << norm << endl;
  gg_->multiplyFourVect(B4vect,1./norm);

  double Chi=getChi(B4vect, coord_ph, vel); // this is EVPA
  //cout << "At r,x,y,z= " << coord_ph[1] << " " << coord_ph[1]*sin(coord_ph[2])*cos(coord_ph[3]) << " " << coord_ph[1]*sin(coord_ph[2])*sin(coord_ph[3]) << " " << coord_ph[1]*cos(coord_ph[2]) << " ; Chi=" << Chi << endl;

  // Computing the angle theta_mag between the magnetic field vector and photon tgt vector in the rest frame of the emitter
  gg_->projectFourVect(&coord_ph[0],B4vect,vel); //Projection of the 4-vector B to 4-velocity to be in the rest frame of the emitter
  double photon_emframe[4]; // photon tgt vector projected in comoving frame
  for (int ii=0;ii<4;ii++){
    photon_emframe[ii]=coord_ph[ii+4];
  }
  gg_->projectFourVect(&coord_ph[0],photon_emframe,vel);
  double bnorm = gg_->norm(&coord_ph[0],B4vect);
  double lnorm = gg_->norm(&coord_ph[0],photon_emframe);
  double lscalb = gg_->ScalarProd(&coord_ph[0],photon_emframe,B4vect);
  double theta_mag = acos(lscalb/(lnorm*bnorm));
  
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

  int use_kappa_synch=1; // put 0 to switch to thermal synch
  if (use_kappa_synch==0){
    // THERMAL SYNCHROTRON
    double besselK2 = bessk(2, 1./Theta);
    //cout << "In Blob: ne, temperature, BB, nu0, besselK2, theta_mag: " << number_density << " " << temperature << " " << BB << " " << nu0 << " " << besselK2 << " " << theta_mag << endl;
    spectrumThermalSynch_->temperature(temperature);
    spectrumThermalSynch_->numberdensityCGS(number_density);
    spectrumThermalSynch_->angle_averaged(0); //  no angle avg of course
    spectrumThermalSynch_->angle_B_pem(theta_mag);   
    spectrumThermalSynch_->cyclotron_freq(nu0);
    spectrumThermalSynch_->besselK2(besselK2);
    //cout << "for anu jnu: " << coord_ph[1] << " " << zz << " " << temperature << " " << number_density << " " << nu0 << " " << thetae << " " << besselK2 << endl;
    spectrumThermalSynch_->radiativeQ(jInu, jQnu, jUnu, jVnu,
				      aInu, aQnu, aUnu, aVnu,
				      rotQnu, rotUnu, rotVnu, nuem, nbnu);
  }else{
    // KAPPA SYNCHRO
    spectrumKappaSynch_->kappaindex(kappaIndex_);
    spectrumKappaSynch_->numberdensityCGS(number_density);
    spectrumKappaSynch_->angle_averaged(0);
    //cout << "In Blob: sin theta (k,B) = " << sin(theta_mag) << endl;
    spectrumKappaSynch_->angle_B_pem(theta_mag);
    //spectrumKappaSynch_->angle_B_pem(0.5); // TEST!!
    spectrumKappaSynch_->cyclotron_freq(nu0);
    spectrumKappaSynch_->thetae(Theta);
    spectrumKappaSynch_->hypergeometric(hypergeom);
    
    spectrumKappaSynch_->radiativeQ(jInu, jQnu, jUnu, jVnu,
				    aInu, aQnu, aUnu, aVnu,
				    rotQnu, rotUnu, rotVnu, nuem, nbnu);
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


void Blob::magneticConfiguration(string config){
  magneticConfig_=config;
}

string Blob::magneticConfiguration() const{
  return magneticConfig_;
}
