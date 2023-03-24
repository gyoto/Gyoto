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
  double modulation = exp(-pow((tcur-timeRef_M_)/timeSigma_M_,2));
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
    double jnu_tot = jnu_synch_kappa[ii],
      anu_tot = anu_synch_kappa[ii];
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
           state_t const &cph,
           double const *co) const {
  // polarized radiativeQ
  double rr, theta, phi, xx, yy, zz=0.;
  switch (gg_->coordKind()) {
  case GYOTO_COORDKIND_SPHERICAL:
    rr = cph[1];
    theta = cph[2];
    phi = cph[3];
    xx = rr*sin(theta)*cos(phi);
    yy = rr*sin(theta)*sin(phi);
    zz = rr*cos(theta);
    break;
  case GYOTO_COORDKIND_CARTESIAN:
    rr = sqrt(cph[1]*cph[1]+cph[2]*cph[2]+cph[3]*cph[3]);
    theta = acos(cph[3]/rr);
    phi = atan(cph[2]/cph[1]);
    xx = cph[1];
    yy = cph[2];
    zz = cph[3];
    break;
  default:
    GYOTO_ERROR("In ThickDisk::radiativeQ(): Unknown coordinate system kind");
  }

  double vel[4], pos[4]; // 4-velocity of emitter
  for (int ii=0;ii<4;ii++){
    vel[ii]=co[ii+4];
    pos[ii]=co[ii];
  }
  //cout << vel[0] << " " << vel[1] << " " << vel[2] << " " << vel[3] << endl;
  
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

  // CHOOSE BFIELD GEOMETRY
  double B4vect[4]={0.,0.,0.,0.};
  double B0=0.,B1=0.,B2=0.,B3=0;
  double spin = 0.;
  if (magneticConfig_=="None")
    GYOTO_ERROR("Specify the magnetic field configuration");
  if (magneticConfig_=="Vertical"){
    double dx1=0.025,
      dx2=0.025;
    double gtt=gg_->gmunu(&cph[0],0,0),
      grr = gg_->gmunu(&cph[0],1,1),
      gthth = gg_->gmunu(&cph[0],2,2),
      gpp=gg_->gmunu(&cph[0],3,3);
    double g_det = gtt+grr+gthth+gpp;
    double F11 = exp(log(rr)-dx1)*sin(theta-dx2*M_PI),
           F12 = exp(log(rr)-dx1)*sin(theta+dx2*M_PI),
           F21 = exp(log(rr)+dx1)*sin(theta-dx2*M_PI),
           F22 = exp(log(rr)+dx1)*sin(theta+dx2*M_PI);
    double B_1 = -(F11-F12+F21-F22)/(2.*dx2*g_det),
           B_2 =  (F11+F12-F21-F22)/(2.*dx1*g_det),
           B_3 = 0.;

    
    B0=B_1*grr*vel[1]+B_2*gthth*vel[2]+B_3*gpp*vel[3];
    B1=(B_1+B0*vel[1])/vel[0];
    B2=(B_2+B0*vel[2])/vel[0];
    B3=(B_3+B0*vel[3])/vel[0];

    double Delta = pow(rr,2)-2.*rr+pow(spin,2.);

    B4vect[0]=B0-B1*2.*pow(rr,2)/Delta;
    B4vect[1]=B1*rr;
    B4vect[2]=B2*M_PI;
    B4vect[3]=B3-B1*spin*rr/Delta;
  }
  else if (magneticConfig_=="Radial"){
    double dx1=0.025,
      dx2=0.025;
    double gtt=gg_->gmunu(&cph[0],0,0),
      grr = gg_->gmunu(&cph[0],1,1),
      gthth = gg_->gmunu(&cph[0],2,2),
      gpp=gg_->gmunu(&cph[0],3,3);
    double g_det = gtt+grr+gthth+gpp;
    double F11 = 1.-cos(theta-dx2*M_PI),
           F12 = 1.-cos(theta+dx2*M_PI),
           F21 = 1.-cos(theta-dx2*M_PI),
           F22 = 1.-cos(theta+dx2*M_PI);
    double B_1 = -(F11-F12+F21-F22)/(2.*dx2*g_det),
           B_2 =  (F11+F12-F21-F22)/(2.*dx1*g_det),
           B_3 = 0.;

    B0=B_1*grr*vel[1]+B_2*gthth*vel[2]+B_3*gpp*vel[3];
    B1=(B_1+B0*vel[1])/vel[0];
    B2=(B_2+B0*vel[2])/vel[0];
    B3=(B_3+B0*vel[3])/vel[0];

    double Delta = pow(rr,2)-2.*rr+pow(spin,2.);

    B4vect[0]=B0-B1*2.*pow(rr,2)/Delta;
    B4vect[1]=B1*rr;
    B4vect[2]=B2*M_PI;
    B4vect[3]=B3-B1*spin*rr/Delta;
  }
  else if (magneticConfig_=="Toroidal"){
    double gtt=gg_->gmunu(&cph[0],0,0),
      gpp=gg_->gmunu(&cph[0],3,3),
      Bp=1.,
      Bt=(-gpp*Bp*vel[3])/(gtt*vel[0]);
    B4vect[0]=Bt;
    B4vect[3]=Bp;
  }
  else
    GYOTO_ERROR("Unknown magnetic field configuration");
  
  double BB0 = 100.; // Gauss
  double norm=gg_->norm(&cph[0], B4vect);
  gg_->multiplyFourVect(B4vect,1./norm);
  //cout << "B4vect: " << B4vect[0] << " " << B4vect[1] << " " << B4vect[2] << " " << B4vect[3] << endl;


  double Chi=getChi(B4vect, cph, vel); // this is EVPA

  // Computing the angle theta_mag between the magnetic field vector and photon tgt vector in the rest frame of the emitter
  gg_->projectFourVect(&cph[0],B4vect,vel); //Projection of the 4-vector B to 4-velocity to be in the rest frame of the emitter
  double photon_emframe[4]; // photon tgt vector projected in comoving frame
  for (int ii=0;ii<4;ii++){
    photon_emframe[ii]=cph[ii+4];
  }
  //cout << "photon_emframe: " << photon_emframe[0] << " " << photon_emframe[1] << " " << photon_emframe[2] << " " << photon_emframe[3] << endl;
  gg_->projectFourVect(&cph[0],photon_emframe,vel);
  double bnorm = gg_->norm(&cph[0],B4vect);
  double lnorm = gg_->norm(&cph[0],photon_emframe);
  double lscalb = gg_->ScalarProd(&cph[0],photon_emframe,B4vect);
  double theta_mag = acos(lscalb/(lnorm*bnorm));
  //cout << "r, theta, phi: " << rr << " " << theta << " " << phi << endl;
  //cout << "B4vect: " << B4vect[0] << " " << B4vect[1] << " " << B4vect[2] << " " << B4vect[3] << endl;
  //cout << "photon_emframe: " << photon_emframe[0] << " " << photon_emframe[1] << " " << photon_emframe[2] << " " << photon_emframe[3] << endl;
  //cout << "bnorm, lnorm, lscalb, lscalb/(lnorm*bnorm): " << bnorm << " " << lnorm << " " << lscalb << " " << lscalb/(lnorm*bnorm) << endl;
  //cout << " " << endl;
  

  double n0 = 6e6; // cm-3
  double Theta0=200.; // Dimensionless temperature

  double coord_spot[4]={co[0]};
  const_cast<Blob*>(this)
    ->getCartesian(coord_spot, 1, coord_spot+1, coord_spot+2, coord_spot+3);
  //above: nasty trick to deal with constness of emission
  double xspot=coord_spot[1], yspot=coord_spot[2], zspot=coord_spot[3];

  double difx=(xx-xspot), dify=(yy-yspot), difz=(zz-zspot);
  double d2 = difx*difx+dify*dify+difz*difz;
  double ds2=1.5*1.5;
  //cout << "d2, ds2:" << ", " << d2 << ", " << ds2 << endl;

  double number_density=n0/3.*exp(-d2/(2.*ds2)),
  Theta=Theta0*pow(rr,-0.84);
  double BB = BB0*pow(rr,-1.);

  double nu0 = GYOTO_ELEMENTARY_CHARGE_CGS*BB
    /(2.*M_PI*GYOTO_ELECTRON_MASS_CGS*GYOTO_C_CGS); // cyclotron freq

  double temperature=GYOTO_ELECTRON_MASS_CGS*GYOTO_C2_CGS*Theta/GYOTO_BOLTZMANN_CGS;

  
  double bessK2 = bessk(2, 1./Theta);
  //cout << "In Blob: ne, temperature, BB, nu0, besselK2, theta_mag: " << number_density << " " << temperature << " " << B0 << " " << nu0 << " " << bessK2 << " " << theta_mag << endl;

  // THERMAL SYNCHRO
  spectrumThermalSynch_->numberdensityCGS(number_density);
  spectrumThermalSynch_->temperature(temperature);
  spectrumThermalSynch_->besselK2(bessK2);
  spectrumThermalSynch_->angle_averaged(0);
  spectrumThermalSynch_->angle_B_pem(theta_mag);
  spectrumThermalSynch_->cyclotron_freq(nu0);

  if (number_density<1.e5) {
    for (size_t ii=0; ii<nbnu; ++ii){
      jInu[ii]=0.;
      jQnu[ii]=0.;
      jUnu[ii]=0.;
      jVnu[ii]=0.;
      aInu[ii]=0.;
      aQnu[ii]=0.;
      aUnu[ii]=0.;
      aVnu[ii]=0.;
      rotQnu[ii]=0.;
      rotUnu[ii]=0.;
      rotVnu[ii]=0.;
    }
  }else{
  spectrumThermalSynch_->radiativeQ(jInu, jQnu, jUnu, jVnu,
                                    aInu, aQnu, aUnu, aVnu,
                                    rotQnu, rotUnu, rotVnu, nuem, nbnu);
  }

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

    if (Inu[ii]<0.)
      GYOTO_ERROR("In Blob::radiativeQ: Inu<0");
    if (Inu[ii]!=Inu[ii] or Onu[ii](0,0)!=Onu[ii](0,0))
      GYOTO_ERROR("In Blob::radiativeQ: Inu or Taunu is nan");
    if (Inu[ii]==Inu[ii]+1. or Onu[ii](0,0)==Onu[ii](0,0)+1.)
      GYOTO_ERROR("In Blob::radiativeQ: Inu or Taunu is infinite");

    //cout << "In Blob: Inu, Qnu, Unu, Vnu, dsem: " << Inu[ii] << ", " << Qnu[ii] << ", " << Unu[ii] << ", " << Vnu[ii] << ", " << dsem << endl;
    //cout << "Onu :" << endl;
    //cout << Omat << endl;
  }
}


void Blob::magneticConfiguration(string config){
  magneticConfig_=config;
}

string Blob::magneticConfiguration() const{
  return magneticConfig_;
}