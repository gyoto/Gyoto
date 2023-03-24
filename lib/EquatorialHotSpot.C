/*
    Copyright 2016, 2018-2020 Frederic Vincent, Thibaut Paumard

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
#include "GyotoEquatorialHotSpot.h"
#include "GyotoPhoton.h"
#include "GyotoPageThorneDisk.h"
#include "GyotoUtils.h"
#include "GyotoFactoryMessenger.h"
#include "GyotoKerrBL.h"
#include "GyotoKerrKS.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <cmath>
#include <limits>
#include <string>
#include <cstring>
#include <time.h> 

using namespace std;
using namespace Gyoto;
using namespace Gyoto::Astrobj;

/// Properties
#include "GyotoProperty.h"
GYOTO_PROPERTY_START(EquatorialHotSpot, "Equatorial hot spot with beaming")
GYOTO_PROPERTY_DOUBLE(EquatorialHotSpot, SpotRadSize, spotRadSize)
GYOTO_PROPERTY_STRING(EquatorialHotSpot, BeamingKind, beaming,
		      "One of: IsotropicBeaming, NormalBeaming, RadialBeaming, "
		      "IsotropicConstant (emission is isotropic and constant"
		      "equals to 1)")
GYOTO_PROPERTY_DOUBLE(EquatorialHotSpot, BeamAngle, beamAngle)
GYOTO_WORLDLINE_PROPERTY_END(EquatorialHotSpot, ThinDisk::properties)

// accessors
void EquatorialHotSpot::spotRadSize(double t) {sizespot_=t;}
double EquatorialHotSpot::spotRadSize() const {return sizespot_;}

void EquatorialHotSpot::beaming(std::string const &b) {
  if (b=="IsotropicBeaming") beaming_=IsotropicBeaming;
  else if (b=="NormalBeaming") beaming_=NormalBeaming;
  else if (b=="RadialBeaming") beaming_=RadialBeaming;
  else if (b=="IsotropicConstant") beaming_=IsotropicConstant;
  else GYOTO_ERROR("Unknown beaming kind");
}
std::string EquatorialHotSpot::beaming() const {
  string b;
  switch (beaming_) {
  case IsotropicBeaming: b="IsotropicBeaming"; break;
  case NormalBeaming:    b="NormalBeaming";    break;
  case RadialBeaming:    b="RadialBeaming";    break;
  case IsotropicConstant: b="IsotropicConstant"; break;
  default: GYOTO_ERROR("Unknown beaming kind");
  }
  return b;
}

void EquatorialHotSpot::beamAngle(double t) {beamangle_=t;}
double EquatorialHotSpot::beamAngle() const {return beamangle_;}

// Needed for legacy XML files
int EquatorialHotSpot::setParameter(string name, string content, string unit) {
  double coord[8];
  char* tc = const_cast<char*>(content.c_str());
  if (name=="InitialCoordinate") {
    name="InitCoord";
    return ThinDisk::setParameter(name, content, unit);
  } else if (name=="Position") {
    if (FactoryMessenger::parseArray(content, coord, 4) != 4)
      GYOTO_ERROR("Worldline \"Position\" requires exactly 4 tokens");
    if (init_vel_) {
      setInitCoord(coord, init_vel_);
      delete[] init_vel_; init_vel_=NULL;
    } else setPosition(coord);
    wait_pos_ = 0;
  } else if (name=="Velocity") {
    if (FactoryMessenger::parseArray(content, coord, 3) != 3)
      GYOTO_ERROR("Worldline \"Velocity\" requires exactly 3 tokens");
    if (wait_pos_) {
      if (init_vel_) delete [] init_vel_;
      init_vel_ = new double[3];
      memcpy(init_vel_, coord, 3*sizeof(double));
    } else setVelocity(coord);
  } else if (name=="NormalBeaming") {
    GYOTO_WARNING << "<" << name << "/> is deprecated, please use "
      "<BeamingKind> " << name << " </BeamingKind> instead";
    beaming(name);
  } else if (name=="NormalBeaming" || name=="RadialBeaming") {
    GYOTO_WARNING << "<" << name << "/> is deprecated, please use \n";
    GYOTO_WARNING << "<BeamingKind> " << name << " </BeamingKind>" << endl;
    GYOTO_WARNING << "<BeamAngle> " << content << "</BeamAngle>" << endl;
    GYOTO_WARNING <<" instead";
    beaming(name);
    setParameter("BeamAngle", content, unit);
  } else return ThinDisk::setParameter(name, content, unit);
  return 0;
}

// Needed for wait_pos_
#ifdef GYOTO_USE_XERCES
void EquatorialHotSpot::fillProperty(Gyoto::FactoryMessenger *fmp, Property const &p) const {
  if (p.name == "InitCoord") {
    if (imin_ <= imax_) {
      state_t coord;
      getInitialCoord(coord);
      // For massive particule, express initial condition with 3-velocity
      double vel[3] = {coord[5]/coord[4], coord[6]/coord[4], coord[7]/coord[4]};
      fmp -> setParameter ("Position", &coord[0], 4);
      fmp -> setParameter ("Velocity", vel, 3);
    }
    return;
  }
  ThinDisk::fillProperty(fmp, p);
}

void EquatorialHotSpot::setParameters(FactoryMessenger* fmp) {
  wait_pos_ = 1;
  ThinDisk::setParameters(fmp);
  wait_pos_ = 0;
  if (init_vel_) {
    delete[] init_vel_; init_vel_=NULL;
    GYOTO_ERROR("Worldline::setParameters(): "
	       "Velocity was found but not Position");
  }
}
#endif
///

Gyoto::Astrobj::EquatorialHotSpot::EquatorialHotSpot()
  : ThinDisk("EquatorialHotSpot"), Worldline(), 
    sizespot_(0.), beaming_(IsotropicBeaming), beamangle_(0.), spectrumThermalSynch_(NULL), magneticConfig_("None")
{
  GYOTO_DEBUG  << "Building EquatorialHotSpot";
  spectrumThermalSynch_ = new Spectrum::ThermalSynchrotron();
}

Gyoto::Astrobj::EquatorialHotSpot::EquatorialHotSpot(const EquatorialHotSpot &o)
  : ThinDisk(o), Worldline(o),
    sizespot_(o.sizespot_), beaming_(o.beaming_), beamangle_(o.beamangle_), spectrumThermalSynch_(NULL), magneticConfig_(o.magneticConfig_)
{
  GYOTO_DEBUG  << "Copying EquatorialHotSpot";
  if (o.spectrumThermalSynch_()) spectrumThermalSynch_=o.spectrumThermalSynch_->clone();
}
EquatorialHotSpot * EquatorialHotSpot::clone() const { 
  return new EquatorialHotSpot(*this); }

Gyoto::Astrobj::EquatorialHotSpot::~EquatorialHotSpot()
{
  GYOTO_DEBUG  << "Destroying EquatorialHotSpot";
}

double EquatorialHotSpot::getMass() const {return 1. ;}

void EquatorialHotSpot::metric(SmartPointer<Metric::Generic> gg) {
  ThinDisk::metric(gg);
  Worldline::metric(gg);
}

void EquatorialHotSpot::setInitialCondition(double coord[8]) {
  if (!metric_) GYOTO_ERROR("Please set metric before calling "
			   "EquatorialHotSpot::setInitialCondition(double*)");
  Worldline::setInitialCondition(metric_, coord, 1);
}

void EquatorialHotSpot::getVelocity(double const pos[4], double vel[4]) {
  double coord_spot[4]={pos[0]};
  const_cast<EquatorialHotSpot*>(this)
    ->getCoord(coord_spot, 1, coord_spot+1, coord_spot+2, coord_spot+3);
  gg_ -> circularVelocity(coord_spot, vel, dir_);
  //cout << "equat Omega= " << vel[3]/vel[0] << endl;
}

double EquatorialHotSpot::emission(double nu_em, double dsem,
				   state_t const &coord_ph,
				   double const coord_obj[8]) const{
  double coord_spot[4]={coord_obj[0]};
  const_cast<EquatorialHotSpot*>(this)
    ->getCartesian(coord_spot, 1, coord_spot+1, coord_spot+2, coord_spot+3);
  //above: nasty trick to deal with constness of emission
  double xspot=coord_spot[1], yspot=coord_spot[2];
  //cout << "spot is at xy= " << xspot << " " << yspot << endl;
  double rr=coord_obj[1], phi=coord_obj[3];
  double difx=(rr*cos(phi)-xspot),
    dify=(rr*sin(phi)-yspot);
  double d2 = difx*difx+dify*dify;
  double ds2=sizespot_*sizespot_;
  if (d2 < 16*ds2){ // we are within 4*rspot,
                    // same as in Schnittman & Bertschinger 2004

    // computing the angle (normal,photon tangent)
    double cosalpha=0.;
    if (beaming_ == NormalBeaming or beaming_ == RadialBeaming){
      double gthth=gg_->gmunu(&coord_ph[0],2,2);
      double pth=coord_ph[6];
      double uemitter[4];
      const_cast<EquatorialHotSpot*>(this)
	->getVelocity(&coord_ph[0],uemitter);
      double pscalu=fabs(gg_->ScalarProd(&coord_ph[0],&coord_ph[4],
					 uemitter));
      if (pscalu==0.) GYOTO_ERROR("Undefined cosalpha!");
      cosalpha = 1./pscalu*sqrt(gthth)*fabs(pth); // = |cos(alpha)|

      if (fabs(cosalpha)>1.)
	GYOTO_ERROR("cosalpha>1!");
    }
      
    // emission Gaussian width
    double sigma2=ds2 ; // following choice of Schnittman & Bertschinger 2004:
                        // sigma = Rspot

    switch (beaming_) {
    case IsotropicBeaming:
      return exp(-d2/(2*sigma2));
    case IsotropicConstant:
      return 1.;
    case NormalBeaming:
      return cosalpha*cosalpha*exp(-d2/(2*sigma2));
    case RadialBeaming:
      return (1.-cosalpha)*(1.-cosalpha)*exp(-d2/(2*sigma2));
    default:
      GYOTO_ERROR("In EquatorialHotSpot::emission:"
		 " incorrect beaming argument");
    }
  }
  // else
  return 0.;
}

void EquatorialHotSpot::radiativeQ(double *Inu, double *Qnu, double *Unu,
				   double *Vnu,
				   Eigen::Matrix4d *Onu,
				   double const *nuem , size_t nbnu,
				   double dsem,
				   state_t const &cph,
				   double const *co) const {
  // polarized radiativeQ

  double vel[4]; // 4-velocity of emitter
  for (int ii=0;ii<4;ii++){
    vel[ii]=co[ii+4];
  }
  
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
  if (magneticConfig_=="None")
    GYOTO_ERROR("Specify the magnetic field configuration");
  if (magneticConfig_=="Vertical"){
    B4vect[2]=-1;
  }
  else if (magneticConfig_=="Radial"){
    B4vect[1]=1;
  }
  else if (magneticConfig_=="Toroidal"){
    double gtt=gg_->gmunu(&cph[0],0,0),
      gpp=gg_->gmunu(&cph[0],3,3),
      Bp=1.,
      Bt=(-gpp*Bp*vel[3])/(gtt*vel[0]);
    B4vect[0]=Bt;
    B4vect[3]=Bp;
  }
  else if (magneticConfig_=="Radial-Azimuthal"){
    double gtt=gg_->gmunu(&cph[0],0,0),
      gpp=gg_->gmunu(&cph[0],3,3),
      Bp=1.,
      Bt=(-gpp*Bp*vel[3])/(gtt*vel[0]);
    B4vect[0]=Bt;
    B4vect[1]=1.;
    B4vect[3]=Bp;
  }
  else
    GYOTO_ERROR("Unknown magnetic field configuration");
  
  double B0 = 100.; // Gauss
  double norm=sqrt(gg_->ScalarProd(&cph[0], B4vect, B4vect));
  gg_->multiplyFourVect(B4vect,1./norm);

  double Chi=getChi(B4vect, cph, vel); // this is EVPA

  double nu0 = GYOTO_ELEMENTARY_CHARGE_CGS*B0
    /(2.*M_PI*GYOTO_ELECTRON_MASS_CGS*GYOTO_C_CGS); // cyclotron freq

  // Computing the angle theta_mag between the magnetic field vector and photon tgt vector in the rest frame of the emitter
  gg_->projectFourVect(&cph[0],B4vect,vel); //Projection of the 4-vector B to 4-velocity to be in the rest frame of the emitter
  double photon_emframe[4]; // photon tgt vector projected in comoving frame
  for (int ii=0;ii<4;ii++){
    photon_emframe[ii]=cph[ii+4]+vel[ii]*gg_->ScalarProd(&cph[0],&cph[4],vel);
  }
  double bnorm = gg_->norm(&cph[0],B4vect);
  double lnorm = gg_->norm(&cph[0],photon_emframe);
  double lscalb = gg_->ScalarProd(&cph[0],photon_emframe,B4vect);
  double theta_mag = acos(lscalb/(lnorm*bnorm));
  

  double n0 = 6e6; // cm-3
  double Theta0=200; // Dimensionless temperature

  double coord_spot[4]={co[0]};
  const_cast<EquatorialHotSpot*>(this)
    ->getCartesian(coord_spot, 1, coord_spot+1, coord_spot+2, coord_spot+3);
  //above: nasty trick to deal with constness of emission
  double xspot=coord_spot[1], yspot=coord_spot[2];
  //cout << "spot is at xy= " << xspot << " " << yspot << endl;
  double rr=co[1], phi=co[3];
  double difx=(rr*cos(phi)-xspot),
    dify=(rr*sin(phi)-yspot);
  double d2 = difx*difx+dify*dify;
  double ds2=sizespot_*sizespot_;

  double number_density=n0/3.*exp(-d2/(2.*ds2)),
  Theta=Theta0*exp(-d2/(2.*ds2))*pow(rr,-0.84);

  if (number_density<1e5) number_density=0.;
  double temperature=GYOTO_ELECTRON_MASS_CGS*GYOTO_C2_CGS*Theta/GYOTO_BOLTZMANN_CGS;

  
  double bessK2 = bessk(2, 1./Theta);
  //cout << "In EquatorialHotSpot: ne, temperature, BB, nu0, besselK2, theta_mag: " << number_density << " " << temperature << " " << B0 << " " << nu0 << " " << bessK2 << " " << theta_mag << endl;

  // THERMAL SYNCHRO
  spectrumThermalSynch_->numberdensityCGS(number_density);
  spectrumThermalSynch_->temperature(temperature);
  spectrumThermalSynch_->besselK2(bessK2);
  spectrumThermalSynch_->angle_averaged(0);
  spectrumThermalSynch_->angle_B_pem(theta_mag);
  spectrumThermalSynch_->cyclotron_freq(nu0);

  spectrumThermalSynch_->radiativeQ(jInu, jQnu, jUnu, jVnu,
                                    aInu, aQnu, aUnu, aVnu,
                                    rotQnu, rotUnu, rotVnu, nuem, nbnu);


  for (size_t ii=0; ii<nbnu; ++ii) {
    //cout << "d2, ds2:" << ", " << d2 << ", " << ds2 << endl;
    //cout << "In EquatorialHotSpot: jInu, jQnu, jUnu, jVnu: " << jInu[ii] << ", " << jQnu[ii] << ", " << jUnu[ii] << ", " << jVnu[ii] << endl;
    //cout << "In EquatorialHotSpot: aInu, aQnu, aUnu, aVnu: " << aInu[ii] << ", " << aQnu[ii] << ", " << aUnu[ii] << ", " << aVnu[ii] << endl;
    //cout << "In EquatorialHotSpot: rQnu, rUnu, rVnu: " << rotQnu[ii] << ", " << rotUnu[ii] << ", " << rotVnu[ii] << endl;
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

    //cout << "In EquatorialHotSpot: Inu, Qnu, Unu, Vnu, dsem: " << Inu[ii] << ", " << Qnu[ii] << ", " << Unu[ii] << ", " << Vnu[ii] << ", " << dsem << endl;
    //cout << "Onu :" << endl;
    //cout << Omat << endl;
  }
}


void EquatorialHotSpot::radiativeQ(double Inu[], // output
				   double Taunu[], // output
				   double const nu_ems[], size_t nbnu, // input
				   double dsem,
				   state_t const &coord_ph,
				   double const coord_obj[8]) const {
  // unpolarized radiativeQ
  for (size_t ii=0; ii<nbnu; ++ii) {
    Inu[ii]=emission(nu_ems[ii], dsem, coord_ph, coord_obj);
    Taunu[ii]=1.;
  }
}

void EquatorialHotSpot::magneticConfiguration(string config){
  magneticConfig_=config;
}

string EquatorialHotSpot::magneticConfiguration() const{
  return magneticConfig_;
}