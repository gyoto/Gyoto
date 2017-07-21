/*
    Copyright 2017 Frederic Vincent

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
#include "GyotoJet.h"
#include "GyotoProperty.h"
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

using namespace std;
using namespace Gyoto;
using namespace Gyoto::Astrobj;

/*
  The formalism in this class is taken from the ZSS paper:
  Zdziarski, Stawarz, Sikora, MNRAS 2017 (in prep at time of coding)
 */

GYOTO_PROPERTY_START(Jet)
GYOTO_PROPERTY_DOUBLE(Jet, BaseJetHeight, baseJetHeight)
GYOTO_PROPERTY_DOUBLE(Jet, BaseJetRadiusOverHeight, baseJetRadiusOverHeight)
GYOTO_PROPERTY_DOUBLE(Jet, GammaMax, gammaMax)
GYOTO_PROPERTY_DOUBLE(Jet, MdotJet, mdotJet)
GYOTO_PROPERTY_DOUBLE(Jet, AlfvenRadiusCoef, alfvenRadiusCoef)
GYOTO_PROPERTY_DOUBLE(Jet, ExpoPL, expoPL)
GYOTO_PROPERTY_END(Jet, Standard::properties)

#define nstep_angint 10 // for angle-averaging integration

// ACCESSORS
void Jet::baseJetHeight(double hh) {baseJetHeight_=hh;}
double Jet::baseJetHeight()const{return baseJetHeight_;}
void Jet::baseJetRadiusOverHeight(double par) {baseJetRadiusOverHeight_=par;}
double Jet::baseJetRadiusOverHeight()const{return baseJetRadiusOverHeight_;}
void Jet::gammaMax(double gam) {gammaMax_=gam;}
double Jet::gammaMax()const{return gammaMax_;}
void Jet::mdotJet(double mdot) {mdotJet_=mdot;}
double Jet::mdotJet()const{return mdotJet_;}
void Jet::alfvenRadiusCoef(double coef) {alfvenRadiusCoef_=coef;}
double Jet::alfvenRadiusCoef()const{return alfvenRadiusCoef_;}
void Jet::expoPL(double index) {expoPL_=index;}
double Jet::expoPL()const{return expoPL_;}

//

Jet::Jet() :
  Standard("Jet"), aa_(0.), baseJetHeight_(1.), baseJetRadiusOverHeight_(1.),
  gammaMax_(1.), mdotJet_(1.), alfvenRadiusCoef_(1.), expoPL_(1.)
{
  GYOTO_DEBUG << endl;
}

Jet::Jet(const Jet& o) :
  Standard(o), aa_(o.aa_), baseJetHeight_(o.baseJetHeight_),
  baseJetRadiusOverHeight_(o.baseJetRadiusOverHeight_),
  gammaMax_(o.gammaMax_), mdotJet_(o.mdotJet_),
  alfvenRadiusCoef_(o.alfvenRadiusCoef_), expoPL_(o.expoPL_)
{
  GYOTO_DEBUG << endl;
  gg_->hook(this);
}
Jet* Jet::clone() const
{ return new Jet(*this); }

Jet::~Jet() {
  GYOTO_DEBUG << endl;
  if (gg_) gg_->unhook(this);
}

double Jet::emission(double nu, double,
		     double *,
		     double coord_obj[8]) const{
  return 1.;
}

void Jet::radiativeQ(double Inu[], // output
		     double Taunu[], // output
		     double nu_ems[], size_t nbnu, // input
		     double dsem,
		     double coord_ph[8],
		     double coord_obj[8]) const {
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
    throwError("In Jet::radiativeQ: Unknown coordinate system kind");
  }

  double Mbh = gg_->mass()*1e3; // cgs BH mass
  double MdotEdd = 1.26e38*Mbh/(GYOTO_SUN_MASS_CGS*GYOTO_C_CGS*GYOTO_C_CGS);
       // cgs Eddington accretion rate
  double jetQ[3];
  JetQuantitiesFromZ(zz,jetQ);
  double rcyljet = jetQ[0],
    rcyljetcgs = jetQ[0]*GYOTO_G_OVER_C_SQUARE_CGS*Mbh;
  double Gamma = jetQ[1];
  double number_density = mdotJet_*MdotEdd
    /(2*M_PI*rcyljetcgs*rcyljetcgs*GYOTO_PROTON_MASS_CGS*GYOTO_C_CGS
      *sqrt(Gamma*Gamma-1.));
  //cout << "jet nb dens= " << number_density << endl;
  
  double aGamma = pow(baseJetHeight_,-0.5);
  double ar = baseJetRadiusOverHeight_*pow(baseJetHeight_,0.5);
  double sigma = 0.5*ar*aGamma;
  double Bphi2 = 4.*M_PI * GYOTO_PROTON_MASS_CGS * GYOTO_C_CGS * GYOTO_C_CGS
    *number_density
    *(gammaMax_/Gamma*(1.+sigma*sigma)-1.);

  double rcylA = 5.; // TEST!! to be changed, with light cylinder exp
  double jetQr[2];
  JetQuantitiesFromR(rcylA,jetQr);
  double zA=jetQr[0], GammaA=jetQr[1];
  double Bphi2zA = 4.*M_PI * GYOTO_PROTON_MASS_CGS * GYOTO_C_CGS * GYOTO_C_CGS
    *number_density
    *(gammaMax_/GammaA*(1.+sigma*sigma)-1.),
    BphizA = sqrt(Bphi2zA);
  double Bp2 = BphizA*GammaA*rcylA*rcylA/(rcyljet*rcyljet);

  double BB = sqrt(Bphi2+Bp2);

  double nu0 = GYOTO_ELEMENTARY_CHARGE_CGS*BB
    /(2.*M_PI*GYOTO_ELECTRON_MASS_CGS*GYOTO_C_CGS); // cyclotron freq

  // Synchrotron emission / absorption
  for (size_t ii=0; ii<nbnu; ++ii){
    double nuem = nu_ems[ii];
    double emis_synch_PL=emissionSynchro_PL_averaged(number_density,
						     nuem,nu0),
      abs_synch_PL=absorptionSynchro_PL_averaged(number_density,
						 nuem,nu0);
    
    double delta_s = dsem*GYOTO_G_OVER_C_SQUARE_CGS*Mbh;
    Inu[ii]=
      emis_synch_PL*GYOTO_INU_CGS_TO_SI*delta_s*exp(-abs_synch_PL*delta_s);
    Taunu[ii]=exp(-abs_synch_PL*delta_s);
    //cout << "emis, abs, Inu, Taunu= " << emis_synch_PL << " " << abs_synch_PL << " " << Inu[ii] << " " << Taunu[ii] << endl;
    // NB: abs_synch is in cgs (cm^-1) as well as delta_s (cm)
  }

}

double Jet::operator()(double const coord[4]) {
  double rcyl=0.; // cylindrical radius
  double zz=0.; // height, z coord
  switch (gg_->coordKind()) {
  case GYOTO_COORDKIND_SPHERICAL:
    rcyl = coord[1]*sin(coord[2]);
    zz   = coord[1]*cos(coord[2]);
    break;
  case GYOTO_COORDKIND_CARTESIAN:
    rcyl = pow(coord[1]*coord[1]+coord[2]*coord[2], 0.5);
    zz   = coord[3];
    break;
  default:
    throwError("In Jet::operator(): Unknown coordinate system kind");
  }

  if  (fabs(zz) < baseJetHeight_) return 1.; // outside jet

  double jetQ[3];
  JetQuantitiesFromZ(zz, jetQ);
  double rcyljet=jetQ[0];
  //cout << "r, rjet, z, theta0, ht= " << rcyl << " " << rcyljet << " " << zz << " " << theta0 << " " << ht << endl;
  
  return rcyl-rcyljet; // inside jet when <=0
}

void Jet::getVelocity(double const pos[4], double vel[4])
{
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
    throwError("In Jet::getVelocity: Unknown coordinate system kind");
  }

  double jetQ[3];
  JetQuantitiesFromZ(zz, jetQ);

  double rcyljet=jetQ[0];
  double Gamma=jetQ[1];
  double dzOverdrho=jetQ[2];

  double beta0 = M_PI/2.,
    beta1 = (atan(dzOverdrho)-beta0)/rcyljet,
    beta = beta0+beta1*rcyl;

  double Gamma2=Gamma*Gamma, v2 = (Gamma2-1.)/Gamma2;
  double grr = gg_->gmunu(pos,1,1), gthth=gg_->gmunu(pos,2,2);
  double mycos = cos(pos[2]+beta), mycos2=mycos*mycos,
    mysin=sin(pos[2]+beta), mysin2=mysin*mysin;
  double Vtilde = sqrt(v2/(grr*mysin2 + gthth*mycos2));
  double Vr = Vtilde*mysin, Vth = Vtilde*mycos;

  double gpp = gg_->gmunu(pos,3,3), gtt = gg_->gmunu(pos,0,0),
    gtp = gg_->gmunu(pos,0,3);
  double utZAMO = sqrt(-gpp/(gtt*gpp-gtp*gtp)),
    uphiZAMO = -utZAMO*gtp/gpp;
  
  vel[0] = Gamma*utZAMO;
  vel[1] = Gamma*Vr;
  vel[2] = Gamma*Vth;
  vel[3] = Gamma*uphiZAMO;

  //cout << "jet stuff= " << zz << " " << ht << " " << ar << " " << 1./theta0 << endl;
  //cout <<"beta stuff= " << dzOverdrho << " " << rcyl << " " << rcyljet << " " << beta0 << " " << beta1 << " " << beta << endl;
  //cout << "velo stuff= " << pos[1] << " " << pos[2] << " " << Vr << " " << Vth << " " << beta << endl;
  //cout << "u2 = " << gg_->ScalarProd(pos,vel,vel) << endl;
  //cout << "4-vel= " << vel[0] << " " << vel[1] << " " << vel[2] << " " << vel[3] << endl;
}

bool Jet::isThreadSafe() const {
  return Standard::isThreadSafe();
}

void Jet::metric(SmartPointer<Metric::Generic> gg) {
  if (gg_) gg_->unhook(this);
  string kin = gg->kind();
  if (kin != "KerrBL" && kin != "KerrKS")
    throwError
      ("Jet::metric(): metric must be KerrBL or KerrKS");
  Generic::metric(gg);
  updateSpin();
  gg->hook(this);
}

void Jet::updateSpin() {
  if (!gg_) return;
  switch (gg_->coordKind()) {
  case GYOTO_COORDKIND_SPHERICAL:
    aa_ = static_cast<SmartPointer<Metric::KerrBL> >(gg_) -> spin();
    break;
  case GYOTO_COORDKIND_CARTESIAN:
    aa_ = static_cast<SmartPointer<Metric::KerrKS> >(gg_) -> spin();
    break;
  default:
    throwError("Jet::updateSpin(): unknown COORDKIND");
  }
}

void Jet::tell(Hook::Teller* msg) {
  if (msg==gg_) updateSpin();
}

void Jet::JetQuantitiesFromZ(const double zz, double qty[3]) const{
  // Computes the jet quantities as a function z,
  // ie rcyl(z), Gamma(z), dz/drcyl(z)

  double aGamma = pow(baseJetHeight_,-0.5); // Gamma(z) = aGamma * sqrt{z}
  double ar = baseJetRadiusOverHeight_*pow(baseJetHeight_,0.5);
            // rcyl(z) = ar * sqrt{z} ;
            // assuming rcyl(zbase) = baseJetRadiusOverHeight_*zbase
  // [aGamma] = 1/sqrt{R} ; [ar] = sqrt{R} ; [ar*aGamma] = 1
  double theta0 = 0.5*ar*aGamma/gammaMax_;
  double ht = gammaMax_*gammaMax_/(aGamma*aGamma);
  // NB: ht is such that baseJetHeight_ <= ht
  // so the jet base is always in the parabola part of the jet shape
  double rcyljet=0.;
  double dzOverdrho=0.;
  double Gamma=0.;
  if (fabs(zz) <= ht) {
    rcyljet = ar*pow(fabs(zz),0.5);
    Gamma = aGamma*pow(fabs(zz),0.5);
    dzOverdrho = 2*rcyljet*rcyljet/(ar*ar);
    if (zz<0.) dzOverdrho*=-1.; // dz/drho<0 for z<0
  }else{
    rcyljet = theta0*(fabs(zz)+ht);
    Gamma = gammaMax_;
    dzOverdrho = 1./theta0;
    if (zz<0.) dzOverdrho*=-1.;
  }

  qty[0]=rcyljet;
  qty[1]=Gamma;
  qty[2]=dzOverdrho;
}

void Jet::JetQuantitiesFromR(const double rr, double qty[2]) const{
  // Computes the jet quantities as a function rcyl,
  // ie z(rcyl), Gamma(rcyl)

  double aGamma = pow(baseJetHeight_,-0.5);
  double ar = baseJetRadiusOverHeight_*pow(baseJetHeight_,0.5);

  double rt = ar*gammaMax_/aGamma;
  double ht = gammaMax_*gammaMax_/(aGamma*aGamma);
  double theta0 = 0.5*ar*aGamma/gammaMax_;

  double zjet=0.;
  double Gamma=0.;
  if (rr <= rt) {
    zjet = rr*rr/(ar*ar);
    Gamma = aGamma*sqrt(zjet);
  }else{
    zjet = rr/theta0 - ht;
    Gamma = gammaMax_;
  }

  qty[0]=zjet;
  qty[1]=Gamma;
}

////// ******* %%%%%%%% ********* /////////
////// ******* %%%%%%%% ********* /////////
////// ******* %%%%%%%% ********* /////////

//// THESE FUNCTIONS ARE COPIED FROM POLISH DOUGHNUTS: TO PUT IN
//// ANOTHER CLASS:

double Jet::emissionSynchro_PL_direction(double number_density_PL,
					 double nuem, double nuc,
					 double theta_mag)
  const {
  // From Petrosian & McTiernan 1983, Phys. Fluids 26 (10), eq. 32
  // Putting g(mu)=1 and using (Y+ + Y_)=2 to get jnu and alphanu.
  // NB: putting g(mu)=1 or 1/2 is not important, it boils down
  // to redefining the % amount delta of PL energy wrt THER energy
  double emis_synch =
    sqrt(3.)*M_PI*GYOTO_ELEMENTARY_CHARGE_CGS*GYOTO_ELEMENTARY_CHARGE_CGS
    *nuc*sin(theta_mag)/(2.*GYOTO_C_CGS)
    *number_density_PL*(expoPL_-1.)
    *pow(3.*nuc*(expoPL_+1.)*sin(theta_mag)/(4.*nuem),0.5*(expoPL_-1.))
    *exp(-0.5*(expoPL_+1.));
  if (emis_synch!=emis_synch) {
    //cout << "stuff= " << nuc << " " << theta_mag << " " << number_density_PL << endl;
    throwError("In Jet::emissionSynchro_PL_direction: "
	       "emissivity is nan");
  }
  if (emis_synch==emis_synch+1.)
    throwError("In Jet::emissionSynchro_PL_direction: "
	       "emissivity is infinite");
  return emis_synch;
}
double Jet::absorptionSynchro_PL_direction(double number_density_PL,
					   double nuem, double nuc,
					   double theta_mag)
  const {
  // From Petrosian & McTiernan 1983, Phys. Fluids 26 (10), eq. 32
  double abs_synch =
    sqrt(3.)*M_PI*GYOTO_ELEMENTARY_CHARGE_CGS*GYOTO_ELEMENTARY_CHARGE_CGS
    *nuc*sin(theta_mag)/(2.*GYOTO_C_CGS)
    *number_density_PL*(expoPL_-1.)
    *pow(3.*nuc*(expoPL_+2.)*sin(theta_mag)/(4.*nuem),0.5*expoPL_)
    *exp(-0.5*(expoPL_+2.))
    *(expoPL_+2.)
    /(GYOTO_ELECTRON_MASS_CGS*nuem*nuem);
  if (abs_synch!=abs_synch) {
    throwError("In Jet::absorptionSynchro_PL_direction: "
	       "abs is nan");
  }
  if (abs_synch==abs_synch+1.)
    throwError("In Jet::absorptionSynchro_PL_direction: "
	       "abs is infinite");
  return abs_synch;
}
double Jet::emissionSynchro_PL_averaged(double number_density_PL,
					double nuem, double nuc)
  const {
  double th0=0., thNm1=M_PI;
  double hh=(thNm1-th0)/double(nstep_angint);
  double emis_synch=0.;
  for (int ii=1;ii<=2*nstep_angint-3;ii+=2){
    double theta=th0+double(ii)/2.*hh;
    emis_synch+=hh*emissionSynchro_PL_direction(number_density_PL,
						nuem,nuc,theta)
      *sin(theta);
  }
  if (emis_synch!=emis_synch) {
    throwError("In Jet::emissionSynchro_PL_averaged: "
	       "emissivity is nan");
  }
  if (emis_synch==emis_synch+1.)
    throwError("In Jet::emissionSynchro_PL_averaged: "
	       "emissivity is infinite");
  return emis_synch/2.;
  //NB: averaged jnu is: \int jnu dOmega = 1/2 * \int jnu*sinth dth
}
double Jet::absorptionSynchro_PL_averaged(double number_density_PL,
					  double nuem, double nuc)
  const {
  double th0=0., thNm1=M_PI;
  double hh=(thNm1-th0)/double(nstep_angint);
  double abs_synch=0.;
  for (int ii=1;ii<=2*nstep_angint-3;ii+=2){
    double theta=th0+double(ii)/2.*hh;
    abs_synch+=hh*absorptionSynchro_PL_direction(number_density_PL,
						 nuem,nuc,theta)
      *sin(theta);
  }
  if (abs_synch!=abs_synch) {
    throwError("In Jet::absorptionSynchro_PL_averaged: "
	       "abs is nan");
  }
  if (abs_synch==abs_synch+1.)
    throwError("In Jet::absorptionSynchro_PL_averaged: "
	       "abs is infinite");
  return abs_synch/2.;
}
