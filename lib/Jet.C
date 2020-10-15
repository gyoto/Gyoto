/*
    Copyright 2017-2019 Frederic Vincent & Thibaut Paumard

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

GYOTO_PROPERTY_START(Jet)
GYOTO_PROPERTY_DOUBLE(Jet, JetOuterOpeningAngle, jetOuterOpeningAngle)
GYOTO_PROPERTY_DOUBLE(Jet, JetInnerOpeningAngle, jetInnerOpeningAngle)
GYOTO_PROPERTY_DOUBLE(Jet, JetBaseHeight, jetBaseHeight)
GYOTO_PROPERTY_DOUBLE(Jet, GammaJet, gammaJet)
GYOTO_PROPERTY_DOUBLE(Jet, JetVphiOverVr, jetVphiOverVr,"this is (r*Vphi/Vr) where V is the jet velocity measured by the ZAMO")
GYOTO_PROPERTY_DOUBLE_UNIT(Jet, BaseNumberDensity, baseNumberDensity)
GYOTO_PROPERTY_DOUBLE(Jet, BaseTemperature, baseTemperature)
GYOTO_PROPERTY_DOUBLE(Jet, TemperatureSlope, temperatureSlope)
GYOTO_PROPERTY_DOUBLE(Jet, MagnetizationParameter,
		      magnetizationParameter)
GYOTO_PROPERTY_DOUBLE(Jet, KappaIndex, kappaIndex, "Index of kappa-distribution synchrotron; leave non-specified to use thermal synchrotron")
GYOTO_PROPERTY_END(Jet, Standard::properties)

#define nstep_angint 10 // for angle-averaging integration
#define default_kappaindex -1 // default (absurd value) kappa index

// ACCESSORS
void Jet::jetOuterOpeningAngle(double ang) {jetOuterOpeningAngle_=ang;}
double Jet::jetOuterOpeningAngle()const{return jetOuterOpeningAngle_;}
void Jet::jetInnerOpeningAngle(double ang) {jetInnerOpeningAngle_=ang;}
double Jet::jetInnerOpeningAngle()const{return jetInnerOpeningAngle_;}
void Jet::jetBaseHeight(double hh) {jetBaseHeight_=hh;}
double Jet::jetBaseHeight()const{return jetBaseHeight_;}
void Jet::gammaJet(double gam) {gammaJet_=gam;}
double Jet::gammaJet()const{return gammaJet_;}
void Jet::jetVphiOverVr(double alpha) {jetVphiOverVr_=alpha;}
double Jet::jetVphiOverVr()const{return jetVphiOverVr_;}
double Jet::baseNumberDensity() const {
  // Converts internal cgs central enthalpy to SI
  double dens=baseNumberDensity_cgs_;
# ifdef HAVE_UDUNITS
  dens = Units::Converter("cm-3", "m-3")(dens);
# else
  GYOTO_WARNING << "Units ignored, please recompile Gyoto with --with-udunits"
		<< endl ;
# endif
  return dens; }
double Jet::baseNumberDensity(string const &unit) const
{
  double dens = baseNumberDensity();
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
void Jet::baseNumberDensity(double dens) {
# ifdef HAVE_UDUNITS
  dens = Units::Converter("m-3", "cm-3")(dens);
# else
  GYOTO_WARNING << "Units ignored, please recompile Gyoto with --with-udunits"
		<< endl ;
# endif
  baseNumberDensity_cgs_=dens;
}
void Jet::baseNumberDensity(double dens, string const &unit) {
  if (unit != "") {
# ifdef HAVE_UDUNITS
    dens = Units::Converter(unit, "m-3")(dens);
# else
    GYOTO_WARNING << "Units ignored, please recompile Gyoto with --with-udunits"
		  << endl ;
# endif
  }
  baseNumberDensity(dens);
}
void Jet::baseTemperature(double tt) {baseTemperature_=tt;}
double Jet::baseTemperature()const{return baseTemperature_;}
void Jet::temperatureSlope(double ss) {temperatureSlope_=ss;}
double Jet::temperatureSlope()const{return temperatureSlope_;}
void Jet::magnetizationParameter(double rr) {
  magnetizationParameter_=rr;}
double Jet::magnetizationParameter()const{
  return magnetizationParameter_;}
void Jet::kappaIndex(double index) {
  spectrumKappaSynch_->kappaindex(index);
}
double Jet::kappaIndex()const{
  return spectrumKappaSynch_->kappaindex();
}

//

Jet::Jet() :
  Standard("Jet"), jetOuterOpeningAngle_(0.785),
  jetInnerOpeningAngle_(0.5), jetBaseHeight_(2.),
  gammaJet_(1.), jetVphiOverVr_(0.),
  baseNumberDensity_cgs_(1.), baseTemperature_(1e10),
  temperatureSlope_(1.),
  magnetizationParameter_(1.)
{
  GYOTO_DEBUG << endl;
  spectrumKappaSynch_   = new Spectrum::KappaDistributionSynchrotron();
  spectrumKappaSynch_->kappaindex(default_kappaindex);
  spectrumThermalSynch_ = new Spectrum::ThermalSynchrotron();
}

Jet::Jet(const Jet& o) :
  Standard(o), jetOuterOpeningAngle_(o.jetOuterOpeningAngle_),
  jetInnerOpeningAngle_(o.jetInnerOpeningAngle_),
  jetBaseHeight_(o.jetBaseHeight_),
  gammaJet_(o.gammaJet_), jetVphiOverVr_(o.jetVphiOverVr_),
  baseNumberDensity_cgs_(o.baseNumberDensity_cgs_),
  baseTemperature_(o.baseTemperature_),
  temperatureSlope_(o.temperatureSlope_),
  magnetizationParameter_(o.magnetizationParameter_),
  spectrumKappaSynch_(NULL),
  spectrumThermalSynch_(NULL)
{
  GYOTO_DEBUG << endl;
  if (gg_) gg_->hook(this);
  if (o.spectrumKappaSynch_()) spectrumKappaSynch_=o.spectrumKappaSynch_->clone();
  if (o.spectrumThermalSynch_()) spectrumThermalSynch_=o.spectrumThermalSynch_->clone();

}
Jet* Jet::clone() const
{ return new Jet(*this); }

Jet::~Jet() {
  GYOTO_DEBUG << endl;
  if (gg_) gg_->unhook(this);
}

void Jet::radiativeQ(double Inu[], // output
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
    GYOTO_ERROR("In Jet::radiativeQ: Unknown coordinate system kind");
  }

  double rcyljetbase = jetBaseHeight_*tan(jetOuterOpeningAngle_);

  //cout << "rcyl, rcylB, zz= " << rcyl << " " << rcyljetbase << " " << zz << endl;

  //rcyl=rcyljetbase;
  //zz=2.; // TEST!!!
  
  double number_density = baseNumberDensity_cgs_
    *(rcyljetbase*rcyljetbase)/(rcyl*rcyl);

  double temperature = baseTemperature_*pow(jetBaseHeight_/fabs(zz),
					    temperatureSlope_);

  double thetae = GYOTO_BOLTZMANN_CGS*temperature
    /(GYOTO_ELECTRON_MASS_CGS*GYOTO_C2_CGS);

  double BB = sqrt(4.*M_PI*magnetizationParameter_
		   *GYOTO_PROTON_MASS_CGS * GYOTO_C_CGS * GYOTO_C_CGS
		   *number_density);
  //cout << "r, z, dens, T, B= " << rcyl << " " << zz << " " << number_density << " " << temperature << " " << BB << endl;
  //cout << "r, z, ne, nebase, B, Bbase= " << coord_ph[1] << " " << zz << " " << number_density << " " << baseNumberDensity_cgs_ << " " << BB << " " << sqrt(8.*M_PI*magnetizationParameter_*GYOTO_PROTON_MASS_CGS * GYOTO_C_CGS * GYOTO_C_CGS*baseNumberDensity_cgs_) << endl;
  //GYOTO_ERROR("testjet");

  double nu0 = GYOTO_ELEMENTARY_CHARGE_CGS*BB
    /(2.*M_PI*GYOTO_ELECTRON_MASS_CGS*GYOTO_C_CGS); // cyclotron freq

  //cout << "jet stuff= " << coord_ph[1] << " " << coord_ph[2] << " " << zz << " " << rcyljetbase << " " << rcyl << " " << number_density << " " << thetae << " " << temperatureSlope_ << " " << nu0 << endl;
  //cout << "jet zz,rcyl,th,ph,ne,Te= " <<  zz << " " << rcyl << " " << coord_ph[2] << " " << coord_ph[3] << " " << number_density << " " << temperature << endl;
  // Use that line for Compton study:
  //cout <<  "jet emis: " << zz << " " << rcyl << " " << number_density << " " << temperature << endl;

  // Emission and absorption synchrotron coefs
  double jnu_synch[nbnu], anu_synch[nbnu];
  for (size_t ii=0; ii<nbnu; ++ii){
    // Initializing to <0 value to create errors if not updated
    jnu_synch[ii]=-1.;
    anu_synch[ii]=-1.;
  }

  if (kappaIndex()!=default_kappaindex){
    // KAPPA-DISTRIB SYNCHROTRON
    spectrumKappaSynch_->numberdensityCGS(number_density);
    spectrumKappaSynch_->angle_averaged(1); // impose angle-averaging
    spectrumKappaSynch_->angle_B_pem(0.); // so we don't care about angle
    spectrumKappaSynch_->cyclotron_freq(nu0);
    spectrumKappaSynch_->thetae(thetae);
    double hypergeom = Gyoto::hypergeom(kappaIndex(), thetae);
    spectrumKappaSynch_->hypergeometric(hypergeom);
    //cout << "jet stuff for kappa: " << nu_ems[0] << " " << number_density << " " << nu0 << " " << thetae << " " << BB << " " << temperature << " " << hypergeom << endl;
    spectrumKappaSynch_->radiativeQ(jnu_synch,anu_synch,
				    nu_ems,nbnu);
  }else{
    // THERMAL SYNCHROTRON
    spectrumThermalSynch_->temperature(temperature);
    spectrumThermalSynch_->numberdensityCGS(number_density);
    spectrumThermalSynch_->angle_averaged(1); // impose angle-averaging
    spectrumThermalSynch_->angle_B_pem(0.);   // so we don't care about angle
    spectrumThermalSynch_->cyclotron_freq(nu0);
    double besselK2 = bessk(2, 1./thetae);
    spectrumThermalSynch_->besselK2(besselK2);
    //cout << "for anu jnu: " << coord_ph[1] << " " << zz << " " << temperature << " " << number_density << " " << nu0 << " " << thetae << " " << besselK2 << endl;
    //cout << "nu passed to synchro= " << nu_ems[0] << endl;
    spectrumThermalSynch_->radiativeQ(jnu_synch,anu_synch,
				      nu_ems,nbnu);
  }

  // RETURNING TOTAL INTENSITY AND TRANSMISSION
  for (size_t ii=0; ii<nbnu; ++ii){

    double jnu_tot = jnu_synch[ii],
      anu_tot = anu_synch[ii];

    //cout << "in jet stuff: " << zz << " " << rcyl << " " << nu_ems[0]  << " " << number_density << " " << nu0 << " " << temperature << " " << thetae << " " << jnu_tot << " " << anu_tot << " " << dsem << endl;

    //cout << "at r,th= " << coord_ph[1] << " " << coord_ph[2] << endl;
    //cout << "at rcyl,z= " << rcyl << " " << zz << endl;
    //cout << "jet jnu anu kappa= " << jnu_tot << " " << anu_tot << endl; //x" " << jnu_tot/anu_tot << " " << dsem << endl;

    // expm1 is a precise implementation of exp(x)-1
    double em1=std::expm1(-anu_tot * dsem * gg_->unitLength());
    Taunu[ii] = em1+1.;
    Inu[ii] = anu_tot == 0. ? jnu_tot * dsem * gg_->unitLength() :
      -jnu_tot / anu_tot * em1;

    if (Inu[ii]<0.)
      GYOTO_ERROR("In Jet::radiativeQ: Inu<0");
    if (Inu[ii]!=Inu[ii] or Taunu[ii]!=Taunu[ii])
      GYOTO_ERROR("In Jet::radiativeQ: Inu or Taunu is nan");
    if (Inu[ii]==Inu[ii]+1. or Taunu[ii]==Taunu[ii]+1.)
      GYOTO_ERROR("In Jet::radiativeQ: Inu or Taunu is infinite");

  }
}

double Jet::operator()(double const coord[4]) {
  //cout << "photon at r,z= " << coord[1] << " " << coord[1]*cos(coord[2]) << endl;
  double rcyl=0.; // cylindrical radius
  double zz=0.; // height, z coord, positive
  switch (gg_->coordKind()) {
  case GYOTO_COORDKIND_SPHERICAL:
    rcyl = coord[1]*sin(coord[2]);
    zz   = fabs(coord[1]*cos(coord[2]));
    break;
  case GYOTO_COORDKIND_CARTESIAN:
    rcyl = pow(coord[1]*coord[1]+coord[2]*coord[2], 0.5);
    zz   = fabs(coord[3]);
    break;
  default:
    GYOTO_ERROR("In Jet::operator(): Unknown coordinate system kind");
  }

  // if (fabs(zz) < jetBaseHeight_) return 1.; // outside jet
  
  // double rcyljetout = fabs(zz*tan(jetOuterOpeningAngle_)),
  //   rcyljetin = fabs(zz*tan(jetInnerOpeningAngle_));

  // //double radbase=10.; // TO BE PARAM
  // //double rcyljetout = radbase
  // //  + fabs((zz-jetBaseHeight_)*tan(jetOuterOpeningAngle_)),
  // //  rcyljetin = radbase + fabs((zz-jetBaseHeight_)*tan(jetInnerOpeningAngle_));
  
  // if  ((rcyl <  rcyljetout) and (rcyl >  rcyljetin)) return -1.; // inside jet
  // else return 1.; // outside jet
  // //cout << "r, rjet, z, theta0, ht= " << rcyl << " " << rcyljet << " " << zz << " " << theta0 << " " << ht << endl;

  double hjetin = rcyl/tan(jetInnerOpeningAngle_),
    hjetout = rcyl/tan(jetOuterOpeningAngle_);

  // double distance = (zz - hjetin)*(zz - hjetout); // <0 inside jet sheath

  // if (zz < jetBaseHeight_) // remove part below jet basis
  //   distance = fabs(distance) + (jetBaseHeight_ - zz);

  //double distance = (zz - jetBaseHeight_ - hjetin)
  //  *(zz - jetBaseHeight_ - hjetout); // <0 inside jet sheath

  double rcyljetin = zz*tan(jetInnerOpeningAngle_),
    rcyljetout = zz*tan(jetOuterOpeningAngle_);

  double distance = (rcyl - rcyljetin)*(rcyl - rcyljetout);
  if (zz < jetBaseHeight_) // remove part below jet basis
    distance = fabs(distance) + (jetBaseHeight_ - zz);

  return distance;
  
}

void Jet::getVelocity(double const pos[4], double vel[4])
{
  double rr = pos[1];
  double Vjet = sqrt(gammaJet_*gammaJet_-1.)/gammaJet_;
  //double Vr = Vjet/(sqrt(gg_->gmunu(pos,1,1)
  //			 +jetVphiOverVr_*jetVphiOverVr_/(rr*rr)*gg_->gmunu(pos,3,3)));
  //double Vphi = jetVphiOverVr_/rr*Vr;
  //cout << "NEW STUFF" << endl;
  //cout << "V2= " << gg_->gmunu(pos,1,1) * Vr*Vr + gg_->gmunu(pos,3,3) * Vphi*Vphi<< " " << (gammaJet_*gammaJet_-1.)/(gammaJet_*gammaJet_) << endl;

  // KerrBL-specific part -- to generalize if possible
  double gpp = gg_->gmunu(pos,3,3), gtt = gg_->gmunu(pos,0,0),
    grr = gg_->gmunu(pos,1,1),
    gtp = gg_->gmunu(pos,0,3);
  double utZAMO = sqrt(-gpp/(gtt*gpp-gtp*gtp)),
    uphiZAMO = -utZAMO*gtp/gpp;
  double Vphi = jetVphiOverVr_*Vjet/sqrt(gpp),
    Vr = sqrt(1-jetVphiOverVr_*jetVphiOverVr_)*Vjet / sqrt(grr);
  //cout << "ZAMO=" << gtt*utZAMO*utZAMO + 2*gtp*utZAMO*uphiZAMO + gpp*uphiZAMO*uphiZAMO << endl;

  // Paper def:
  vel[0] = gammaJet_*utZAMO;
  vel[1] = -gammaJet_*Vr;
  vel[2] = 0.;
  //vel[3] = gammaJet_*uphiZAMO;
  vel[3] = gammaJet_*(uphiZAMO+Vphi);

  double u2 = gg_->ScalarProd(pos,vel,vel);
  double tol = 1e-6;
  if (fabs(u2+1)>tol) throwError("In Jett::getVelo: bad jet velocity");

  // TEST
  //vel[0] = 1.;
  //vel[1] = 0.;
  //vel[2] = 0.;
  //vel[3] = 0.;
  // June 2019 def:
  /*double gammar = 1., gammaphi = 1.5,
    ur = sqrt(gammar*gammar-1.)/gammar,
    uphi = 1./pos[1]*sqrt(gammaphi*gammaphi-1.)/gammaphi,
    grr = gg_->gmunu(pos,1,1),
    aa = gtt, bb = 2.*gtp*uphi, cc = 1.+gpp*uphi*uphi+grr*ur*ur,
    Delta = bb*bb-4.*aa*cc,
    ut1 = (-bb - sqrt(Delta))/(2.*aa),
    ut2 = (-bb + sqrt(Delta))/(2.*aa);
  vel[0]=ut1;
  vel[1]=ur;
  vel[2]=0.;
  vel[3]=uphi;*/

  //cout << "V2= " << gg_->gmunu(pos,1,1)*Vr*Vr + gg_->gmunu(pos,3,3)*Vphi*Vphi << endl;
  //cout << "u2 = " << gg_->ScalarProd(pos,vel,vel) << endl;
}

bool Jet::isThreadSafe() const {
  return Standard::isThreadSafe()
    && (!spectrumKappaSynch_ || spectrumKappaSynch_->isThreadSafe())
    && (!spectrumThermalSynch_ || spectrumThermalSynch_->isThreadSafe());
}

void Jet::metric(SmartPointer<Metric::Generic> gg) {
  if (gg_) gg_->unhook(this);
  string kin = gg->kind();
  //if (kin != "KerrBL" or kin!="NumericalMetricLorene")
  //  GYOTO_ERROR
  //    ("Jet::metric(): metric must be KerrBL");
  // NB: KerrBL needed for ZAMO velocity in getVelocity,
  // could be generalized if needed
  Generic::metric(gg);
}
