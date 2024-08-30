/*
    Copyright 2017-2024 Frederic Vincent, Thibaut Paumard, Paloma Thévenet

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
GYOTO_PROPERTY_BOOL(Jet, Parabolic, NonParabolic, parabolic,
		    "Chose whether the jet sheath has a parabolic shape or if not it will be a conical shape.")
GYOTO_PROPERTY_BOOL(Jet, Outflowing, NonOutflowing, outflowing,
		    "Is the jet outflowing or inflowing?")
GYOTO_PROPERTY_DOUBLE(Jet, JetStagnationRadius, jetStagnationRadius, "The jet is outflowing above this radius, inflowing below")
GYOTO_PROPERTY_DOUBLE(Jet, JetShapeInnerParabolaParam, jetShapeInnerParabolaParam, "Such that the jet (sheath) inner boundary's shape follows: z = jetShapeInnerParabolaParam_ * rcyl^2, where rcyl is the cylindrical radius; the smaller this param, the more open the jet; for parabolic jet")
GYOTO_PROPERTY_DOUBLE(Jet, JetShapeOuterParabolaParam, jetShapeOuterParabolaParam, "Such that the jet (sheath) outer boundary's shape follows: z = jetShapeOuterParabolaParam_ * rcyl^2, where rcyl is the cylindrical radius; the smaller this param, the more open the jet; for parabolic jet")
GYOTO_PROPERTY_DOUBLE(Jet, JetInnerOpeningAngle, jetInnerOpeningAngle,"Jet sheath wall inner BL-theta opening angle; for conical jet")
GYOTO_PROPERTY_DOUBLE(Jet, JetOuterOpeningAngle, jetOuterOpeningAngle,"Jet sheath wall outer BL-theta opening angle; for conical jet")
GYOTO_PROPERTY_DOUBLE(Jet, JetInnerRadius, jetInnerRadius, "Jet basis")
GYOTO_PROPERTY_DOUBLE(Jet, GammaJet, gammaJet, "Constant jet Lorentz factor")
GYOTO_PROPERTY_DOUBLE(Jet, JetVphiOverVr, jetVphiOverVr,"this is V^(phi)/V^(r), ratio of phi to r component of jet velocity as observed by ZAMO, expressed in an orthonormal basis")
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
// chose parabolic or not (=conical)
void Jet::parabolic(bool parabol) {parabolic_=parabol;}
bool Jet::parabolic() const { return parabolic_; }
// chose outflowing or not (=inflowing)
void Jet::outflowing(bool out) {outflowing_=out;}
bool Jet::outflowing() const { return outflowing_; }
// for parabolic:
void Jet::jetShapeInnerParabolaParam(double param) {
  jetShapeInnerParabolaParam_=param;
  if (jetShapeOuterParabolaParam_ != -1.){ // default parabola value is -1
    if (jetShapeOuterParabolaParam_ > jetShapeInnerParabolaParam_)
      GYOTO_ERROR("The outer parabola parameter cannot be bigger "
		  "than the inner one");
  }  
}
double Jet::jetShapeInnerParabolaParam() const {return jetShapeInnerParabolaParam_;}
void Jet::jetShapeOuterParabolaParam(double param) {
  jetShapeOuterParabolaParam_=param;
  if (jetShapeInnerParabolaParam_ != -1.){ // default parabola value is -1
    if (jetShapeOuterParabolaParam_ > jetShapeInnerParabolaParam_)
      GYOTO_ERROR("The outer parabola parameter cannot be bigger "
		  "than the inner one");
  }  
}
double Jet::jetShapeOuterParabolaParam() const {return jetShapeOuterParabolaParam_;}
// for conical:
void Jet::jetOuterOpeningAngle(double ang) {
  jetOuterOpeningAngle_=ang;
  if (jetInnerOpeningAngle_ != -1.){ // default angle value is -1
    if (jetOuterOpeningAngle_ < jetInnerOpeningAngle_)
    GYOTO_ERROR("The outer opening angle cannot be smaller "
		"than the inner opening angle");
  }
}
double Jet::jetOuterOpeningAngle()const{return jetOuterOpeningAngle_;}
void Jet::jetInnerOpeningAngle(double ang) {
  jetInnerOpeningAngle_=ang;
  if (jetOuterOpeningAngle_ != -1.){ // default angle value is -1
    if (jetOuterOpeningAngle_ < jetInnerOpeningAngle_)
      GYOTO_ERROR("The outer opening angle cannot be smaller "
		  "than the inner opening angle");
  }
}
double Jet::jetInnerOpeningAngle()const{return jetInnerOpeningAngle_;}
// for both:
void Jet::jetStagnationRadius(double stag) {jetStagnationRadius_=stag;}
double Jet::jetStagnationRadius() const {return jetStagnationRadius_;}
void Jet::jetInnerRadius(double hh) {jetInnerRadius_=hh;}
double Jet::jetInnerRadius()const{return jetInnerRadius_;}
void Jet::jetVphiOverVr(double alpha) {jetVphiOverVr_=alpha;}
double Jet::jetVphiOverVr()const{return jetVphiOverVr_;}
void Jet::gammaJet(double gam) {gammaJet_=gam;}
double Jet::gammaJet()const{return gammaJet_;}
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
  Standard("Jet"),
  parabolic_(1), outflowing_(1), jetShapeInnerParabolaParam_(-1.),
  jetShapeOuterParabolaParam_(-1.),
  jetOuterOpeningAngle_(-1.),
  jetInnerOpeningAngle_(-1.), jetInnerRadius_(2.),
  jetStagnationRadius_(0.),
  gammaJet_(1.), jetVphiOverVr_(0.),
  baseNumberDensity_cgs_(1.), baseTemperature_(1e10),
  temperatureSlope_(1.), magneticConfig_("None"),
  magnetizationParameter_(1.)
{
  GYOTO_DEBUG << endl;
  spectrumKappaSynch_   = new Spectrum::KappaDistributionSynchrotron();
  spectrumKappaSynch_->kappaindex(default_kappaindex);
  spectrumThermalSynch_ = new Spectrum::ThermalSynchrotron();
}

Jet::Jet(const Jet& o) :
  Standard(o),
  parabolic_(o.parabolic_),
  outflowing_(o.outflowing_),
  jetShapeInnerParabolaParam_(o.jetShapeInnerParabolaParam_),
  jetShapeOuterParabolaParam_(o.jetShapeOuterParabolaParam_),
  jetOuterOpeningAngle_(o.jetOuterOpeningAngle_),
  jetInnerOpeningAngle_(o.jetInnerOpeningAngle_),
  jetInnerRadius_(o.jetInnerRadius_),
  jetStagnationRadius_(o.jetStagnationRadius_),
  gammaJet_(o.gammaJet_), jetVphiOverVr_(o.jetVphiOverVr_),
  baseNumberDensity_cgs_(o.baseNumberDensity_cgs_),
  baseTemperature_(o.baseTemperature_),
  temperatureSlope_(o.temperatureSlope_),
  magnetizationParameter_(o.magnetizationParameter_),
  magneticConfig_(o.magneticConfig_),
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
  double rr=0.; // spherical radius
  switch (gg_->coordKind()) {
  case GYOTO_COORDKIND_SPHERICAL:
    rr = coord_ph[1];
    rcyl = coord_ph[1]*sin(coord_ph[2]);
    zz   = coord_ph[1]*cos(coord_ph[2]);
    break;
  case GYOTO_COORDKIND_CARTESIAN:
    rcyl = pow(coord_ph[1]*coord_ph[1]+coord_ph[2]*coord_ph[2], 0.5);
    zz   = coord_ph[3];
    rr = pow(coord_ph[1]*coord_ph[1]+coord_ph[2]*coord_ph[2]
	     +coord_ph[3]*coord_ph[3], 0.5);
    break;
  default:
    GYOTO_ERROR("In Jet::radiativeQ: Unknown coordinate system kind");
  }

  double rcyljetbase = 0.;
  if (parabolic_==0){ // conical jet
    rcyljetbase = jetInnerRadius_*tan(jetOuterOpeningAngle_);
    //cout << "rcyl stuff: " << jetInnerRadius_ << " " << tan(jetOuterOpeningAngle_) << " " << rcyljetbase << endl;
  }else{ // parabolic jet
    // We want the cylindrical radius of the jet base, knowing that the
    // jet base is at spherical radius jetInnerRadius_, and that
    // abs(z) = jetShapeOuterParabolaParam_ * rcyl^2.
    // We have: jetInnerRadius_^2 = rcyljetbase^2 + zjetbase^2
    //                     = zjetbase/jetShapeOuterParabolaParam_ + zjetbase^2
    // So we get a second order equation with discriminant:
    // Delta = 1/jetShapeOuterParabolaParam_^2 + 4*jetInnerRadius_^2
    double discri = 1./(jetShapeOuterParabolaParam_*jetShapeOuterParabolaParam_) + 4*jetInnerRadius_*jetInnerRadius_;
    if (discri<0) GYOTO_ERROR("Bad discriminant!");
    double  zjetbase = 0.5*(-1/jetShapeOuterParabolaParam_ + sqrt(discri));
    rcyljetbase = sqrt(zjetbase/jetShapeOuterParabolaParam_);
  }

  //cout << "In emission t, rcyl, zz= " << coord_ph[0]<< " " << rcyl << " " << zz << endl;

  //rcyl=rcyljetbase;
  //zz=2.; // TEST!!!
  
  double number_density = baseNumberDensity_cgs_
    *(rcyljetbase*rcyljetbase)/(rcyl*rcyl); // using cylindrical radius here by conservation of mass in jet

  double temperature = baseTemperature_*pow(jetInnerRadius_/rr,
					    temperatureSlope_);
  // NB: Vincent+19 torus+jet Sgr model considers a T(z) rather than
  // a T(r), see below. I don't see why (in 2024), so I prefer to consider
  // a T(r) power law. The difference on the image is small anyway.
  //baseTemperature_*pow(jetInnerRadius_/fabs(zz),
  //		    temperatureSlope_); // 2019 version

  //cout << "ne T= " << jetInnerRadius_ << " " << jetOuterOpeningAngle_ << " " << baseNumberDensity_cgs_ << " " << rcyljetbase << " " << rcyl << " " << number_density << " " << temperature << endl;

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
  double rr=0.; // spherical radius
  switch (gg_->coordKind()) {
  case GYOTO_COORDKIND_SPHERICAL:
    rr = coord[1];
    rcyl = coord[1]*sin(coord[2]);
    zz   = fabs(coord[1]*cos(coord[2]));
    break;
  case GYOTO_COORDKIND_CARTESIAN:
    rr = pow(coord[1]*coord[1]+coord[2]*coord[2]+coord[3]*coord[3], 0.5);
    rcyl = pow(coord[1]*coord[1]+coord[2]*coord[2], 0.5);
    zz   = fabs(coord[3]);
    break;
  default:
    GYOTO_ERROR("In Jet::operator(): Unknown coordinate system kind");
  }

  /*
         *** Jet sheath geometry ***

	 *** 1. Conical

	 In the 2D plane (cylindrical radius rcyl, altitude above equat plane z)
	 we consider two lines: one with BL angle theta=jetInnerOpeningAngle_,
	 and one with BL angle theta=jetOuterOpeningAngle_. Ther oigin is at
	 (0,0) in this plane. The jet sheath is in between these two lines,
	 removing the region with BL r<jetInnerRadius_ (typically chosen
	 at the horizon). The distance function implemented here isas always
	 such that distance<0 inside the sheath. So considering some point 
	 P(rcyl,z) along a geodesic, we return:
	       distance = (rcyl - rcyljetin)*(rcyl - rcyljetout),
	 where rcyljetin is the cylindrical radius of the inner jet sheath 
	 boundary at the altiude z of point P, and idem for the outer one.
	 This distance is negative iff inside the sheath, 
	 and below jetInnerRadius_ we switch to simply
	       distance(r<jetInnerRadius_) = jetInnerRadius_ - r
	 which is obviously positive and thus excluded.
	 The distance function is not continuous at r=jetInnerRadius_
	 but we dont care.

	 *** 2. Parabolic

	 Exactly similar but we consider parabolas 
	      z = jetShapeInnerParabolaParam_ * rcyl²
	 and 
	      z = jetShapeOuterParabolaParam_ * rcyl²
	 so the jet sheath is in between the two parabolas.
	 The distance function is defined in the exact same way.
   */

  double rcyljetin=-1., rcyljetout=-1.; // stupid initialization for debugging
  if (parabolic_==true){
    rcyljetin = sqrt(zz/jetShapeInnerParabolaParam_);
    rcyljetout = sqrt(zz/jetShapeOuterParabolaParam_);
    //cout << "rcyl in out= " << zz << " " << rcyljetin << " " << rcyljetout << endl;
  }else{
    rcyljetin = zz*tan(jetInnerOpeningAngle_);
    rcyljetout = zz*tan(jetOuterOpeningAngle_);
  }

  //if (rr<5.)
  //  cout << "at point " << rcyl << " " << zz << " the jet boundaries are " << rcyljetin << " " << rcyljetout << endl;

  double distance = (rcyl - rcyljetin)*(rcyl - rcyljetout);
 
  if (rr < jetInnerRadius_) // remove part below jet basis
    distance = jetInnerRadius_ - rr; 

  return distance;
  
}

void Jet::getVelocity(double const pos[4], double vel[4])
{
  double rr = pos[1], theta = pos[2], rcyl = rr*sin(theta); // MAKE IT INDEP OF COORD KIND
  double Vjet = sqrt(gammaJet_*gammaJet_-1.)/gammaJet_;
  //double Vr = Vjet/(sqrt(gg_->gmunu(pos,1,1)
  //			 +jetVphiOverVr_*jetVphiOverVr_/(rr*rr)*gg_->gmunu(pos,3,3)));
  //double Vphi = jetVphiOverVr_/rr*Vr;
  //cout << "NEW STUFF" << endl;
  //cout << "V2= " << gg_->gmunu(pos,1,1) * Vr*Vr + gg_->gmunu(pos,3,3) * Vphi*Vphi<< " " << (gammaJet_*gammaJet_-1.)/(gammaJet_*gammaJet_) << endl;

  // KerrBL-specific part -- to generalize if possible
  double gpp = gg_->gmunu(pos,3,3), gtt = gg_->gmunu(pos,0,0),
    gthth = gg_->gmunu(pos,2,2),
    grr = gg_->gmunu(pos,1,1),
    gtp = gg_->gmunu(pos,0,3);
  double utZAMO = sqrt(-gpp/(gtt*gpp-gtp*gtp)),
    uphiZAMO = -utZAMO*gtp/gpp;

  // Lets now define the velo of the jet V as measured by the ZAMO
  double Vphi=0., Vr=0., Vth=0.;
  if (parabolic_==true){ // See FV 2024 notes for details on equations.
    double sz = 1., zz = rr*cos(theta);
    if (zz<0.) sz = -1.; // sz is the sign of z.
    // rvel and thvel are the components of the velocity vector along the
    // parabolic jet sheath in the absence of phi motion, in the orthonormal
    // (e_r = \partial_r/sqrt(grr), e_theta = \partial_th/sqrt(gthth)) frame
    // associated to BL coordinates
    double rvel = sin(theta) + 2.*sz*jetShapeInnerParabolaParam_*rcyl*cos(theta),
      thvel = cos(theta) - 2.*sz*jetShapeInnerParabolaParam_*rcyl*sin(theta);
    // See FV 2024 notes for the sz trick that allows to deal with both
    // z>0 and z<0 parts of jet
    // Unitary vector u along parabolic sheath in the absence of phi motion
    double ur = rvel/sqrt(rvel*rvel + thvel*thvel),
      uth = thvel/sqrt(rvel*rvel + thvel*thvel);
    // so: u = u^(r) e_r + u^(th) e_theta is a unit vector along the velocity
    // field lines. The parentheses (r) and (theta) remind that we are
    // in the orthonormal basis.

    // Let's add a phi degree of liberty to the motion. So we go from
    // the unit vector u to a unit vector u':
    // u' = (u + alpha*e_phi) / sqrt(1+alpha^2)
    //    = u^(r)/sqrt(1+alpha^2) e_r + u^th/sqrt(1+alpha^2) e_th
    //                                + alpha/sqrt(1+alpha^2) e_phi
    // We want to parametrize this by the ratio:
    // V^(phi) / V^(r) = alpha / u^(r)
    // Below I call u^(phi) = alpha.
    double uph = jetVphiOverVr_ * ur;
    // Then the unit vector u' reads in components:
    double tmp = sqrt(1+uph*uph), uprimer = ur/tmp, uprimeth = uth/tmp,
      uprimeph = uph/tmp;
    // Now the velocity reads, in the coordinate basis (\partial_mu)
    
    Vr = Vjet*uprimer/sqrt(grr);
    Vth = Vjet*uprimeth/sqrt(gthth);
    Vphi = Vjet*uprimeph/sqrt(gpp);
    
    // It is obvious to see that if jetVphiOverVr_=0, we find the
    // correct velocity in the orthonormal basis: Vjet * (ur e_r + uth e_th)
    //cout << "V sum= " << Vr*Vr/grr + Vth*Vth/gthth + Vphi*Vphi/gpp << " " << Vjet*Vjet << endl;
  }else{
    // Here in the absence of phi motion, we simply want a velo along e_r,
    // so: V = Vjet e_r, and adding a phi dependence leads to
    // V = V^(r) e_r + V^(phi) e_phi, and we parametrize using
    // alpha = V^(phi)/V^(r), with (V^(r))^2 + (V^(phi))^2 = 1;
    // These relations immediately lead to (keeping in mind that
    // V^r = V^(r)/sqrt(grr) and similarly for phi):
    
    Vr = Vjet / (sqrt(grr) * sqrt(1+jetVphiOverVr_*jetVphiOverVr_));
    Vphi = jetVphiOverVr_/sqrt(1+jetVphiOverVr_*jetVphiOverVr_) \
      * Vjet/sqrt(gpp);
    
    // If jetVphiOverVr_=0, V is purely radial as it should for conical case
  }

  // The 4-velo is then defined by the standard formula:
  // u = gamma * (u_zamo + V):
  
  vel[0] = gammaJet_*utZAMO;
  if (outflowing_==1){ // outflow
    if (rr > jetStagnationRadius_){ // outflow above stagnation
                                    //(at zero radius by default,
                                    //so pure outflow if not set by user)
      vel[1] = gammaJet_*Vr;
      vel[2] = gammaJet_*Vth;
    }else{ // inflow below stagnation
      vel[1] = -gammaJet_*Vr;
      vel[2] = -gammaJet_*Vth;
    }
  }else{ // pure inflow
    vel[1] = -gammaJet_*Vr;
    vel[2] = -gammaJet_*Vth;
  }
  vel[3] = gammaJet_*(uphiZAMO+Vphi);

  double u2 = gg_->ScalarProd(pos,vel,vel);
  double tol = 1e-6;
  //cout << "vel norm= " << fabs(u2+1) << endl;
  if (fabs(u2+1)>tol or u2!=u2) throwError("In Jett::getVelo: bad jet velocity");

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

void Jet::radiativeQ(double *Inu, double *Qnu, double *Unu,
           double *Vnu,
           Eigen::Matrix4d *Onu,
           double const *nuem , size_t nbnu,
           double dsem,
           state_t const &coord_ph,
           double const *co) const {
  
  // polarized radiativeQ
  double rcyl=0.; // cylindrical radius
  double zz=0.; // height, z coord
  double rr=0.; // spherical radius  
  switch (gg_->coordKind()) {
  case GYOTO_COORDKIND_SPHERICAL:
    rr = coord_ph[1];
    rcyl = coord_ph[1]*sin(coord_ph[2]);
    zz   = coord_ph[1]*cos(coord_ph[2]);
    break;
  case GYOTO_COORDKIND_CARTESIAN:
    rcyl = pow(coord_ph[1]*coord_ph[1]+coord_ph[2]*coord_ph[2], 0.5);
    zz   = coord_ph[3];
    rr = pow(coord_ph[1]*coord_ph[1]+coord_ph[2]*coord_ph[2]
	     +coord_ph[3]*coord_ph[3], 0.5);
    break;
  default:
    GYOTO_ERROR("In Jet::radiativeQ: Unknown coordinate system kind");
  }

  double vel[4]; // 4-velocity of emitter
  for (int ii=0;ii<4;ii++){
    vel[ii]=co[ii+4];
  }

  double rcyljetbase = 0.;
  if (parabolic_==0){ // conical jet
    rcyljetbase = jetInnerRadius_*tan(jetOuterOpeningAngle_);
    //cout << "rcyl stuff: " << jetInnerRadius_ << " " << tan(jetOuterOpeningAngle_) << " " << rcyljetbase << endl;
  }else{ // parabolic jet
    // We want the cylindrical radius of the jet base, knowing that the
    // jet base is at spherical radius jetInnerRadius_, and that
    // abs(z) = jetShapeOuterParabolaParam_ * rcyl^2.
    // We have: jetInnerRadius_^2 = rcyljetbase^2 + zjetbase^2
    //                     = zjetbase/jetShapeOuterParabolaParam_ + zjetbase^2
    // So we get a second order equation with discriminant:
    // Delta = 1/jetShapeOuterParabolaParam_^2 + 4*jetInnerRadius_^2
    double discri = 1./(jetShapeOuterParabolaParam_*jetShapeOuterParabolaParam_) + 4*jetInnerRadius_*jetInnerRadius_;
    if (discri<0) GYOTO_ERROR("Bad discriminant!");
    double  zjetbase = 0.5*(-1/jetShapeOuterParabolaParam_ + sqrt(discri));
    rcyljetbase = sqrt(zjetbase/jetShapeOuterParabolaParam_);
  }

  //cout << "rcyl, rcylB, zz= " << rcyl << " " << rcyljetbase << " " << zz << endl;

  //rcyl=rcyljetbase;
  //zz=2.; // TEST!!!
  
  double number_density = baseNumberDensity_cgs_
    *(rcyljetbase*rcyljetbase)/(rcyl*rcyl);

  double temperature = baseTemperature_*pow(jetInnerRadius_/rr,
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

  // CHOOSE BFIELD GEOMETRY
  // Note: Bfield is simply a normalized spacelike vector, we only
  // need its direction, the norm is computed independently.

  double B4vect[4]={0.,0.,0.,0.};
  computeB4vect(B4vect, magneticConfig_, co, coord_ph);
    
  //cout << "***BL mf in orthonorm frame= "<< B4vect[0]  << " " << B4vect[1]  << " " << coord_ph[1]*B4vect[2]  << " " << coord_ph[1]*abs(sin(coord_ph[2]))*B4vect[3]  << endl;

  //cout << "***k photon in orthonorm frame= " << coord_ph[4] << " " << coord_ph[5] << " " << coord_ph[1]*coord_ph[6] << " " << coord_ph[1]*abs(sin(coord_ph[2]))*coord_ph[7] << endl;
  
  double norm=sqrt(gg_->ScalarProd(&coord_ph[0], B4vect, B4vect));
  //cout << "norm= " << norm << endl;
  
  if (fabs(norm-1.)>GYOTO_DEFAULT_ABSTOL) GYOTO_ERROR("Bad mf normalization");

  double Chi=getChi(B4vect, coord_ph, vel); // this is EVPA
  //cout << "At r,x,y,z= " << coord_ph[1] << " " << coord_ph[1]*sin(coord_ph[2])*cos(coord_ph[3]) << " " << coord_ph[1]*sin(coord_ph[2])*sin(coord_ph[3]) << " " << coord_ph[1]*cos(coord_ph[2]) << " ; Chi=" << Chi << endl;

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
  double theta_mag = acos(lscalb/(lnorm*bnorm));

  if (theta_mag<0. or theta_mag>M_PI) throwError("Jet: bad B angle");

  
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

  if (kappaIndex()!=default_kappaindex){
    // KAPPA-DISTRIB SYNCHROTRON
    spectrumKappaSynch_->kappaindex(kappaIndex());
    spectrumKappaSynch_->numberdensityCGS(number_density);
    spectrumKappaSynch_->angle_averaged(0);
    spectrumKappaSynch_->angle_B_pem(theta_mag);
    spectrumKappaSynch_->cyclotron_freq(nu0);
    spectrumKappaSynch_->thetae(thetae);
    double hypergeom = Gyoto::hypergeom(kappaIndex(), thetae);
    spectrumKappaSynch_->hypergeometric(hypergeom);

    spectrumKappaSynch_->radiativeQ(jInu, jQnu, jUnu, jVnu,
				    aInu, aQnu, aUnu, aVnu,
				    rotQnu, rotUnu, rotVnu, nuem, nbnu);
  }else{
    // THERMAL SYNCHROTRON
    double besselK2 = bessk(2, 1./thetae);
    spectrumThermalSynch_->temperature(temperature);
    spectrumThermalSynch_->numberdensityCGS(number_density);
    spectrumThermalSynch_->angle_averaged(0); //  no angle avg of course
    spectrumThermalSynch_->angle_B_pem(theta_mag);
    spectrumThermalSynch_->cyclotron_freq(nu0);
    spectrumThermalSynch_->besselK2(besselK2);
    spectrumThermalSynch_->radiativeQ(jInu, jQnu, jUnu, jVnu,
				      aInu, aQnu, aUnu, aVnu,
				      rotQnu, rotUnu, rotVnu, nuem, nbnu);
  }

  // RETURNING TOTAL INTENSITY AND TRANSMISSION
  for (size_t ii=0; ii<nbnu; ++ii) {

    //cout << "at r,th= " << coord_ph[1] << " " << coord_ph[2] << endl;
    //cout << "at rcyl,z= " << rcyl << " " << zz << endl;
    //cout << "jet jnu anu kappa= " << jnu_tot << " " << anu_tot << endl; //x" " << jnu_tot/anu_tot << " " << dsem << endl;
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

    //cout << "in jet stuff: " << zz << " " << rcyl << " " << nu_ems[0]  << " " << number_density << " " << nu0 << " " << temperature << " " << thetae << " " << jnu_tot << " " << anu_tot << " " << dsem << endl;

    

    if (Inu[ii]<0.)
      GYOTO_ERROR("In Jet::radiativeQ: Inu<0");
    if (Inu[ii]!=Inu[ii] or Onu[ii](0,0)!=Onu[ii](0,0))
      GYOTO_ERROR("In Jet::radiativeQ: Inu or Taunu is nan");
    if (Inu[ii]==Inu[ii]+1. or Onu[ii](0,0)==Onu[ii](0,0)+1.)
      GYOTO_ERROR("In Jet::radiativeQ: Inu or Taunu is infinite");

  }
}

void Jet::magneticConfiguration(string config){
  magneticConfig_=config;
}

string Jet::magneticConfiguration() const{
  return magneticConfig_;
}
