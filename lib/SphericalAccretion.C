/*
    Copyright 2021 Frederic Vincent & Thibaut Paumard

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
#include "GyotoSphericalAccretion.h"
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
#include <random>

using namespace std;
using namespace Gyoto;
using namespace Gyoto::Astrobj;

GYOTO_PROPERTY_START(SphericalAccretion)
GYOTO_PROPERTY_DOUBLE(SphericalAccretion, SphericalAccretionInnerRadius, sphericalAccretionInnerRadius)
GYOTO_PROPERTY_DOUBLE_UNIT(SphericalAccretion,
			   NumberDensityAtInnerRadius,
			   numberDensityAtInnerRadius)
GYOTO_PROPERTY_DOUBLE(SphericalAccretion, DensitySlope, densitySlope)
GYOTO_PROPERTY_DOUBLE(SphericalAccretion,
		      TemperatureAtInnerRadius,
		      temperatureAtInnerRadius)
GYOTO_PROPERTY_DOUBLE(SphericalAccretion, TemperatureSlope, temperatureSlope)
GYOTO_PROPERTY_DOUBLE(SphericalAccretion, MagnetizationParameter,
		      magnetizationParameter)
GYOTO_PROPERTY_BOOL(SphericalAccretion, UseSelfAbsorption, NoUseSelfAbsorption,
		    useSelfAbsorption)
GYOTO_PROPERTY_END(SphericalAccretion, Standard::properties)

// ACCESSORS
bool SphericalAccretion::useSelfAbsorption() const {return use_selfabsorption_;}
void SphericalAccretion::useSelfAbsorption(bool abs) {use_selfabsorption_=abs;}
void SphericalAccretion::sphericalAccretionInnerRadius(double hh) {sphericalAccretionInnerRadius_=hh;}
double SphericalAccretion::sphericalAccretionInnerRadius()const{return sphericalAccretionInnerRadius_;}
double SphericalAccretion::numberDensityAtInnerRadius() const {
  // Converts internal cgs central enthalpy to SI
  double dens=numberDensityAtInnerRadius_cgs_;
# ifdef HAVE_UDUNITS
  dens = Units::Converter("cm-3", "m-3")(dens);
# else
  GYOTO_WARNING << "Units ignored, please recompile Gyoto with --with-udunits"
		<< endl ;
# endif
  return dens; }
double SphericalAccretion::numberDensityAtInnerRadius(string const &unit) const
{
  double dens = numberDensityAtInnerRadius();
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
void SphericalAccretion::numberDensityAtInnerRadius(double dens) {
# ifdef HAVE_UDUNITS
  dens = Units::Converter("m-3", "cm-3")(dens);
# else
  GYOTO_WARNING << "Units ignored, please recompile Gyoto with --with-udunits"
		<< endl ;
# endif
  numberDensityAtInnerRadius_cgs_=dens;
}
void SphericalAccretion::numberDensityAtInnerRadius(double dens, string const &unit) {
  if (unit != "") {
# ifdef HAVE_UDUNITS
    dens = Units::Converter(unit, "m-3")(dens);
# else
    GYOTO_WARNING << "Units ignored, please recompile Gyoto with --with-udunits"
		  << endl ;
# endif
  }
  numberDensityAtInnerRadius(dens);
}
void SphericalAccretion::densitySlope(double ss) {densitySlope_=ss;}
double SphericalAccretion::densitySlope()const{return densitySlope_;}
void SphericalAccretion::temperatureAtInnerRadius(double tt) {temperatureAtInnerRadius_=tt;}
double SphericalAccretion::temperatureAtInnerRadius()const{return temperatureAtInnerRadius_;}
void SphericalAccretion::temperatureSlope(double ss) {temperatureSlope_=ss;}
double SphericalAccretion::temperatureSlope()const{return temperatureSlope_;}
void SphericalAccretion::magnetizationParameter(double rr) {
  magnetizationParameter_=rr;}
double SphericalAccretion::magnetizationParameter()const{
  return magnetizationParameter_;}
//

SphericalAccretion::SphericalAccretion() :
  Standard("SphericalAccretion"),
  sphericalAccretionInnerRadius_(2.),
  numberDensityAtInnerRadius_cgs_(1.), temperatureAtInnerRadius_(1e10),
  temperatureSlope_(1.),
  densitySlope_(2.),
  magnetizationParameter_(1.),
  use_selfabsorption_(1)
{
  GYOTO_DEBUG << endl;
  spectrumThermalSynch_ = new Spectrum::ThermalSynchrotron();
}

SphericalAccretion::SphericalAccretion(const SphericalAccretion& o) :
  Standard(o),
  sphericalAccretionInnerRadius_(o.sphericalAccretionInnerRadius_),
  numberDensityAtInnerRadius_cgs_(o.numberDensityAtInnerRadius_cgs_),
  temperatureAtInnerRadius_(o.temperatureAtInnerRadius_),
  temperatureSlope_(o.temperatureSlope_),
  densitySlope_(o.densitySlope_),
  magnetizationParameter_(o.magnetizationParameter_),
  spectrumThermalSynch_(NULL),
  use_selfabsorption_(o.use_selfabsorption_)
{
  GYOTO_DEBUG << endl;
  if (gg_) gg_->hook(this);
  if (o.spectrumThermalSynch_()) spectrumThermalSynch_=o.spectrumThermalSynch_->clone();

}
SphericalAccretion* SphericalAccretion::clone() const
{ return new SphericalAccretion(*this); }

SphericalAccretion::~SphericalAccretion() {
  GYOTO_DEBUG << endl;
  if (gg_) gg_->unhook(this);
}

void SphericalAccretion::radiativeQ(double Inu[], // output
		     double Taunu[], // output
		     double const nu_ems[], size_t nbnu, // input
		     double dsem,
		     state_t const &coord_ph,
		     double const coord_obj[8]) const {
  
  double rr, rcyl, theta, zz=0.;
  switch (gg_->coordKind()) {
  case GYOTO_COORDKIND_SPHERICAL:
    rr = coord_ph[1];
    rcyl = coord_ph[1]*sin(coord_ph[2]);
    theta = coord_ph[2];
    zz   = coord_ph[1]*cos(coord_ph[2]);
    break;
  case GYOTO_COORDKIND_CARTESIAN:
    rcyl = pow(coord_ph[1]*coord_ph[1]+coord_ph[2]*coord_ph[2], 0.5);
    rr = sqrt(coord_ph[1]*coord_ph[1]+coord_ph[2]*coord_ph[2]
	      +coord_ph[3]*coord_ph[3]);
    theta   = acos(coord_ph[3]/rr);
    zz   = coord_ph[3];
    break;
  default:
    GYOTO_ERROR("In SphericalAccretion::radiativeQ(): Unknown coordinate system kind");
  }

  double number_density = numberDensityAtInnerRadius_cgs_
    *pow(sphericalAccretionInnerRadius_/rr, densitySlope_);
    //*(sphericalAccretionInnerRadius_*sphericalAccretionInnerRadius_)/(rr*rr);

  //cout << "Spherical r, z, rho= " << rcyl << " " << zz << " " << number_density << endl;

  double temperature = temperatureAtInnerRadius_
    *pow(sphericalAccretionInnerRadius_/rr, temperatureSlope_);

  double thetae = GYOTO_BOLTZMANN_CGS*temperature
    /(GYOTO_ELECTRON_MASS_CGS*GYOTO_C2_CGS);

  double BB = sqrt(4.*M_PI*magnetizationParameter_
		   *GYOTO_PROTON_MASS_CGS * GYOTO_C_CGS * GYOTO_C_CGS
		   *number_density);

  double nu0 = GYOTO_ELEMENTARY_CHARGE_CGS*BB
    /(2.*M_PI*GYOTO_ELECTRON_MASS_CGS*GYOTO_C_CGS); // cyclotron freq

  //cout << "emission stuff: " << sphericalAccretionInnerRadius_ << " " << rr << " " << number_density << " " << temperature << " " << BB << endl;

 
  // Emission and absorption synchrotron coefs
  double jnu_synch[nbnu], anu_synch[nbnu];
  for (size_t ii=0; ii<nbnu; ++ii){
    // Initializing to <0 value to create errors if not updated
    jnu_synch[ii]=-1.;
    anu_synch[ii]=-1.;
  }

  // THERMAL SYNCHROTRON
  spectrumThermalSynch_->temperature(temperature);
  spectrumThermalSynch_->numberdensityCGS(number_density);
  spectrumThermalSynch_->angle_averaged(1); // impose angle-averaging
  spectrumThermalSynch_->angle_B_pem(0.);   // so we don't care about angle
  spectrumThermalSynch_->cyclotron_freq(nu0);
  double besselK2 = bessk(2, 1./thetae);
  spectrumThermalSynch_->besselK2(besselK2);
  //cout << "besselK: " << besselK2 << endl;
  //cout << "nu passed to synchro= " << nu_ems[0] << endl;
  spectrumThermalSynch_->radiativeQ(jnu_synch,anu_synch,
				    nu_ems,nbnu);

  // RETURNING TOTAL INTENSITY AND TRANSMISSION
  for (size_t ii=0; ii<nbnu; ++ii){

    double jnu_tot = jnu_synch[ii], // 1./(rr*rr),
      anu_tot=0.; // will stay 0 if use_selfabsorption_=False
    if (use_selfabsorption_)
      anu_tot = anu_synch[ii]; // else is updated

    //if (nu_ems[ii]>1e9){
    //  cout << "anu, jnu= " << anu_tot <<  " " << jnu_tot << endl;
    //}

    double em1=std::expm1(-anu_tot * dsem * gg_->unitLength());
    Taunu[ii] = em1+1.;
    Inu[ii] = anu_tot == 0. ? jnu_tot * dsem * gg_->unitLength() :
      -jnu_tot / anu_tot * em1;

    if (Inu[ii]<0.)
      GYOTO_ERROR("In SphericalAccretion::radiativeQ: Inu<0");
    if (Inu[ii]!=Inu[ii] or Taunu[ii]!=Taunu[ii])
      GYOTO_ERROR("In SphericalAccretion::radiativeQ: Inu or Taunu is nan");
    if (Inu[ii]==Inu[ii]+1. or Taunu[ii]==Taunu[ii]+1.)
      GYOTO_ERROR("In SphericalAccretion::radiativeQ: Inu or Taunu is infinite");

  }
}

double SphericalAccretion::operator()(double const coord[4]) {
  // zpos: modulus of altitude above equatorial plane
  // rproj: radius projected in the equatorial plane
  double rr;
  
  switch (gg_ -> coordKind()) {
  case GYOTO_COORDKIND_SPHERICAL:
    rr  = coord[1];
    break;
  case GYOTO_COORDKIND_CARTESIAN:
    rr  = sqrt(coord[1]*coord[1]+coord[2]*coord[2]+coord[3]*coord[3]);
    break;
  default:
    GYOTO_ERROR("SphericalAccretion::operator(): unknown COORDKIND");
  }
  
  return sphericalAccretionInnerRadius_ - rr; // >0 outside, <0 inside 
}

void SphericalAccretion::getVelocity(double const pos[4], double vel[4])
{

  double rr = pos[1];
  double gtt = gg_->gmunu(pos,0,0),
    grr = gg_->gmunu(pos,1,1),
    guptt = gg_->gmunu_up(pos,0,0),
    guptp = gg_->gmunu_up(pos,0,3),
    guprr = gg_->gmunu_up(pos,1,1);

  // 4-vel obtained by imposing: u_t=-1, u_phi=0, u^theta=0
  // see FV notes SphericalVelocity.pdf for details
  vel[0] = -guptt;
  vel[1] = -sqrt((-1.-guptt)*guprr);
  vel[2] = 0;
  vel[3] = -guptp;

  // CHECKS
  
  double tol=1e-4;
  
  // DEBUG: compare to Falcke+00 shadow paper formulas //////////////////
  // double spin=0.8;
  // double delta=rr*rr - 2.*rr + spin*spin;
  // double theta=pos[2];
  // double AA = (rr*rr + spin*spin)*(rr*rr + spin*spin) - spin*spin*delta*sin(theta)*sin(theta);
  // double vr_F00 = -delta*sqrt(2.*rr*(rr*rr + spin*spin))/AA,
  //   Omega_F00 = 2.*spin*rr/AA;
  // if (fabs(vr_F00-vel[1]/vel[0])>tol or fabs(Omega_F00-vel[3]/vel[0])>tol)
  //   throwError("In SphericalAccretion::getVelo different from Falcke+00");
  // --> perfect agreement 210211
  ////////////////////////////////////////////////////////////////////////

  // Check 4vel normalization
  double u2 = gg_->ScalarProd(pos,vel,vel);
  //cout << "4vel,u2= " << rr << " " << pos[2] << " " << gtt << " " << grr << " " << vel[0] << " " << vel[1] << " " << vel[2] << " " << vel[3] << " " << u2 << endl;

  if (fabs(u2+1.)>tol or u2!=u2) {
    cerr << " *** 4-velocity squared norm= " << u2 << endl;
    throwError("In SphericalAccretion: 4vel "
	       "is not properly normalized!");
  }
  
}

bool SphericalAccretion::isThreadSafe() const {
  return Standard::isThreadSafe()
    && (!spectrumThermalSynch_ || spectrumThermalSynch_->isThreadSafe());
}

void SphericalAccretion::metric(SmartPointer<Metric::Generic> gg) {
  if (gg_) gg_->unhook(this);
  string kin = gg->kind();
  //if (kin != "KerrBL" or kin!="NumericalMetricLorene")
  //  GYOTO_ERROR
  //    ("SphericalAccretion::metric(): metric must be KerrBL");
  // NB: KerrBL needed for ZAMO velocity in getVelocity,
  // could be generalized if needed
  Generic::metric(gg);
}
