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
#include "GyotoThickDisk.h"
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

GYOTO_PROPERTY_START(ThickDisk)
GYOTO_PROPERTY_DOUBLE(ThickDisk, ThickDiskOpeningAngle, thickDiskOpeningAngle,
		      "Angle between spin axis and disk surface, so "
		      "if it is pi/2, the disk is razor-thin")
GYOTO_PROPERTY_DOUBLE(ThickDisk, ThickDiskInnerRadius, thickDiskInnerRadius)
GYOTO_PROPERTY_DOUBLE_UNIT(ThickDisk,
			   NumberDensityAtInnerRadius,
			   numberDensityAtInnerRadius)
GYOTO_PROPERTY_DOUBLE(ThickDisk,
		      TemperatureAtInnerRadius,
		      temperatureAtInnerRadius)
GYOTO_PROPERTY_DOUBLE(ThickDisk, TemperatureSlope, temperatureSlope)
GYOTO_PROPERTY_DOUBLE(ThickDisk, MagnetizationParameter,
		      magnetizationParameter)
GYOTO_PROPERTY_END(ThickDisk, Standard::properties)

// ACCESSORS
void ThickDisk::thickDiskOpeningAngle(double ang) {
  if (ang > M_PI/2.) throwError("In ThickDisk::thickDiskOpeningAngle "
				"opening angle should be <pi/2 rad"); 
  thickDiskOpeningAngle_=ang;}
double ThickDisk::thickDiskOpeningAngle()const{return thickDiskOpeningAngle_;}
void ThickDisk::thickDiskInnerRadius(double hh) {thickDiskInnerRadius_=hh;}
double ThickDisk::thickDiskInnerRadius()const{return thickDiskInnerRadius_;}
double ThickDisk::numberDensityAtInnerRadius() const {
  // Converts internal cgs central enthalpy to SI
  double dens=numberDensityAtInnerRadius_cgs_;
# ifdef HAVE_UDUNITS
  dens = Units::Converter("cm-3", "m-3")(dens);
# else
  GYOTO_WARNING << "Units ignored, please recompile Gyoto with --with-udunits"
		<< endl ;
# endif
  return dens; }
double ThickDisk::numberDensityAtInnerRadius(string const &unit) const
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
void ThickDisk::numberDensityAtInnerRadius(double dens) {
# ifdef HAVE_UDUNITS
  dens = Units::Converter("m-3", "cm-3")(dens);
# else
  GYOTO_WARNING << "Units ignored, please recompile Gyoto with --with-udunits"
		<< endl ;
# endif
  numberDensityAtInnerRadius_cgs_=dens;
}
void ThickDisk::numberDensityAtInnerRadius(double dens, string const &unit) {
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
void ThickDisk::temperatureAtInnerRadius(double tt) {temperatureAtInnerRadius_=tt;}
double ThickDisk::temperatureAtInnerRadius()const{return temperatureAtInnerRadius_;}
void ThickDisk::temperatureSlope(double ss) {temperatureSlope_=ss;}
double ThickDisk::temperatureSlope()const{return temperatureSlope_;}
void ThickDisk::magnetizationParameter(double rr) {
  magnetizationParameter_=rr;}
double ThickDisk::magnetizationParameter()const{
  return magnetizationParameter_;}
//

ThickDisk::ThickDisk() :
  Standard("ThickDisk"), thickDiskOpeningAngle_(0.785),
  thickDiskInnerRadius_(2.),
  numberDensityAtInnerRadius_cgs_(1.), temperatureAtInnerRadius_(1e10),
  temperatureSlope_(1.),
  magnetizationParameter_(1.)
{
  GYOTO_DEBUG << endl;
  spectrumThermalSynch_ = new Spectrum::ThermalSynchrotron();
}

ThickDisk::ThickDisk(const ThickDisk& o) :
  Standard(o), thickDiskOpeningAngle_(o.thickDiskOpeningAngle_),
  thickDiskInnerRadius_(o.thickDiskInnerRadius_),
  numberDensityAtInnerRadius_cgs_(o.numberDensityAtInnerRadius_cgs_),
  temperatureAtInnerRadius_(o.temperatureAtInnerRadius_),
  temperatureSlope_(o.temperatureSlope_),
  magnetizationParameter_(o.magnetizationParameter_),
  spectrumThermalSynch_(NULL)
{
  GYOTO_DEBUG << endl;
  if (gg_) gg_->hook(this);
  if (o.spectrumThermalSynch_()) spectrumThermalSynch_=o.spectrumThermalSynch_->clone();

}
ThickDisk* ThickDisk::clone() const
{ return new ThickDisk(*this); }

ThickDisk::~ThickDisk() {
  GYOTO_DEBUG << endl;
  if (gg_) gg_->unhook(this);
}

void ThickDisk::radiativeQ(double Inu[], // output
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
    GYOTO_ERROR("In ThickDisk::radiativeQ(): Unknown coordinate system kind");
  }

  double number_density = numberDensityAtInnerRadius_cgs_
    *(thickDiskInnerRadius_*thickDiskInnerRadius_)/(rr*rr);

  double temperature = temperatureAtInnerRadius_
    *pow(thickDiskInnerRadius_/rr, temperatureSlope_);

  double r0 = 4., phi0 = 0., phi = coord_ph[3],
    sigr = 2., sigp = M_PI/4.; // spin0: r0=9; spin08: r0=4
  double gaussr = 1./(sqrt(2.*M_PI)*sigr)
    * exp(-0.5*(rcyl-r0)*(rcyl-r0)/(sigr*sigr));
  double dphi = fabs(phi-phi0), dphibis = fabs(phi-2.*M_PI-phi0);
  if (dphi > dphibis){
    dphi = dphibis;
  }
  double gaussp = 1./(sqrt(2.*M_PI)*sigp)
    * exp(-0.5*dphi*dphi/(sigp*sigp));
  double gauss2d = gaussr*gaussp;
  double T0 = 1.6e11; // this is to be tuned: 4e11 too small, 6e11 looks good
  double DeltaTemperature = T0*gauss2d;

  //temperature+=DeltaTemperature;

  double thetae = GYOTO_BOLTZMANN_CGS*temperature
    /(GYOTO_ELECTRON_MASS_CGS*GYOTO_C2_CGS);

  double BB = sqrt(4.*M_PI*magnetizationParameter_
		   *GYOTO_PROTON_MASS_CGS * GYOTO_C_CGS * GYOTO_C_CGS
		   *number_density);

  double B0 = 100.; // for ne_inner=5.5e4, B_inner=10.2G; B0=50 too small, B0=100 looks good 
  double DeltaB = B0*gauss2d;

  //BB += DeltaB;

  // // Random generator: mersenne_twister_engine seeded with rd()
  // std::random_device rd;
  // std::mt19937 generator(rd());
  // // Define a real uniform distribution within some bounds
  // std::uniform_real_distribution<double> distribution(0.9,1.1);
  // double randnb = distribution(generator); // draw a random number

  
  //cout << "r ne B= " << coord_ph[1] << " " << number_density << " " << BB << endl;
  //cout << "IN DISK r, z, ne, B= " << coord_ph[1] << " " << zz << " " << number_density << " " << BB << endl;
  //GYOTO_ERROR("testjet");

  double nu0 = GYOTO_ELEMENTARY_CHARGE_CGS*BB
    /(2.*M_PI*GYOTO_ELECTRON_MASS_CGS*GYOTO_C_CGS); // cyclotron freq

  //cout << "jet stuff= " << coord_ph[1] << " " << coord_ph[2] << " " << zz << " " << rcyljetbase << " " << rcyl << " " << number_density << " " << thetae << " " << temperatureSlope_ << " " << nu0 << endl;
  //cout << "jet zz,rcyl,th,ph,ne,Te= " <<  zz << " " << rcyl << " " << coord_ph[2] << " " << coord_ph[3] << " " << number_density << " " << temperature << endl;
  // Use that line for Compton study:
  //cout <<  zz << " " << rcyl << " " << number_density << " " << temperature << endl;

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
  //cout << "for anu jnu: " << coord_ph[1] << " " << zz << " " << temperature << " " << number_density << " " << nu0 << " " << thetae << " " << besselK2 << endl;
  //cout << "nu passed to synchro= " << nu_ems[0] << endl;
  spectrumThermalSynch_->radiativeQ(jnu_synch,anu_synch,
				    nu_ems,nbnu);

  // RETURNING TOTAL INTENSITY AND TRANSMISSION
  for (size_t ii=0; ii<nbnu; ++ii){

    double jnu_tot = jnu_synch[ii],
      anu_tot = anu_synch[ii];

    //cout << "in disk stuff: " << zz << " " << rcyl << " " << nu_ems[0]  << " " << number_density << " " << nu0 << " " << temperature << " " << thetae << " " << jnu_tot << " " << anu_tot << " " << dsem << endl;

    //cout << "in jet stuff: " << number_density << " " << nu0 << " " << thetae << " " << hypergeom << " " << jnu_tot << " " << anu_tot << " " << dsem << endl;

    //cout << "at r,th= " << coord_ph[1] << " " << coord_ph[2] << endl;
    //cout << "at rcyl,z= " << rcyl << " " << zz << endl;
    //cout << "jet jnu anu kappa= " << jnu_tot << " " << anu_tot << endl; //x" " << jnu_tot/anu_tot << " " << dsem << endl;

    // expm1 is a precise implementation of exp(x)-1

    // TEST compare to Nalewajko
    //anu_tot=0.;
    //jnu_tot=1./sqrt(nu_ems[0]);
    
    double em1=std::expm1(-anu_tot * dsem * gg_->unitLength());
    Taunu[ii] = em1+1.;
    Inu[ii] = anu_tot == 0. ? jnu_tot * dsem * gg_->unitLength() :
      -jnu_tot / anu_tot * em1;

    if (Inu[ii]<0.)
      GYOTO_ERROR("In ThickDisk::radiativeQ: Inu<0");
    if (Inu[ii]!=Inu[ii] or Taunu[ii]!=Taunu[ii])
      GYOTO_ERROR("In ThickDisk::radiativeQ: Inu or Taunu is nan");
    if (Inu[ii]==Inu[ii]+1. or Taunu[ii]==Taunu[ii]+1.)
      GYOTO_ERROR("In ThickDisk::radiativeQ: Inu or Taunu is infinite");

  }
}

double ThickDisk::operator()(double const coord[4]) {
  // zpos: modulus of altitude above equatorial plane
  // rproj: radius projected in the equatorial plane
  double zpos=0., rproj=0.;
  
  switch (gg_ -> coordKind()) {
  case GYOTO_COORDKIND_SPHERICAL:
    rproj  = coord[1]*sin(coord[2]);
    zpos  = fabs(coord[1]*cos(coord[2]));
    break;
  case GYOTO_COORDKIND_CARTESIAN:
    zpos  = fabs(coord[3]);
    rproj  = sqrt(coord[1]*coord[1]+coord[2]*coord[2]);
    break;
  default:
    GYOTO_ERROR("ThickDisk::operator(): unknown COORDKIND");
  }
  
  //if (rproj < thickDiskInnerRadius_)
  //  return 1.; // outside disk
  // I think its ok for rproj < rin also TBC
  
  double zdisk = 0.;
  if (rproj > thickDiskInnerRadius_)
    zdisk = (rproj - thickDiskInnerRadius_)
      * tan(M_PI/2. - thickDiskOpeningAngle_) ; // altitude of disk at rproj
  // zdisk is fixed at zero rproj <= rinner,
  // then the distance to the disk is always positive

  //  cout << "in disk distance, zpos, zdisk= " << zpos << " " << zdisk << endl;
  
  return zpos - zdisk; // >0 outside, <0 inside flared disk
}

void ThickDisk::getVelocity(double const pos[4], double vel[4])
{
  
  double risco = 0.;
  if (gg_->kind()!="Minkowski" && gg_->kind()!="Hayward")
    risco=gg_->getRms(); // prograde Kerr ISCO
  // let risco=0 if metric is Minko; then ISCO not defined
  // let also risco=0 for Hayward as we would need to
  // compute it numerically and give it in xml Metric field,
  // not implemented so far
  
  //cout << "risco= " << risco << endl;

  // First, project the current location onto the eqplane
  double sinth = sin(pos[2]);
  double coord[4] = {pos[0], pos[1]*sinth, M_PI*0.5, pos[3]};
  double rcyl = coord[1]; // cylindrical radius of current location
  
  // Important Remark: the profile below, both above and below ISCO,
  // is such that u \cdot u is *not* equal to -1 if one uses the gmunu of the
  // actual point of emission ("pos", here). The 4vel is normalized when using
  // the gmunu of the projected point onto the eqplane ("coord", here).
  // This is Falcke+00 treatment, and RIAF people treatment, as much as I
  // understand.
  // Finding a prescription satisfying u \cdot u = -1 with the gmunu
  // of the point of emission is not trivial, see FV notes BelowISCO.pdf.
  
  if (rcyl > risco){
    gg_ -> circularVelocity(coord, vel, 1); // Keplerian 4vel at projected point onto eqplane, normalized with gmunu(coord), *not* with gmunu(pos)
    
    //gg_ -> circularVelocity(pos, vel, 1); // --> this call would lead to a
    // properly normalized 4vel, with the gmunu(pos). But then the 4vel
    // would be discontinuous at risco if one just copy-pastes C75 prescription
    // below ISCO.
  }else{
    // Impose KerrBL here, not clear yet for other metrics
    string kin = gg_->kind();
    if (kin != "KerrBL")
      GYOTO_ERROR("ThickDisk: KerrBL needed below ISCO!");
    double SPIN = static_cast<SmartPointer<Metric::KerrBL> >(gg_) -> spin();

    // Copy-paste Cunnigham 1975 equatorial profile at all z,
    // using rcyl for the radial coord
    double lambda_ms = (risco*risco - 2.*SPIN*sqrt(risco) + SPIN*SPIN)/(pow(risco,1.5) - 2.*sqrt(risco) + SPIN),
      gamma_ms = sqrt(1.-2./(3.*risco)),
      delta = rcyl*rcyl - 2.*rcyl + SPIN*SPIN,
      hh = (2.*rcyl - SPIN*lambda_ms)/delta;

    vel[0] = gamma_ms*(1.+2./rcyl*(1.+hh)); // this is: -Ems*g^{tt} + Lms*g^{tp}
    vel[1] = -sqrt(2./(3.*risco))*pow(risco/rcyl-1.,1.5); // this is: -sqrt{(-1 - g_{tt}*u^t - g_{pp}*u^p - 2*g_{tp}*u^t*u^p)/grr}
    vel[2] = 0.;
    vel[3] = gamma_ms/(rcyl*rcyl)*(lambda_ms+SPIN*hh);

    // Some earlier unconverged tests:
    
    // // metric coef at the actual point of emission
    // double guptt = gg_->gmunu_up(pos,0,0), guptp = gg_->gmunu_up(pos,0,3),
    //   guppp = gg_->gmunu_up(pos,3,3), gtt = gg_->gmunu(pos,0,0),
    //   gtp = gg_->gmunu(pos,0,3), gpp = gg_->gmunu(pos,3,3),
    //   grr=gg_->gmunu(pos,1,1);
    // // compute geodesic constants E, L at the point P with same z as
    // // the point of emission, but with rcyl=risco
    // double zz = pos[1]*cos(pos[2]);
    // double theta_P = atan2(risco,zz);
    // double r_P = sqrt(risco*risco + zz*zz);
    // double pos_P[4] = {pos[0],r_P,theta_P,pos[3]};
    // double vel_P[4];
    // // Compute Keplerian 4vel at point P
    // gg_ -> circularVelocity(pos_P, vel_P, 1);
    // // Metric coef at point P
    // double gtt_P = gg_->gmunu(pos_P,0,0), gtp_P = gg_->gmunu(pos_P,0,3),
    //   gpp_P = gg_->gmunu(pos_P,3,3);
    // // E and L at this point P
    // double Ems = -gtt_P*vel_P[0]-gtp_P*vel_P[3],
    //   Lms = gtp_P*vel_P[0] + gpp_P*vel_P[3];
    // // We assume that E keeps the same value Ems for all rcyl<risco
    // // and that L is linearly interpolated between 0 at horizon
    // // and Lms ar risco
    // double rhor = 1. + sqrt(1. - SPIN*SPIN), Lhor = 0.; 
    // double Linterp = Lhor + (rcyl - rhor)/(risco - rhor)*(Lms-Lhor);
    // // u^t and u^phi at the emission point assuming same E=Ems and L=Linterp
    // double ut = -Ems*guptt + Linterp*guptp, up = -Ems*guptp + Linterp*guppp;
    // vel[0] = ut;
    // vel[1] = -sqrt((-1.-gtt*ut*ut-gpp*up*up-2.*gtp*ut*up)/grr);
    // vel[2] = 0.;
    // vel[3] = up;
    // double u2 = gg_->ScalarProd(pos,vel,vel);
    // double tol=1e-5;
    // if (fabs(u2+1.)>tol) {
    //   cerr << setprecision(10) << "u2= " << u2 << " " << fabs(u2+1.) << endl;
    //   throwError("In ThickDisk: 4vel below "
    // 		 "ISCO not properly normalized!");
    // }
  }
  
  //cout << "at rcyl= " << rcyl << " u2 = " << gg_->ScalarProd(pos,vel,vel) << endl;
  
}

bool ThickDisk::isThreadSafe() const {
  return Standard::isThreadSafe()
    && (!spectrumThermalSynch_ || spectrumThermalSynch_->isThreadSafe());
}

void ThickDisk::metric(SmartPointer<Metric::Generic> gg) {
  if (gg_) gg_->unhook(this);
  string kin = gg->kind();
  //if (kin != "KerrBL" or kin!="NumericalMetricLorene")
  //  GYOTO_ERROR
  //    ("ThickDisk::metric(): metric must be KerrBL");
  // NB: KerrBL needed for ZAMO velocity in getVelocity,
  // could be generalized if needed
  Generic::metric(gg);
}
