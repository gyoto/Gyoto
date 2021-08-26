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
GYOTO_PROPERTY_DOUBLE(ThickDisk, ThickDiskZGaussianSigma,
		      thickDiskZGaussianSigma,"The standard deviation "
		      "of the Gaussian G modulating the density "
		      "with altitude z out of the equatorial plane, "
		      "divided by the cylindrical radius. So "
		      "G(z) = exp(-z^2 / 2*(rcyl*thickDiskZGaussianSigma_)^2)")
GYOTO_PROPERTY_DOUBLE(ThickDisk, ThickDiskInnerRadius, thickDiskInnerRadius)
GYOTO_PROPERTY_BOOL(ThickDisk, UseSelfAbsorption, NoUseSelfAbsorption,
		    useSelfAbsorption)
// GYOTO_PROPERTY_BOOL(ThickDisk,
// 		    AngleAveraged, NoAngleAveraged,
// 		    angleAveraged)
GYOTO_PROPERTY_VECTOR_DOUBLE(ThickDisk, VeloParam, veloParam,
			     "The two coef alpha and beta such that "
			     "u^r = u^r_circ + (1-alpha)*(u^r_rad - u^r_circ)"
			     " and similarly for Omega and beta.")
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
void ThickDisk::thickDiskZGaussianSigma(double sig) {
  thickDiskZGaussianSigma_=sig;}
double ThickDisk::thickDiskZGaussianSigma() const {
  return thickDiskZGaussianSigma_;}
void ThickDisk::thickDiskInnerRadius(double hh) {thickDiskInnerRadius_=hh;}
double ThickDisk::thickDiskInnerRadius()const{return thickDiskInnerRadius_;}
bool ThickDisk::useSelfAbsorption() const {return use_selfabsorption_;}
void ThickDisk::useSelfAbsorption(bool abs) {use_selfabsorption_=abs;}
void ThickDisk::veloParam(std::vector<double> const &v) {
  size_t nn = v.size();
  if (nn!=2)
    GYOTO_ERROR("In ThickDisk: choose exactly 2 velocity parameters");
  alpha_veloparam_ = v[0];
  beta_veloparam_  = v[1];
  if (alpha_veloparam_<0. || alpha_veloparam_>1.
      || beta_veloparam_<0. || beta_veloparam_>1.){
    GYOTO_ERROR("In ThickDisk: velocity parameters should be "
		"between 0 and 1!");
  }
}
std::vector<double> ThickDisk::veloParam() const {
  std::vector<double> v(2, 0.);
  v[0] = alpha_veloparam_;
  v[1] = beta_veloparam_;
  return v;
}
// bool ThickDisk::angleAveraged() const
// {return angle_averaged_;}
// void ThickDisk::angleAveraged(bool ang)
// {
//   angle_averaged_=ang;
//   spectrumThermalSynch_->angle_averaged(ang);
// }
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
  Standard("ThickDisk"),
  thickDiskInnerRadius_(2.),
  thickDiskZGaussianSigma_(1.),
  use_selfabsorption_(1),
  alpha_veloparam_(1.),
  beta_veloparam_(1.),
  numberDensityAtInnerRadius_cgs_(1.), temperatureAtInnerRadius_(1e10),
  temperatureSlope_(1.),
  magnetizationParameter_(1.)
{
  GYOTO_DEBUG << endl;
  spectrumThermalSynch_ = new Spectrum::ThermalSynchrotron();
}

ThickDisk::ThickDisk(const ThickDisk& o) :
  Standard(o), 
  thickDiskInnerRadius_(o.thickDiskInnerRadius_),
  thickDiskZGaussianSigma_(o.thickDiskZGaussianSigma_),
  use_selfabsorption_(o.use_selfabsorption_),
  alpha_veloparam_(o.alpha_veloparam_),
  beta_veloparam_(o.beta_veloparam_),
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
  double zsigma = thickDiskZGaussianSigma_*rcyl;
  double  expofact_zscaling = exp(-zz*zz/(2.*zsigma*zsigma)); // simple Gaussian modulation around z=0; RIAF model (Broderick+11) use zsigma=rcyl
  number_density *= expofact_zscaling;
  //cout << "z/r fact, dens= " << expofact_zscaling << " " << number_density << endl;
  //cout << "ThickDisk r, z, rho= " << rcyl << " " << zz << " " << number_density << endl;
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

  /*if (zz<0.1){
    cout << "at r rcyl z= " << rr << " " << rcyl << " " << zz << endl;
    cout << "ne,T,B= " << number_density << " " << temperature << " " << BB << endl;
    cout << "beta= " << 8.*M_PI*number_density*GYOTO_BOLTZMANN_CGS*temperature/(BB*BB) << endl;
    }*/

  double B0 = 100.; // for ne_inner=5.5e4, B_inner=10.2G; B0=50 too small, B0=100 looks good 
  double DeltaB = B0*gauss2d;

  //BB += DeltaB;

  // // Random generator: mersenne_twister_engine seeded with rd()
  // std::random_device rd;
  // std::mt19937 generator(rd());
  // // Define a real uniform distribution within some bounds
  // std::uniform_real_distribution<double> distribution(0.9,1.1);
  // double randnb = distribution(generator); // draw a random number

  
  //cout << "IN DISK r, z, ne, B= " << coord_ph[1] << " " << zz << " " << number_density << " " << BB << endl;
  //GYOTO_ERROR("testjet");

  double nu0 = GYOTO_ELEMENTARY_CHARGE_CGS*BB
    /(2.*M_PI*GYOTO_ELECTRON_MASS_CGS*GYOTO_C_CGS); // cyclotron freq

  //cout << "jet stuff= " << coord_ph[1] << " " << coord_ph[2] << " " << zz << " " << rcyljetbase << " " << rcyl << " " << number_density << " " << thetae << " " << temperatureSlope_ << " " << nu0 << endl;

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
  spectrumThermalSynch_->angle_averaged(1.); // impose angle-averaging
  spectrumThermalSynch_->angle_B_pem(0.);   // so we don't care about angle
  spectrumThermalSynch_->cyclotron_freq(nu0);
  double besselK2 = bessk(2, 1./thetae);
  spectrumThermalSynch_->besselK2(besselK2);
  //cout << "for anu jnu: " << coord_ph[1] << " " << zz << " " << temperature << " " << number_density << " " << nu0 << " " << thetae << " " << besselK2 << endl;
  //cout << "nu passed to synchro= " << nu_ems[0] << endl;

  //if (number_density==0.) {
  if (number_density<numberDensityAtInnerRadius_cgs_/1e10) { // CHECK THAT
    // Can happen due to strongly-killing z-expo factor
    // if zsigma is small. Then leads to 0/0 in synchro stuff. TBC
    for (size_t ii=0; ii<nbnu; ++ii){
      jnu_synch[ii]=0.;
      anu_synch[ii]=0.;
    }
  }else{
    spectrumThermalSynch_->radiativeQ(jnu_synch,anu_synch,
				      nu_ems,nbnu);
  }

  // RETURNING TOTAL INTENSITY AND TRANSMISSION
  for (size_t ii=0; ii<nbnu; ++ii){

    double jnu_tot = jnu_synch[ii],
      anu_tot = anu_synch[ii];

    //cout << "in disk stuff: " << zz << " " << rcyl << " " << nu_ems[0]  << " " << number_density << " " << nu0 << " " << temperature << " " << thetae << " " << jnu_tot << " " << anu_tot << " " << dsem << endl;

    // expm1 is a precise implementation of exp(x)-1
    double em1=std::expm1(-anu_tot * dsem * gg_->unitLength());
    if (use_selfabsorption_){
      // with self absorption
      Taunu[ii] = em1+1.;
      Inu[ii] = anu_tot == 0. ? jnu_tot * dsem * gg_->unitLength() :
	-jnu_tot / anu_tot * em1;
    }else{
      // no absorption
      Taunu[ii] = 1.;
      Inu[ii] = jnu_tot * dsem * gg_->unitLength() ;
    }
    //cout << "unit length= " << gg_->unitLength() << " " << jnu_tot * gg_->unitLength() <<  endl;

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

  // // 2019 paper version:
  // double zdisk = 0.;   // zdisk is fixed at zero rproj <= rinner,
  // // then the distance to the disk is always positive
  // double rproj_lim=thickDiskInnerRadius_;
  // if (rproj > thickDiskInnerRadius_){
  //   if (rproj > rproj_lim) // usual linear surface above rproj_lim
  //     zdisk = (rproj - thickDiskInnerRadius_)
  // 	* tan(M_PI/2. - thickDiskOpeningAngle_) ; 
  //   else // parabola surface below, connecting continuously
  //     zdisk = (rproj - thickDiskInnerRadius_)
  // 	*(rproj - thickDiskInnerRadius_)/(rproj_lim - thickDiskInnerRadius_)
  // 	* tan(M_PI/2. - thickDiskOpeningAngle_) ;
  // }  
  // return zpos - zdisk; // >0 outside, <0 inside flared disk

  // 2021 version without surface
  return -1.; // matter is everywhere

}

void ThickDisk::getVelocity(double const pos[4], double vel[4])
{
  double rcyl = pos[1]*sin(pos[2]);// cylindrical radius of current location

  double vel_circ[4], vel_rad[4];

  double rr = pos[1];
  double gtt = gg_->gmunu(pos,0,0),
    grr = gg_->gmunu(pos,1,1),
    gpp = gg_->gmunu(pos,3,3),
    gtp = gg_->gmunu(pos,0,3),
    guptt = gg_->gmunu_up(pos,0,0),
    guptp = gg_->gmunu_up(pos,0,3),
    guppp = gg_->gmunu_up(pos,3,3),
    guprr = gg_->gmunu_up(pos,1,1);
  
  // CIRCULAR VELOCITY
  
  // u_\mu = (u_t,0,0,u_phi) = -u_t (-1,0,0,ll)
  // with ll = rcyl^{3/2}/(rcyl+1.)
  
  double mycst=1; // Gold+20 choice, Keplerian is mycst=-2, 
  // see Appendix of 2021 M87 paper.
  double ll=pow(rcyl,1.5)/(rcyl+mycst);
  double u_t_minus=sqrt(-1./(guptt - 2.*guptp*ll + guppp*ll*ll));
  double u_t = -u_t_minus, u_phi = u_t_minus*ll;
  
  vel_circ[0] = guptt*u_t + guptp*u_phi;
  vel_circ[1] = 0.;
  vel_circ[2] = 0.;
  vel_circ[3] = guptp*u_t + guppp*u_phi;
  
  double Omega_circ = vel_circ[3]/vel_circ[0];
  
  // RADIAL VELOCITY
  
  // 4-vel obtained by imposing: u_t=-1, u_phi=0, u^theta=0
  // see FV notes SphericalVelocity.pdf for details
  vel_rad[0] = -guptt;
  vel_rad[1] = -sqrt((-1.-guptt)*guprr);
  vel_rad[2] = 0;
  vel_rad[3] = -guptp;

  double Omega_rad = vel_rad[3]/vel_rad[0];

  // MIXED VELOCITY

  double alpha=alpha_veloparam_, beta=beta_veloparam_;

  vel[1] = vel_circ[1] + (1-alpha)*(vel_rad[1]-vel_circ[1]);
  vel[2] = 0.;
  double Omega = Omega_circ + (1-beta)*(Omega_rad-Omega_circ);
  double normfact = gtt + 2*Omega*gtp + Omega*Omega*gpp;
  if (normfact>0) throwError("In ThickDisk::getVelocity: velocity "
			     "prescription non physical.");
  vel[0] = sqrt(-(1. + grr*vel[1]*vel[1])/normfact);
  vel[3] = Omega*vel[0];

  //cout << "at rcyl, th-pi/2= " << rcyl << " " << fabs(pos[2]-M_PI/2.) << " u2 = " << gg_->ScalarProd(pos,vel,vel) << endl;
  double tol=0.03, normcur=gg_->ScalarProd(pos,vel,vel) ;
  //cout << "4vel at r z= " << pos[1]*sin(pos[2]) << " " << pos[1]*cos(pos[2]) << " " << vel[0] << " " << vel[1] << " " << vel[2] << " " << vel[3] << " " << normcur << endl; // gg_->ScalarProd(pos,vel,vel)
  
  if ((fabs(normcur+1.)>tol) ||
      (normcur!=normcur) ||
      (normcur==normcur+1)) {
    cerr << setprecision(10) << "at rcyl th= " << rcyl << " " << pos[2] << ", u2= " << normcur << endl;
    throwError("In ThickDisk: 4vel not properly normalized!");
  }
  
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
