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
GYOTO_PROPERTY_DOUBLE(ThickDisk, DensitySlope, densitySlope)
GYOTO_PROPERTY_DOUBLE(ThickDisk,
		      TemperatureAtInnerRadius,
		      temperatureAtInnerRadius)
GYOTO_PROPERTY_DOUBLE(ThickDisk, TemperatureSlope, temperatureSlope)
GYOTO_PROPERTY_DOUBLE(ThickDisk, MagnetizationParameter,
		      magnetizationParameter)
GYOTO_PROPERTY_END(ThickDisk, Standard::properties)

// Global variable to put to 1 to use the formalism
// of Vos+22 to which we compare Gyoto in the polarization paper.
// For debugging/code comparison only.
// For standard use, put it to 0 to use Gyoto formalism.
#define USE_IPOLE_FORMALISM 0

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
void ThickDisk::densitySlope(double ss) {densitySlope_=ss;}
double ThickDisk::densitySlope()const{return densitySlope_;}
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
  densitySlope_(2.),
  magnetizationParameter_(1.),
  magneticConfig_("None")
{
  GYOTO_DEBUG << endl;
  spectrumThermalSynch_ = new Spectrum::ThermalSynchrotron();
  if (USE_IPOLE_FORMALISM!=0 and USE_IPOLE_FORMALISM!=1){
    GYOTO_ERROR("Bad value of USE_IPOLE_FORMALISM should be 0 or 1!");
  }
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
  densitySlope_(o.densitySlope_),
  magnetizationParameter_(o.magnetizationParameter_),
  spectrumThermalSynch_(NULL),
  magneticConfig_(o.magneticConfig_)
{
  GYOTO_DEBUG << endl;
  if (gg_) gg_->hook(this);
  if (o.spectrumThermalSynch_()) spectrumThermalSynch_=o.spectrumThermalSynch_->clone();
  if (USE_IPOLE_FORMALISM!=0 and USE_IPOLE_FORMALISM!=1){
    GYOTO_ERROR("Bad value of USE_IPOLE_FORMALISM should be 0 or 1!");
  }

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

  if (rr<thickDiskInnerRadius_) {
    GYOTO_WARNING << "ThickDisk should typically have an inner radius at the "
      "horizon, here r<rin!" << endl;
  }

  double number_density = numberDensityAtInnerRadius_cgs_
    *pow(thickDiskInnerRadius_/rr, densitySlope_);
    //*(thickDiskInnerRadius_*thickDiskInnerRadius_)/(rr*rr);
    //*pow(thickDiskInnerRadius_/rr,3.);
  //cout << "nb density before expo= " << number_density << endl;
  double zsigma = thickDiskZGaussianSigma_*rcyl;
  double  expofact_zscaling = exp(-zz*zz/(2.*zsigma*zsigma)); // simple Gaussian modulation around z=0; RIAF model (Broderick+11) use zsigma=rcyl
  //cout << "ne before expo= " << rr << " " <<thickDiskInnerRadius_ << " " <<  numberDensityAtInnerRadius_cgs_ << " " << number_density << endl;
  number_density *= expofact_zscaling;
  
  //cout << "z/r fact, dens= " << expofact_zscaling << " " << number_density << endl;
  //cout << "ThickDisk r, z, rho= " << rcyl << " " << zz << " " << number_density << endl;
  double temperature = temperatureAtInnerRadius_
    *pow(thickDiskInnerRadius_/rr, temperatureSlope_);
    //*pow(thickDiskInnerRadius_/rr, 2.);
  //cout << "params in disk= " << thickDiskInnerRadius_ << " " << densitySlope_<< " " << temperatureSlope_ << " " << temperatureAtInnerRadius_ << " " <<  thickDiskZGaussianSigma_ << " " << magneticConfig_ << " " << endl;
  //throwError("test disk");

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
  //if (fabs(coord_ph[1]*cos(coord_ph[2]))<0.5)
  //cout << "thickdisk stuff= " << coord_ph[1] << " " << zz << " " << number_density << " " << BB << " " << thetae << endl;

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
  spectrumThermalSynch_->angle_B_pem(M_PI/2.);   // so we don't care about angle
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

    //if (rcyl<1.5 && fabs(zz)<0.01) jnu_tot=3e-15; // TEST!!

    //cout << "in disk stuff: " << zz << " " << rcyl << " " << nu_ems[0]  << " " << number_density << " " << nu0 << " " << temperature << " " << thetae << " " << jnu_tot << " " << anu_tot << " " << dsem << endl;
    //if (nu_ems[ii]>1e9)
    //  cout << dsem << " " << nu_ems[ii] << " " << jnu_tot << " " << anu_tot << endl;
    //if (nu_ems[ii]>1e9 && fabs(coord_ph[1]*cos(coord_ph[2]))<0.5)
    //  cout << "disk nu jnu anu in cgs= " <<nu_ems[ii] << " " << jnu_tot*10. << " " << anu_tot*0.01 << endl;

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
  if (USE_IPOLE_FORMALISM==0){
    // Vincent+22 formalism for circular/radial/mixed velocity profile
    
    double vel_circ[4], vel_rad[4];
    double rcyl = pos[1]*sin(pos[2]);// cylindrical radius of current location

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

    // TEST!!!////// no spatial velocity/////
    //vel[0] = sqrt(-1./gtt);
    //vel[1] = 0.;
    //vel[2] = 0.;
    //vel[3] = 0.;
    /////////////////////////////////////////

    //cout << "at rcyl, th-pi/2= " << rcyl << " " << fabs(pos[2]-M_PI/2.) << " u2 = " << gg_->ScalarProd(pos,vel,vel) << endl;
    double tol=0.03, normcur=gg_->ScalarProd(pos,vel,vel) ;
    //cout << "4vel at r z= " << pos[1]*sin(pos[2]) << " " << pos[1]*cos(pos[2]) << " " << vel[0] << " " << vel[1] << " " << vel[2] << " " << vel[3] << " " << normcur << endl; // gg_->ScalarProd(pos,vel,vel)
    
    if ((fabs(normcur+1.)>tol) ||
        (normcur!=normcur) ||
        (normcur==normcur+1)) {
      cerr << setprecision(10) << "at rcyl th= " << rcyl << " " << pos[2] << ", u2= " << normcur << endl;
      throwError("In ThickDisk: 4vel not properly normalized!");
    }
  }else{
    // Ipole formalism for comparison to Vos+22
    
    double rr = pos[1], theta=pos[2];
    double Risco = 6.;
    double a = 0.;
    if (rr>Risco){
      double sth2=sin(theta)*sin(theta),
      cth2=cos(theta)*cos(theta),
      rho2=rr*rr+a*a*cth2,
      a2=a*a,
      r2=rr*rr,
      DD=1.-2./rr+a2/r2,
      mu= 1. +a2*cth2/r2;

      double g_tt = -(1.-2./(rr*mu)),
      g_tp = -2.*a*sth2/(rr*mu),
      g_pp = r2*sth2*(1.+a2/r2+2.*a2*sth2/(r2*rr*mu));

      /*if (g_tt!=gg_->gmunu(pos,0,0)
          or g_tp!=gg_->gmunu(pos,0,3)
          or g_pp!=gg_->gmunu(pos,3,3)){
        cout << "g_tt : " << g_tt << ", " << gg_->gmunu(pos,0,0) << endl;
        cout << "g_tp : " << g_tp << ", " << gg_->gmunu(pos,0,3) << endl;
        cout << "g_pp : " << g_pp << ", " << gg_->gmunu(pos,3,3) << endl;
        GYOTO_ERROR("metric set by hand not equal to real metric.");
      }*/

      double omega = 1./(pow(rr,1.5)+a),
      m1oA2=g_tt+2.*omega*g_tp+omega*omega*g_pp,
      AA=sqrt(-1./m1oA2);

      vel[0]=AA;
      vel[1]=0.;
      vel[2]=0.;
      vel[3]=AA*omega;
    }else{
      double vel_cov[4]={0.,0.,0.,0.};
      vel_cov[0]=-1./sqrt(-gg_->gmunu_up(pos,0,0));
      vel_cov[1]=0.;
      vel_cov[2]=0.;
      vel_cov[3]=0.;

      vel[0]=gg_->gmunu_up(pos,0,0)*vel_cov[0];
      vel[1]=0.;
      vel[2]=0.;
      vel[3]=gg_->gmunu_up(pos,0,3)*vel_cov[3];
    }
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

void ThickDisk::radiativeQ(double *Inu, double *Qnu, double *Unu,
           double *Vnu,
           Eigen::Matrix4d *Onu,
           double const *nuem , size_t nbnu,
           double dsem,
           state_t const &coord_ph,
           double const *co) const {
  // polarized radiativeQ
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

  double vel[4]; // 4-velocity of emitter
  for (int ii=0;ii<4;ii++){
    vel[ii]=co[ii+4];
  }
  //cout << vel[0] << " " << vel[1] << " " << vel[2] << " " << vel[3] << endl;

  double number_density, zsigma;
  if (USE_IPOLE_FORMALISM==1){
    number_density = numberDensityAtInnerRadius_cgs_
      *pow(rr, densitySlope_); // scaling with (r/M)^{densitySlope_}
    zsigma = thickDiskZGaussianSigma_*rr;
  }else{
    number_density = numberDensityAtInnerRadius_cgs_
      *pow(thickDiskInnerRadius_/rr, densitySlope_); // scaling with (r/rin)^{-densitySlope_} ; careful with sign
    zsigma = thickDiskZGaussianSigma_*rcyl;
  }
  //cout << "nb density before expo= " << number_density << endl;

  double expofact_zscaling = exp(-zz*zz/(2.*zsigma*zsigma)); // simple Gaussian modulation around z=0; RIAF model (Broderick+11) use zsigma=rcyl
  //cout << "ne before expo= " << rr << " " <<thickDiskInnerRadius_ << " " <<  numberDensityAtInnerRadius_cgs_ << " " << number_density << endl;
  number_density *= expofact_zscaling;

  //number_density=6e5; // TEST
  
  if (rr<thickDiskInnerRadius_) number_density=0.; //1e-3; //0.; //NB: actually Vos+22 use density=1e-3 below ISCO; their rmin is 1.05*rhor; for the Aimar+23 paper we took a zero density below ISCO in ipole as well.
  double thetae, temperature, BB;
  if (USE_IPOLE_FORMALISM==1){
    thetae = 200.*pow(rr, temperatureSlope_);
    temperature = GYOTO_ELECTRON_MASS_CGS*GYOTO_C2_CGS*thetae/GYOTO_BOLTZMANN_CGS;

    double BB0 = 100.; // Gauss
    BB = BB0*pow(rr,-1.);
  }else{
    temperature = temperatureAtInnerRadius_
      *pow(thickDiskInnerRadius_/rr, temperatureSlope_);
    thetae = GYOTO_BOLTZMANN_CGS*temperature
      /(GYOTO_ELECTRON_MASS_CGS*GYOTO_C2_CGS);
    //double BB0 = 100.; // Gauss
    BB = sqrt(4.*M_PI*magnetizationParameter_
	      *GYOTO_PROTON_MASS_CGS * GYOTO_C_CGS * GYOTO_C_CGS
	      *number_density);
    //BB = BB0*pow(rr,-1.);
    //BB = sqrt(4.*M_PI*magnetizationParameter_
    //		   *GYOTO_PROTON_MASS_CGS * GYOTO_C_CGS * GYOTO_C_CGS
    //		   *number_density);
  }
  
  double nu0 = GYOTO_ELEMENTARY_CHARGE_CGS*BB
    /(2.*M_PI*GYOTO_ELECTRON_MASS_CGS*GYOTO_C_CGS); // cyclotron freq

  //cout << "jet stuff= " << coord_ph[1] << " " << coord_ph[2] << " " << zz << " " << rcyljetbase << " " << rcyl << " " << number_density << " " << thetae << " " << temperatureSlope_ << " " << nu0 << endl;

  /**************************/
  /* CHOOSE BFIELD GEOMETRY */
  /**************************/
  
  double B4vect[4]={0.,0.,0.,0.};
  if (USE_IPOLE_FORMALISM==1){ // or USE_IPOLE_FORMALISM==0){ // TEST
    string kin = gg_->kind();
    if (kin != "KerrBL") GYOTO_ERROR("ThickDisk in Ipole formalism should be in Kerr!");
    double spin = static_cast<SmartPointer<Metric::KerrBL> >(gg_) -> spin();
    if (spin!=0.) GYOTO_ERROR("ThickDisk in Ipole formalism should be in Schwarzschild!");
    
    computeB4vect_ipole(B4vect, magneticConfig_, co, coord_ph, spin);
    
  }else{
    computeB4vect(B4vect, magneticConfig_, co, coord_ph);
  }
  //cout << "B squared norm:" << gg_->ScalarProd(&coord_ph[0], B4vect, B4vect) << endl;
  double norm=sqrt(gg_->ScalarProd(&coord_ph[0], B4vect, B4vect));

  if (USE_IPOLE_FORMALISM==0){
    // Gyoto B formalism gives a unit vector
    if (fabs(norm-1.)>GYOTO_DEFAULT_ABSTOL) GYOTO_ERROR("Bad mf normalization");
  }else{
    // Ipole B formalism does not
    gg_->multiplyFourVect(B4vect,1./norm);
  }
  //cout << "B norm= " << norm << endl;

  double Chi=getChi(B4vect, coord_ph, vel); // this is EVPA
  //cout << "At r,th,ph Chi[deg]= " << coord_ph[1] << " " << coord_ph[2] << " " << coord_ph[3] << " " << Chi*180./M_PI << endl;

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

  if (theta_mag<0. or theta_mag>M_PI) throwError("ThickDisk: bad B angle");

  Eigen::Matrix4d Omat, Pmat;
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

  double besselK2 = bessk(2, 1./thetae);
  //if (fabs(zz)<10.) cout << "In ThickDisk: ne, temperature, BB, nu0, besselK2, theta_mag: " << number_density << " " << temperature << " " << BB << " " << nu0 << " " << besselK2 << " " << theta_mag << endl;

  // THERMAL SYNCHROTRON
  spectrumThermalSynch_->temperature(temperature);
  spectrumThermalSynch_->numberdensityCGS(number_density);
  spectrumThermalSynch_->angle_averaged(0); // no angle avg of course
  spectrumThermalSynch_->angle_B_pem(theta_mag); 
  spectrumThermalSynch_->cyclotron_freq(nu0);
  spectrumThermalSynch_->besselK2(besselK2);
  //cout << "for anu jnu: " << coord_ph[1] << " " << zz << " " << temperature << " " << number_density << " " << nu0 << " " << thetae << " " << besselK2 << endl;
  //cout << "nu passed to synchro= " << nuem[0] << endl;

  //if (number_density==0.) {
  //if (number_density<1.e4) {
  if (number_density<numberDensityAtInnerRadius_cgs_/1e10) {
    
    // Can happen due to strongly-killing z-expo factor
    // if zsigma is small. Then leads to 0/0 in synchro stuff. TBC
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

  // RETURNING TOTAL INTENSITY AND TRANSMISSION
  for (size_t ii=0; ii<nbnu; ++ii) {
    //if (fabs(zz)<10.) cout << "In ThickDisk: rr, jInu, jQnu, jUnu, jVnu: " << rr <<  " " << jInu[ii] << ", " << jQnu[ii] << ", " << jUnu[ii] << ", " << jVnu[ii] << endl;
    //cout << "In ThickDisk: aInu, aQnu, aUnu, aVnu: " << aInu[ii] << ", " << aQnu[ii] << ", " << aUnu[ii] << ", " << aVnu[ii] << endl;
    //cout << "In ThickDisk: rQnu, rUnu, rVnu: " << rotQnu[ii] << ", " << rotUnu[ii] << ", " << rotVnu[ii] << endl;
    //cout << "RADSTUFF: " << "r= " << coord_ph[1] << " " << ", th= " << coord_ph[2]*180./M_PI << " Rad transf stuff: " << jInu[ii]/(nuem[ii]*nuem[ii])*10. << ", " << jQnu[ii]/(nuem[ii]*nuem[ii])*10. << ", " << jUnu[ii]/(nuem[ii]*nuem[ii])*10. << ", " << jVnu[ii]/(nuem[ii]*nuem[ii])*10. << " " << nuem[ii]*aInu[ii]*0.01 << ", " << nuem[ii]*aQnu[ii]*0.01 << ", " << nuem[ii]*aUnu[ii]*0.01 << ", " << nuem[ii]*aVnu[ii]*0.01 << " " << nuem[ii]*rotQnu[ii]*0.01 << ", " << nuem[ii]*rotUnu[ii]*0.01 << ", " << nuem[ii]*rotVnu[ii]*0.01 << endl;
    Eigen::Vector4d JstokesDs=rotateJs(jInu[ii], jQnu[ii], jUnu[ii], jVnu[ii], Chi)*dsem*gg_->unitLength(), Jstokes=rotateJs(jInu[ii], jQnu[ii], jUnu[ii], jVnu[ii], Chi);
    //cout << Jstokes << endl;
    Omat = Omatrix(aInu[ii], aQnu[ii], aUnu[ii], aVnu[ii], rotQnu[ii], rotUnu[ii], rotVnu[ii], Chi, dsem);
    Pmat = Pmatrix(aInu[ii], aQnu[ii], aUnu[ii], aVnu[ii], rotQnu[ii], rotUnu[ii], rotVnu[ii], sin(2.*Chi), cos(2.*Chi), dsem);
    //cout << Omat << endl;
    // Computing the increment of the Stokes parameters. Equivalent to dInu=exp(-anu*dsem)*jnu*dsem in the non-polarised case.
    Eigen::Vector4d StokesFirst=Omat*JstokesDs,
      StokesMonika=Pmat*Jstokes, // Monika's version
      Stokes=StokesFirst; //StokesMonika; //StokesFirst;
    //cout << "StokesFirst= " << StokesFirst << endl;
    //cout << "StokesMonika= " << StokesMonika << endl;
    //cout << Stokes << endl;
    Inu[ii] = Stokes(0);
    Qnu[ii] = Stokes(1);
    Unu[ii] = Stokes(2);
    Vnu[ii] = Stokes(3);
    Onu[ii] = Omat;

    //cout << "In ThickDisk: r,th,ph, Inu, Qnu, Unu, Vnu, dsem, LP: " << coord_ph[1] << " " << coord_ph[2] << " " << coord_ph[3] << " " << Inu[ii] << ", " << Qnu[ii] << ", " << Unu[ii] << ", " << Vnu[ii] << ", " << dsem << ", " << pow(Qnu[ii]*Qnu[ii]+Unu[ii]*Unu[ii],0.5)/Inu[ii] << endl;

    if (Inu[ii]<0.)
      GYOTO_ERROR("In ThickDisk::radiativeQ(): Inu<0");
    if (Inu[ii]!=Inu[ii] or Onu[ii](0,0)!=Onu[ii](0,0))
      GYOTO_ERROR("In ThickDisk::radiativeQ(): Inu or Taunu is nan");
    if (Inu[ii]==Inu[ii]+1. or Onu[ii](0,0)==Onu[ii](0,0)+1.)
      GYOTO_ERROR("In ThickDisk::radiativeQ(): Inu or Taunu is infinite");
  }
}

void ThickDisk::magneticConfiguration(string config){
  magneticConfig_=config;
}

string ThickDisk::magneticConfiguration() const{
  return magneticConfig_;
}
