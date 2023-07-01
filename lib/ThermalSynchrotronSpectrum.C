/*
  Copyright 2018-2020 Frederic Vincent, Thibaut Paumard

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

#include "GyotoThermalSynchrotronSpectrum.h"
#include "GyotoDefs.h"
#include "GyotoUtils.h"
#include <cmath>
#ifdef GYOTO_USE_XERCES
#include "GyotoFactory.h"
#include "GyotoFactoryMessenger.h"
#endif
using namespace Gyoto;

#include "GyotoProperty.h"
GYOTO_PROPERTY_START(Spectrum::ThermalSynchrotron,
		     "Thermal synchrotron emission")
GYOTO_PROPERTY_END(Spectrum::ThermalSynchrotron, Generic::properties)

// Global variable to put to 1 to use the formalism
// of Vos+22 to which we compare Gyoto in the polarization paper.
// For debugging/code comparison only.
// For standard use, put it to 0 to use Gyoto formalism (Marszewski+21).
#define USE_IPOLE_FORMALISM 0

#define nstep_angint 100 // for angle-averaging integration
Spectrum::ThermalSynchrotron::ThermalSynchrotron()
: Spectrum::Generic("ThermalSynchrotron"),
  spectrumBB_(NULL), T_(10000.), numberdensityCGS_(0.),
  angle_B_pem_(0.), cyclotron_freq_(1.),
  angle_averaged_(0), bessel_K2_(1.)
{
  // A BB spectrum is needed to compute alpha_nu=j_nu/BB
  spectrumBB_ = new Spectrum::BlackBody(); 
}
Spectrum::ThermalSynchrotron::ThermalSynchrotron(const ThermalSynchrotron &o)
: Spectrum::Generic(o),
  T_(o.T_),
  spectrumBB_(NULL),
  numberdensityCGS_(o.numberdensityCGS_),
  angle_B_pem_(o.angle_B_pem_),
  cyclotron_freq_(o.cyclotron_freq_),
  angle_averaged_(o.angle_averaged_),
  bessel_K2_(o.bessel_K2_)
{
  if (o.spectrumBB_()) spectrumBB_=o.spectrumBB_->clone();
}

double Spectrum::ThermalSynchrotron::temperature() const { return T_; }
void Spectrum::ThermalSynchrotron::temperature(double tt) {
  T_ = tt; 
  spectrumBB_->temperature(T_);
}
double Spectrum::ThermalSynchrotron::numberdensityCGS() const { 
  return numberdensityCGS_; }
void Spectrum::ThermalSynchrotron::numberdensityCGS(double rho) { 
  numberdensityCGS_ = rho; }
double Spectrum::ThermalSynchrotron::angle_B_pem() const { 
  return angle_B_pem_; }
void Spectrum::ThermalSynchrotron::angle_B_pem(double angle) { 
  angle_B_pem_ = angle; }
double Spectrum::ThermalSynchrotron::cyclotron_freq() const { 
  return cyclotron_freq_; }
void Spectrum::ThermalSynchrotron::cyclotron_freq(double freq) { 
  cyclotron_freq_ = freq; }
bool Spectrum::ThermalSynchrotron::angle_averaged() const { 
  return angle_averaged_; }
void Spectrum::ThermalSynchrotron::angle_averaged(bool ang) { 
  angle_averaged_ = ang; }
double Spectrum::ThermalSynchrotron::besselK2() const { 
  return bessel_K2_; }
void Spectrum::ThermalSynchrotron::besselK2(double bessel) { 
  bessel_K2_ = bessel; }

  
Spectrum::ThermalSynchrotron * Spectrum::ThermalSynchrotron::clone() const
{ return new Spectrum::ThermalSynchrotron(*this); }

double Spectrum::ThermalSynchrotron::operator()(double nu) const {
  GYOTO_ERROR("In ThermalSynch: "
	     "Synchrotron emission not defined for optically thick case");
  return 0.;
}
double Spectrum::ThermalSynchrotron::operator()(double nu, 
						double , 
						double ds) const{
  double dsCGS = ds*100.; // ds should be given in SI
  // Returns intensity increment in SI:
  return jnuCGS(nu)*dsCGS*exp(-alphanuCGS(nu)*dsCGS)*GYOTO_INU_CGS_TO_SI;
}

double Spectrum::ThermalSynchrotron::jnuCGS(double nu) const{
  //std::cout << "in jnu brems= " << cst_ << " " <<  numberdensityCGS_ << " " << T_ << std::endl;
  double Theta_elec
    = T_*GYOTO_BOLTZMANN_CGS/(GYOTO_ELECTRON_MASS_CGS*GYOTO_C2_CGS);

  //std::cout << "in synch ther thetate ne nu0= " << Theta_elec << " " << numberdensityCGS_ << " " << cyclotron_freq_ << std::endl;

  int useWZ00=0; // 1 to use WZ00, 0 to use Pandya+16
  double emis_synch=0.;
  
  if (useWZ00==1){
    // Wardzinski & Zdziarski 2000
    double gamma0=0., chi0=0.;
    double sth=sin(angle_B_pem_), cth=cos(angle_B_pem_);
    if (Theta_elec<=0.08){
      gamma0 = sqrt(1+2.*nu*Theta_elec/cyclotron_freq_
		    *pow(1.+9.*nu*Theta_elec*sth*sth/(2.*cyclotron_freq_)
			 ,-0.3333333333));
      chi0 = sqrt((2.*Theta_elec*(gamma0*gamma0-1.))
		  /(gamma0*(3.*gamma0*gamma0-1.)));
    }else{
      gamma0 = sqrt(1.+pow(4.*nu*Theta_elec/(3.*cyclotron_freq_*sth)
			   ,0.6666666666));
      chi0 = sqrt(2.*Theta_elec/(3.*gamma0));
    }
    double tt = sqrt(gamma0*gamma0-1.)*sth,
      nn = nu*(1.+tt*tt)/(cyclotron_freq_*gamma0);
    double Z0 = pow((tt*exp(1./sqrt(1.+tt*tt)))/(1.+sqrt(1.+tt*tt)),2.*nn);
    double K2 = bessel_K2_;
    //std::cout << "bessel= " << K2 << std::endl;
    double ne0 = numberdensityCGS_/Theta_elec*gamma0*sqrt(gamma0*gamma0-1.)/K2
      *exp(-gamma0/Theta_elec);
    // this is j_nu synchro:
    emis_synch =
      M_PI*GYOTO_ELEMENTARY_CHARGE_CGS*GYOTO_ELEMENTARY_CHARGE_CGS
      /(2.*GYOTO_C_CGS)*sqrt(cyclotron_freq_*nu)*chi0*ne0
      *(1.+2.*cth*cth/(sth*sth*gamma0*gamma0))
      *pow(1.-(1.-1./(gamma0*gamma0))*cth*cth,0.25)
      *Z0;
    //std::cout << "stuff in emis synch ther= " << cyclotron_freq_ << " " << nu << " " << chi0 << " " << ne0 << " " << gamma0 << " " << Z0 << " " << angle_B_pem_ << " " << emis_synch << std::endl;
  }else if (USE_IPOLE_FORMALISM==1){
    // ipole https://github.com/moscibrodzka/ipole/blob/ipole-v2.0/model_radiation.c
    if (angle_B_pem_<=0 or angle_B_pem_>=M_PI)
      emis_synch=0;
    else{
      double EE=GYOTO_ELEMENTARY_CHARGE_CGS,S3=1.73205090765888,
        ME=GYOTO_ELECTRON_MASS_CGS, CL=GYOTO_C_CGS;

      double nuc=3.*cyclotron_freq_*sin(angle_B_pem_)/2.*Theta_elec*Theta_elec+1.,
        xx=nu/nuc;
      double I_I=2.5651*(1 + 1.92*pow(xx,-1./3.) + 0.9977*pow(xx,-2./3.)) * exp(-1.8899 * pow(xx,1./3.));
      emis_synch=numberdensityCGS_*EE*EE*nu/2./S3/CL/Theta_elec/Theta_elec*I_I;
    }
  }else{
    // Pandya, Zhang, Chandra, Gammie, 2016
    // Marszewski, Prather, Joshi, Pandya, Gammie 2021
    double nus = 2./9.*cyclotron_freq_*Theta_elec*Theta_elec*sin(angle_B_pem_),
      xx = nu/nus,
      Js = exp(-pow(xx,1./3.))*sqrt(2.)*M_PI/27.*sin(angle_B_pem_)*	\
      pow(pow(xx,1./2.)+pow(2.,11./12.)*pow(xx,1./6.),2.);
    emis_synch = numberdensityCGS_*					\
      GYOTO_ELEMENTARY_CHARGE_CGS*GYOTO_ELEMENTARY_CHARGE_CGS*cyclotron_freq_/ \
      GYOTO_C_CGS*\
      Js;
  }

  return emis_synch;
}

double Spectrum::ThermalSynchrotron::jQnuCGS(double nu) const{
  //std::cout << "in jnu brems= " << cst_ << " " <<  numberdensityCGS_ << " " << T_ << std::endl;
  double Theta_elec
    = T_*GYOTO_BOLTZMANN_CGS/(GYOTO_ELECTRON_MASS_CGS*GYOTO_C2_CGS);

  //std::cout << "in synch ther thetate ne nu0= " << Theta_elec << " " << numberdensityCGS_ << " " << cyclotron_freq_ << std::endl;

  double emis_synch=0.;
  if (USE_IPOLE_FORMALISM==1){
    // ipole https://github.com/moscibrodzka/ipole/blob/ipole-v2.0/model_radiation.c
    if (angle_B_pem_<=0 or angle_B_pem_>=M_PI)
      emis_synch=0;
    else{
      double EE=GYOTO_ELEMENTARY_CHARGE_CGS,S3=1.73205090765888,
        ME=GYOTO_ELECTRON_MASS_CGS, CL=GYOTO_C_CGS;

      double nuc=3.*cyclotron_freq_*sin(angle_B_pem_)/2.*Theta_elec*Theta_elec+1.,
        xx=nu/nuc;
      double I_Q=2.5651*(1 + 0.93193*pow(xx,-1./3.) + 0.499873*pow(xx,-2./3.)) * exp(-1.8899 * pow(xx,1./3.));
      emis_synch=numberdensityCGS_*EE*EE*nu/2./S3/CL/Theta_elec/Theta_elec*I_Q;
    }
  }else{
    // Marszewski, Prather, Joshi, Pandya, Gammie 2021
    double nus = 2./9.*cyclotron_freq_*Theta_elec*Theta_elec*sin(angle_B_pem_),
      xx = nu/nus,
      Js = exp(-pow(xx,1./3.))*sqrt(2.)*M_PI/27.*sin(angle_B_pem_)* \
      pow(pow(xx,1./2.)+((7.*pow(Theta_elec,24./25.)+35.)/(10.*pow(Theta_elec,24./25.)+75.))*pow(2.,11./12.)*pow(xx,1./6.),2.);
    emis_synch = numberdensityCGS_*         \
      GYOTO_ELEMENTARY_CHARGE_CGS*GYOTO_ELEMENTARY_CHARGE_CGS*cyclotron_freq_/ \
      GYOTO_C_CGS*\
      Js;
  }

  return emis_synch;
}

double Spectrum::ThermalSynchrotron::jUnuCGS(double nu) const{
  // Marszewski, Prather, Joshi, Pandya, Gammie 2021
  return 0.;
}

double Spectrum::ThermalSynchrotron::jVnuCGS(double nu) const{
  //std::cout << "in jnu brems= " << cst_ << " " <<  numberdensityCGS_ << " " << T_ << std::endl;
  double Theta_elec
    = T_*GYOTO_BOLTZMANN_CGS/(GYOTO_ELECTRON_MASS_CGS*GYOTO_C2_CGS);

  //std::cout << "in synch ther thetate ne nu0= " << Theta_elec << " " << numberdensityCGS_ << " " << cyclotron_freq_ << std::endl;

  double emis_synch=0.;
  if (USE_IPOLE_FORMALISM==1){
    // ipole https://github.com/moscibrodzka/ipole/blob/ipole-v2.0/model_radiation.c
    if (angle_B_pem_<=0 or angle_B_pem_>=M_PI)
      emis_synch=0;
    else{
      double EE=GYOTO_ELEMENTARY_CHARGE_CGS,S3=1.73205090765888,
        ME=GYOTO_ELECTRON_MASS_CGS, CL=GYOTO_C_CGS;
      
      double nuc=3.*cyclotron_freq_*sin(angle_B_pem_)/2.*Theta_elec*Theta_elec+1.,
        xx=nu/nuc;
      double I_V=(1.81348/xx+3.42319*pow(xx,-2./3.)+0.0292545*pow(xx,-0.5)+2.03773*pow(xx,-1./3.)) * exp(-1.8899 * pow(xx,1./3.));
      emis_synch=2.*numberdensityCGS_*EE*EE*nu/tan(angle_B_pem_)/3./S3/CL/Theta_elec/Theta_elec/Theta_elec*I_V;
    }
  }else{
    // Marszewski, Prather, Joshi, Pandya, Gammie 2021
    double nus = 2./9.*cyclotron_freq_*Theta_elec*Theta_elec*sin(angle_B_pem_),
      xx = nu/nus,
      Js = exp(-pow(xx,1./3.))*cos(angle_B_pem_)/Theta_elec* \
      (M_PI/3.+M_PI/3.*pow(xx,1./3.)+2./300.*sqrt(xx)+2.*M_PI/19.*pow(xx,2./3.));
    emis_synch = numberdensityCGS_*         \
      GYOTO_ELEMENTARY_CHARGE_CGS*GYOTO_ELEMENTARY_CHARGE_CGS*cyclotron_freq_/ \
      GYOTO_C_CGS*\
      Js;
  }

  return emis_synch;
}

double Spectrum::ThermalSynchrotron::alphanuCGS(double nu) const{
  double BB  = (*spectrumBB_)(nu)/GYOTO_INU_CGS_TO_SI; // B_nu in cgs
  double jnu = jnuCGS(nu);
  if (BB==0.){
    if (jnu==0.) return 0.;
    else GYOTO_ERROR("In ThermalSynch: alphanu undefined!");
  }
  // Kirchhoff's law:
  return jnuCGS(nu)/BB;
}

double Spectrum::ThermalSynchrotron::alphaQnuCGS(double nu) const{
  double BB  = (*spectrumBB_)(nu)/GYOTO_INU_CGS_TO_SI; // B_nu in cgs
  double jnu = jQnuCGS(nu);
  if (BB==0.){
    if (jnu==0.) return 0.;
    else GYOTO_ERROR("In ThermalSynch: alphanu undefined!");
  }
  // Kirchhoff's law:
  return jQnuCGS(nu)/BB;
}

double Spectrum::ThermalSynchrotron::alphaUnuCGS(double nu) const{
  double BB  = (*spectrumBB_)(nu)/GYOTO_INU_CGS_TO_SI; // B_nu in cgs
  double jnu = jUnuCGS(nu);
  if (BB==0.){
    if (jnu==0.) return 0.;
    else GYOTO_ERROR("In ThermalSynch: alphanu undefined!");
  }
  // Kirchhoff's law:
  return jUnuCGS(nu)/BB;
}

double Spectrum::ThermalSynchrotron::alphaVnuCGS(double nu) const{
  double BB  = (*spectrumBB_)(nu)/GYOTO_INU_CGS_TO_SI; // B_nu in cgs
  double jnu = jVnuCGS(nu);
  if (BB==0.){
    if (jnu==0.) return 0.;
    else GYOTO_ERROR("In ThermalSynch: alphanu undefined!");
  }
  // Kirchhoff's law:
  return jVnuCGS(nu)/BB;
}

double Spectrum::ThermalSynchrotron::rQnuCGS(double nu) const{
  double Theta_elec
    = T_*GYOTO_BOLTZMANN_CGS/(GYOTO_ELECTRON_MASS_CGS*GYOTO_C2_CGS);

  double rho_Q=0;
  if (USE_IPOLE_FORMALISM==1){
    if (angle_B_pem_<=0 or angle_B_pem_>=M_PI)
      rho_Q=0;
    else{
      double EE=GYOTO_ELEMENTARY_CHARGE_CGS,S3=1.73205090765888,
        ME=GYOTO_ELECTRON_MASS_CGS, CL=GYOTO_C_CGS, S2=1.41421356237310;

      double wp2=4.*M_PI*numberdensityCGS_*EE*EE/ME,
        omega0=2.*M_PI*cyclotron_freq_,
        Xe=Theta_elec*sqrt(S2*sin(angle_B_pem_)*(1.e3*omega0/2./M_PI/nu));
      double extraterm=(0.011*exp(-Xe/47.2)-pow(2.,-1./3.)/pow(3.,23./6.)*M_PI*1.e4*pow(Xe+1.e-16,-8./3.))*(0.5+0.5*tanh((log(Xe)-log(120.))/0.1)),
        jffunc=2.011*exp(-pow(Xe,1.035)/4.7)-cos(Xe*0.5)*exp(-pow(Xe,1.2)/2.73)-0.011*exp(-Xe/47.2)+extraterm;
      rho_Q=2.*M_PI*nu/2./CL*wp2*omega0*omega0/pow(2.*M_PI*nu,4)*jffunc*(Theta_elec/(2.*Theta_elec*Theta_elec)+6.*Theta_elec)*sin(angle_B_pem_)*sin(angle_B_pem_);
    }
  }else{
    // Marszewski, Prather, Joshi, Pandya, Gammie 2021
    double nus = 2./9.*cyclotron_freq_*Theta_elec*Theta_elec*sin(angle_B_pem_),
      xx = nu/nus,
      f0=2.011*exp(-19.78*pow(xx,-0.5175))-cos(39.89*pow(xx,-0.5))*exp(-70.16*pow(xx,-0.6))-0.011*exp(-1.69*pow(xx,-0.5)),
      fm=f0+(0.011*exp(-1.69*pow(xx,-0.5))-0.003135*pow(xx,4./3.))*1./2.*(1.+tanh(10*log(0.6648*pow(xx,-0.5))));
      
      rho_Q=numberdensityCGS_*pow(GYOTO_ELEMENTARY_CHARGE_CGS,2.)*pow(cyclotron_freq_,2.)*pow(sin(angle_B_pem_),2.)/ \
      (GYOTO_ELECTRON_MASS_CGS*GYOTO_C_CGS*pow(nu,3.))* \
      fm*(bessk1(pow(Theta_elec,-1))/bessk(2., pow(Theta_elec,-1))+6.*Theta_elec);
  }
  //return 0.;
  return rho_Q;
}

double Spectrum::ThermalSynchrotron::rUnuCGS(double nu) const{
  // Marszewski, Prather, Joshi, Pandya, Gammie 2021
  return 0.;
}

double Spectrum::ThermalSynchrotron::rVnuCGS(double nu) const{
  double Theta_elec
    = T_*GYOTO_BOLTZMANN_CGS/(GYOTO_ELECTRON_MASS_CGS*GYOTO_C2_CGS);

  double rho_V=0;
  if (USE_IPOLE_FORMALISM==1){
    double EE=GYOTO_ELEMENTARY_CHARGE_CGS,S3=1.73205090765888,
        ME=GYOTO_ELECTRON_MASS_CGS, CL=GYOTO_C_CGS, S2=1.41421356237310;

    double wp2=4.*M_PI*numberdensityCGS_*EE*EE/ME,
      omega0=2.*M_PI*cyclotron_freq_,
      Xe=Theta_elec*sqrt(S2*sin(angle_B_pem_)*(1.e3*omega0/2./M_PI/nu)),
      Je=0.43793091*log(1.+0.00185777*pow(Xe,1.50316886));

    if (Theta_elec>3.0){
      rho_V=2.*M_PI*nu/CL*wp2*omega0/pow(2.*M_PI*nu,3)*((-log(1./Theta_elec/2.)-0.5772)-Je)/(2*Theta_elec*Theta_elec)*cos(angle_B_pem_);
    }else if (0.2<Theta_elec and Theta_elec<=3.0){
      rho_V=2.*M_PI*nu/CL*wp2*omega0/pow(2.*M_PI*nu,3)*(bessk0(1./Theta_elec)-Je)/bessk(2.,1./Theta_elec)*cos(angle_B_pem_);
    }else{
      rho_V=2.*M_PI*nu/CL*wp2*omega0/pow(2.*M_PI*nu,3)*cos(angle_B_pem_);
    }
  }else{
    // Marszewski, Prather, Joshi, Pandya, Gammie 2021
    double nus = 2./9.*cyclotron_freq_*Theta_elec*Theta_elec*sin(angle_B_pem_),
      xx = nu/nus,
      DeltaJ5=0.4379*log(1.+1.3414*pow(xx,-0.7515));

    rho_V=2.*numberdensityCGS_*pow(GYOTO_ELEMENTARY_CHARGE_CGS,2.)*cyclotron_freq_/ \
    (GYOTO_ELECTRON_MASS_CGS*GYOTO_C_CGS*pow(nu,2.))* \
    (bessk0(pow(Theta_elec,-1.))-DeltaJ5)/bessk(2., pow(Theta_elec,-1.))*cos(angle_B_pem_);
  }
  //return 0.;
  return rho_V;
}

void Spectrum::ThermalSynchrotron::radiativeQ(double jnu[], // output
						double alphanu[], // output
						double const nu_ems[],
						size_t nbnu
						) {
  double Theta_elec
    = T_*GYOTO_BOLTZMANN_CGS/(GYOTO_ELECTRON_MASS_CGS*GYOTO_C2_CGS);
  
  double thetae_min_ther = 0.01;
  if (Theta_elec < thetae_min_ther) {
    // Below this value, 0/0 problems arise. From mma it is clear
    // that jnu goes quickly to 0 for thetae<0.01
    for (size_t ii=0; ii< nbnu; ++ii){
      jnu[ii]=0.;
      alphanu[ii]=0.;
    }
    return;
  }
  
  for (size_t ii=0; ii< nbnu; ++ii){
    double nu = nu_ems[ii];
    double BB  = (*spectrumBB_)(nu) ;
    double jnucur=0., anucur=0.;
    if (!angle_averaged_){
      jnucur = jnuCGS(nu);
    }else{
      double th0=0.01, thNm1=M_PI-0.01; // sin(theta) must never be 0
      double hh=(thNm1-th0)/double(nstep_angint);
      double theta=th0;
      angle_B_pem(theta);
      double jnusinprev=jnuCGS(nu)*sin(theta), jnusinnext=jnusinprev;
      for (int jj=1;jj<=nstep_angint;jj++){
      	theta=th0+double(jj)*hh;
      	angle_B_pem(theta);
      	jnusinnext=jnuCGS(nu)*sin(theta);
      	jnucur+=0.5*0.5*hh*(jnusinprev+jnusinnext);
      	jnusinprev=jnusinnext;
      	//NB: averaged jnu is: 1/4pi * \int jnu dOmega = 1/2 * \int jnu*sinth dth
      }
    }
    
    // OUTPUTS
    jnu[ii]= jnucur * GYOTO_JNU_CGS_TO_SI ;
    if (BB==0.){
      if (jnucur==0.) alphanu[ii]=0.;
      else GYOTO_ERROR("In ThermalSynch: alphanu undefined!");
    }else
      alphanu[ii]=jnu[ii]/BB;
    
  }
}


void Spectrum::ThermalSynchrotron::radiativeQ(double jInu[], double jQnu[], double jUnu[], double jVnu[], // Output
        double aInu[], double aQnu[], double aUnu[], double aVnu[], // Output
        double rQnu[], double rUnu[], double rVnu[], // Output
        double const nu_ems[],
        size_t nbnu ){

  double Theta_elec
    = T_*GYOTO_BOLTZMANN_CGS/(GYOTO_ELECTRON_MASS_CGS*GYOTO_C2_CGS);
  
  double thetae_min_ther = 0.01;
  if (Theta_elec < thetae_min_ther) {
    // Below this value, 0/0 problems arise. From mma it is clear
    // that jnu goes quickly to 0 for thetae<0.01
    for (size_t ii=0; ii< nbnu; ++ii){
      jInu[ii]=0.;
      jQnu[ii]=0.;
      jUnu[ii]=0.;
      jVnu[ii]=0.;
      aInu[ii]=0.;
      aQnu[ii]=0.;
      aUnu[ii]=0.;
      aVnu[ii]=0.;
      rQnu[ii]=0.;
      rUnu[ii]=0.;
      rVnu[ii]=0.;
    }
    return;
  }
  
  for (size_t ii=0; ii< nbnu; ++ii){
    double nu = nu_ems[ii];
    double jInucur=0., jQnucur=0.,jUnucur=0.,jVnucur=0.;
    double aInucur=0., aQnucur=0., aUnucur=0., aVnucur=0.;
    double rQnucur=0., rUnucur=0., rVnucur=0.; 
    if (!angle_averaged_){
      jInucur = jnuCGS(nu);
      jQnucur = jQnuCGS(nu);
      jUnucur = jUnuCGS(nu);
      jVnucur = jVnuCGS(nu);
      aInucur = alphanuCGS(nu);
      aQnucur = alphaQnuCGS(nu);
      aUnucur = alphaUnuCGS(nu);
      aVnucur = alphaVnuCGS(nu);
      rQnucur = rQnuCGS(nu);
      rUnucur = rUnuCGS(nu);
      rVnucur = rVnuCGS(nu);
      //std::cout << "jQ/jI :" << jQnucur/jInucur << std::endl; 
    }else{
      double th0=0.01, thNm1=M_PI-0.01; // sin(theta) must never be 0
      double hh=(thNm1-th0)/double(nstep_angint);
      double theta=th0;
      angle_B_pem(theta);

      double jInusinprev=jnuCGS(nu)*sin(theta), jInusinnext=jInusinprev;
      double jQnusinprev=jQnuCGS(nu)*sin(theta), jQnusinnext=jQnusinprev;
      double jUnusinprev=jUnuCGS(nu)*sin(theta), jUnusinnext=jUnusinprev;
      double jVnusinprev=jVnuCGS(nu)*sin(theta), jVnusinnext=jVnusinprev;
      double aInusinprev=alphanuCGS(nu)*sin(theta), aInusinnext=aInusinprev;
      double aQnusinprev=alphaQnuCGS(nu)*sin(theta), aQnusinnext=aQnusinprev;
      double aUnusinprev=alphaUnuCGS(nu)*sin(theta), aUnusinnext=aUnusinprev;
      double aVnusinprev=alphaVnuCGS(nu)*sin(theta), aVnusinnext=aVnusinprev;
      double rQnusinprev=rQnuCGS(nu)*sin(theta), rQnusinnext=rQnusinprev;
      double rUnusinprev=rUnuCGS(nu)*sin(theta), rUnusinnext=rUnusinprev;
      double rVnusinprev=rVnuCGS(nu)*sin(theta), rVnusinnext=rVnusinprev;

      for (int jj=1;jj<=nstep_angint;jj++){
        theta=th0+double(jj)*hh;
        angle_B_pem(theta);

        jInusinnext=jnuCGS(nu)*sin(theta);
        jQnusinnext=jQnuCGS(nu)*sin(theta);
        jUnusinnext=jUnuCGS(nu)*sin(theta);
        jVnusinnext=jVnuCGS(nu)*sin(theta);
        aInusinnext=alphanuCGS(nu)*sin(theta);
        aQnusinnext=alphaQnuCGS(nu)*sin(theta);
        aUnusinnext=alphaUnuCGS(nu)*sin(theta);
        aVnusinnext=alphaVnuCGS(nu)*sin(theta);
        rQnusinnext=rQnuCGS(nu)*sin(theta);
        rUnusinnext=rUnuCGS(nu)*sin(theta);
        rVnusinnext=rVnuCGS(nu)*sin(theta);

        jInucur+=0.5*0.5*hh*(jInusinprev+jInusinnext);
        jQnucur+=0.5*0.5*hh*(jQnusinprev+jQnusinnext);
        jUnucur+=0.5*0.5*hh*(jUnusinprev+jUnusinnext);
        jVnucur+=0.5*0.5*hh*(jVnusinprev+jVnusinnext);
        aInucur+=0.5*0.5*hh*(aInusinprev+aInusinnext);
        aQnucur+=0.5*0.5*hh*(aQnusinprev+aQnusinnext);
        aUnucur+=0.5*0.5*hh*(aUnusinprev+aUnusinnext);
        aVnucur+=0.5*0.5*hh*(aVnusinprev+aVnusinnext);
        rQnucur+=0.5*0.5*hh*(rQnusinprev+rQnusinnext);
        rUnucur+=0.5*0.5*hh*(rUnusinprev+rUnusinnext);
        rVnucur+=0.5*0.5*hh*(rVnusinprev+rVnusinnext);

        jInusinprev=jInusinnext;
        jQnusinprev=jQnusinnext;
        jUnusinprev=jUnusinnext;
        jVnusinprev=jVnusinnext;
        aInusinprev=aInusinnext;
        aQnusinprev=aQnusinnext;
        aUnusinprev=aUnusinnext;
        aVnusinprev=aVnusinnext;
        rQnusinprev=rQnusinnext;
        rUnusinprev=rUnusinnext;
        rVnusinprev=rVnusinnext;
        //NB: averaged jnu is: 1/4pi * \int jnu dOmega = 1/2 * \int jnu*sinth dth
      }
    }
    
    // OUTPUTS
    jInu[ii]=jInucur * GYOTO_JNU_CGS_TO_SI;
    jQnu[ii]=jQnucur * GYOTO_JNU_CGS_TO_SI;
    jUnu[ii]=jUnucur * GYOTO_JNU_CGS_TO_SI;
    jVnu[ii]=jVnucur * GYOTO_JNU_CGS_TO_SI;
    aInu[ii]=aInucur * GYOTO_ANU_CGS_TO_SI;
    aQnu[ii]=aQnucur * GYOTO_ANU_CGS_TO_SI;
    aUnu[ii]=aUnucur * GYOTO_ANU_CGS_TO_SI;
    aVnu[ii]=aVnucur * GYOTO_ANU_CGS_TO_SI;
    rQnu[ii]=rQnucur * GYOTO_ANU_CGS_TO_SI;
    rUnu[ii]=rUnucur * GYOTO_ANU_CGS_TO_SI;
    rVnu[ii]=rVnucur * GYOTO_ANU_CGS_TO_SI;
    
  }
}
