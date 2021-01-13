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

double Spectrum::ThermalSynchrotron::jnuCGS(double nu, double* stokesQ, double* stokesU, double* stokesV) const{
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
  }else{
    // Pandya, Zhang, Chandra, Gammie, 2016
    double nus = 2./9.*cyclotron_freq_*Theta_elec*Theta_elec*sin(angle_B_pem_),
      xx = nu/nus,
      Js_I = exp(-pow(xx,1./3.))*sqrt(2.)*M_PI/27.*sin(angle_B_pem_)*	\
      pow(pow(xx,1./2.)+pow(2.,11./12.)*pow(xx,1./6.),2.);
    emis_synch = numberdensityCGS_*					\
      GYOTO_ELEMENTARY_CHARGE_CGS*GYOTO_ELEMENTARY_CHARGE_CGS*cyclotron_freq_/ \
      GYOTO_C_CGS*\
      Js_I;

    if (stokesQ!=NULL){
        double Js_Q = -exp(-pow(xx,1./3.))*sqrt(2.)*M_PI/27.*sin(angle_B_pem_)*  \
      pow(pow(xx,1./2.)+((7.*pow(Theta_elec,24./25.)+35.)/(10.*pow(Theta_elec,24./25.)+75.))*pow(2.,11./12.)*pow(xx,1./6.),2.);
        *stokesQ = numberdensityCGS_*         \
      GYOTO_ELEMENTARY_CHARGE_CGS*GYOTO_ELEMENTARY_CHARGE_CGS*cyclotron_freq_/ \
      GYOTO_C_CGS*\
      Js_Q;
      }
    if (stokesU!=NULL){*stokesU = 0.;}

    if (stokesV!=NULL){
        double Js_V = -exp(-pow(xx,1./3.))*(37.-87.*sin(angle_B_pem_-28./25.))/(100.*(Theta_elec+1.))*\
      pow(1.+(pow(Theta_elec,3./5.)/25.+7./10.)*pow(xx,9./25.),5./3.);
        *stokesV = numberdensityCGS_*         \
      GYOTO_ELEMENTARY_CHARGE_CGS*GYOTO_ELEMENTARY_CHARGE_CGS*cyclotron_freq_/ \
      GYOTO_C_CGS*\
      Js_V;
      }
  }

  return emis_synch;
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
  }else{
    // Pandya, Zhang, Chandra, Gammie, 2016
    double nus = 2./9.*cyclotron_freq_*Theta_elec*Theta_elec*sin(angle_B_pem_),
      xx = nu/nus,
      Js = exp(-pow(xx,1./3.))*sqrt(2.)*M_PI/27.*sin(angle_B_pem_)* \
      pow(pow(xx,1./2.)+pow(2.,11./12.)*pow(xx,1./6.),2.);
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
	//NB: averaged jnu is: \int jnu dOmega = 1/2 * \int jnu*sinth dth
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
