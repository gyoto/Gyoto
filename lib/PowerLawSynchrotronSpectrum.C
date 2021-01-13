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

#include "GyotoPowerLawSynchrotronSpectrum.h"
#include "GyotoDefs.h"
#include <cmath>
#ifdef GYOTO_USE_XERCES
#include "GyotoFactory.h"
#include "GyotoFactoryMessenger.h"
#endif
using namespace Gyoto;

#include "GyotoProperty.h"
GYOTO_PROPERTY_START(Spectrum::PowerLawSynchrotron,
		     "Powerlaw synchrotron emission")
GYOTO_PROPERTY_END(Spectrum::PowerLawSynchrotron, Generic::properties)

#define nstep_angint 10 // for angle-averaging integration
#define usePMT83 0 // 1 to use PMT83 jnu and alphanu, 0 to use Pandya+16
//#define gamma_min 1.
//#define gamma_max DBL_MAX

Spectrum::PowerLawSynchrotron::PowerLawSynchrotron()
: Spectrum::Generic("PowerLawSynchrotron"),
  numberdensityCGS_(0.),
  angle_B_pem_(0.), cyclotron_freq_(1.),
  PLindex_(0.), angle_averaged_(0), 
  gamma_min_(1.), gamma_max_(DBL_MAX)
{}
Spectrum::PowerLawSynchrotron::PowerLawSynchrotron(const PowerLawSynchrotron &o)
: Spectrum::Generic(o),
  spectrumBB_(NULL),
  numberdensityCGS_(o.numberdensityCGS_),
  angle_B_pem_(o.angle_B_pem_),
  cyclotron_freq_(o.cyclotron_freq_),
  PLindex_(o.PLindex_),
  angle_averaged_(o.angle_averaged_),
  gamma_min_(o.gamma_min_),
  gamma_max_(o.gamma_max_)
{
  if (o.spectrumBB_()) spectrumBB_=o.spectrumBB_->clone();
}

double Spectrum::PowerLawSynchrotron::numberdensityCGS() const { 
  return numberdensityCGS_; }
void Spectrum::PowerLawSynchrotron::numberdensityCGS(double rho) { 
  numberdensityCGS_ = rho; }
double Spectrum::PowerLawSynchrotron::angle_B_pem() const { 
  return angle_B_pem_; }
void Spectrum::PowerLawSynchrotron::angle_B_pem(double angle) { 
  angle_B_pem_ = angle; }
double Spectrum::PowerLawSynchrotron::cyclotron_freq() const { 
  return cyclotron_freq_; }
void Spectrum::PowerLawSynchrotron::cyclotron_freq(double freq) { 
  cyclotron_freq_ = freq; }
double Spectrum::PowerLawSynchrotron::PLindex() const { 
  return PLindex_; }
void Spectrum::PowerLawSynchrotron::PLindex(double ind) { 
  PLindex_ = ind; }
bool Spectrum::PowerLawSynchrotron::angle_averaged() const { 
  return angle_averaged_; }
void Spectrum::PowerLawSynchrotron::angle_averaged(bool ang) { 
  angle_averaged_ = ang; }
double Spectrum::PowerLawSynchrotron::gamma_min() const { 
  return gamma_min_; }
void Spectrum::PowerLawSynchrotron::gamma_min(double gmin) { 
  gamma_min_ = gmin; }
double Spectrum::PowerLawSynchrotron::gamma_max() const { 
  return gamma_max_; }
void Spectrum::PowerLawSynchrotron::gamma_max(double gmax) { 
  gamma_max_ = gmax; }
  
Spectrum::PowerLawSynchrotron * Spectrum::PowerLawSynchrotron::clone() const
{ return new Spectrum::PowerLawSynchrotron(*this); }

double Spectrum::PowerLawSynchrotron::operator()(double nu) const {
  GYOTO_ERROR("In PLSynch: "
	     "Synchrotron emission not defined for optically thick case");
  return 0.;
}
double Spectrum::PowerLawSynchrotron::operator()(double nu, 
						double , 
						double ds) const{
  double dsCGS = ds*100.; // ds should be given in SI
  // Returns intensity increment in SI:
  return jnuCGS(nu)*dsCGS*exp(-alphanuCGS(nu)*dsCGS)*GYOTO_INU_CGS_TO_SI;
}

double Spectrum::PowerLawSynchrotron::jnuCGS(double nu, double* stokesQ, double* stokesU, double* stokesV) const{
  double emis_synch = 0.;

  if (usePMT83==1){
    /* 
       From Petrosian & McTiernan 1983, Phys. Fluids 26 (10), eq. 32
       Putting g(mu)=1 and using (Y+ + Y_)=2 to get jnu and alphanu.
       NB: putting g(mu)=1 or 1/2 is not important, it boils down
       to redefining the % amount delta of PL energy wrt THER energy
    */
    //std::cout << "PL synch stuff= " << cyclotron_freq_ << " " << angle_B_pem_ << " " << PLindex_ << " " << numberdensityCGS_ << " " << angle_averaged_ << std::endl;
    emis_synch =
      sqrt(3.)*M_PI*GYOTO_ELEMENTARY_CHARGE_CGS*GYOTO_ELEMENTARY_CHARGE_CGS
      *cyclotron_freq_*sin(angle_B_pem_)/(2.*GYOTO_C_CGS)
      *numberdensityCGS_*(PLindex_-1.)
      *pow(3.*cyclotron_freq_*(PLindex_+1.)
	   *sin(angle_B_pem_)/(4.*nu),0.5*(PLindex_-1.))
      *exp(-0.5*(PLindex_+1.));
  }else{
    // Pandya, Zhang, Chandra, Gammie, 2016
    if (gamma_max_<sqrt(nu/cyclotron_freq_))
      GYOTO_ERROR("In PLSynchro: increase gamma_max_");
    // Ensure gamma_min_^2 < nu/nu0 < gamma_max_^2
    
    double sinth = sin(angle_B_pem_);
    double Js_I = pow(3.,PLindex_/2.)*(PLindex_-1.)*sinth/	\
      (2.*(PLindex_+1.)*(pow(gamma_min_,1.-PLindex_) -
			 pow(gamma_max_,1.-PLindex_))) * \
      tgamma((3.*PLindex_-1.)/12.)*tgamma((3.*PLindex_+19.)/12.) *	\
      pow(nu/(cyclotron_freq_*sinth),(1.-PLindex_)/2.);
    emis_synch = numberdensityCGS_*					\
      GYOTO_ELEMENTARY_CHARGE_CGS*GYOTO_ELEMENTARY_CHARGE_CGS*cyclotron_freq_/ \
      GYOTO_C_CGS*\
      Js_I;

      if (stokesQ!=NULL){
        double Js_Q = -Js_I*(PLindex_+1.)/(PLindex_+7./3.);
        *stokesQ = numberdensityCGS_*         \
      GYOTO_ELEMENTARY_CHARGE_CGS*GYOTO_ELEMENTARY_CHARGE_CGS*cyclotron_freq_/ \
      GYOTO_C_CGS*\
      Js_Q;
      }
      if (stokesU!=NULL){
        *stokesU = 0.;
      }
      if (stokesV!=NULL){
        double Js_V = -Js_I*171./250.*pow(PLindex_,0.5)/tan(angle_B_pem_)*pow(nu/(3.*cyclotron_freq_*sinth),-0.5);
        *stokesV = numberdensityCGS_*         \
      GYOTO_ELEMENTARY_CHARGE_CGS*GYOTO_ELEMENTARY_CHARGE_CGS*cyclotron_freq_/ \
      GYOTO_C_CGS*\
      Js_V;
      }
  }
  
  return emis_synch;
}

double Spectrum::PowerLawSynchrotron::alphanuCGS(double nu, double* stokesQ, double* stokesU, double* stokesV) const{
  double abs_synch = 0.;
    
  if (usePMT83==1){
    // From Petrosian & McTiernan 1983, Phys. Fluids 26 (10), eq. 32
    abs_synch =
      sqrt(3.)*M_PI*GYOTO_ELEMENTARY_CHARGE_CGS*GYOTO_ELEMENTARY_CHARGE_CGS
      *cyclotron_freq_*sin(angle_B_pem_)/(2.*GYOTO_C_CGS)
      *numberdensityCGS_*(PLindex_-1.)
      *pow(3.*cyclotron_freq_*(PLindex_+2.)*sin(angle_B_pem_)
	   /(4.*nu),0.5*PLindex_)
      *exp(-0.5*(PLindex_+2.))
      *(PLindex_+2.)
      /(GYOTO_ELECTRON_MASS_CGS*nu*nu);
  }else{
    if (gamma_max_<sqrt(nu/cyclotron_freq_))
      GYOTO_ERROR("In PLSynchro: increase gamma_max_");
    // Ensure gamma_min_^2 < nu/nu0 < gamma_max_^2

    double sinth = sin(angle_B_pem_);
    double As = pow(3.,(PLindex_+1.)/2.)*(PLindex_-1.)/	\
      (4.*(pow(gamma_min_,1.-PLindex_) -
			 pow(gamma_max_,1.-PLindex_))) *			\
      tgamma((3.*PLindex_+12.)/12.)*tgamma((3.*PLindex_+22.)/12.) *	\
      pow(nu/(cyclotron_freq_*sinth),-(PLindex_+2.)/2.);
    abs_synch = numberdensityCGS_*					\
      GYOTO_ELEMENTARY_CHARGE_CGS*GYOTO_ELEMENTARY_CHARGE_CGS/		\
      (nu*GYOTO_ELECTRON_MASS_CGS*GYOTO_C_CGS)*				\
      As;

    if (stokesQ!=NULL){
        double Js_Q = -Js_I*3./4.*pow(PLindex_-1.,43./500.);
        *stokesQ = numberdensityCGS_*         \
      GYOTO_ELEMENTARY_CHARGE_CGS*GYOTO_ELEMENTARY_CHARGE_CGS*cyclotron_freq_/ \
      GYOTO_C_CGS*\
      Js_Q;
      }
      if (stokesU!=NULL){
        *stokesU = 0.;
      }
      if (stokesV!=NULL){
        double Js_V = -Js_I*7./4.*pow(0.71*PLindex_+22./625.,197./500.)*\
          pow(pow(sinth,-48./25.)-1.,64./125.)*pow(nu/(cyclotron_freq_*sinth),-0.5);
        *stokesV = numberdensityCGS_*         \
      GYOTO_ELEMENTARY_CHARGE_CGS*GYOTO_ELEMENTARY_CHARGE_CGS*cyclotron_freq_/ \
      GYOTO_C_CGS*\
      Js_V;
      }
  }
  
  return abs_synch;
}

void Spectrum::PowerLawSynchrotron::radiativeQ(double jnu[], // output
						double alphanu[], // output
						double const nu_ems[],
						size_t nbnu
						) {
  for (size_t ii=0; ii< nbnu; ++ii){
    double nu = nu_ems[ii];
    double jnucur=0., anucur=0.;
    if (!angle_averaged_){
      jnucur = jnuCGS(nu);
      anucur = alphanuCGS(nu);
    }else{
      double th0=0.01, thNm1=M_PI-0.01; // avoiding sinth=0.
      double hh=(thNm1-th0)/double(nstep_angint);
      double theta=th0;
      angle_B_pem(theta);
      double jnusinprev=jnuCGS(nu)*sin(theta), jnusinnext=jnusinprev;
      double anusinprev=alphanuCGS(nu)*sin(theta), anusinnext=anusinprev;
      for (int jj=1;jj<=nstep_angint;jj++){
	theta=th0+double(jj)/2.*hh;
	angle_B_pem(theta);
	jnusinnext=jnuCGS(nu)*sin(theta);
	anusinnext=alphanuCGS(nu)*sin(theta);
	jnucur+=0.5*0.5*hh*(jnusinprev+jnusinnext);
	anucur+=0.5*0.5*hh*(anusinprev+anusinnext);
	//NB: averaged jnu is: \int jnu dOmega = 1/2 * \int jnu*sinth dth
      }
    }
    
    // OUTPUTS
    jnu[ii]= jnucur * GYOTO_JNU_CGS_TO_SI;
    alphanu[ii]= anucur * GYOTO_ANU_CGS_TO_SI;
    
  }
}
