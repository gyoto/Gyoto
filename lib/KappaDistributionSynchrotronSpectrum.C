/*
  Copyright 2018 Frederic Vincent, Thibaut Paumard

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

#include "GyotoKappaDistributionSynchrotronSpectrum.h"
#include "GyotoDefs.h"
#include <cmath>
#include <cstdlib> /* atof */
#ifdef GYOTO_USE_XERCES
#include "GyotoFactory.h"
#include "GyotoFactoryMessenger.h"
#endif
using namespace Gyoto;
using namespace std;

#include "GyotoProperty.h"
GYOTO_PROPERTY_START(Spectrum::KappaDistributionSynchrotron,
		     "Powerlaw synchrotron emission")
GYOTO_PROPERTY_END(Spectrum::KappaDistributionSynchrotron, Generic::properties)



#define nstep_angint 10 // for angle-averaging integration

Spectrum::KappaDistributionSynchrotron::KappaDistributionSynchrotron()
: Spectrum::Generic("KappaDistributionSynchrotron"),
  numberdensityCGS_(0.),
  angle_B_pem_(0.), cyclotron_freq_(1.), thetae_(1.),
  kappaindex_(0.), angle_averaged_(0), hypergeometric_(1)
{}
Spectrum::KappaDistributionSynchrotron::KappaDistributionSynchrotron(const KappaDistributionSynchrotron &o)
: Spectrum::Generic(o),
  spectrumBB_(NULL),
  numberdensityCGS_(o.numberdensityCGS_),
  angle_B_pem_(o.angle_B_pem_),
  cyclotron_freq_(o.cyclotron_freq_),
  thetae_(o.thetae_),
  kappaindex_(o.kappaindex_),
  hypergeometric_(o.hypergeometric_),
  angle_averaged_(o.angle_averaged_)
{
  if (o.spectrumBB_()) spectrumBB_=o.spectrumBB_->clone();
}

double Spectrum::KappaDistributionSynchrotron::numberdensityCGS() const { 
  return numberdensityCGS_; }
void Spectrum::KappaDistributionSynchrotron::numberdensityCGS(double rho) { 
  numberdensityCGS_ = rho; }
double Spectrum::KappaDistributionSynchrotron::angle_B_pem() const { 
  return angle_B_pem_; }
void Spectrum::KappaDistributionSynchrotron::angle_B_pem(double angle) { 
  angle_B_pem_ = angle; }
double Spectrum::KappaDistributionSynchrotron::cyclotron_freq() const { 
  return cyclotron_freq_; }
void Spectrum::KappaDistributionSynchrotron::cyclotron_freq(double freq) { 
  cyclotron_freq_ = freq; }
double Spectrum::KappaDistributionSynchrotron::thetae() const { 
  return thetae_; }
void Spectrum::KappaDistributionSynchrotron::thetae(double th) { 
  thetae_ = th; }
double Spectrum::KappaDistributionSynchrotron::kappaindex() const { 
  return kappaindex_; }
void Spectrum::KappaDistributionSynchrotron::kappaindex(double ind) { 
  kappaindex_ = ind; }
double Spectrum::KappaDistributionSynchrotron::hypergeometric() const { 
  return hypergeometric_; }
void Spectrum::KappaDistributionSynchrotron::hypergeometric(double hh) { 
  hypergeometric_ = hh; }
bool Spectrum::KappaDistributionSynchrotron::angle_averaged() const { 
  return angle_averaged_; }
void Spectrum::KappaDistributionSynchrotron::angle_averaged(bool ang) { 
  angle_averaged_ = ang; }
  
Spectrum::KappaDistributionSynchrotron * Spectrum::KappaDistributionSynchrotron::clone() const
{ return new Spectrum::KappaDistributionSynchrotron(*this); }

double Spectrum::KappaDistributionSynchrotron::operator()(double nu) const {
  throwError("In PLSynch: "
	     "Synchrotron emission not defined for optically thick case");
  return 0.;
}
double Spectrum::KappaDistributionSynchrotron::operator()(double nu, 
						double , 
						double ds) const{
  double dsCGS = ds*100.; // ds should be given in SI
  // Returns intensity increment in SI:
  return jnuCGS(nu)*dsCGS*exp(-alphanuCGS(nu)*dsCGS)*GYOTO_INU_CGS_TO_SI;
}

double Spectrum::KappaDistributionSynchrotron::jnuCGS(double nu) const{  

  // Pandya, Zhang, Chandra, Gammie, 2016

  //cout << "in kappa jnu stuff: " << cyclotron_freq_ << " " << thetae_ << " " << kappaindex_ << endl;
  
  double sinth = sin(angle_B_pem_),
    nuk = cyclotron_freq_*pow(thetae_*kappaindex_,2.)*sinth,
    Xk = nu/nuk,
    Js_low = pow(Xk,1./3.)*sinth*4.*M_PI*tgamma(kappaindex_-4./3.)/ \
    (pow(3.,7./3.)*tgamma(kappaindex_-2.)),
    Js_high = pow(Xk,-(kappaindex_-2.)/2.)*sinth*pow(3.,(kappaindex_-1.)/2.)* \
    (kappaindex_-1.)*(kappaindex_-2.)/4.*tgamma(kappaindex_/4.-1./3.)* \
    tgamma(kappaindex_/4.+4./3.),
    expo = 3.*pow(kappaindex_,-3./2.),
    Js = pow(pow(Js_low,-expo) + pow(Js_high,-expo),-1./expo);
  
  double emis_synch = numberdensityCGS_*				\
    GYOTO_ELEMENTARY_CHARGE_CGS*GYOTO_ELEMENTARY_CHARGE_CGS*cyclotron_freq_/ \
    GYOTO_C_CGS*							\
    Js;
  
  
  return emis_synch;
}

double Spectrum::KappaDistributionSynchrotron::alphanuCGS(double nu) const{

  // Pandya, Zhang, Chandra, Gammie, 2016

  //cout << "in kappa anu stuff: " << cyclotron_freq_ << " " << thetae_ << " " << kappaindex_ << " " << hypergeometric_ << endl;

  double sinth = sin(angle_B_pem_),
    nuk = cyclotron_freq_*pow(thetae_*kappaindex_,2.)*sinth,
    Xk = nu/nuk,
    As_low = pow(Xk,-2./3.)*pow(3.,1./6.)*10./41.* \
    2.*M_PI/pow(thetae_*kappaindex_,10./3.-kappaindex_) * \
    (kappaindex_-1.)*(kappaindex_-2.)*kappaindex_/(3.*kappaindex_-1.) * \
    tgamma(5./3.)*hypergeometric_,
    As_high = pow(Xk,-(1.+kappaindex_)/2.)*pow(M_PI,3./2.)/3. * \
    (kappaindex_-1.)*(kappaindex_-2.)*kappaindex_/pow(thetae_*kappaindex_,3.)* \
    (2.*tgamma(2.+kappaindex_/2.)/(2.+kappaindex_)-1.)* \
    (pow(3./kappaindex_,19./4.)+3./5.),
    expo = pow(-7./4.+8./5.*kappaindex_,-43./50.),
    As = pow(pow(As_low,-expo) + pow(As_high,-expo),-1./expo);
  
  double abs_synch = numberdensityCGS_*					\
    GYOTO_ELEMENTARY_CHARGE_CGS*GYOTO_ELEMENTARY_CHARGE_CGS/		\
    (nu*GYOTO_ELECTRON_MASS_CGS*GYOTO_C_CGS)*				\
    As;
  
  return abs_synch;
}

void Spectrum::KappaDistributionSynchrotron::radiativeQ(double jnu[], // output
						double alphanu[], // output
						double nu_ems[],
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
