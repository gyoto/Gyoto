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

#include "GyotoKappaDistributionSynchrotronSpectrum.h"
#include "GyotoDefs.h"
#include <cmath>
#ifdef GYOTO_USE_XERCES
#include "GyotoFactory.h"
#include "GyotoFactoryMessenger.h"
#endif
using namespace Gyoto;
using namespace std;

#include "GyotoProperty.h"
GYOTO_PROPERTY_START(Spectrum::KappaDistributionSynchrotron,
         "Kappa synchrotron emission")
GYOTO_PROPERTY_END(Spectrum::KappaDistributionSynchrotron, Generic::properties)



#define nstep_angint 100 // for angle-averaging integration

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
  GYOTO_ERROR("In PLSynch: "
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
  
  double emis_synch = numberdensityCGS_*        \
    GYOTO_ELEMENTARY_CHARGE_CGS*GYOTO_ELEMENTARY_CHARGE_CGS*cyclotron_freq_/ \
    GYOTO_C_CGS*              \
    Js;
  
  //cout << "in kappa spec angleB jnu= " << angle_B_pem_ << " " << emis_synch << endl;
  return emis_synch;
}

double Spectrum::KappaDistributionSynchrotron::jQnuCGS(double nu) const{  

  // Marszewski, Prather, Joshi, Pandya, Gammie 2021

  //cout << "in kappa jnu stuff: " << cyclotron_freq_ << " " << thetae_ << " " << kappaindex_ << endl;
  
  double sinth = sin(angle_B_pem_),
    nuk = cyclotron_freq_*pow(thetae_*kappaindex_,2.)*sinth,
    Xk = nu/nuk,
    Js_low = 0.5*pow(Xk,1./3.)*sinth*4.*M_PI*tgamma(kappaindex_-4./3.)/ \
    (pow(3.,7./3.)*tgamma(kappaindex_-2.)),
    Js_high = (pow(4./5.,2.)+kappaindex_/50.)*pow(Xk,-(kappaindex_-2.)/2.)*sinth*pow(3.,(kappaindex_-1.)/2.)* \
    (kappaindex_-1.)*(kappaindex_-2.)/4.*tgamma(kappaindex_/4.-1./3.)* \
    tgamma(kappaindex_/4.+4./3.),
    expo = 3.7*pow(kappaindex_,-8./5.),
    Js = pow(pow(Js_low,-expo) + pow(Js_high,-expo),-1./expo);
  
  double emis_synch = numberdensityCGS_*        \
    GYOTO_ELEMENTARY_CHARGE_CGS*GYOTO_ELEMENTARY_CHARGE_CGS*cyclotron_freq_/ \
    GYOTO_C_CGS*              \
    Js;
  
  //cout << "in kappa spec angleB jnu= " << angle_B_pem_ << " " << emis_synch << endl;
  return emis_synch;
}

double Spectrum::KappaDistributionSynchrotron::jUnuCGS(double nu) const{  

  // Marszewski, Prather, Joshi, Pandya, Gammie 2021
  return 0.;
}

double Spectrum::KappaDistributionSynchrotron::jVnuCGS(double nu) const{  

  // Marszewski, Prather, Joshi, Pandya, Gammie 2021

  //cout << "in kappa jnu stuff: " << cyclotron_freq_ << " " << thetae_ << " " << kappaindex_ << endl;
  
  double sinth = sin(angle_B_pem_),
    nuk = cyclotron_freq_*pow(thetae_*kappaindex_,2.)*sinth,
    Xk = nu/nuk,
    Js_low = pow(3./4.,2.)*pow(pow(sinth,-12./5.)-1.,12./25.)*pow(kappaindex_, -66./125.)/thetae_*\
    pow(Xk,-7./20.)*pow(Xk,1./3.)*sinth*4.*M_PI*tgamma(kappaindex_-4./3.)/ \
    (pow(3.,7./3.)*tgamma(kappaindex_-2.)),
    Js_high = pow(7./8.,2.)*pow(pow(sinth,-5./2.)-1.,11./25.)*pow(kappaindex_, -11./25.)/thetae_*\
    pow(Xk,-1./2.)*pow(Xk,-(kappaindex_-2.)/2.)*sinth*pow(3.,(kappaindex_-1.)/2.)* \
    (kappaindex_-1.)*(kappaindex_-2.)/4.*tgamma(kappaindex_/4.-1./3.)* \
    tgamma(kappaindex_/4.+4./3.),
    expo = 3.*pow(kappaindex_,-3./2.),
    Js = pow(pow(Js_low,-expo) + pow(Js_high,-expo),-1./expo)*cos(angle_B_pem_)/abs(cos(angle_B_pem_));
  
  double emis_synch = numberdensityCGS_*        \
    GYOTO_ELEMENTARY_CHARGE_CGS*GYOTO_ELEMENTARY_CHARGE_CGS*cyclotron_freq_/ \
    GYOTO_C_CGS*              \
    Js;
  
  //cout << "in kappa spec angleB jnu= " << angle_B_pem_ << " " << emis_synch << endl;
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
  
  double abs_synch = numberdensityCGS_*         \
    GYOTO_ELEMENTARY_CHARGE_CGS*GYOTO_ELEMENTARY_CHARGE_CGS/    \
    (nu*GYOTO_ELECTRON_MASS_CGS*GYOTO_C_CGS)*       \
    As;
  
  //cout << "in kappa spec angleB anu= " << angle_B_pem_ << " " << abs_synch << endl;

  return abs_synch;
}

double Spectrum::KappaDistributionSynchrotron::alphaQnuCGS(double nu) const{

  // Marszewski, Prather, Joshi, Pandya, Gammie 2021

  //cout << "in kappa anu stuff: " << cyclotron_freq_ << " " << thetae_ << " " << kappaindex_ << " " << hypergeometric_ << endl;

  double sinth = sin(angle_B_pem_),
    nuk = cyclotron_freq_*pow(thetae_*kappaindex_,2.)*sinth,
    Xk = nu/nuk,
    As_low = 25./48.*pow(Xk,-2./3.)*pow(3.,1./6.)*10./41.* \
    2.*M_PI/pow(thetae_*kappaindex_,10./3.-kappaindex_) * \
    (kappaindex_-1.)*(kappaindex_-2.)*kappaindex_/(3.*kappaindex_-1.) * \
    tgamma(5./3.)*hypergeometric_,
    As_high = pow(Xk,-(1.+kappaindex_)/2.)*pow(M_PI,3./2.)/3. * \
    (kappaindex_-1.)*(kappaindex_-2.)*kappaindex_/pow(thetae_*kappaindex_,3.)* \
    (2.*tgamma(2.+kappaindex_/2.)/(2.+kappaindex_)-1.)* \
    (pow(21.,2.)*pow(kappaindex_,-pow(12./5.,2.))+11./20.),
    expo = 7./5.*pow(kappaindex_,-23./20.),
    As = pow(pow(As_low,-expo) + pow(As_high,-expo),-1./expo);
  
  double abs_synch = numberdensityCGS_*         \
    GYOTO_ELEMENTARY_CHARGE_CGS*GYOTO_ELEMENTARY_CHARGE_CGS/    \
    (nu*GYOTO_ELECTRON_MASS_CGS*GYOTO_C_CGS)*       \
    As;
  
  //cout << "in kappa spec angleB anu= " << angle_B_pem_ << " " << abs_synch << endl;

  return abs_synch;
}

double Spectrum::KappaDistributionSynchrotron::alphaUnuCGS(double nu) const{

  // Marszewski, Prather, Joshi, Pandya, Gammie 2021
  return 0.;
}

double Spectrum::KappaDistributionSynchrotron::alphaVnuCGS(double nu) const{

  // Marszewski, Prather, Joshi, Pandya, Gammie 2021

  //cout << "in kappa anu stuff: " << cyclotron_freq_ << " " << thetae_ << " " << kappaindex_ << " " << hypergeometric_ << endl;

  double sinth = sin(angle_B_pem_),
    nuk = cyclotron_freq_*pow(thetae_*kappaindex_,2.)*sinth,
    Xk = nu/nuk,
    As_low = 77./100./thetae_*pow(pow(sinth,-114./50.)-1.,223./500.)*pow(kappaindex_,-7./10.) *\
    pow(Xk,-7./20.)*pow(Xk,-2./3.)*pow(3.,1./6.)*10./41.* \
    2.*M_PI/pow(thetae_*kappaindex_,10./3.-kappaindex_) * \
    (kappaindex_-1.)*(kappaindex_-2.)*kappaindex_/(3.*kappaindex_-1.) * \
    tgamma(5./3.)*hypergeometric_,
    As_high = pow(Xk,-(1.+kappaindex_)/2.)*pow(M_PI,3./2.)/3. * \
    (kappaindex_-1.)*(kappaindex_-2.)*kappaindex_/pow(thetae_*kappaindex_,3.)* \
    (2.*tgamma(2.+kappaindex_/2.)/(2.+kappaindex_)-1.)* \
    143./10.*pow(thetae_,-116./125.)*pow(pow(sinth,-41./20.)-1.,1./2.) *\
    (pow(13.,2.)*pow(kappaindex_,-8.)+13./2500.*kappaindex_-263./5000.+47./200./kappaindex_)*pow(Xk,-0.5),
    expo = 61./50.*pow(kappaindex_,-142./125.)+7./1000.,
    As = pow(pow(As_low,-expo) + pow(As_high,-expo),-1./expo)*cos(angle_B_pem_)/abs(cos(angle_B_pem_));
  
  double abs_synch = numberdensityCGS_*         \
    GYOTO_ELEMENTARY_CHARGE_CGS*GYOTO_ELEMENTARY_CHARGE_CGS/    \
    (nu*GYOTO_ELECTRON_MASS_CGS*GYOTO_C_CGS)*       \
    As;
  
  //cout << "in kappa spec angleB anu= " << angle_B_pem_ << " " << abs_synch << endl;

  return abs_synch;
}

double Spectrum::KappaDistributionSynchrotron::rQnuCGS(double nu) const{
  
  // Marszewski, Prather, Joshi, Pandya, Gammie 2021
  double rho_Q=0;
  double sinth = sin(angle_B_pem_),
    nuk = cyclotron_freq_*pow(thetae_*kappaindex_,2.)*sinth,
    Xk = nu/nuk,
    rho_kappa=0;
    if (Xk<0.1)
      GYOTO_ERROR("Xk too low to compute rhoQ with these formula");
    switch(int(kappaindex_*10))
    {
      case 35:
        rho_kappa=17.*thetae_-3.*sqrt(thetae_)+7.*sqrt(thetae_)*exp(-5.*thetae_)*\
          (1-exp(-pow(Xk,0.84)/30.)-sin(Xk/10.)*exp(-3./2.*pow(Xk,0.471)));
          break;
      case 40:
        rho_kappa=46./3.*thetae_-5./3.*sqrt(thetae_)+17./3.*sqrt(thetae_)*exp(-5.*thetae_)*\
          (1-exp(-pow(Xk,0.84)/18.)-sin(Xk/6.)*exp(-7./4.*pow(Xk,0.5)));
          break;
      case 45:
        rho_kappa=14.*thetae_-13./8.*sqrt(thetae_)+9./2.*sqrt(thetae_)*exp(-5.*thetae_)*\
          (1-exp(-pow(Xk,0.84)/12.)-sin(Xk/4.)*exp(-2.*pow(Xk,0.525)));
          break;
      case 50:
        rho_kappa=25./2.*thetae_-sqrt(thetae_)+5.*sqrt(thetae_)*exp(-5.*thetae_)*\
          (1-exp(-pow(Xk,0.84)/8.)-sin(3.*Xk/8.)*exp(-9./4.*pow(Xk,0.541)));
          break;
      default:
        GYOTO_ERROR("Faraday coefficients not defined for values of kappa different of 3.5, 4., 4.5, 5.");
    }
  rho_Q=numberdensityCGS_*pow(GYOTO_ELEMENTARY_CHARGE_CGS*cyclotron_freq_*sinth,2.)/ \
    (GYOTO_ELECTRON_MASS_CGS*GYOTO_C_CGS*pow(nu,3.))* \
    rho_kappa;

  return rho_Q;
}

double Spectrum::KappaDistributionSynchrotron::rUnuCGS(double nu) const{
  // Marszewski, Prather, Joshi, Pandya, Gammie 2021
  return 0.;
}

double Spectrum::KappaDistributionSynchrotron::rVnuCGS(double nu) const{
  
  // Marszewski, Prather, Joshi, Pandya, Gammie 2021
  double rho_V=0;
  double sinth = sin(angle_B_pem_),
    nuk = cyclotron_freq_*pow(thetae_*kappaindex_,2.)*sinth,
    Xk = nu/nuk,
    rho_kappa=0;
    if (Xk<0.1)
      GYOTO_ERROR("Xk too low to compute rhoQ with these formula");
    switch(int(kappaindex_*10))
    {
      case 35:
        rho_kappa=(pow(thetae_,2.)+2.*thetae_+1.)/ \
          (25./8.*pow(thetae_,2.)+4.*thetae_+1.)* \
          (1-0.17*log(1.+0.447*pow(Xk,-0.5)));
          break;
      case 40:
        rho_kappa=(pow(thetae_,2.)+54.*thetae_+50.)/ \
          (30./11.*pow(thetae_,2.)+134.*thetae_+50.)* \
          (1-0.17*log(1.+0.391*pow(Xk,-0.5)));
          break;
      case 45:
        rho_kappa=(pow(thetae_,2.)+43.*thetae_+38.)/ \
          (7./3.*pow(thetae_,2.)+185./2.*thetae_+38.)* \
          (1-0.17*log(1.+0.348*pow(Xk,-0.5)));
          break;
      case 50:
        rho_kappa=(thetae_+13./14.)/ \
          (2.*thetae_+13./14.)* \
          (1-0.17*log(1.+0.313*pow(Xk,-0.5)));
          break;
      default:
        GYOTO_ERROR("Faraday coefficients not defined for values of kappa different of 3.5, 4., 4.5, 5.");
    }
  rho_V=2.*numberdensityCGS_*pow(GYOTO_ELEMENTARY_CHARGE_CGS,2.)*cyclotron_freq_*cos(angle_B_pem_)*bessk0(1./thetae_)/ \
    (GYOTO_ELECTRON_MASS_CGS*GYOTO_C_CGS*pow(nu,2.)*bessk(2.,1./thetae_))* \
    rho_kappa;

  return rho_V;
}

void Spectrum::KappaDistributionSynchrotron::radiativeQ(double jnu[], // output
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
        theta=th0+double(jj)*hh;
        angle_B_pem(theta);
        jnusinnext=jnuCGS(nu)*sin(theta);
        anusinnext=alphanuCGS(nu)*sin(theta);
        jnucur+=0.5*0.5*hh*(jnusinprev+jnusinnext);
        anucur+=0.5*0.5*hh*(anusinprev+anusinnext);
        jnusinprev=jnusinnext;
        anusinprev=anusinnext;
        //NB: averaged jnu is: \int jnu dOmega = 1/2 * \int jnu*sinth dth
      }
    }
    
    // OUTPUTS
    jnu[ii]= jnucur * GYOTO_JNU_CGS_TO_SI;
    alphanu[ii]= anucur * GYOTO_ANU_CGS_TO_SI;
    
  }
}

void Spectrum::KappaDistributionSynchrotron::radiativeQ(double jInu[], double jQnu[], double jUnu[], double jVnu[], // Output
        double aInu[], double aQnu[], double aUnu[], double aVnu[], // Output
        double rQnu[], double rUnu[], double rVnu[], // Output
        double const nu_ems[],
        size_t nbnu ){  
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