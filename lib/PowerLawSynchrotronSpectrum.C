/*
  Copyright 2018 Frederic Vincent

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
#include <cstdlib> /* atof */
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

Spectrum::PowerLawSynchrotron::PowerLawSynchrotron()
: Spectrum::Generic("PowerLawSynchrotron"),
  numberdensityCGS_(0.),
  angle_B_pem_(0.), cyclotron_freq_(1.),
  PLindex_(0.), angle_averaged_(0)
{}

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
  
Spectrum::PowerLawSynchrotron * Spectrum::PowerLawSynchrotron::clone() const
{ return new Spectrum::PowerLawSynchrotron(*this); }

double Spectrum::PowerLawSynchrotron::operator()(double nu) const {
  throwError("In PLSynch: "
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

double Spectrum::PowerLawSynchrotron::jnuCGS(double nu) const{
  /* 
     From Petrosian & McTiernan 1983, Phys. Fluids 26 (10), eq. 32
     Putting g(mu)=1 and using (Y+ + Y_)=2 to get jnu and alphanu.
     NB: putting g(mu)=1 or 1/2 is not important, it boils down
     to redefining the % amount delta of PL energy wrt THER energy
  */
  //std::cout << "PL synch stuff= " << cyclotron_freq_ << " " << angle_B_pem_ << " " << PLindex_ << " " << numberdensityCGS_ << " " << angle_averaged_ << std::endl;
  double emis_synch =
    sqrt(3.)*M_PI*GYOTO_ELEMENTARY_CHARGE_CGS*GYOTO_ELEMENTARY_CHARGE_CGS
    *cyclotron_freq_*sin(angle_B_pem_)/(2.*GYOTO_C_CGS)
    *numberdensityCGS_*(PLindex_-1.)
    *pow(3.*cyclotron_freq_*(PLindex_+1.)
	 *sin(angle_B_pem_)/(4.*nu),0.5*(PLindex_-1.))
    *exp(-0.5*(PLindex_+1.));
  
  return emis_synch;
}

double Spectrum::PowerLawSynchrotron::alphanuCGS(double nu) const{
  // From Petrosian & McTiernan 1983, Phys. Fluids 26 (10), eq. 32
  double abs_synch =
    sqrt(3.)*M_PI*GYOTO_ELEMENTARY_CHARGE_CGS*GYOTO_ELEMENTARY_CHARGE_CGS
    *cyclotron_freq_*sin(angle_B_pem_)/(2.*GYOTO_C_CGS)
    *numberdensityCGS_*(PLindex_-1.)
    *pow(3.*cyclotron_freq_*(PLindex_+2.)*sin(angle_B_pem_)
	 /(4.*nu),0.5*PLindex_)
    *exp(-0.5*(PLindex_+2.))
    *(PLindex_+2.)
    /(GYOTO_ELECTRON_MASS_CGS*nu*nu);
  
  return abs_synch;
}

void Spectrum::PowerLawSynchrotron::radiativeQ(double jnu[], // output
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
      double th0=0., thNm1=M_PI;
      double hh=(thNm1-th0)/double(nstep_angint);
      for (int jj=1;jj<=2*nstep_angint-3;jj+=2){
	double theta=th0+double(jj)/2.*hh;
	angle_B_pem(theta);
	jnucur+=0.5*hh*jnuCGS(nu)*sin(theta);
	//NB: averaged jnu is: \int jnu dOmega = 1/2 * \int jnu*sinth dth
	anucur+=0.5*hh*alphanuCGS(nu)*sin(theta);
      }
    }
    
    // OUTPUTS
    jnu[ii]= jnucur * GYOTO_JNU_CGS_TO_SI;
    alphanu[ii]= anucur * GYOTO_ANU_CGS_TO_SI;
    
  }
}
