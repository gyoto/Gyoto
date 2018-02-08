/*
  Copyright 2014-2018 Frederic Vincent

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

#include "GyotoThermalBremsstrahlungSpectrum.h"
#include "GyotoDefs.h"
#include <cmath>
#include <cstdlib> /* atof */
#ifdef GYOTO_USE_XERCES
#include "GyotoFactory.h"
#include "GyotoFactoryMessenger.h"
#endif
using namespace Gyoto;

#include "GyotoProperty.h"
GYOTO_PROPERTY_START(Spectrum::ThermalBremsstrahlung,
		     "Thermal bremsstrahlung emission")
GYOTO_PROPERTY_END(Spectrum::ThermalBremsstrahlung, Generic::properties)

// This awful constant is the constant part of the thermal brems j_nu
static double const cst_ = 1/(4.*M_PI)
    *(pow(2.,5)*M_PI*pow(GYOTO_ELEMENTARY_CHARGE_CGS,6))
    /(3.*GYOTO_ELECTRON_MASS_CGS*pow(GYOTO_C_CGS,3))
    *sqrt(2*M_PI/(3.*GYOTO_BOLTZMANN_CGS*GYOTO_ELECTRON_MASS_CGS));

Spectrum::ThermalBremsstrahlung::ThermalBremsstrahlung()
: Spectrum::Generic("ThermalBremsstrahlung"),
  spectrumBB_(NULL), T_(10000.), numberdensityCGS_(0.)
{
  Tm1_=1./T_; Tm05_=sqrt(Tm1_);
  // A BB spectrum is needed to compute alpha_nu=j_nu/BB
  spectrumBB_ = new Spectrum::BlackBody(); 
}

double Spectrum::ThermalBremsstrahlung::temperature() const { return T_; }
void Spectrum::ThermalBremsstrahlung::temperature(double tt) {
  T_ = tt; Tm1_=1./T_; Tm05_=sqrt(Tm1_);
  spectrumBB_->temperature(T_);
}
double Spectrum::ThermalBremsstrahlung::numberdensityCGS() const { 
  return numberdensityCGS_; }
void Spectrum::ThermalBremsstrahlung::numberdensityCGS(double rho) { 
  numberdensityCGS_ = rho;
}
  
Spectrum::ThermalBremsstrahlung * Spectrum::ThermalBremsstrahlung::clone() const
{ return new Spectrum::ThermalBremsstrahlung(*this); }

double Spectrum::ThermalBremsstrahlung::operator()(double nu) const {
  throwError("In ThermalBrems: "
	     "Bremsstrahlung emission not defined for optically thick case");
  return 0.;
}
double Spectrum::ThermalBremsstrahlung::operator()(double nu, 
						   double , 
						   double ds) const{
  double dsCGS = ds*100.; // ds should be given in SI
  // Returns intensity increment in SI:
  return jnuCGS(nu)*dsCGS*exp(-alphanuCGS(nu)*dsCGS)*GYOTO_INU_CGS_TO_SI;
}

double Spectrum::ThermalBremsstrahlung::jnuCGS(double nu) const{
  //std::cout << "in brems cst,ne,Te= " << cst_ << " " <<  numberdensityCGS_ << " " << T_ << std::endl;
  double temp_e = T_*GYOTO_BOLTZMANN_CGS/(GYOTO_ELECTRON_MASS_CGS*GYOTO_C2_CGS);
  double jnu=0.;
  if (temp_e < 0.01){
    // Brems for non-relat emitter
    jnu = cst_*Tm05_*numberdensityCGS_*numberdensityCGS_
      *exp(-GYOTO_PLANCK_OVER_BOLTZMANN*nu*Tm1_);
  } else {
    // Brems for relat emitter
    // Formula from Stepney&Guilbert (1983) eq.2-3; Svensson (1982) eqs.19-20
    // TO BE DOUBLE-CHECKED
    double Fee=0., Fei=0.;
    if (temp_e < 1.) {
      Fee = 20./(9.*sqrt(M_PI))*(44.- 3.*M_PI*M_PI)*pow(temp_e,3./2.)
	*(1.+1.1*temp_e+temp_e*temp_e - 1.25*pow(temp_e,2.5));
      Fei = 4.*sqrt(2.*temp_e / (M_PI*M_PI*M_PI))*(1.+1.781*pow(temp_e,1.34));
    }else{
      Fee = 24.*temp_e*(log(2.*exp(-GYOTO_EULER_MASCHERONI)*temp_e)+1.25);
      Fei = 9.*temp_e / (2.*M_PI)*(log(1.123*temp_e + 0.42) + 1.5);
    }
    
    double fee = numberdensityCGS_*numberdensityCGS_*GYOTO_C_CGS
      *GYOTO_ELECTRON_CLASSICAL_RADIUS_CGS
      *GYOTO_ELECTRON_CLASSICAL_RADIUS_CGS
      *GYOTO_ELECTRON_MASS_CGS*GYOTO_C2_CGS*GYOTO_ALPHA_F*Fee;
    double fei = numberdensityCGS_*numberdensityCGS_*GYOTO_THOMSON_CGS
      *GYOTO_C_CGS*GYOTO_ALPHA_F*GYOTO_ELECTRON_MASS_CGS*GYOTO_C2_CGS*Fei;
    
    double gaunt=0.;
    double fact=T_*GYOTO_BOLTZMANN_CGS/(GYOTO_PLANCK_CGS*nu);
    if (fact<=1.) gaunt=sqrt(3./M_PI*fact);
    else gaunt=sqrt(3.)/M_PI*log(4./GYOTO_EULER_MASCHERONI*fact) ;
    jnu = 1./(4.*M_PI)				
      *GYOTO_PLANCK_OVER_BOLTZMANN*1./T_		
      *exp(-GYOTO_PLANCK_OVER_BOLTZMANN*nu/T_)	
      *(fei+fee)*gaunt;
  }
  return jnu;
}

double Spectrum::ThermalBremsstrahlung::alphanuCGS(double nu) const{
  double BB  = (*spectrumBB_)(nu)/GYOTO_INU_CGS_TO_SI; // B_nu in cgs
  double jnu = jnuCGS(nu);
  if (BB==0.){
    if (jnu==0.) return 0.;
    else throwError("In ThermalBrems: alphanu undefined!");
  }
  // Kirchhoff's law:
  return jnuCGS(nu)/BB;
}

void Spectrum::ThermalBremsstrahlung::radiativeQ(double jnu[], // output
						double alphanu[], // output
						double nu_ems[],
						size_t nbnu
						) {
  for (size_t ii=0; ii< nbnu; ++ii){
    
    double nu = nu_ems[ii];
    double BB  = (*spectrumBB_)(nu);
    
    jnu[ii]=this->jnuCGS(nu)*GYOTO_JNU_CGS_TO_SI;
    if (BB==0.){
      if (jnu[ii]==0.) alphanu[ii]=0.;
      else throwError("In ThermalBrems: alphanu undefined!");
    }else
      alphanu[ii]=jnu[ii]/BB;
    
  }
  
}
