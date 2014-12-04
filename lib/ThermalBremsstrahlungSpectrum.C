/*
  Copyright 2014 Frederic Vincent

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

/// Properties

#include "GyotoProperty.h"
GYOTO_PROPERTY_START(Spectrum::ThermalBremsstrahlung)
GYOTO_PROPERTY_DOUBLE(Spectrum::ThermalBremsstrahlung, Temperature, temperature)
GYOTO_PROPERTY_END(Spectrum::ThermalBremsstrahlung, Generic::properties)

///



Spectrum::ThermalBremsstrahlung::ThermalBremsstrahlung() :
  T_(10000.), massdensityCGS_(0.),
  spectrumBB_(NULL),Spectrum::Generic("ThermalBremsstrahlung") {
  Tm1_=1./T_; Tm05_=sqrt(Tm1_);
  // This awful constant is the constant part of the thermal brems j_nu
  cst_ = 1/(4.*M_PI)
    *(pow(2.,5)*M_PI*pow(GYOTO_ELEMENTARY_CHARGE_CGS,6))
    /(3.*GYOTO_ELECTRON_MASS_CGS*pow(GYOTO_C_CGS,3))
    *sqrt(2*M_PI/(3.*GYOTO_BOLTZMANN_CGS*GYOTO_ELECTRON_MASS_CGS))
    *1./pow(GYOTO_ATOMIC_MASS_UNIT_CGS,2);
  // A BB spectrum is needed to compute alpha_nu=j_nu/BB
  spectrumBB_ = new Spectrum::BlackBody(); 
}

double Spectrum::ThermalBremsstrahlung::temperature() const { return T_; }
void Spectrum::ThermalBremsstrahlung::temperature(double tt) {
  T_ = tt; Tm1_=1./T_; Tm05_=sqrt(Tm1_);
  spectrumBB_->temperature(T_);
}
double Spectrum::ThermalBremsstrahlung::massdensityCGS() const { 
  return massdensityCGS_; }
void Spectrum::ThermalBremsstrahlung::massdensityCGS(double rho) { 
  massdensityCGS_ = rho; }
  
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
  return  cst_*Tm05_*massdensityCGS_*massdensityCGS_
    *exp(-GYOTO_PLANCK_OVER_BOLTZMANN*nu*Tm1_);
}

double Spectrum::ThermalBremsstrahlung::alphanuCGS(double nu) const{
  double BB=(*spectrumBB_)(nu)/GYOTO_INU_CGS_TO_SI; // B_nu in cgs
  if (BB==0.){
    throwError("In ThermalBrems: "
	       "bad temperature");
  }
  // Kirchhoff's law:
  return jnuCGS(nu)/BB;
}
