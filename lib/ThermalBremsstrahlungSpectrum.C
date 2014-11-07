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

Spectrum::ThermalBremsstrahlung::ThermalBremsstrahlung() :
  spectrumBB_(NULL),Spectrum::Generic("ThermalBremsstrahlung") {
  // This awful constant is the constant part of the thermal brems j_nu
  cst_ = 1/(4.*M_PI)
    *(pow(2.,5)*M_PI*pow(GYOTO_ELEMENTARY_CHARGE_CGS,6))
    /(3.*GYOTO_ELECTRON_MASS_CGS*pow(GYOTO_C_CGS,3))
    *sqrt(2*M_PI/(3.*GYOTO_BOLTZMANN_CGS*GYOTO_ELECTRON_MASS_CGS))
    *1./pow(GYOTO_ATOMIC_MASS_UNIT_CGS,2);
  // A BB spectrum is needed to compute alpha_nu=j_nu/BB
  spectrumBB_ = new Spectrum::BlackBody(); 
}
  
Spectrum::ThermalBremsstrahlung * Spectrum::ThermalBremsstrahlung::clone() const
{ return new Spectrum::ThermalBremsstrahlung(*this); }

double Spectrum::ThermalBremsstrahlung::operator()(double nu) const {
  return 0.; // what should this return in opt thin case?
}

double Spectrum::ThermalBremsstrahlung::jnu(double nu, double temp, 
					    double massdensity) {
  return  cst_*1./sqrt(temp)*massdensity*massdensity
    *exp(-GYOTO_PLANCK_CGS*nu/(GYOTO_BOLTZMANN_CGS*temp));
}

double Spectrum::ThermalBremsstrahlung::alphanu(double nu, double temp, 
					    double massdensity) {
  spectrumBB_->temperature(temp);
  double BB=(*spectrumBB_)(nu)/GYOTO_INU_CGS_TO_SI; // B_nu in cgs
  if (BB==0.){
    throwError("In ThermalBrems: "
	       "bad temperature");
  }
  // Kirchhoff's law:
  return jnu(nu,temp,massdensity)/BB;
}

#ifdef GYOTO_USE_XERCES
void Spectrum::ThermalBremsstrahlung::fillElement(FactoryMessenger *fmp) const {
  Spectrum::Generic::fillElement(fmp);
}

#endif
