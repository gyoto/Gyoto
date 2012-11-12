/*
    Copyright 2011 Thibaut Paumard

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

#include "GyotoBlackBodySpectrum.h"
#include "GyotoDefs.h"
#include <cmath>
#ifdef GYOTO_USE_XERCES
#include "GyotoFactory.h"
#include "GyotoFactoryMessenger.h"
#endif
using namespace Gyoto;

Spectrum::BlackBody::BlackBody() :
  Spectrum::Generic("BlackBody"), T_(10000.),
  cst_(2.*GYOTO_PLANCK_OVER_C_SQUARE) {Tm1_=1./T_;}
Spectrum::BlackBody::BlackBody(double T, double c) :
  Spectrum::Generic("BlackBody"), T_(T), cst_(c) {Tm1_=1./T_;}
Spectrum::BlackBody * Spectrum::BlackBody::clone() const
{ return new Spectrum::BlackBody(*this); }

double Spectrum::BlackBody::getTemperature() const { return T_; }
void Spectrum::BlackBody::setTemperature(double c) { T_ = c; Tm1_=1./T_; }
double Spectrum::BlackBody::getScaling() const { return cst_; }
void Spectrum::BlackBody::setScaling(double c) { cst_ = c; }

double Spectrum::BlackBody::operator()(double nu) const {
  return  cst_*nu*nu*nu
    /(exp(GYOTO_PLANCK_OVER_BOLTZMANN*nu*Tm1_)-1.);
}

void Spectrum::BlackBody::setParameter(std::string name,
				       std::string content,
				       std::string unit) {
  char * tc=const_cast<char*>(content.c_str());
  if (name=="Temperature") setTemperature(atof(tc));
  else if (name=="Scaling") setScaling(atof(tc));
  else Spectrum::Generic::setParameter(name, content, unit);
}

#ifdef GYOTO_USE_XERCES
void Spectrum::BlackBody::fillElement(FactoryMessenger *fmp) const {
  fmp->setParameter("Temperature", T_);
  fmp->setParameter("Scaling", cst_);
  Spectrum::Generic::fillElement(fmp);
}

#endif
