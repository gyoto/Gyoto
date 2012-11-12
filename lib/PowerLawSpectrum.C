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

#include "GyotoPowerLawSpectrum.h"
#include "GyotoFactoryMessenger.h"

#include <cmath>
#ifdef GYOTO_USE_XERCES
#include "GyotoFactory.h"
#endif
using namespace Gyoto;

Spectrum::PowerLaw::PowerLaw() :
  Spectrum::Generic("PowerLaw"), constant_(1.), exponent_(0.) {}
Spectrum::PowerLaw::PowerLaw(double p, double c) :
  Spectrum::Generic("PowerLaw"), constant_(c), exponent_(p) {}
Spectrum::PowerLaw * Spectrum::PowerLaw::clone() const
{ return new Spectrum::PowerLaw(*this); }

double Spectrum::PowerLaw::getConstant() const { return constant_; }
void Spectrum::PowerLaw::setConstant(double c) { constant_ = c; }
double Spectrum::PowerLaw::getExponent() const { return exponent_; }
void Spectrum::PowerLaw::setExponent(double c) { exponent_ = c; }

double Spectrum::PowerLaw::operator()(double nu) const {
  return constant_ * pow(nu, exponent_);
}

#ifdef GYOTO_USE_XERCES
void Spectrum::PowerLaw::fillElement(FactoryMessenger *fmp) const {
  fmp->setParameter("Exponent", exponent_);
  fmp->setParameter("Constant", constant_);
  Spectrum::Generic::fillElement(fmp);
}

void Gyoto::Spectrum::PowerLaw::setParameter(std::string name,
			      std::string content,
			      std::string unit) {
  char * tc=const_cast<char*>(content.c_str());
  if (name=="Exponent") setExponent(atof(tc));
  else if (name=="Constant") setConstant(atof(tc));
  else Spectrum::Generic::setParameter(name, content, unit);
}

#endif

