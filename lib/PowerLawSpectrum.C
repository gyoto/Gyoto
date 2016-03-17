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
#include "GyotoUtils.h"

#include <cstdlib>

#include <cmath>
#ifdef GYOTO_USE_XERCES
#include "GyotoFactory.h"
#endif
using namespace Gyoto;

/// Properties

#include "GyotoProperty.h"
GYOTO_PROPERTY_START(Spectrum::PowerLaw)
GYOTO_PROPERTY_DOUBLE(Spectrum::PowerLaw, Exponent, exponent)
GYOTO_PROPERTY_DOUBLE(Spectrum::PowerLaw, Constant, constant)
GYOTO_PROPERTY_VECTOR_DOUBLE(Spectrum::PowerLaw, CutOffIneV, cutoffinev)
GYOTO_PROPERTY_END(Spectrum::PowerLaw, Generic::properties)

///

Spectrum::PowerLaw::PowerLaw() :
Spectrum::Generic("PowerLaw"), constant_(1.), exponent_(0.),
  minfreq_(DBL_MIN), maxfreq_(DBL_MAX){}
Spectrum::PowerLaw::PowerLaw(double p, double c) :
  Spectrum::Generic("PowerLaw"), constant_(c), exponent_(p),
  minfreq_(DBL_MIN), maxfreq_(DBL_MAX){}
Spectrum::PowerLaw * Spectrum::PowerLaw::clone() const
{ return new Spectrum::PowerLaw(*this); }

double Spectrum::PowerLaw::constant() const { return constant_; }
void Spectrum::PowerLaw::constant(double c) { constant_ = c; }
double Spectrum::PowerLaw::exponent() const { return exponent_; }
void Spectrum::PowerLaw::exponent(double c) { exponent_ = c; }
void Spectrum::PowerLaw::cutoffinev(std::vector<double> const &v) {
  if (v.size() != 2)
    throwError("In PowerLawSpectrum: Only 2 arguments to define"
	       " cutoffs");
  minfreq_ = v[0]*GYOTO_eV2Hz;
  maxfreq_ = v[1]*GYOTO_eV2Hz;
}
std::vector<double> Spectrum::PowerLaw::cutoffinev() const {
  std::vector<double> v (2, 0.);
  v[0]=minfreq_; v[1]=maxfreq_;
  return v;
}

double Spectrum::PowerLaw::operator()(double nu) const {
if (nu<minfreq_ || nu>maxfreq_) return 0.; // cutoffs
  return constant_ * pow(nu, exponent_);
}
