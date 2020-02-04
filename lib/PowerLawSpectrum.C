/*
    Copyright 2011-2012, 2014, 2016, 2019-2020 Thibaut Paumard & Frédéric Vincent

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

#include <cmath>
#ifdef GYOTO_USE_XERCES
#include "GyotoFactory.h"
#endif
using namespace Gyoto;

/// Properties

#include "GyotoProperty.h"
GYOTO_PROPERTY_START(Spectrum::PowerLaw,
		     "'Constant * nu[Hz]^Exponent' between CutOff[0] and CutOff[1]")
GYOTO_PROPERTY_DOUBLE(Spectrum::PowerLaw, Exponent, exponent, "Exponent of power law")
GYOTO_PROPERTY_DOUBLE(Spectrum::PowerLaw, Constant, constant, "Constant in front of power law")
GYOTO_PROPERTY_VECTOR_DOUBLE_UNIT(Spectrum::PowerLaw, CutOff, cutoff,
				  "Cut-off frequencies in any unit convertible to Hz, m or eV "
				  "(default: '0 DBL_MAX'; default unit: Hz).")
GYOTO_PROPERTY_END(Spectrum::PowerLaw, Generic::properties)

///

Spectrum::PowerLaw::PowerLaw() :
Spectrum::Generic("PowerLaw"), constant_(1.), exponent_(0.),
  minfreq_(0.), maxfreq_(DBL_MAX){}
Spectrum::PowerLaw::PowerLaw(double p, double c) :
  Spectrum::Generic("PowerLaw"), constant_(c), exponent_(p),
  minfreq_(0.), maxfreq_(DBL_MAX){}
Spectrum::PowerLaw * Spectrum::PowerLaw::clone() const
{ return new Spectrum::PowerLaw(*this); }

double Spectrum::PowerLaw::constant() const { return constant_; }
void Spectrum::PowerLaw::constant(double c) { constant_ = c; }
double Spectrum::PowerLaw::exponent() const { return exponent_; }
void Spectrum::PowerLaw::exponent(double c) { exponent_ = c; }
void Spectrum::PowerLaw::cutoff(std::vector<double> const &v, std::string const &u) {
  cutoff({Units::ToHerz(v[0], u), Units::ToHerz(v[1], u)});
}
void Spectrum::PowerLaw::cutoff(std::vector<double> const &v) {
  if (v.size() != 2)
    GYOTO_ERROR("CutOff needs exactly two cut-off frequencies");
  minfreq_ = v[0];
  maxfreq_ = v[1];
  if (minfreq_ > maxfreq_) {
    double tmp = minfreq_;
    minfreq_ = maxfreq_;
    maxfreq_ = tmp;
  }
}

std::vector<double> Spectrum::PowerLaw::cutoff(std::string const &u) const {
  return {Units::FromHerz(minfreq_, u), Units::FromHerz(maxfreq_, u)};
}
std::vector<double> Spectrum::PowerLaw::cutoff() const {
  return {minfreq_, maxfreq_};
}

double Spectrum::PowerLaw::operator()(double nu) const {
if (nu<minfreq_ || nu>maxfreq_) return 0.; // cutoffs
  return constant_ * pow(nu, exponent_);
}
