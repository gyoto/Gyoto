/*
    Copyright 2011-2012, 2014, 2017, 2019-2020 Thibaut Paumard & Frederic Vincent

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

/// Properties

#include "GyotoProperty.h"
GYOTO_PROPERTY_START(Spectrum::BlackBody)
GYOTO_PROPERTY_DOUBLE(Spectrum::BlackBody, Temperature, temperature)
GYOTO_PROPERTY_DOUBLE(Spectrum::BlackBody, Scaling, scaling)
GYOTO_PROPERTY_DOUBLE(Spectrum::BlackBody, ColorCorrection, colorCorrection)
GYOTO_PROPERTY_END(Spectrum::BlackBody, Generic::properties)

///

Spectrum::BlackBody::BlackBody() :
Spectrum::Generic("BlackBody"), T_(10000.), colorcor_(1.), colorcorm4_(1.),
  cst_(2.*GYOTO_PLANCK_OVER_C_SQUARE) {Tm1_=1./T_;}
Spectrum::BlackBody::BlackBody(double T, double c) :
  Spectrum::Generic("BlackBody"), T_(T), cst_(c), colorcor_(1.), colorcorm4_(1.)
{Tm1_=1./T_;}
Spectrum::BlackBody * Spectrum::BlackBody::clone() const
{ return new Spectrum::BlackBody(*this); }

double Spectrum::BlackBody::temperature() const { return T_; }
void Spectrum::BlackBody::temperature(double c) { T_ = c; Tm1_=1./T_; }
double Spectrum::BlackBody::scaling() const { return cst_; }
void Spectrum::BlackBody::scaling(double c) { cst_ = c; }
double Spectrum::BlackBody::colorCorrection() const { return colorcor_; }
void Spectrum::BlackBody::colorCorrection(double c) {
  colorcor_ = c; colorcorm4_ = 1./(c*c*c*c); T_*=c; Tm1_/=c;}

double Spectrum::BlackBody::operator()(double nu) const {
  return  colorcorm4_*cst_*nu*nu*nu
    /(expm1(GYOTO_PLANCK_OVER_BOLTZMANN*nu*Tm1_));
}

