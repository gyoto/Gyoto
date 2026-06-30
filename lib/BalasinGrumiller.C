/*
    Copyright 2025 Filipe Costa, Frédéric Vincent, Thibaut Paumard

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


#include <iostream>
#include <string>
#include <sstream>


#include "GyotoUtils.h"
#include "GyotoWorldline.h"
#include "GyotoBalasinGrumiller.h"
#include "GyotoError.h"
#include "GyotoProperty.h"
#include <cmath>

using namespace std ; 
using namespace Gyoto ; 
using namespace Gyoto::Metric ;

// **** Section I - Contructor, metric and Christoffel symbols ******
GYOTO_PROPERTY_START(BalasinGrumiller, "your custom description")
GYOTO_PROPERTY_DOUBLE(BalasinGrumiller, V0, V0value,
		      "V0 (adimensioned, 0.000733333).")
GYOTO_PROPERTY_DOUBLE(BalasinGrumiller, R, Rvalue, "R (kpc, 100).")
GYOTO_PROPERTY_DOUBLE(BalasinGrumiller, r0, r0value, "r0 (kpc, 1).")
GYOTO_PROPERTY_END(BalasinGrumiller, Generic::properties)

BalasinGrumiller::BalasinGrumiller()
  : Gyoto::Metric::Generic(GYOTO_COORDKIND_SPHERICAL, "BalasinGrumiller"),
V0value_(0.000733333), Rvalue_(100.), r0value_(1.) {
  } 

//cloner
BalasinGrumiller* BalasinGrumiller::clone() const { return new BalasinGrumiller(*this); }

// Mutators
void BalasinGrumiller::V0value(const double a) {
  V0value_=a;
 tellListeners();
}

void BalasinGrumiller::Rvalue(const double a) {
  Rvalue_=a;
 tellListeners();
}

void BalasinGrumiller::r0value(const double a) {
  r0value_=a;
 tellListeners();
}

// Accessors
double BalasinGrumiller::V0value() const { return V0value_ ; }
double BalasinGrumiller::Rvalue() const { return Rvalue_ ; }
double BalasinGrumiller::r0value() const { return r0value_ ; }

void BalasinGrumiller::gmunu(double g[4][4], const double * pos) const
{
double r = pos[1]; // Extract radial coordinate
double theta = pos[2]; // Extract theta coordinate
  double W = (Rvalue_ - r0value_) * V0value_ + (V0value_ * (-sqrt(Rvalue_*Rvalue_ + r*r - 2.0*Rvalue_*r*cos(theta)) 
         - sqrt(Rvalue_*Rvalue_ + r*r + 2.0*Rvalue_*r*cos(theta)) 
         + sqrt(r0value_*r0value_ + r*r - 2.0*r0value_*r*cos(theta)) 
         + sqrt(r0value_*r0value_ + r*r + 2.0*r0value_*r*cos(theta)))) / 2.0;

  // double DWDr = (V0value_ * (-(r - Rvalue_*cos(theta)) / sqrt(Rvalue_*Rvalue_ + r*r - 2.0*Rvalue_*r*cos(theta))
  //        - (r + Rvalue_*cos(theta)) / sqrt(Rvalue_*Rvalue_ + r*r + 2.0*Rvalue_*r*cos(theta))
  //        + (r - r0value_*cos(theta)) / sqrt(r0value_*r0value_ + r*r - 2.0*r0value_*r*cos(theta))
  //        + (r + r0value_*cos(theta)) / sqrt(r0value_*r0value_ + r*r + 2.0*r0value_*r*cos(theta)))) / 2.0;

  // double DWDth = (V0value_ * Rvalue_ * r * sin(theta) * (
  //        1.0 / sqrt(Rvalue_*Rvalue_ + r*r - 2.0*Rvalue_*r*cos(theta))
  //        - 1.0 / sqrt(Rvalue_*Rvalue_ + r*r + 2.0*Rvalue_*r*cos(theta))
  //        - 1.0 / sqrt(r0value_*r0value_ + r*r - 2.0*r0value_*r*cos(theta))
  //        + 1.0 / sqrt(r0value_*r0value_ + r*r + 2.0*r0value_*r*cos(theta)))) / 2.0;
  size_t mu, nu;
  for (mu=0; mu<4; ++mu)
    for (nu=mu+1; nu<4; ++nu)
      g[mu][nu]=g[nu][mu]=0;
	
g[0][0] = - 1.0;
g[0][3] = g[3][0]= W;
g[1][1] = 1.0;
g[2][2] = r * r;
g[3][3] = r * r * sin(theta) * sin(theta)-W*W; // sin(theta)^2

}

int BalasinGrumiller::christoffel(double dst[4][4][4], const double pos[8]) const {
//  GYOTO_DEBUG<<endl;
  size_t alpha, mu, nu;
  for (alpha=0; alpha<4; ++alpha)
    for (mu=0; mu<4; ++mu)
      for (nu=0; nu<4; ++nu)
	dst[alpha][mu][nu]=0.; // Initialize to zero

  double r=pos[1], theta=pos[2];
double W = (Rvalue_ - r0value_) * V0value_ + (V0value_ * (-sqrt(Rvalue_*Rvalue_ + r*r - 2.0*Rvalue_*r*cos(theta)) 
         - sqrt(Rvalue_*Rvalue_ + r*r + 2.0*Rvalue_*r*cos(theta)) 
         + sqrt(r0value_*r0value_ + r*r - 2.0*r0value_*r*cos(theta)) 
         + sqrt(r0value_*r0value_ + r*r + 2.0*r0value_*r*cos(theta)))) / 2.0;

  double DWDr = (V0value_ * (-(r - Rvalue_*cos(theta)) / sqrt(Rvalue_*Rvalue_ + r*r - 2.0*Rvalue_*r*cos(theta))
         - (r + Rvalue_*cos(theta)) / sqrt(Rvalue_*Rvalue_ + r*r + 2.0*Rvalue_*r*cos(theta))
         + (r - r0value_*cos(theta)) / sqrt(r0value_*r0value_ + r*r - 2.0*r0value_*r*cos(theta))
         + (r + r0value_*cos(theta)) / sqrt(r0value_*r0value_ + r*r + 2.0*r0value_*r*cos(theta)))) / 2.0;

  double DWDth = (V0value_ * Rvalue_ * r * sin(theta) * (
         1.0 / sqrt(Rvalue_*Rvalue_ + r*r - 2.0*Rvalue_*r*cos(theta))
         - 1.0 / sqrt(Rvalue_*Rvalue_ + r*r + 2.0*Rvalue_*r*cos(theta))
         - 1.0 / sqrt(r0value_*r0value_ + r*r - 2.0*r0value_*r*cos(theta))
         + 1.0 / sqrt(r0value_*r0value_ + r*r + 2.0*r0value_*r*cos(theta)))) / 2.0;


dst[0][0][1] = dst[0][1][0] = (DWDr*W)/(2.0*r*r*sin(theta)*sin(theta));
dst[0][0][2] = dst[0][2][0] = (DWDth*W)/(2.0*r*r*sin(theta)*sin(theta));
dst[0][1][3] = dst[0][3][1] = -0.5*(-2.0*W*r + DWDr*r*r + DWDr*W*W/(sin(theta)*sin(theta)))/(r*r);
dst[0][2][3] = dst[0][3][2] = -0.5*DWDth + W*cos(theta)/sin(theta) - (DWDth*W*W)/(2.0*r*r*sin(theta)*sin(theta));
dst[1][0][3] = dst[1][3][0] = -0.5*DWDr;
dst[1][2][2] = -r;
dst[1][3][3] = DWDr*W - r*sin(theta)*sin(theta);
dst[2][0][3] = dst[2][3][0] = -0.5*DWDth/(r*r);
dst[2][1][2] = dst[2][2][1] = 1.0/r;
dst[2][3][3] = (DWDth*W)/(r*r) - cos(theta)*sin(theta);
dst[3][0][1] = dst[3][1][0] = DWDr/(2.0*r*r*sin(theta)*sin(theta));
dst[3][0][2] = dst[3][2][0] = DWDth/(2.0*r*r*sin(theta)*sin(theta));
dst[3][1][3] = dst[3][3][1] = 1.0/r - (DWDr*W)/(2.0*r*r*sin(theta)*sin(theta));
dst[3][2][3] = dst[3][3][2] = cos(theta)/sin(theta) - (DWDth*W)/(2.0*r*r*sin(theta)*sin(theta));


  return 0;
}

