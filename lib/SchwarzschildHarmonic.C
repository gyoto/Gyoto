/*
    Copyright 2021 Frederic Vincent & Thibaut Paumard

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

#include "GyotoSchwarzschildHarmonic.h"

#include <iostream>
#include <cstdlib>
#include <cmath>

using namespace Gyoto;
using namespace Gyoto::Metric;
using namespace std;

/// Properties
#include "GyotoProperty.h"
GYOTO_PROPERTY_START(SchwarzschildHarmonic, "Schwarzschild in harmonic coordinates")
GYOTO_PROPERTY_END(SchwarzschildHarmonic, Generic::properties)

// accessors

Gyoto::Metric::SchwarzschildHarmonic::SchwarzschildHarmonic()
  : Generic(GYOTO_COORDKIND_SPHERICAL, "SchwarzschildHarmonic")
{
  GYOTO_DEBUG << endl;
}

Gyoto::Metric::SchwarzschildHarmonic::SchwarzschildHarmonic(const SchwarzschildHarmonic & orig)
  : Generic(orig)
{
  GYOTO_DEBUG << endl;
}

// default copy constructor should be fine 
SchwarzschildHarmonic * SchwarzschildHarmonic::clone() const {
  return new SchwarzschildHarmonic(*this); }

Gyoto::Metric::SchwarzschildHarmonic::~SchwarzschildHarmonic()
{
  GYOTO_DEBUG << endl;
}

double SchwarzschildHarmonic::gmunu(const double * pos, int mu, int nu) const {
  double rr = pos[1];
  if (rr<=0.) GYOTO_ERROR("In SchwarzschildHarmonic::gmunu: r<0!");

  double sth2, cth2;
  sincos(pos[2], &sth2, &cth2);
  sth2*=sth2; cth2*=cth2;
  double r2=rr*rr;

  if ((mu==0) && (nu==0)) return -(rr-1.)/(rr+1.);
  if ((mu==1) && (nu==1)) return (rr+1.)/(rr-1.);
  if ((mu==2) && (nu==2)) return (rr+1)*(rr+1.);
  if ((mu==3) && (nu==3)) return (rr+1)*(rr+1.)*sth2;

  return 0.;
} 

int SchwarzschildHarmonic::christoffel(double dst[4][4][4], double const pos[4]) const
{
  int a, mu, nu;
  for (a=0; a<4; ++a)
    for(mu=0; mu<4; ++mu)
      for(nu=0; nu<4; ++nu)
	dst[a][mu][nu]=0.;

  double rr = pos[1], r2=rr*rr;
  double sth, cth;
  sincos(pos[2], &sth, &cth);
  if (rr==0. || sth==0.) GYOTO_ERROR("In SchwarzschildHarmonic::christoffel: "
				    "bad coord");

  dst[0][0][1]=dst[0][1][0]=1./(r2-1.);
  dst[1][0][0]=(rr-1.)/(r2*rr+3.*r2+3.*rr+1.);
  dst[2][1][2]=dst[2][2][1]=1./(rr+1.);
  dst[3][1][3]=dst[3][3][1]=1./(rr+1.);
  dst[1][1][1]=-1./(r2-1.);
  dst[2][3][3]=-cth*sth;
  dst[3][2][3]=dst[3][3][2]=cth/sth;
  dst[1][2][2]=-rr+1.;
  dst[1][3][3]=-(rr-1.)*sth*sth;

  return 0;
}

int SchwarzschildHarmonic::isStopCondition(double const * const coord) const {
  double rsink = 1. + GYOTO_KERR_HORIZON_SECURITY;
  return coord[1] < rsink ;
}

void SchwarzschildHarmonic::circularVelocity(double const * coor, double* vel,
					     double dir) const {
  double sinth = sin(coor[2]);
  double rBL = coor[1]+1.;
  double coord[4] = {coor[0], rBL*sinth, M_PI*0.5, coor[3]};
  
  vel[1] = vel[2] = 0.;
  vel[3] = 1./(dir*pow(coord[1], 1.5));
  
  vel[0] = SysPrimeToTdot(coor, vel+1);
  vel[3] *= vel[0];
}
