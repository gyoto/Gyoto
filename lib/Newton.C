/*
    Copyright 2020 Thibaut Paumard

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

#include "GyotoUtils.h"
#include "GyotoFactoryMessenger.h"
#include "GyotoNewton.h"
#include "GyotoError.h"
#include "GyotoProperty.h"
#include <cmath>

using namespace std ; 
using namespace Gyoto ; 
using namespace Gyoto::Metric ; 

//// Property list:
//
// Note that none of those lines ends with punctation. "," and ";" are
// added by the macros where needed. Three steps:
//  1- GYOTO_PROPERTY_START(<classname>)
//  2- For each Property we want to support, a line such as:
//       GYOTO_PROPERTY_<type>(<classname>, <propertyname>, <accessorname>)
//     Note that the BOOL type is a bit special: the argument
//     <propertyname> is replaced by two arguments: <name_if_true> and
//     <name_if_false>.
//  3- GYOTO_PROPERTY_END(<classname>, <pointer to parent's property list>)
//
////
GYOTO_PROPERTY_START(Newton,
		     "Classical gravitation.")
GYOTO_PROPERTY_BOOL(Newton, Spherical, Cartesian, spherical,
		    "Whether to use spherical or Cartesian coordinates.")
GYOTO_PROPERTY_END(Newton, Generic::properties)

// This is the minimal constructor: it just sets the coordinate kind and
// the metric kind name.
Newton::Newton() :
  Generic(GYOTO_COORDKIND_CARTESIAN, "Newton")
{}

// The cloner is necessary. If the metric class is not trivial (e.g. contains
// arrays), it may be necessary to implement the copy constructor as well. 
Newton* Newton::clone() const { return new Newton(*this); }

void Newton::spherical(bool t) {
  coordKind(t?GYOTO_COORDKIND_SPHERICAL:GYOTO_COORDKIND_CARTESIAN);
}

bool Newton::spherical() const {
  return coordKind() == GYOTO_COORDKIND_SPHERICAL;
}

void Newton::gmunu(double g[4][4], const double * pos) const
{
  GYOTO_ERROR("Newton metric is degenerate");
  g[0][0]=pos[0];
}

int Newton::christoffel(double dst[4][4][4], const double pos[8]) const {
  GYOTO_ERROR("Newton metric is degenerate");
  dst[0][0][0]=pos[0];
  return 0;
}

/*
Let : Y=[x0,x1,x2,x3,x0_dot,x1_dot,x2_dot,x3_dot] (dot=d/dtau, tau=proper time)
diff is such as : Y_dot=diff(Y)
The general equation of geodesics is used.
 */
int Metric::Generic::diff(const state_t &x,
			  state_t &dxdt) const {
  if (x.size()<8) GYOTO_ERROR("x should have at least 8 elements");
  if (x.size() != dxdt.size()) GYOTO_ERROR("x.size() should be the same as dxdt.size()");
  if (x[4]<1e-6) return 1;
  int nvec = (x.size()-4)/4;
  dxdt[0]=x[4];
  dxdt[1]=x[5];
  dxdt[2]=x[6];
  dxdt[3]=x[7];
  double dst[4][4][4];
  int retval=christoffel(dst, x.data());
  if (retval) return retval;
  for(int alpha=0; alpha<4; ++alpha) {
    for (int v=1; v<=nvec; ++v)
      dxdt[alpha+4*v]=0.;
    for (int i=0;i<4;i++)
      for (int j=0;j<4;j++)
	for (int v=1; v<=nvec; ++v) 
	  dxdt[alpha+v*4] -= dst[alpha][i][j]*x[4+i]*x[v*4+j];
  }
  return 0;
}
