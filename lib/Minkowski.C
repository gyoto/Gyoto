/*
    Copyright 2014 Thibaut Paumard

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
#include "GyotoMinkowski.h"
#include "GyotoError.h"
#include <cmath>

using namespace std ; 
using namespace Gyoto ; 
using namespace Gyoto::Metric ; 

// This is the minimal constructor: it just sets the coordinate kind and
// the metric kind name.
Minkowski::Minkowski() :
  Generic(GYOTO_COORDKIND_CARTESIAN, "Minkowski")
{}

// The cloner is necessary. If the metric class is not trivial (e.g. contains
// arrays), it may be necessary to implement the copy constructor as well. 
Minkowski* Minkowski::clone() const { return new Minkowski(*this); }

void Minkowski::gmunu(double g[4][4], const double * pos) const
{
  GYOTO_DEBUG<<endl;
  size_t mu, nu;
  for (mu=0; mu<4; ++mu)
    for (nu=mu+1; nu<4; ++nu)
      g[mu][nu]=g[nu][mu]=0;

  g[0][0]=-1;
  if (coordKind()==GYOTO_COORDKIND_CARTESIAN) {
    for (mu=1; mu<4; ++mu) g[mu][mu]=1.;
    GYOTO_DEBUG<<"done"<<endl;
    return;
  }

  double r=pos[1], theta=pos[2];
  double tmp=r*sin(theta);

  g[1][1]=1.;
  g[2][2]=r*r;
  g[3][3]=tmp*tmp;
  GYOTO_DEBUG<<"done"<<endl;

}

int Minkowski::christoffel(double dst[4][4][4], const double pos[8]) const {
  GYOTO_DEBUG<<endl;
  size_t alpha, mu, nu;
  for (alpha=0; alpha<4; ++alpha)
    for (mu=0; mu<4; ++mu)
      for (nu=0; nu<4; ++nu)
	dst[alpha][mu][nu]=0.;
  if (coordKind()==GYOTO_COORDKIND_CARTESIAN) return 0;

  double r=pos[1], theta=pos[2], sth, cth;
  sincos(theta, &sth, &cth);

  dst[1][2][2]=-r;                         // Gamma^r_th_th  = -r
  dst[1][3][3]=-r*sth*sth;                 // Gamma^r_ph_ph  = -r*sin²(th)
  dst[2][1][2]=dst[2][2][1]= 1./r;         // Gamma^th_r_th  = 1/r
  dst[2][3][3]=-sth*cth;                   // Gamma^th_ph_ph = -sin(th)*cos(th)
  dst[3][1][3]=dst[3][3][1]= dst[2][1][2]; // Gamma^ph_r_ph  = 1/r
  dst[3][2][3]=dst[3][3][2]= tan(M_PI_2 - pos[2]); // Gamma^ph_th_ph = cotan(th)

  return 0;
}

// It's only necessary to provide one of the two forms for gmunu and
// Christoffel. The preferred, most efficient form is given above. The
// second form is given below, as an example.
double Minkowski::gmunu(const double * pos, int mu, int nu) const {
  if (mu<0 || nu<0 || mu>3 || nu>3)
    throwError ("KerrKS::gmunu: incorrect value for mu or nu");

  if (mu!=nu) return 0.;
  if (mu==0)  return -1.;

  double tmp;

  switch (coordKind()) {
  case GYOTO_COORDKIND_CARTESIAN:
    return 1.;
  case GYOTO_COORDKIND_SPHERICAL:
    switch (mu) {
    case 1:
      return 1.;
    case 2:
      return pos[1]*pos[1];
    case 3:
      tmp=pos[1]*sin(pos[2]);
      return tmp*tmp;
    }
  }

  throwError("BUG: this point should not be reached.");
  return 0.;
  
} 

double Minkowski::christoffel(const double pos[8], const int alpha, const int mmu, const int nnu) const
{
  if (coordKind()==GYOTO_COORDKIND_CARTESIAN) return 0.;

  if (alpha==0) return 0;
  
  double tmp, tmp2; int mu, nu;
  if (nnu<mmu) {nu=mmu; mu=nnu; }
  else { nu=nnu; mu=mmu; }

  switch (alpha) {

  case 1:
    if (mu!=nu) return 0.;
    switch (mu) {
    case 2:
      return -pos[1];                       // Gamma^r_th_th  = -r
    case 3:
      tmp=sin(pos[2]);
      return -pos[1]*tmp*tmp;               // Gamma^r_ph_ph  = -r*sin²(th)
    default:
      return 0.;
    }
  case 2:
    if (mu==1 && nu==2) return 1./pos[1];   // Gamma^th_r_th  = 1/r
    if (mu==3 && nu==3) {
      sincos(pos[2], &tmp, &tmp2);
      return -tmp*tmp2;                     // Gamma^th_ph_ph = -sin(th)*cos(th)
    }
    return 0.;
  case 3:
    if (nu!=3) return 0.;
    if (mu==1) return 1./pos[1];            // Gamma^ph_r_ph  = 1/r
    if (mu==2) return tan(M_PI_2 - pos[2]); // Gamma^ph_th_ph = cotan(th)
    return 0;
  }
 
  throwError("BUG: this point should not be reached.");
  return 0.;

}


// Fillelement is required to be able to export the Metric to an XML file.
#ifdef GYOTO_USE_XERCES
void Minkowski::fillElement(Gyoto::FactoryMessenger *fmp) {
  if (coordKind()==GYOTO_COORDKIND_SPHERICAL)
    fmp -> setParameter ( "Spherical");
  Generic::fillElement(fmp);
}
#endif

// setParameter is the minimal interface to read data from the XML.
void Minkowski::setParameter(string name, string content, string unit){
  if      (name=="Spherical")     coordKind(GYOTO_COORDKIND_SPHERICAL);
  else if (name=="Cartesian")     coordKind(GYOTO_COORDKIND_CARTESIAN);
  else Generic::setParameter(name, content, unit);
}
