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

Minkowski::Minkowski() :
  Generic(GYOTO_COORDKIND_CARTESIAN),
  WIP("Metric::Minkowski")
{kind("Minkowski");}

Minkowski* Minkowski::clone() const { return new Minkowski(*this); };

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
      
} 

double Minkowski::christoffel(const double pos[8], const int alpha, const int mmu, const int nnu) const
{
  if (coordKind()==GYOTO_COORDKIND_CARTESIAN) return 0.;

  if (alpha==0) return 0;
  
  double tmp, tmp2; int mu, nu;
  if (nu<mu) {nu=mmu; mu=nnu; }
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
 
}

#ifdef GYOTO_USE_XERCES
void Minkowski::fillElement(Gyoto::FactoryMessenger *fmp) {
  if (coordKind()==GYOTO_COORDKIND_SPHERICAL)
    fmp -> setParameter ( "Spherical");
  Generic::fillElement(fmp);
}
#endif

void Minkowski::setParameter(string name, string content, string unit){
  if      (name=="Spherical")     coordKind(GYOTO_COORDKIND_SPHERICAL);
  else if (name=="Cartesian")     coordKind(GYOTO_COORDKIND_CARTESIAN);
  else Generic::setParameter(name, content, unit);
}
