/*
    Copyright 2014, 2015, 2019-2020 Thibaut Paumard & Frédéric Vincent

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
GYOTO_PROPERTY_START(Minkowski,
		     "Flat space-time.")
GYOTO_PROPERTY_BOOL(Minkowski, Spherical, Cartesian, spherical,
		    "Whether to use spherical or Cartesian coordinates.")
GYOTO_PROPERTY_END(Minkowski, Generic::properties)

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
    GYOTO_ERROR ("Minkowski::gmunu: incorrect value for mu or nu");

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

  GYOTO_ERROR("BUG: this point should not be reached.");
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
 
  GYOTO_ERROR("BUG: this point should not be reached.");
  return 0.;

}

void Minkowski::observerTetrad(obskind_t obskind,
			       double const pos[4], double fourvel[4],
			       double screen1[4], double screen2[4], 
			       double screen3[4]) const{
  if (coordKind()!=GYOTO_COORDKIND_SPHERICAL){
    GYOTO_ERROR("In Minkowski::observerTetrad: "
	       "coordinates should be spherical-like");
  }
  if (obskind==GYOTO_OBSKIND_KEPLERIAN){
    double gtt = gmunu(pos,0,0),
      grr      = gmunu(pos,1,1),
      gthth    = gmunu(pos,2,2),
      gpp      = gmunu(pos,3,3);
    double omega = 1./(pow(pos[1],1.5));
    double ut2 = -1/(gtt+gpp*omega*omega);
    if (ut2 <= 0. || grr<=0. || gthth <=0.) {
      GYOTO_ERROR("In Minkowski::observerTetrad: "
		 "bad values");
    }
    double ut = sqrt(ut2);
    double fourv[4]={ut,0.,0.,omega*ut};
    double e3[4] = {0.,-1./sqrt(grr),0.,0.};
    double e2[4] = {0.,0.,-1./sqrt(gthth),0.};
    
    double fact1 = gpp*omega/gtt,
      fact2 = gtt*fact1*fact1+gpp;
    if (fact2 <= 0.) GYOTO_ERROR("In Minkowski::observerTetrad: "
				"bad values");
    double a2 = 1./sqrt(fact2), a1 = -a2*fact1;
    double e1[4] = {-a1,0.,0.,-a2};
    
    for (int ii=0;ii<4;ii++){
      fourvel[ii]=fourv[ii];
      screen1[ii]=e1[ii];
      screen2[ii]=e2[ii];
      screen3[ii]=e3[ii];
    }
  }
  Generic::observerTetrad(obskind,pos,fourvel,screen1,screen2,screen3);
}

void Minkowski::spherical(bool t) {
  coordKind(t?GYOTO_COORDKIND_SPHERICAL:GYOTO_COORDKIND_CARTESIAN);
}

bool Minkowski::spherical() const {
  return coordKind() == GYOTO_COORDKIND_SPHERICAL;
}

int Minkowski::diff(const state_t &xi,
		    state_t &dxdt,
		    double mass) const {
  // Check input
  if (xi.size()<8) GYOTO_ERROR("x should have at least 8 elements");
  if (xi.size() != dxdt.size())
    GYOTO_ERROR("x.size() should be the same as dxdt.size()");

  // If not Keplerian or if null geodesic, use generic implementation
  if (!keplerian_ || !mass)
    return Generic::diff(xi, dxdt, mass);

  // We are computing a Keplerian, time-like geodesic.

  // There is no parallel transport for time-like geodesics!
  if (xi.size() > 8) GYOTO_ERROR("No parallel transport for time-like geodesics");

  // x[4:8] is actually dx[0:4]/dt
  dxdt[0]=xi[4];
  dxdt[1]=xi[5];
  dxdt[2]=xi[6];
  dxdt[3]=xi[7];

  // Now the actual equation of motion: d²x/dt²=-ur/r²

  double t, x, y, z, tdot, xdot, ydot, zdot,
    r, theta, phi, rdot, thetadot, phidot,
    tdotdot, xdotdot, ydotdot, zdotdot,
    rdotdot, thetadotdot, phidotdot, r3, tdot3;
  double sth, cth, sph, cph;

  t=xi[0]; tdot=xi[4];

  // Convert to Cartesian

  switch (coordKind()) {
  case GYOTO_COORDKIND_CARTESIAN:
    x=xi[1];
    y=xi[2];
    z=xi[3];
    xdot=xi[5];
    ydot=xi[6];
    zdot=xi[7];
    r3=pow(x*x+y*y+z*z, 1.5);
    break;
  case GYOTO_COORDKIND_SPHERICAL:
    r=xi[1];
    theta=xi[2];
    phi=xi[3];
    rdot=xi[5];
    thetadot=xi[6];
    phidot=xi[7];
    r3=r*r*r;
    sincos(theta, &sth, &cth);
    sincos(phi, &sph, &cph);
    x=r*sth*cph;
    y=r*sth*sph;
    z=r*cth;
    xdot=rdot*sth*cph+r*thetadot*cth*cph-r*phidot*sth*sph;
    ydot=rdot*sth*sph+r*thetadot*cth*sph+r*phidot*sth*cph;
    zdot=rdot*cth-r*thetadot*sth;
    break;
  default:
    GYOTO_ERROR("unimplemented COORDKIND");
  }

  if (r3==0) return 1;

  // Compute second derivatives according to Newton

  // First tdotdot
  /*
    Newton's law is yields:,
     xi'' = -xi / r³ for 1 <= i <= 3 (1)
    We also have
     xi'  = xidot/tdot
     xi'' = d(xidot/tdot)/dt = d(xidot/tdot)/dtau / tdot
     xi'' = (xidotdot*tdot - tdotdot*xidot) / tdot³
    thus we can rewrite Newton's law (1) as
     xidotdot*tdot - tdotdot*xidot = -xi * tdot³ / r³ for i in 1..3 (2)
    We also know that the norm of the quadri velocity is -1:
     -tdot² + sum(xidot²) = -1    (3)
    Derivating (3), it comes:
     tdot*tdotdot = sum(xidot*xidotdot) (4)
    Multiplying (2) by xidot and summing over 1..3:
     tdot*sum(xidot*xidotdot)-tdotdot*sum(xidot²)=-(tdot³/r³)*sum(xi*xidot)
    Thanks to (4), we can replace sum(xidot*xidotdot) by tdot*tdotdot,
    and then factorize tdotdot on the left-hand side:
     tdotdot*(tdot²-sum(xidot²))=-(tdot³/r³)*sum(xi*xidot)
    We recognize (4) and finally:
     tdotdot=-(tdot³/r³)*sum(xi*xidot)
     */
  tdot3=tdot*tdot*tdot;
  tdotdot=-tdot3*(x*xdot+y*ydot+z*zdot)/r3;

  // Then the rest
  /*
    We use (2) again to get xidotdot:
    xidotdot=(-xi*tdot³/r³+tdotdot*xidot)/tdot
  */
  xdotdot=(-x*tdot3/r3+tdotdot*xdot)/tdot;
  ydotdot=(-y*tdot3/r3+tdotdot*ydot)/tdot;
  zdotdot=(-z*tdot3/r3+tdotdot*zdot)/tdot;


  // Convert back to COORDKIND

  dxdt[4]=tdotdot;
  switch (coordKind()) {
  case GYOTO_COORDKIND_CARTESIAN:
    dxdt[5]=xdotdot;
    dxdt[6]=ydotdot;
    dxdt[7]=zdotdot;
    break;
  case GYOTO_COORDKIND_SPHERICAL:
    rdotdot=(xdotdot*x+ydotdot*y+zdotdot*z
	     +xdot*xdot+ydot*ydot+zdot*zdot
	     -rdot*rdot)/r;
    thetadotdot=(
		 (z*rdotdot-zdotdot*r)
		 -((z*rdot-zdot*r)*(2.*rdot*r*r*r-zdot*z*r*r-z*z*rdot*r)
		   /(r*r*r*r-z*z*r*r))
		 )
      *pow(r*r*r*r-z*z*r*r, -0.5);
    phidotdot=((ydotdot*x-y*xdotdot)*(x*x+y*y)
	       -2.*(ydot*x-y*xdot)*(x*xdot+y*ydot))
      /((x*x+y*y)*(x*x+y*y));
    dxdt[5]=rdotdot;
    dxdt[6]=thetadotdot;
    dxdt[7]=phidotdot;
    break;
  default:
    GYOTO_ERROR("unimplemented COORDKIND");
  }
  return 0;
}
