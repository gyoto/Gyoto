/*
    Copyright 2011-2015 Frederic Vincent, Thibaut Paumard

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
#include "GyotoKerrKS.h"
#include "GyotoWorldline.h"
#include "GyotoError.h"
#include "GyotoProperty.h"

#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <string>
#include <cstring>

using namespace std ; 
using namespace Gyoto ; 
using namespace Gyoto::Metric ; 

GYOTO_PROPERTY_START(KerrKS,
		     "Metric around a rotating black-hole, in Cartesian Kerr-Schild coordinates.")
GYOTO_PROPERTY_DOUBLE(KerrKS, Spin, spin,
		      "Spin parameter (adimensioned, 0).")
GYOTO_PROPERTY_DOUBLE(KerrKS, HorizonSecurity, horizonSecurity,
		      "Thickness of sink layer around horizon (geometrical units, 0.01).")
GYOTO_PROPERTY_END(KerrKS, Generic::properties)

KerrKS::KerrKS():
  Generic(GYOTO_COORDKIND_CARTESIAN, "KerrKS"),
  spin_(0.), a2_(0.), rsink_(2.+GYOTO_KERR_HORIZON_SECURITY),
  drhor_(GYOTO_KERR_HORIZON_SECURITY)
{}

KerrKS * KerrKS::clone () const { return new KerrKS(*this); }

// Mutators
void KerrKS::spin(const double a) {
  spin_=a;
  a2_=spin_*spin_;
  rsink_=1.+sqrt(1.-a2_)+drhor_;
  tellListeners();
}

void KerrKS::horizonSecurity(const double drhor) {
  drhor_=drhor;
  rsink_=1.+sqrt(1.-a2_)+drhor_;
  tellListeners();
}

// Accessors
double KerrKS::spin() const { return spin_ ; }
double KerrKS::horizonSecurity() const {return drhor_; }

int KerrKS::setParameter(std::string name,
			    std::string content,
			    std::string unit) {
  if (name=="GenericIntegrator") {
    GYOTO_WARNING << "Specifying GenericIntegrator is useless and obsolete\n";
  } else if (name=="SpecificIntegrator") {
    GYOTO_SEVERE << "SpecificIntegrator is not supported anymore\n";
  } else return Metric::Generic::setParameter(name, content, unit);
  return 0;
}

void KerrKS::gmunu(double g[4][4], const double * pos) const {
  double
    x=pos[1], y=pos[2], z=pos[3],
    x2=x*x, y2=y*y, z2=z*z,
    tau=x2+y2+z2-a2_,
    r2=0.5*(tau+sqrt(tau*tau+4*a2_*z2)),
    r=sqrt(r2),
    r3=r2*r, r4=r2*r2, r2_a2=r2+a2_,
    f=2.*r3/(r4+a2_*z2);
  double k[4]=
    {
      1.,
      (r*x+spin_*y)/r2_a2,
      (r*y-spin_*x)/r2_a2,
      z/r
    };
  for (int mu=0; mu<4; ++mu)
    for (int nu=0; nu<=mu;++nu)
      g[mu][nu]=g[nu][mu]=f*k[mu]*k[nu];
  g[0][0] -= 1.;
  for (int mu=1; mu<4;  ++mu) g[mu][mu]+=1.;
}

// The thee macros below are used two or three times each to optimally implement
// gmunu_up, jacobian, and gmunu_up_and_jacobian

// intermediate results
#define _COMMON_TERMS \
  double						\
    x=pos[1], y=pos[2], z=pos[3],				\
    x2=x*x, y2=y*y, z2=z*z, a2z2=a2_*z2,		\
    x2_y2_z2=x2+y2+z2,					\
    tau=x2_y2_z2-a2_,					\
    rho2=tau*tau+4.*a2z2, rho=sqrt(rho2),		\
    r2=0.5*(tau+rho),					\
    r=sqrt(r2), r3=r2*r, r4=r2*r2, r2_a2=r2+a2_,	\
    rx_ay=r*x+spin_*y, ry_ax=r*y-spin_*x,		\
    f=2.*r3/(r4+a2_*z2);//, fr2=f*r2;

// gmunu_up
#define _FILL_GMUNU_UP							\
  {									\
    double frac = f/							\
      (-fr2*(rx_ay*rx_ay + ry_ax*ry_ax)+ r2_a2*r2_a2 *(-r2 +fr2 -f*z2));\
    double kup[4]=							\
    {									\
      -r*r2_a2,								\
      r*rx_ay,								\
      r*ry_ax,								\
      r2_a2*z								\
    };									\
    for (mu=0; mu<4; ++mu) {						\
      for (nu=0; nu<=mu;++nu) {						\
        gup[mu][nu]=gup[nu][mu]=frac*kup[mu]*kup[nu];			\
      }									\
    }									\
    gup[0][0] -= 1.;							\
    for (mu=1; mu<4; ++mu) gup[mu][mu] += 1.;				\
  }

// jacobian
#define _FILL_JACOBIAN					\
  {							\
  double k[4]=						\
    {							\
     1.,						\
     rx_ay/r2_a2,					\
     ry_ax/r2_a2,					\
     z/r						\
    };							\
  double								\
  a4=a2_*a2_,								\
    r4_a2z2=r4+a2z2,							\
    temp=-(2.*r3*(r4-3.*a2z2))/(r4_a2z2*r4_a2z2*rho),			\
    temp2=(a4+2.*r2*x2_y2_z2 - a2_* (x2_y2_z2 - 4.* z2 + rho));		\
  double df[4]=								\
    {									\
     0.,								\
     x*temp,								\
     y*temp,								\
     -((4.*r*z*(2.* a4*a2_ + (a2_ + 2.*r2)*x2_y2_z2*x2_y2_z2 +		\
		a4*(-3.*x2 - 3.*y2 + z2 - 2.*rho) +			\
		a2_*(x2 + y2 - z2)*rho))/(rho*temp2*temp2))		\
    };									\
  double								\
  frac1=1./(r2_a2*r2_a2*rho),						\
    frac2=z/(r2_a2*r*rho),						\
    frac3=-z/(r*rho);							\
  double dk[4][4]=							\
    {									\
      /* d/dt */								\
      {0., 0., 0., 0.},							\
      /* d/dx */							\
      {									\
       0.,								\
       (r3*(x2+rho)-rx_ay*x*(x2+y2+z2+rho)+a2_*(rx_ay*x+r*(x2+rho)))*frac1, \
       (x*(r3*y+a2_*(ry_ax+r*y)-ry_ax*(x2+y2+z2))-(spin_*r2_a2+ry_ax*x)*rho)*frac1, \
       x*frac3								\
      },								\
      /* d/dy */								\
      {									\
       0.,								\
       (a2_*(rx_ay+r*x)*y+r2_a2*spin_*rho-y*(-r3*x+rx_ay*(x2+y2+z2+rho)))*frac1, \
       (r3*(y2+rho)-ry_ax*y*(x2+y2+z2+rho)+a2_*(ry_ax*y+r*(y2+rho)))*frac1, \
       y*frac3								\
      },								\
      /* d/dz */							\
      {									\
       0.,								\
       ((a2_-r2)*x-2*spin_*r*y)*frac2,					\
       ((a2_-r2)*y+2*spin_*r*x)*frac2,					\
       (2.*r2- (z2*(a2_ + x2 + y2 + z2 + rho))/rho)/(2.*r3)		\
      }									\
    };									\
      for(a=0; a<4; ++a)						\
	for (mu=0; mu<4; ++mu)						\
	  for (nu=0; nu<=mu;++nu)					\
	    jac[a][mu][nu]=jac[a][nu][mu]=df[a]*k[mu]*k[nu]+f*dk[a][mu]*k[nu]+f*k[mu]*dk[a][nu]; \
  }


void KerrKS::gmunu_up(double gup[4][4], const double * pos) const {
  size_t mu, nu;
  _COMMON_TERMS;
  double fr2=f*r2;
  _FILL_GMUNU_UP;
}

void KerrKS::jacobian(double jac[4][4][4], const double * pos) const {
  size_t a, mu, nu, i;
  _COMMON_TERMS;
  _FILL_JACOBIAN;
}

void KerrKS::gmunu_up_and_jacobian(double gup[4][4], double jac[4][4][4], const double * pos) const {
  size_t a, mu, nu, i;
  _COMMON_TERMS;
  double fr2=f*r2;
  _FILL_GMUNU_UP;
  _FILL_JACOBIAN;
}

double KerrKS::gmunu(const double * pos, int mu, int nu) const {
  if (mu<0 || nu<0 || mu>3 || nu>3) GYOTO_ERROR ("KerrKS::gmunu: incorrect value for mu or nu");
  //double x=pos[0], y=pos[1], z=pos[2];
  double x=pos[1], y=pos[2], z=pos[3];
  double x2=x*x;
  double y2=y*y;
  double z2=z*z;
  double temp=x2+y2+z2-a2_;
  double rr=sqrt(0.5*(temp+sqrt(temp*temp+4*a2_*z2)));
  double r2=rr*rr;
  double r3=rr*r2;
  double r4=rr*r3;
  double fact=2.*r3/(r4+a2_*z2);

  double res=0.;
  if (mu==nu) {
    if ((mu==0) && (nu==0)) res=fact-1.;
    if ((mu==1) && (nu==1)) res= 1.+fact*pow((rr*x+spin_*y)/(r2+a2_),2);
    if ((mu==2) && (nu==2)) res= 1.+fact*pow((rr*y-spin_*x)/(r2+a2_),2);
    if ((mu==3) && (nu==3)) res= 1.+fact*z2/r2;
  }
  if (nu<mu) {int vu=nu; nu=mu; mu=vu;}
  if (mu==0) {
    if (nu==1) res= fact/(r2+a2_)*(rr*x+spin_*y);
    if (nu==2) res= fact/(r2+a2_)*(rr*y-spin_*x);
    if (nu==3) res= fact*z/rr;
  }
  if (mu==1) {
    if (nu==2) res= fact/pow(r2+a2_,2)*(x*y*(r2-a2_)+spin_*rr*(y2-x2));
    if (nu==3) res= fact/(r2+a2_)*(rr*x+spin_*y)*z/rr;
  }
  if ((mu==2) && (nu==3)) res= fact/(r2+a2_)*(rr*y-spin_*x)*z/rr;

  return res;
  
} 

void KerrKS::circularVelocity(double const coor[4], double vel[4],
			      double dir) const {

  if (keplerian_) {
    // If keplerian_ is true, let the generic implementation return
    // the Keplerian velocity instead of the true circular velocity
    Generic::circularVelocity(coor, vel, dir);
    return;
  }

  double rcross=sqrt ( coor[1]*coor[1] + coor[2]*coor[2] - spin_*spin_);
  double Omega=dir*pow(rcross*rcross*rcross, -0.5);//angular Keplerian velocity
  
  vel[1] = -coor[2]*Omega;
  vel[2] =  coor[1]*Omega;
  vel[3] = 0.;
  vel[0] = SysPrimeToTdot(coor, vel+1);
  vel[1] *= vel[0];
  vel[2] *= vel[0];

}

int KerrKS::isStopCondition(double const * const coord) const {
  double
    x=coord[1], y=coord[2], z=coord[3],
    Tdot=coord[4], xdot=coord[5], ydot=coord[6], zdot=coord[7],
    x2=x*x, y2=y*y, z2=z*z, a2z2=a2_*z2,
    x2_y2_z2=x2+y2+z2,
    tau=x2_y2_z2-a2_,
    rho2=tau*tau+4.*a2z2, rho=sqrt(rho2),
    r2=0.5*(tau+rho),
    r=sqrt(r2);
  double rdot=(x*xdot+y*ydot+z*zdot+a2_*z*zdot/r2)/(r+a2z2/(r*r2));

  //  return (r<rsink_ && rdot >0 && Tdot>0);
  return (r<rsink_);
}
