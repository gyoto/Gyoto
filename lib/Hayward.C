/*
  Copyright 2018-2019 Frederic Lamy, Frederic Vincent, Thibaut Paumard
  
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
#include "GyotoHayward.h"
#include "GyotoWorldline.h"
#include "GyotoError.h"
#include "GyotoProperty.h"

#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <string>
#include <sstream>

using namespace std ;
using namespace Gyoto ;
using namespace Gyoto::Metric ;

GYOTO_PROPERTY_START(Hayward,
                     "Metric around a rotating black-hole, in spherical Boyer-Lindquist coordinates. Cf. Lamy et al. 2018, https://arxiv.org/abs/1802.01635")
GYOTO_PROPERTY_DOUBLE(Hayward, Spin, spin,
                      "Spin parameter (adimensioned, 0).")
GYOTO_PROPERTY_DOUBLE(Hayward, Charge, charge,
                      "Charge parameter (adimensioned, 0).")
GYOTO_PROPERTY_END(Hayward, Generic::properties)


Hayward::Hayward() :
Generic(GYOTO_COORDKIND_SPHERICAL, "Hayward"),
  spin_(0.), a2_(0.), a3_(0.), a4_(0.), charge_(0.), b2_(0.)
{}

// default copy constructor should be fine
Hayward * Hayward::clone () const { return new Hayward(*this); }

// Mutators
void Hayward::spin(const double a) {
  spin_=a;
  a2_=spin_*spin_;
  a3_=a2_*spin_;
  a4_=a2_*a2_;
  tellListeners();
}

void Hayward::charge(const double b) {
  charge_=b;
  b2_=charge_*charge_;
  tellListeners();
}

// Accessors
double Hayward::spin() const { return spin_ ; }
double Hayward::charge() const { return charge_ ; }

double Hayward::getSpecificAngularMomentum(double rr) const {
  // this is l = -u_phi/u_t for a circular equatorial 4-velocity 
  double r2=rr*rr,r3=r2*rr, r5=r3*r2;
  double aa=spin_, a2=aa*spin_, a3=a2*spin_;
  double bb=charge_, b2=bb*charge_;
  double md=r3+2.*b2_;
  double m=r3/md;
  double mdot=-3.*r5/md/md+3.*r2/md;
  double sqrtr=sqrt(rr);
  double sqrt1=sqrt(m-rr*mdot);
  
  return(((a2*rr+r3+2.*a2*m)*sqrt1*sqrtr-(a3+3.*spin_*r2)*m+(a3*rr+spin_*r3)*mdot)/(a2*rr*mdot+r3+2*sqrt1*spin_*m*sqrtr-(a2+2.*r2)*m));
}

double Hayward::getPotential(double const pos[4], double l_cst) const {
  // this is W = -ln(|u_t|) for a circular equatorial 4-velocity
  double  gtt = gmunu(pos,0,0);
  double  gtp = gmunu(pos,0,3);
  double  gpp = gmunu(pos,3,3);
  double  Omega = -(gtp + l_cst * gtt)/(gpp + l_cst * gtp) ;
  
  double  W = 0.5 * log(abs(gtt + 2. * Omega * gtp + Omega*Omega * gpp))
    - log(abs(gtt + Omega * gtp)) ;
  return  W ;
}

void Hayward::circularVelocity(double const coor[4], double vel[4],
			       double dir) const {
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG<<"coor=["<<coor[0]<<", "<<coor[1]<<", "<<coor[2]<<", "<<coor[3]
	     <<"], dir="<<dir<<endl;
# endif
  
  // project on equatorial plane
  double sinth = sin(coor[2]);
  double coord[4] = {coor[0], coor[1]*sinth, M_PI*0.5, coor[3]};
  
  vel[1] = vel[2] = 0.;
  
  // Omega Keplerian from Eq. 30 of Lamy+18
  double b4 = b2_*b2_;
  double rr = coord[1], r2=rr*rr, r3=r2*rr, r4=r2*r2,
    r5=r4*rr, r6=r5*rr, r7=r6*rr;
  double sqrtarg = -(4.*b2_*r2-r5) / (4.*b4+4.*b2_*r3+r6);
  double nume = 4.*spin_*b2_*rr - spin_*r4
    + dir* (4.*b4+4.*b2_*r3+r6)*sqrt(sqrtarg),
    deno = r7 - (a2_-4.*b2_)*r4+4.*(a2_*b2_+b4)*rr;
  
  vel[3] = nume/deno;
  
  vel[0] = SysPrimeToTdot(coor, vel+1);
  vel[3] *= vel[0];
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG_ARRAY(vel,4);
# endif
}

//Computation of metric coefficients in covariant form

void Hayward::gmunu(double g[4][4], const double * pos) const
{
  double r = pos[1];
  double sth2, cth2;
  sincos(pos[2], &sth2, &cth2);
  sth2*=sth2; cth2*=cth2;
  
  int mu, nu;
  for (mu=0; mu<4; ++mu)
    for (nu=0; nu<4; ++nu)
      g[mu][nu]=0.;
  
  
  if (r>=1.)
    {
      double u=1./r, u2=u*u, u3=u2*u, u4=u3*u, u5=u4*u, u6=u5*u, u7=u6*u;
      double a2b2=a2_*b2_;
      
      g[0][0] =-(2.*a2b2*u5*cth2+a2_*u2*cth2+2.*b2_*u3-2.*u+1.)/(a2_*u2*cth2+1.)/(2.*b2_*u3+1.);
      g[1][1] =(a2_*u2*cth2+1.)*(2.*b2_*u3+1.)/(2.*a2b2*u5+2.*b2_*u3+a2_*u2-2.*u+1.);
      g[2][2] =(a2_*u2*cth2+1.)/(u2); //pb de 1/u2
      g[3][3] =(2.*a4_*b2_*u7*cth2+2.*a2b2*u5*cth2+a4_*u4*cth2+2.*a2b2*u5+2.*a2_*u3*sth2+a2_*u2*cth2+2.*b2_*u3+a2_*u2+1.)*sth2/(a2_*u2*cth2+1.)/(2.*b2_*u3+1.)/(u2);
      g[0][3] = g[3][0] =-2.*spin_*u*sth2/(a2_*u2*cth2+1.)/(2.*b2_*u3+1.);
    }
  
  if (r>=0. && r<1.)
    {
      double r2=r*r, r3=r2*r;
      double m=r3/(r3+2.*b2_);
      double sigma=r2+a2_*cth2;
      double delta=r2-2.*m*r+a2_;
      
      g[0][0] = -1.+2.*m*r/sigma;
      g[1][1] = sigma/delta;
      g[2][2] = sigma;
      g[3][3] = (r2+a2_+2.*m*r*a2_*sth2/sigma)*sth2;
      g[0][3] = g[3][0] = -2.*spin_*m*r*sth2/sigma;
    }
  
  if (r<0.)
    {
      double r2=r*r, r3=r2*r;
      double m=-r3/(-r3+2.*b2_);
      double sigma=r2+a2_*cth2;
      double delta=r2-2.*m*r+a2_;
      
      g[0][0] = -1.+2.*m*r/sigma;
      g[1][1] = sigma/delta;
      g[2][2] = sigma;
      g[3][3] = (r2+a2_+2.*m*r*a2_*sth2/sigma)*sth2;
      g[0][3] = g[3][0] = -2.*spin_*m*r*sth2/sigma;
    }
}




double Hayward::gmunu(const double * pos, int mu, int nu) const
{
  double r = pos[1];
  double sth2, cth2;
  sincos(pos[2], &sth2, &cth2);
  sth2*=sth2; cth2*=cth2;
  
  if (r>=1.)
    {
      double u=1./r, u2=u*u, u3=u2*u, u4=u3*u, u5=u4*u, u6=u5*u, u7=u6*u;
      double a2b2=a2_*b2_;
      
      
      if ((mu==0) && (nu==0)) return -(2.*a2b2*u5*cth2+a2_*u2*cth2+2.*b2_*u3-2.*u+1.)/(a2_*u2*cth2+1.)/(2.*b2_*u3+1.);
      if ((mu==1) && (nu==1)) return (a2_*u2*cth2+1.)*(2.*b2_*u3+1.)/(2.*a2b2*u5+2.*b2_*u3+a2_*u2-2.*u+1.);
      if ((mu==2) && (nu==2)) return (a2_*u2*cth2+1.)/(u2);
      if ((mu==3) && (nu==3)) return (2.*a4_*b2_*u7*cth2+2.*a2b2*u5*cth2+a4_*u4*cth2+2.*a2b2*u5+2.*a2_*u3*sth2+a2_*u2*cth2+2.*b2_*u3+a2_*u2+1.)*sth2/(a2_*u2*cth2+1.)/(2.*b2_*u3+1.)/(u2);
      if (((mu==0) && (nu==3)) || ((mu==3) && (nu==0))) return -2.*spin_*u*sth2/(a2_*u2*cth2+1.)/(2.*b2_*u3+1.);
    }
  
  if (r>=0 && r<1.)
    {
      double r2=r*r, r3=r2*r, r4=r2*r2;
      double d=r3+2.*b2_;
      double m=r3/d;
      double sigma=r2+a2_*cth2;
      double delta=r2-2.*m*r+a2_;
      
      if ((mu==0) && (nu==0)) return -1.+2.*m*r/sigma;
      if ((mu==1) && (nu==1)) return sigma/delta;
      if ((mu==2) && (nu==2)) return sigma;
      if ((mu==3) && (nu==3)) return (r2+a2_+2.*m*r*a2_*sth2/sigma)*sth2;
      if (((mu==0) && (nu==3)) || ((mu==3) && (nu==0))) return  -2.*spin_*m*r*sth2/sigma;
    }
  
  if (r<0.)
    {
      double r2=r*r, r3=r2*r;
      double m=-r3/(-r3+2.*b2_);
      double sigma=r2+a2_*cth2;
      double delta=r2-2.*m*r+a2_;
      
      if ((mu==0) && (nu==0)) return -1.+2.*m*r/sigma;
      if ((mu==1) && (nu==1)) return sigma/delta;
      if ((mu==2) && (nu==2)) return sigma;
      if ((mu==3) && (nu==3)) return (r2+a2_+2.*m*r*a2_*sth2/sigma)*sth2;
      if (((mu==0) && (nu==3)) || ((mu==3) && (nu==0))) return  -2.*spin_*m*r*sth2/sigma;
    }
  
  
  return 0.;
  
}


//Computation of metric coefficients in contravariant form

void Hayward::gmunu_up(double gup[4][4], const double * pos) const
{
  double r = pos[1];
  double sth2, cth2;
  sincos(pos[2], &sth2, &cth2);
  sth2*=sth2; cth2*=cth2;
  double a2b2=a2_*b2_, a4b2=a2_*a2b2;
  
  int mu, nu;
  for (mu=0; mu<4; ++mu) for (nu=0; nu<4; ++nu) gup[mu][nu]=0.;
  
  if (r>=1.)
    {
      double u=1./r, u2=u*u, u3=u2*u, u4=u3*u, u5=u4*u, u6=u5*u, u7=u6*u;
      
      gup[0][0]=-(2.*a4_*b2_*u7*cth2+2.*a2b2*u5*cth2+a4_*u4*cth2+2.*a2b2*u5+2.*a2_*u3*sth2+a2_*u2*cth2+2.*b2_*u3+a2_*u2+1.)/(2.*a2b2*u5+2.*b2_*u3+a2_*u2-2.*u+1.)/(a2_*u2*cth2+1.);
      gup[1][1]=(2.*a2b2*u5+2.*b2_*u3+a2_*u2-2.*u+1.)/(a2_*u2*cth2+1.)/(2.*b2_*u3+1.);
      gup[2][2]=u2/(a2_*u2*cth2+1.);
      gup[3][3]=(2.*a2b2*u5*cth2+2.*b2_*u3+a2_*u2*cth2-2.*u+1.)*u2/(2.*a4_*b2_*u7*cth2+2.*a2b2*u5*cth2+a4_*u4*cth2+2.*a2b2*u5+2.*a2_*u3*sth2+a2_*u2*cth2-2.*a2_*u3+2.*b2_*u3+a2_*u2-2.*u+1.)/(sth2);
      gup[0][3] = gup[3][0] = -2.*spin_*u3/(2.*a2b2*u5+2.*b2_*u3+a2_*u2-2.*u+1.)/(a2_*u2*cth2+1.);
    }
  
  if (r>=0. && r<1.)
    {
      double r2=r*r, r3=r2*r, r4=r2*r2, r5=r4*r, r6=r5*r, r7=r6*r;
      double d=r3+2.*b2_;
      double m=r3/d;
      double sigma=r2+a2_*cth2;
      double delta=r2-2.*m*r+a2_;
      
      gup[0][0] = -(a4_*r3*cth2+a2_*r5*cth2+2.*a4b2*cth2+2.*a2b2*r2*cth2+2.*a2_*r4*sth2+a2_*r5+r7+2.*a2b2*r2+2.*b2_*r4)/sigma/(a2_*r3+r5+2.*a2b2+2.*b2_*r2-2.*r4);
      gup[1][1] = (a2_*r3+r5+2.*a2b2+2.*b2_*r2-2.*r4)/sigma/d;
      gup[2][2] = 1./sigma;
      gup[3][3] = (a2_*r3*cth2+2.*a2b2*cth2+r5+2.*b2_*r2-2.*r4)/(a4_*r3*cth2+a2_*r5*cth2+2.*a4b2*cth2+2.*a2b2*r2*cth2+2.*a2_*r4*sth2+a2_*r5+r7+2.*a2b2*r2-2.*a2_*r4+2.*b2_*r4-2.*r6)/sth2;
      gup[0][3] = gup[3][0] = -2.*spin_*r4/sigma/(a2_*r3+r5+2.*a2b2+2.*b2_*r2-2.*r4);
    }
  
  if (r<0.)
    {
      double r2=r*r, r3=r2*r, r4=r2*r2, r5=r4*r, r6=r5*r, r7=r6*r;
      double m=-r3/(-r3+2.*b2_);
      double sigma=r2+a2_*cth2;
      double delta=r2-2.*m*r+a2_;
      double dm=r3-2.*b2_;
      
      gup[0][0] = -(a4_*r3*cth2+a2_*r5*cth2-2.*a4b2*cth2-2.*a2b2*r2*cth2+2.*a2_*r4*sth2+a2_*r5+r7-2.*a2b2*r2-2.*b2_*r4)/sigma/(a2_*r3+r5-2.*a2b2-2.*b2_*r2-2.*r4);
      gup[1][1] = (a2_*r3+r5-2.*a2b2-2.*b2_*r2-2.*r4)/sigma/dm;
      gup[2][2] = 1./sigma;
      gup[3][3] = (a2_*r3*cth2-2.*a2b2*cth2+r5-2.*b2_*r2-2.*r4)/(a4_*r3*cth2+a2_*r5*cth2-2.*a4b2*cth2-2.*a2b2*r2*cth2+2.*a2_*r4*sth2+a2_*r5+r7-2.*a2b2*r2-2.*a2_*r4-2.*b2_*r4-2.*r6)/sth2;
      gup[0][3] = gup[3][0] = -2.*spin_*r4/sigma/(a2_*r3+r5-2.*a2b2-2.*b2_*r2-2.*r4);
    }
  
  
  //    cout << "g^00=" << gup[0][0] << endl ;
  //    cout << "g^11=" << gup[1][1] << endl ;
  //    cout << "g^22=" << gup[2][2] << endl ;
  //    cout << "g^33=" << gup[3][3] << endl ;
  //    cout << "g^03=" << gup[0][3] << endl ;
  
}


double Hayward::gmunu_up(const double * pos, int mu, int nu) const
{
  double r = pos[1];
  double sth2, cth2;
  sincos(pos[2], &sth2, &cth2);
  sth2*=sth2; cth2*=cth2;
  double a2b2=a2_*b2_, a4b2=a2_*a2b2;
  
  if (r>=1.)
    {
      double u=1./r, u2=u*u, u3=u2*u, u4=u3*u, u5=u4*u, u6=u5*u, u7=u6*u;
      
      if ((mu==0) && (nu==0)) return -(2.*a4_*b2_*u7*cth2+2.*a2b2*u5*cth2+a4_*u4*cth2+2.*a2b2*u5+2.*a2_*u3*sth2+a2_*u2*cth2+2.*b2_*u3+a2_*u2+1.)/(2.*a2b2*u5+2.*b2_*u3+a2_*u2-2.*u+1.)/(a2_*u2*cth2+1.);
      if ((mu==1) && (nu==1)) return (2.*a2b2*u5+2.*b2_*u3+a2_*u2-2.*u+1.)/(a2_*u2*cth2+1.)/(2.*b2_*u3+1.);
      if ((mu==2) && (nu==2)) return u2/(a2_*u2*cth2+1.);
      if ((mu==3) && (nu==3)) return (2.*a2b2*u5*cth2+2.*b2_*u3+a2_*u2*cth2-2.*u+1.)*u2/(2.*a4_*b2_*u7*cth2+2.*a2b2*u5*cth2+a4_*u4*cth2+2.*a2b2*u5+2.*a2_*u3*sth2+a2_*u2*cth2-2.*a2_*u3+2.*b2_*u3+a2_*u2-2.*u+1.)/(sth2);
      if (((mu==0) && (nu==3)) || ((mu==3) && (nu==0)))return -2.*spin_*u3/(2.*a2b2*u5+2.*b2_*u3+a2_*u2-2.*u+1.)/(a2_*u2*cth2+1.);
    }
  
  if (r>=0. && r<1.)
    {
      double r2=r*r, r3=r2*r, r4=r2*r2, r5=r4*r, r6=r5*r, r7=r6*r;
      double d=r3+2.*b2_;
      double m=r3/d;
      double sigma=r2+a2_*cth2;
      double delta=r2-2.*m*r+a2_;
      
      if ((mu==0) && (nu==0)) return -(a4_*r3*cth2+a2_*r5*cth2+2.*a4b2*cth2+2.*a2b2*r2*cth2+2.*a2_*r4*sth2+a2_*r5+r7+2.*a2b2*r2+2.*b2_*r4)/sigma/(a2_*r3+r5+2.*a2b2+2.*b2_*r2-2.*r4);
      if ((mu==1) && (nu==1)) return (a2_*r3+r5+2.*a2b2+2.*b2_*r2-2.*r4)/sigma/d;
      if ((mu==2) && (nu==2)) return 1./sigma;
      if ((mu==3) && (nu==3)) return (a2_*r3*cth2+2.*a2b2*cth2+r5+2.*b2_*r2-2.*r4)/(a4_*r3*cth2+a2_*r5*cth2+2.*a4b2*cth2+2.*a2b2*r2*cth2+2.*a2_*r4*sth2+a2_*r5+r7+2.*a2b2*r2-2.*a2_*r4+2.*b2_*r4-2.*r6)/sth2;
      if (((mu==0) && (nu==3)) || ((mu==3) && (nu==0)))return -2*spin_*r4/sigma/(a2_*r3+r5+2.*a2b2+2.*b2_*r2-2.*r4);
    }
  
  if (r<0.)
    {
      double r2=r*r, r3=r2*r, r4=r2*r2, r5=r4*r, r6=r5*r, r7=r6*r;
      double m=-r3/(-r3+2.*b2_);
      double sigma=r2+a2_*cth2;
      double delta=r2-2.*m*r+a2_;
      double dm=r3-2.*b2_;
      
      if ((mu==0) && (nu==0)) return -(a4_*r3*cth2+a2_*r5*cth2-2.*a4b2*cth2-2.*a2b2*r2*cth2+2.*a2_*r4*sth2+a2_*r5+r7-2.*a2b2*r2-2.*b2_*r4)/sigma/(a2_*r3+r5-2.*a2b2-2.*b2_*r2-2.*r4);
      if ((mu==1) && (nu==1)) return (a2_*r3+r5-2.*a2b2-2.*b2_*r2-2.*r4)/sigma/dm;
      if ((mu==2) && (nu==2)) return 1./sigma;
      if ((mu==3) && (nu==3)) return (a2_*r3*cth2-2.*a2b2*cth2+r5-2.*b2_*r2-2.*r4)/(a4_*r3*cth2+a2_*r5*cth2-2.*a4b2*cth2-2.*a2b2*r2*cth2+2.*a2_*r4*sth2+a2_*r5+r7-2.*a2b2*r2-2.*a2_*r4-2.*b2_*r4-2.*r6)/sth2;
      if (((mu==0) && (nu==3)) || ((mu==3) && (nu==0)))return  -2.*spin_*r4/sigma/(a2_*r3+r5+2.*a2b2+2.*b2_*r2-2.*r4);
    }
  
  return 0.;
}



int Hayward::christoffel(double dst[4][4][4], double const pos[4]) const
{
  int a, mu, nu;
  for (a=0; a<4; ++a)
    for(mu=0; mu<4; ++mu)
      for(nu=0; nu<4; ++nu)
	dst[a][mu][nu]=0.;
  
  double r = pos[1];
  double sth, cth;
  sincos(pos[2], &sth, &cth);
  double sth2 = sth*sth,
    cth2 = cth*cth,
    cth4=cth2*cth2,
    sth4=sth2*sth2,
    a2cthsth=a2_*cth*sth,
    a2b2=a2_*b2_,
    a4b2=a4_*b2_,
    b4_=b2_*b2_,
    a4b4=a4_*b4_,
    a2b4=a2_*b4_;
  double a3=a2_*spin_,
    a6=a4_*a2_;
  
  if (r>=1.)
    {
      
      double u=1./r, u2=u*u, u3=u2*u, u4=u3*u, u5=u4*u, u6=u5*u, u7=u6*u, u8=u7*u, u9=u8*u, u10=u9*u;
      
      dst[1][1][1]=-(4.*a2b2*b2_*u7*cth2-4.*a2b2*b2_*u7-8.*a2b2*u5*cth2+4.*a2b2*u4*cth2-4.*a2b2*u4-a2_*u2*cth2-4.*b2_*u3+a2_*u*cth2-a2_*u+1.)*u2/(2.*a2b2*u5+2.*b2_*u3+a2_*u2-2.*u+1.)/(a2_*u2*cth2+1.)/(2.*b2_*u3+1.);
      dst[1][2][1]=dst[1][1][2]=-a2cthsth*u2/(a2_*u2*cth2+1.);
      dst[1][2][2]=-(2.*a2b2*u5+2.*b2_*u3+a2_*u2-2.*u+1)/(a2_*u2*cth2+1.)/(2.*b2_*u3+1.)/u;
      dst[1][3][3]=-(4.*a4b4*u10*cth4+8.*a4b2*u8*cth2*sth2+4.*a4b2*u7*cth4+8.*a2b4*u8*cth2+a4_*u5*cth2*sth2+a4_*u4*cth4+4.*a2b2*u6*sth2+8.*a2b2*u5*cth2+4.*b4_*u6-a2_*u3*sth2+2.*a2_*u2*cth2+4.*b2_*u3+1.)*(2.*a2b2*u5+2.*b2_*u3+a2_*u2-2.*u+1)*sth2/(a2_*u2*cth2+1.)/(a2_*u2*cth2+1.)/(a2_*u2*cth2+1.)/(2.*b2_*u3+1.)/(2.*b2_*u3+1.)/(2.*b2_*u3+1.)/u;
      dst[1][3][0]=dst[1][0][3]=(8.*a2b2*u5*cth2+a2_*u2*cth2+4.*b2_*u3-1.)*(2.*a2b2*u5+2.*b2_*u3+a2_*u2-2.*u+1.)*spin_*u2*sth2/(a2_*u2*cth2+1.)/(a2_*u2*cth2+1.)/(a2_*u2*cth2+1.)/(2.*b2_*u3+1.)/(2.*b2_*u3+1.)/(2.*b2_*u3+1.);
      dst[1][0][0]=(8.*a2b2*u5*sth2-8.*a2b2*u5+a2_*u2*sth2-4.*b2_*u3-a2_*u2+1.)*(2.*a2b2*u5+2.*b2_*u3+a2_*u2-2.*u+1.)*u2/(a2_*u2*cth2+1.)/(a2_*u2*cth2+1.)/(a2_*u2*cth2+1.)/(2.*b2_*u3+1.)/(2.*b2_*u3+1.)/(2.*b2_*u3+1.);
      dst[2][1][1]=(2.*b2_*u3+1.)*a2cthsth*u4/(2.*a2b2*u5+2.*b2_*u3+a2_*u2-2.*u+1.)/(a2_*u2*cth2+1.);
      dst[2][2][1]=dst[2][1][2]=u/(a2_*u2*cth2+1.);
      dst[2][2][2]=-a2cthsth*u2/(a2_*u2*cth2+1.);
      dst[2][3][3]=-(2.*a4_*a2b2*u9*cth4+2.*a4b2*u7*cth4+a6*u6*cth4+4.*a4b2*u7*cth2-2.*a4_*u5*cth4+a4_*u4*cth4+4.*a2b2*u5*cth2+2.*a4_*u4*cth2+2.*a4_*u5+2.*a2b2*u5-4.*a2_*u3*cth2+2.*a2_*u2*cth2+4.*a2_*u3+2.*b2_*u3+a2_*u2+1.)*cth*sth/(a2_*u2*cth2+1.)/(a2_*u2*cth2+1.)/(a2_*u2*cth2+1.)/(2.*b2_*u3+1.);
      dst[2][0][3]=dst[2][3][0]=2.*(a2_*u2+1.)*spin_*u3*cth*sth/(a2_*u2*cth2+1.)/(a2_*u2*cth2+1.)/(a2_*u2*cth2+1.)/(2.*b2_*u3+1.);
      dst[2][0][0]=-2.*a2_*u5*cth*sth/(a2_*u2*cth2+1.)/(a2_*u2*cth2+1.)/(a2_*u2*cth2+1.)/(2.*b2_*u3+1.);
      dst[3][3][1]=dst[3][1][3]=(4.*a4b4*u10*cth4+8.*a4b2*u8*cth2*sth2+4.*a4b2*u7*cth4+8.*a2b4*u8*cth2+a4_*u5*cth2*sth2-4.*a2b2*u6*cth2+a4_*u4*cth4+4.*a2b2*u6*sth2+8.*a2b2*u5*cth2+4.*b4_*u6-2.*a2_*u3*cth2-a2_*u3*sth2-4.*b2_*u4+2.*a2_*u2*cth2+4.*b2_*u3-2.*u+1.)*u/(2.*a2b2*u5+2.*b2_*u3+a2_*u2-2.*u+1.)/(a2_*u2*cth2+1.)/(a2_*u2*cth2+1.)/(2.*b2_*u3+1.);
      dst[3][3][2]=dst[3][2][3]=(2.*a4b2*u7*cth4+a4_*u4*cth4+4.*a2b2*u5*cth2-2.*a2_*u3*cth2+2.*a2_*u2*cth2+2.*a2_*u3+2.*b2_*u3+1.)*cth/(a2_*u2*cth2+1.)/(a2_*u2*cth2+1.)/(2.*b2_*u3+1.)/sth;
      dst[3][0][1]=dst[3][1][0]=-(8.*a2b2*u5*cth2+a2_*u2*cth2+4.*b2_*u3-1.)*spin_*u4/(2.*a2b2*u5+2.*b2_*u3+a2_*u2-2.*u+1.)/(a2_*u2*cth2+1.)/(a2_*u2*cth2+1.)/(2.*b2_*u3+1.);
      dst[3][0][2]=dst[3][2][0]=-2.*spin_*u3*cth/(a2_*u2*cth2+1.)/(a2_*u2*cth2+1.)/(2.*b2_*u3+1.)/sth;
      dst[0][3][1]=dst[0][1][3]=(8.*a4b2*u7*cth2+4.*a2b2*u5*cth2+a4_*u4*cth2+4.*a2b2*u5-a2_*u2*cth2-a2_*u2-3.)*spin_*u2*sth2/(2.*a2b2*u5+2.*b2_*u3+a2_*u2-2.*u+1.)/(a2_*u2*cth2+1.)/(a2_*u2*cth2+1.)/(2.*b2_*u3+1.);
      dst[0][3][2]=dst[0][2][3]=-2.*(a2_*u2*sth2-a2_*u2-1.)*a3_*u3*cth*sth*sth2/(a2_*u2*cth2+1.)/(a2_*u2*cth2+1.)/(a2_*u2*cth2+1.)/(2.*b2_*u3+1.);
      dst[0][0][1]=dst[0][1][0]=(8.*a2b2*u5*sth2-8.*a2b2*u5+a2_*u2*sth2-4.*b2_*u3-a2_*u2+1.)*(a2_*u2+1.)*u2/(2.*a2b2*u5+2.*b2_*u3+a2_*u2-2.*u+1.)/(a2_*u2*cth2+1.)/(a2_*u2*cth2+1.)/(2.*b2_*u3+1.);
      dst[0][0][2]=dst[0][2][0]=-2.*a2_*u3*cth*sth/(a2_*u2*cth2+1.)/(a2_*u2*cth2+1.)/(2.*b2_*u3+1.);
      
      
      //cout << "t=" << pos[0] << " r=" << pos[1] << " theta=" << pos[2] << " dt=" << pos[4] << " dr=" << pos[5] << " dth=" << pos[6] << " dph=" << pos[7] << endl ;
    }
  
  if (r>=0. && r<1.)
    {
      double r2=r*r, r3=r2*r, r4=r3*r, r5=r4*r, r6=r5*r, r7=r6*r, r8=r7*r, r9=r8*r, r10=r9*r;
      double d=r3+2.*b2_;
      double Sigma=r2+a2_*cth2, Sigma2=Sigma*Sigma, Sigma3=Sigma2*Sigma;
      
      dst[1][1][1]=(a2_*r6*sth2+4.*a2b2*r3*sth2-a2_*r5*sth2+4.*a2b2*b2_*sth2-8.*a2b2*r2*sth2+a2_*r5-r7+8.*a2b2*r2+4.*b2_*r4)*r/(a2_*r3+r5+2.*a2b2+2.*b2_*r2-2.*r4)/Sigma/d;
      dst[1][2][1]=dst[1][1][2]=-a2cthsth/Sigma;
      dst[1][2][2]=-(a2_*r3+r5+2.*a2b2+2.*b2_*r2-2.*r4)*r/Sigma/d;
      dst[1][3][3]=-(a4_*r6*cth4+4.*a4_*b2_*r3*cth4+a4_*r5*cth2*sth2+2.*a2_*r8*cth2+4.*a4b4*cth4+8.*a4b2*r2*cth2*sth2+8.*a2b2*r5*cth2-a2_*r7*sth2+r10+8.*a2b4*r2*cth2+4.*a2b2*r4*sth2+4.*b2_*r7+4.*b4_*r4)*(a2_*r3+r5+2.*a2b2+2.*b2_*r2-2.*r4)*r*sth2/Sigma3/d/d/d;
      dst[1][3][0]=dst[1][0][3]=(a2_*r3*cth2+8.*a2b2*cth2-r5+4.*b2_*r2)*(a2_*r3+r5+2.*a2b2+2.*b2_*r2-2.*r4)*spin_*r3*sth2/Sigma3/d/d/d;
      dst[1][0][0]=-(a2_*r3*cth2+8.*a2b2*cth2-r5+4.*b2_*r2)*(a2_*r3+r5+2.*a2b2+2.*b2_*r2-2.*r4)*r3/Sigma3/d/d/d;
      dst[2][1][1]=d*a2cthsth/(a2_*r3+r5+2.*a2b2+2.*b2_*r2-2.*r4)/Sigma;
      dst[2][2][1]=dst[2][1][2]=r/Sigma;
      dst[2][2][2]=-a2cthsth/Sigma;
      dst[2][3][3]=-(a4_*a2_*r3*cth4+a4_*r5*cth4+2.*a4_*a2b2*cth4+2.*a4b2*r2*cth4-2*a4_*r4*cth4+2.*a4_*r5*cth2+2.*a2_*r7*cth2+4.*a2b2*a2_*r2*cth2+4.*a2b2*r4*cth2-4.*a2_*r6*cth2+a2_*r7+r9+2.*a4_*r4+2.*a2b2*r4+4.*a2_*r6+2.*b2_*r6)*sth*cth/Sigma3/d;
      dst[2][0][3]=dst[2][3][0]=2.*(a2_+r2)*spin_*r4*cth*sth/Sigma3/d;
      dst[2][0][0]=-2.*a2_*r4*cth*sth/Sigma3/d;
      dst[3][3][1]=dst[3][1][3]=(a4_*r6*cth4+4.*a4_*b2_*r3*cth4+a4_*r5*cth2*sth2+2.*a2_*r8*cth2+4.*a4b4*cth4+8.*a4_*b2_*r2*cth2*sth2+8.*a2b2*r5*cth2-2.*a2_*r7*cth2-a2_*r7*sth2+r10+8.*a2b4*r2*cth2-4.*a2b2*r4*cth2+4.*a2b2*r4*sth2+4.*b2_*r7-2.*r9+4.*b4_*r4-4.*b2_*r6)*r/Sigma2/(a2_*r3+r5+2.*a2b2+2.*b2_*r2-2.*r4)/d;
      dst[3][3][2]=dst[3][2][3]=(a4_*r3*cth4+2.*a4_*b2_*cth4+2.*a2_*r5*cth2+4.*a2b2*r2*cth2-2.*a2_*r4*cth2+r7+2.*a2_*r4+2.*b2_*r4)*cth/Sigma2/d/sth;
      dst[3][0][1]=dst[3][1][0]=-(a2_*r3*cth2+8.*a2b2*cth2-r5+4.*b2_*r2)*spin_*r3/Sigma2/(a2_*r3+r5+2.*a2b2+2.*b2_*r2-2.*r4)/d;
      dst[3][0][2]=dst[3][2][0]=-2.*spin_*r4*cth/Sigma2/d/sth;
      dst[0][3][1]=dst[0][1][3]=(a4_*r3*cth2-a2_*r5*cth2+8.*a4b2*cth2+4.*a2b2*r2*cth2-a2_*r5-3.*r7+4.*a2b2*r2)*spin_*r3*sth2/Sigma2/(a2_*r3+r5+2.*a2b2+2.*b2_*r2-2.*r4)/d;
      dst[0][3][2]=dst[0][2][3]=-2.*(a2_*sth2-a2_-r2)*a3_*r4*cth*sth2*sth/Sigma3/d;
      dst[0][0][1]=dst[0][1][0]=(a2_*r3*sth2+8.*a2b2*sth2-a2_*r3+r5-8.*a2b2-4.*b2_*r2)*(a2_+r2)*r3/(a2_*r3+r5+2.*a2b2+2.*b2_*r2-2.*r4)/Sigma2/d;
      dst[0][0][2]=dst[0][2][0]=-2.*a2_*r4*cth*sth/Sigma2/d;
      
      // cout.precision(16);
      // cout << "t=" << pos[0] << " r=" << pos[1] << " theta=" << pos[2] << " dt=" << pos[4] << " dr=" << pos[5] << " dth=" << pos[6] << " dph=" << pos[7] << endl ;
      
    }
  
  if (r<0.)
    {
      double r2=r*r, r3=r2*r, r4=r3*r, r5=r4*r, r6=r5*r, r7=r6*r, r8=r7*r, r9=r8*r, r10=r9*r;
      double sth3=sth2*sth;
      double Sigma=r2+a2_*cth2, Sigma2=Sigma*Sigma, Sigma3=Sigma2*Sigma;
      double dm=r3-2.*b2_, dm3=dm*dm*dm;
      
      
      dst[1][1][1]=(a2_*r6*sth2 - 4.*a2b2*r3*sth2 - a2_*r5*sth2 + 4.*a2b4*sth2 + 8.*a2b2*r2*sth2 + a2_*r5 - r7 - 8.*a2b2*r2 - 4.*b2_*r4)*r/(a2_*r3 + r5 - 2.*a2b2 - 2.*b2_*r2 - 2.*r4)/Sigma/dm;
      dst[1][2][1]=dst[1][1][2]=-a2cthsth/Sigma;
      dst[1][2][2]=-(a2_*r3 + r5 - 2.*a2b2 - 2.*b2_*r2 - 2.*r4)*r/Sigma/dm;
      dst[1][3][3]=-(a4_*r6*cth4 - 4.*a4_*b2_*r3*cth4 + a4_*r5*cth2*sth2 + 2.*a2_*r8*cth2 + 4.*a4_*b4_*cth4 - 8.*a4_*b2_*r2*cth2*sth2 - 8.*a2b2*r5*cth2 - a2_*r7*sth2 + r10 + 8.*a2b4*r2*cth2 - 4.*a2b2*r4*sth2 - 4.*b2_*r7 + 4.*b4_*r4)*(a2_*r3 + r5 - 2.*a2b2 - 2.*b2_*r2 - 2.*r4)*r*sth2/Sigma3/dm3;
      dst[1][3][0]=dst[1][0][3]=(a2_*r3*cth2 - 8.*a2b2*cth2 - r5 - 4.*b2_*r2)*(a2_*r3 + r5 - 2.*a2b2 - 2.*b2_*r2 - 2.*r4)*spin_*r3*sth2/Sigma3/dm3;
      dst[1][0][0]=-(a2_*r3*cth2 - 8.*a2b2*cth2 - r5 - 4.*b2_*r2)*(a2_*r3 + r5 - 2.*a2b2 - 2.*b2_*r2 - 2.*r4)*r3/Sigma3/dm3;
      dst[2][1][1]=dm*a2cthsth/(a2_*r3 + r5 - 2.*a2b2 - 2.*b2_*r2 - 2.*r4)/Sigma;
      dst[2][2][1]=dst[2][1][2]=r/Sigma;
      dst[2][2][2]=-a2cthsth/Sigma;
      dst[2][3][3]=-(a6*r3*cth4 + a4_*r5*cth4 - 2*a6*b2_*cth4 - 2*a4_*b2_*r2*cth4 - 2.*a4_*r4*cth4 + 2.*a4_*r5*cth2 + 2.*a2_*r7*cth2 - 4.*a4_*b2_*r2*cth2 - 4.*a2b2*r4*cth2 - 4.*a2_*r6*cth2 + a2_*r7 + r9 + 2.*a4_*r4 - 2.*a2b2*r4 + 4.*a2_*r6 - 2.*b2_*r6)*cth*sth/Sigma3/dm;
      dst[2][0][3]=dst[2][3][0]=2*(a2_ + r2)*spin_*r4*cth*sth/Sigma3/dm;
      dst[2][0][0]=-2*a2_*r4*cth*sth/Sigma3/dm;
      dst[3][3][1]=dst[3][1][3]=(a4_*r6*cth4 - 4.*a4_*b2_*r3*cth4 + a4_*r5*cth2*sth2 + 2.*a2_*r8*cth2 + 4.*a4_*b4_*cth4 - 8.*a4_*b2_*r2*cth2*sth2 - 8.*a2b2*r5*cth2 - 2.*a2_*r7*cth2 - a2_*r7*sth2 + r10 + 8.*a2b4*r2*cth2 + 4.*a2b2*r4*cth2 - 4.*a2b2*r4*sth2 - 4.*b2_*r7 - 2.*r9 + 4.*b4_*r4 + 4.*b2_*r6)*r/(a2_*r3 + r5 - 2.*a2b2 - 2.*b2_*r2 - 2.*r4)/Sigma2/dm;
      dst[3][3][2]=dst[3][2][3]=(a4_*r3*cth4 - 2.*a4_*b2_*cth4 + 2.*a2_*r5*cth2 - 4.*a2b2* r2*cth2 - 2.*a2_*r4*cth2 + r7 + 2.*a2_*r4 - 2.*b2_*r4)*cth/Sigma2/dm/sth;
      dst[3][0][1]=dst[3][1][0]=-(a2_*r3*cth2 - 8.*a2b2*cth2 - r5 - 4.*b2_*r2)*spin_*r3/(a2_*r3 + r5 - 2.*a2b2 - 2.*b2_*r2 - 2.*r4)/Sigma2/dm;
      dst[3][0][2]=dst[3][2][0]=-2.*spin_*r4*cth/Sigma2/dm/sth;
      dst[0][3][1]=dst[0][1][3]=(a4_*r3*cth2 - a2_*r5*cth2 - 8.*a4_*b2_*cth2 - 4.*a2b2*r2*cth2 - a2_*r5 - 3*r7 - 4.*a2b2*r2)*spin_*r3*sth2/(a2_*r3 + r5 - 2.*a2b2 - 2.*b2_*r2 - 2.*r4)/Sigma2/dm;
      dst[0][3][2]=dst[0][2][3]=-2.*(a2_*sth2 - a2_ - r2)*a3*r4*cth*sth3/Sigma3/dm;
      dst[0][0][1]=dst[0][1][0]=-(a2_*r3*cth2 - 8.*a2b2*cth2 - r5 - 4.*b2_*r2)*(a2_ + r2)*r3/(a2_*r3 + r5 - 2.*a2b2 - 2.*b2_*r2 - 2.*r4)/Sigma2/dm;
      dst[0][0][2]=dst[0][2][0]=-2.*a2_*r4*cth*sth/Sigma2/dm;
      
      //cout.precision(16);
      //cout << "t=" << pos[0] << " r=" << pos[1] << " theta=" << pos[2] << " dt=" << pos[4] << " dr=" << pos[5] << " dth=" << pos[6] << " dph=" << pos[7] << endl ;
    }
  
  return 0;
}



// Optimized version
double Hayward::ScalarProd(const double* pos,
                           const double* u1, const double* u2) const {
  double g[4][4];
  gmunu(g, pos);
  double res=
    g[0][0]*u1[0]*u2[0]+
    g[1][1]*u1[1]*u2[1]+
    g[2][2]*u1[2]*u2[2]+
    g[3][3]*u1[3]*u2[3]+
    g[0][3]*u1[0]*u2[3]+
    g[3][0]*u1[3]*u2[0];
  
# if GYOTO_DEBUG_ENABLED
  GYOTO_IF_DEBUG
    GYOTO_DEBUG_ARRAY(pos, 4);
  GYOTO_DEBUG_ARRAY(u1, 4);
  GYOTO_DEBUG_ARRAY(u2, 4);
  GYOTO_DEBUG   << "ScalarProd(pos, u1, u2)="
		<<  res
		<< endl;
  GYOTO_ENDIF_DEBUG
# endif
    
    return res;
  
}
