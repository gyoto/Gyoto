/*
  Copyright 2026 Frederic Vincent & Thibaut Paumard
  
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

#include "GyotoReissnerNordstrom.h"

#include <iostream>
#include <cstdlib>
#include <cmath>

using namespace Gyoto;
using namespace Gyoto::Metric;
using namespace std;

/// Properties
#include "GyotoProperty.h"
GYOTO_PROPERTY_START(ReissnerNordstrom, "Schwarzschild in harmonic coordinates")
GYOTO_PROPERTY_DOUBLE(ReissnerNordstrom, Charge, charge,
		      "Charge parameter (adimensioned, 0).")
GYOTO_PROPERTY_END(ReissnerNordstrom, Generic::properties)

// accessors

Gyoto::Metric::ReissnerNordstrom::ReissnerNordstrom() :
Generic(GYOTO_COORDKIND_SPHERICAL, "ReissnerNordstrom"),
  charge_(0.) 
{
  GYOTO_DEBUG << endl;
}

Gyoto::Metric::ReissnerNordstrom::ReissnerNordstrom(const ReissnerNordstrom & orig)
  : Generic(orig), charge_(orig.charge_)
{
  GYOTO_DEBUG << endl;
}

// default copy constructor should be fine 
ReissnerNordstrom * ReissnerNordstrom::clone() const {
  return new ReissnerNordstrom(*this); }

Gyoto::Metric::ReissnerNordstrom::~ReissnerNordstrom()
{
  GYOTO_DEBUG << endl;
}

// Mutators
void ReissnerNordstrom::charge(const double charge) {
  charge_=charge;
}

// Accessors
double ReissnerNordstrom::charge() const { return charge_ ; }

double ReissnerNordstrom::gmunu(const double * pos, int mu, int nu) const {
  double rr = pos[1];
  if (rr<=0.) GYOTO_ERROR("In ReissnerNordstrom::gmunu: r<0!");

  double sth2, cth2;
  sincos(pos[2], &sth2, &cth2);
  sth2*=sth2; cth2*=cth2;

  if ((mu==0) && (nu==0)) return -(1 - 2./rr + charge_*charge_/(rr*rr));
  if ((mu==1) && (nu==1)) return 1./(1. - 2./rr + charge_*charge_/(rr*rr));
  if ((mu==2) && (nu==2)) return rr*rr;
  if ((mu==3) && (nu==3)) return rr*rr*sth2;

  return 0.;
}

double ReissnerNordstrom::gmunu_up(const double * pos, int mu, int nu) const {
  double rr = pos[1];
  if (rr<=0.) GYOTO_ERROR("In ReissnerNordstrom::gmunu: r<0!");

  double sth2, cth2;
  sincos(pos[2], &sth2, &cth2);
  sth2*=sth2; cth2*=cth2;

  if ((mu==0) && (nu==0)) return -1./(1 - 2./rr + charge_*charge_/(rr*rr));
  if ((mu==1) && (nu==1)) return (1. - 2./rr + charge_*charge_/(rr*rr));
  if ((mu==2) && (nu==2)) return 1./(rr*rr);
  if ((mu==3) && (nu==3)) return 1./(rr*rr*sth2);

  return 0.;
} 

int ReissnerNordstrom::christoffel(double dst[4][4][4], double const pos[4]) const
{
  int a, mu, nu;
  for (a=0; a<4; ++a)
    for(mu=0; mu<4; ++mu)
      for(nu=0; nu<4; ++nu)
	dst[a][mu][nu]=0.;

  double rr = pos[1], r2=rr*rr, r3 = rr*r2;
  double sth, cth;
  sincos(pos[2], &sth, &cth);
  if (rr==0. || sth==0.) GYOTO_ERROR("In ReissnerNordstrom::christoffel: "
				    "bad coord");

  double q2 = charge_*charge_;
  
  dst[0][0][1]=dst[0][1][0]= -(q2 - rr) / (r3 - 2*r2 + q2*rr);
  dst[1][0][0]=-(q2*q2 - 3.*q2*rr - r3 + (2.+q2)*r2)/(r3*r2);
  dst[2][1][2]=dst[2][2][1]=1./rr;
  dst[3][1][3]=dst[3][3][1]=1./rr;
  dst[1][1][1]=(q2 - rr)/ (r3 - 2*r2 + q2*rr);
  dst[2][3][3]= - cth*sth;
  dst[3][2][3]=dst[3][3][2]= cth/sth;
  dst[1][2][2]= - (q2 - 2*rr + r2)/rr;
  dst[1][3][3]= - (q2 - 2*rr + r2)/rr * sth*sth;

  return 0;
}

int ReissnerNordstrom::isStopCondition(double const * const coord) const {
  double rsink = 0.;
  if (charge_<1.){ // RN black hole solution with event horizon
    double rhor = 1 + sqrt(1 - charge_*charge_);
    rsink = rhor + GYOTO_KERR_HORIZON_SECURITY;
  }  

  return coord[1] < rsink ;
}

void ReissnerNordstrom::circularVelocity(double const * coor, double* vel,
					     double dir) const {
  double sinth = sin(coor[2]);
  double coord[4] = {coor[0], coor[1]*sinth, M_PI*0.5, coor[3]};
  if (charge_*charge_ > coord[1])
    GYOTO_ERROR("RN::circularVelocity: bad radius to evaluate Omega");
  
  vel[1] = vel[2] = 0.;
  vel[3] = 1./(dir*pow(coord[1], 1.5)) * sqrt(1 - charge_*charge_/coord[1]);
  
  vel[0] = SysPrimeToTdot(coor, vel+1);
  vel[3] *= vel[0];
}

// PolishDoughnut specific functions
double ReissnerNordstrom::getPotential(double const pos[4], double l_cst) const {
  // this is W = -ln(|u_t|) for a circular equatorial 4-velocity
  // Careful this is not the same as in publications where W is defined
  // as Wpubli = +ln(|u_t|). Thus, the center of the doughnut, which is
  // the minimum of Wpubli, is the maximum of WGyoto.
  
  double rr = pos[1], r2 = rr*rr, sth = sin(pos[2]), sth2 = sth*sth,
    q2 = charge_*charge_,
    term = r2 - 2.*rr + q2;
  if (r2*sth2 - term * l_cst * l_cst / r2 == 0.){
    cout << "At r,sth2= " << rr << " " << sth2 << endl;
    GYOTO_ERROR("bad values in potential");
  }
  double logarg = term * sth2 / (r2*sth2 - term * l_cst * l_cst / r2);

  double  gtt = gmunu(pos,0,0);
  double  gpp = gmunu(pos,3,3);
  
  if (logarg < 0){
    cout << "At r,sth2= " << rr << " " << sth2 << endl;
    cout << "ut2= " << -gtt*gpp/(gpp+l_cst*l_cst*gtt) << endl;
    GYOTO_ERROR("bad values in potential");
  }
  double W = -1./2. * log(logarg);

  // double  gtt = gmunu(pos,0,0);
  // double  gtp = gmunu(pos,0,3);
  // double  gpp = gmunu(pos,3,3);
  // double  Omega = -(gtp + l_cst * gtt)/(gpp + l_cst * gtp) ;
  
  // double  W = 0.5 * log(abs(gtt + 2. * Omega * gtp + Omega*Omega * gpp)) 
  //   - log(abs(gtt + Omega * gtp)) ;
  
  return  W ;
}
double ReissnerNordstrom::getSpecificAngularMomentum(double rr) const {
  // this is l = -u_phi/u_t for a circular equatorial 4-velocity
  double qq=charge_, q2=qq*qq, r2=rr*rr;
  if (rr<q2 or 1. - 2./rr + q2/r2 == 0.){
    cout << "r, q2, term= " << rr << " " << q2 << " " << 1. - 2./rr + q2/r2 << endl;
    GYOTO_ERROR("bad values in l !");
  }
  return sqrt(rr-q2)/(1. - 2./rr + q2/r2);
}
