/*
  Copyright 2026 Karim Abd El Dayem, Frederic Vincent & Thibaut Paumard
  
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

#include "GyotoKerr2PN.h"

#include <iostream>
#include <cstdlib>
#include <cmath>

using namespace Gyoto;
using namespace Gyoto::Metric;
using namespace std;

/// Properties
#include "GyotoProperty.h"
GYOTO_PROPERTY_START(Kerr2PN, "Kerr spacetime at second post-Newtonian order, "
		     "in Cartesian harmonic coordinates")
GYOTO_PROPERTY_DOUBLE(Kerr2PN, Spin, spin,
		      "Spin parameter (dimensionless, 0).")
GYOTO_PROPERTY_DOUBLE(Kerr2PN, Quadrupole, quadrupole,
		      "Mass quadrupole parameter (dimensionless, -spin^2).")
GYOTO_PROPERTY_BOOL(Kerr2PN, KerrQuadrupole, NoKerrQuadrupole,
		    kerrQuadrupole,
		    "Should we fix the quadrupole to the Kerr value -spin^2? (bool, True)")
GYOTO_PROPERTY_END(Kerr2PN, Generic::properties)

// accessors

Gyoto::Metric::Kerr2PN::Kerr2PN() :
Generic(GYOTO_COORDKIND_CARTESIAN, "Kerr2PN"),
  spin_(0.), Q_(0.), quadrupole_is_kerr_(1)
{
  GYOTO_DEBUG << endl;
}

Gyoto::Metric::Kerr2PN::Kerr2PN(const Kerr2PN & orig)
  : Generic(orig), spin_(orig.spin_), Q_(orig.Q_),
    quadrupole_is_kerr_(orig.quadrupole_is_kerr_)
{
  GYOTO_DEBUG << endl;
}

// default copy constructor should be fine 
Kerr2PN * Kerr2PN::clone() const {
  return new Kerr2PN(*this); }

Gyoto::Metric::Kerr2PN::~Kerr2PN()
{
  GYOTO_DEBUG << endl;
}

// Mutators
void Kerr2PN::spin(const double spin) {
  spin_=spin;

  if (quadrupole_is_kerr_){
    // If quadrupole tied to Kerr spin, update it
    Q_=-spin*spin; // Kerr quadrupole formula Q2 = -G^2 M^3/c^4 * spin^2.
  } // else they are independent parameters
}
void Kerr2PN::quadrupole(const double QQ) {
  Q_=QQ;

  if (quadrupole_is_kerr_){
    // If quadrupole tied to Kerr spin, update spin
    if (Q_>0.)
      GYOTO_ERROR("Kerr quadrupole should be <0!");
    spin_=pow(-Q_,0.5); // Q_ = -spin^2 in geometrized units
  } // else they are independent parameters
}
void Kerr2PN::kerrQuadrupole(bool t) {
  quadrupole_is_kerr_=t;}

// Accessors
double Kerr2PN::spin() const { return spin_ ; }
double Kerr2PN::quadrupole() const { return Q_ ; }
bool Kerr2PN::kerrQuadrupole() const {
  return quadrupole_is_kerr_;}

double Kerr2PN::gmunu(const double * pos, int mu, int nu) const {
  double xx = pos[1], yy = pos[2], zz = pos[3],
    rr = pow(xx*xx+yy*yy+zz*zz,0.5);

  if ((mu==0) && (nu==0)) return -1.;
  if ((mu==1) && (nu==1)) return 1.;
  if ((mu==2) && (nu==2)) return 1.;
  if ((mu==3) && (nu==3)) return 1.;

  return 0.;
}

// double Kerr2PN::gmunu_up(const double * pos, int mu, int nu) const {
//   double rr = pos[1];
//   if (rr<=0.) GYOTO_ERROR("In Kerr2PN::gmunu: r<0!");

//   double sth2, cth2;
//   sincos(pos[2], &sth2, &cth2);
//   sth2*=sth2; cth2*=cth2;

//   //if ((mu==0) && (nu==0)) return -1./(1 - 2./rr + charge_*charge_/(rr*rr));
//   //if ((mu==1) && (nu==1)) return (1. - 2./rr + charge_*charge_/(rr*rr));
//   //if ((mu==2) && (nu==2)) return 1./(rr*rr);
//   //if ((mu==3) && (nu==3)) return 1./(rr*rr*sth2);

//   return 0.;
// } 

int Kerr2PN::christoffel(double dst[4][4][4], double const pos[4]) const
{
  int a, mu, nu;
  for (a=0; a<4; ++a)
    for(mu=0; mu<4; ++mu)
      for(nu=0; nu<4; ++nu)
	dst[a][mu][nu]=0.;

  //cout << "In Kerr2PN" << endl;

  return 0;
}

int Kerr2PN::isStopCondition(double const * const coord) const {
  double rsink = 0.;
  if (spin_*spin_<1.){ // Black hole solution with event horizon
    double rhor = 1 + sqrt(1 - spin_*spin_);
    rsink = rhor + GYOTO_KERR_HORIZON_SECURITY;
  }

  double xx=coord[1], yy=coord[2], zz=coord[3],
    rr = pow(xx*xx+yy*yy+zz*zz,0.5);

  return rr < rsink ;
}

// void Kerr2PN::circularVelocity(double const * coor, double* vel,
// 					     double dir) const {
//   double sinth = sin(coor[2]);
//   double coord[4] = {coor[0], coor[1]*sinth, M_PI*0.5, coor[3]};
//   if (charge_*charge_ > coord[1])
//     GYOTO_ERROR("RN::circularVelocity: bad radius to evaluate Omega");
  
//   vel[1] = vel[2] = 0.;
//   vel[3] = 1./(dir*pow(coord[1], 1.5)) * sqrt(1 - charge_*charge_/coord[1]);
  
//   vel[0] = SysPrimeToTdot(coor, vel+1);
//   vel[3] *= vel[0];
// }

// // PolishDoughnut specific functions
// double Kerr2PN::getPotential(double const pos[4], double l_cst) const {
//   // this is W = -ln(|u_t|) for a circular equatorial 4-velocity
//   // Careful this is not the same as in publications where W is defined
//   // as Wpubli = +ln(|u_t|). Thus, the center of the doughnut, which is
//   // the minimum of Wpubli, is the maximum of WGyoto.
  
//   double rr = pos[1], r2 = rr*rr, sth = sin(pos[2]), sth2 = sth*sth,
//     q2 = charge_*charge_,
//     term = r2 - 2.*rr + q2;
//   if (r2*sth2 - term * l_cst * l_cst / r2 == 0.){
//     cout << "At r,sth2= " << rr << " " << sth2 << endl;
//     GYOTO_ERROR("bad values in potential");
//   }
//   double logarg = term * sth2 / (r2*sth2 - term * l_cst * l_cst / r2);

//   double  gtt = gmunu(pos,0,0);
//   double  gpp = gmunu(pos,3,3);
  
//   if (logarg < 0){
//     cout << "At r,sth2= " << rr << " " << sth2 << endl;
//     cout << "ut2= " << -gtt*gpp/(gpp+l_cst*l_cst*gtt) << endl;
//     GYOTO_ERROR("bad values in potential");
//   }
//   double W = -1./2. * log(logarg);

//   // double  gtt = gmunu(pos,0,0);
//   // double  gtp = gmunu(pos,0,3);
//   // double  gpp = gmunu(pos,3,3);
//   // double  Omega = -(gtp + l_cst * gtt)/(gpp + l_cst * gtp) ;
  
//   // double  W = 0.5 * log(abs(gtt + 2. * Omega * gtp + Omega*Omega * gpp)) 
//   //   - log(abs(gtt + Omega * gtp)) ;
  
//   return  W ;
// }
// double Kerr2PN::getSpecificAngularMomentum(double rr) const {
//   // this is l = -u_phi/u_t for a circular equatorial 4-velocity
//   double qq=charge_, q2=qq*qq, r2=rr*rr;
//   if (rr<q2 or 1. - 2./rr + q2/r2 == 0.){
//     cout << "r, q2, term= " << rr << " " << q2 << " " << 1. - 2./rr + q2/r2 << endl;
//     GYOTO_ERROR("bad values in l !");
//   }
//   return sqrt(rr-q2)/(1. - 2./rr + q2/r2);
// }
