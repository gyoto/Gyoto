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

/*
  ---------------------------------------------------------------------
  Helper functions for the 2PN harmonic-coordinate Kerr metric.

  Conventions (see the document "Christoffel Symbols in the Kerr Metric,
  Harmonic Coordinates, Second-Order Post-Newtonian Approximation", 
  K. Abd El Dayem):

   - Geometrized units with G = c = M = 1 (as is standard in Gyoto:
     the Mass Property only rescales physical <-> geometrical units,
     it never appears explicitly inside gmunu()/christoffel()).
   - Cartesian harmonic coordinates x^mu = (t, x, y, z), pos[0..3].
   - The spin vector direction e is NOT stored by this class (only the
     scalar spin_ = chi is), so we adopt the standard convention that
     the spin axis is aligned with z: e = (0,0,1), J = (0,0,spin_).
   - n = x/r is the unit vector from the SMBH to the field point,
     n.e = z/r.
  ---------------------------------------------------------------------
*/

namespace {  

  // 3D Levi-Civita symbol, indices in {0,1,2} <-> {x,y,z}
  inline double epsilon3(int i, int j, int k) {
    if ( (i==0&&j==1&&k==2) || (i==1&&j==2&&k==0) || (i==2&&j==0&&k==1) )
      return 1.;
    if ( (i==0&&j==2&&k==1) || (i==2&&j==1&&k==0) || (i==1&&j==0&&k==2) )
      return -1.;
    return 0.;
  }

  // d_m g_{0l} = 2 eps_{lmk} J^k / r^3 - 6 eps_{lnk} J^k x^m x^n / r^5
  // m, l in {0,1,2} <-> {x,y,z}; xv = (x,y,z); Jv = (Jx,Jy,Jz)
  inline double dg0(int m, int l,
		     const double xv[3], const double Jv[3],
		     double r3, double r5) {
    double term1 = 0.;
    for (int k=0; k<3; ++k) term1 += 2.*epsilon3(l,m,k)*Jv[k];
    term1 /= r3;

    double term2 = 0.;
    for (int n=0; n<3; ++n)
      for (int k=0; k<3; ++k)
	term2 += -6.*epsilon3(l,n,k)*Jv[k]*xv[m]*xv[n];
    term2 /= r5;

    return term1 + term2;
  }

}

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
GYOTO_PROPERTY_DOUBLE(Kerr2PN, PNOrder, pnOrder,
		      "Post-Newtonian truncation order: 1 or 2 (default 2).")
GYOTO_PROPERTY_END(Kerr2PN, Generic::properties)




// TESTTTTTTTTTTTTTTTTTTTTTTTTTT TO SEE IF GYOTO UPDATESSSSSSSSSSSSSSSSSSSSSSSSSSSSS
#include <iostream>
#include <stdexcept>

// accessors

Gyoto::Metric::Kerr2PN::Kerr2PN() :
Generic(GYOTO_COORDKIND_CARTESIAN, "Kerr2PN"),
  spin_(0.), Q_(0.), quadrupole_is_kerr_(1), pnorder_(2)   // default: full 2PN
{
  // TESTTTTTTTTTTTTTTTTTTTTTTTTTT TO SEE IF GYOTO UPDATESSSSSSSSSSSSSSSSSSSSSSSSSSSSS
  // --- YOUR DEBUG CODE GOES HERE ---
  std::cout << "\n==========================================" << std::endl;
  std::cout << "   >>> KERR2PN MY NEW2 CODE IS RUNNING <<<   " << std::endl;
  std::cout << "==========================================\n" << std::endl;
  // throw std::runtime_error("FORCE STOP: KERR2PN UPDATED SUCCESSFULLY!");
  // // ---------------------------------
  // TESTTTTTTTTTTTTTTTTTTTTTTTTTT TO SEE IF GYOTO UPDATESSSSSSSSSSSSSSSSSSSSSSSSSSSSS
  // --- YOUR DEBUG CODE GOES HERE ---
  std::cout << "\n==========================================" << std::endl;
  std::cout << "   >>> KERRHARMONIC NEW CODE IS RUNNING <<<   " << std::endl;
  std::cout << "==========================================\n" << std::endl;
  // throw std::runtime_error("FORCE STOP: KERR2PN UPDATED SUCCESSFULLY!");
  // // ---------------------------------

  GYOTO_DEBUG << endl;
}


Gyoto::Metric::Kerr2PN::Kerr2PN(const Kerr2PN & orig)
  : Generic(orig), spin_(orig.spin_), Q_(orig.Q_),
    quadrupole_is_kerr_(orig.quadrupole_is_kerr_),
    pnorder_(orig.pnorder_)
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
  quadrupole_is_kerr_=t;
}
void Kerr2PN::pnOrder(double n) {
  int ni = int(n+0.5);
  if (ni!=1 && ni!=2)
    GYOTO_ERROR("Kerr2PN::pnOrder: order must be 1 or 2");
  pnorder_ = (unsigned int)(ni);
}
double Kerr2PN::pnOrder() const { return double(pnorder_); }

// Accessors
double Kerr2PN::spin() const { return spin_ ; }
double Kerr2PN::quadrupole() const { return Q_ ; }
bool Kerr2PN::kerrQuadrupole() const {
  return quadrupole_is_kerr_;}

// double Kerr2PN::gmunu(const double * pos, int mu, int nu) const {
//   double xx = pos[1], yy = pos[2], zz = pos[3];
//   double r2 = xx*xx+yy*yy+zz*zz, rr = sqrt(r2), r3 = r2*rr, r4 = r2*r2;
//   double ne = zz/rr;             // n.e, e=(0,0,1)
//   double chi2 = spin_*spin_;

//   // g_00  (Eq. 4)
//   if (mu==0 && nu==0) {
//     double fQ = 3.*ne*ne - 1.;
//     return -1. + 2./rr - 2./r2 + 2./r3 - chi2/r3*fQ;
//   }

//   // g_0i = g_i0  (Eq. 5), J=(0,0,spin_)
//   if ( (mu==0 && nu>=1 && nu<=3) || (nu==0 && mu>=1 && mu<=3) ) {
//     int ii = (mu==0) ? nu : mu; // spatial index (1..3) of the requested component
//     double xv[3] = {xx, yy, zz};
//     double Jv[3] = {0., 0., spin_};
//     double cross[3];
//     cross[0] = xv[1]*Jv[2] - xv[2]*Jv[1];
//     cross[1] = xv[2]*Jv[0] - xv[0]*Jv[2];
//     cross[2] = xv[0]*Jv[1] - xv[1]*Jv[0];
//     return 2.*cross[ii-1]/r3;
//   }

//   // g_ij  (Eq. 6)
//   if (mu>=1 && mu<=3 && nu>=1 && nu<=3) {
//     double xv[3] = {xx, yy, zz};
//     double delta = (mu==nu) ? 1. : 0.;
//     double A = 1. + 2./rr + 1./r2;
//     return A*delta + xv[mu-1]*xv[nu-1]/r4;
//   }

//   return 0.;
// }

double Kerr2PN::gmunu(const double * pos, int mu, int nu) const {
  double xx = pos[1], yy = pos[2], zz = pos[3];
  double r2 = xx*xx+yy*yy+zz*zz, rr = sqrt(r2), r3 = r2*rr, r4 = r2*r2;
  double ne = zz/rr;
  double chi2 = spin_*spin_;

  // g_00  (Eq. 4)
  if (mu==0 && nu==0) {
    double val = -1.;
    val += 2./rr - 2./r2;                       // Newtonian, always present
    if (pnorder_>=2) {
      val +=  + 2./r3;
      double fQ = 3.*ne*ne - 1.;
      val += -chi2/r3*fQ;
    }
    // // TESTTTTTTTTTTTTTTTTTTTTTTTTTT TO SEE IF GYOTO UPDATESSSSSSSSSSSSSSSSSSSSSSSSSSSSS
    // if (pnorder_<=1) {
    //   val =0;
    // }
    return val;
  }

  // g_0i = g_i0  (Eq. 5) -- 1.5PN, needs the full 2PN treatment to be consistent
  if ( (mu==0 && nu>=1 && nu<=3) || (nu==0 && mu>=1 && mu<=3) ) {
    if (pnorder_<2) return 0.;
    int ii = (mu==0) ? nu : mu;
    double xv[3] = {xx, yy, zz};
    double Jv[3] = {0., 0., spin_};
    double cross[3];
    cross[0] = xv[1]*Jv[2] - xv[2]*Jv[1];
    cross[1] = xv[2]*Jv[0] - xv[0]*Jv[2];
    cross[2] = xv[0]*Jv[1] - xv[1]*Jv[0];
    return 2.*cross[ii-1]/r3;
  }

  // g_ij  (Eq. 6)
  if (mu>=1 && mu<=3 && nu>=1 && nu<=3) {
    double xv[3] = {xx, yy, zz};
    double delta = (mu==nu) ? 1. : 0.;
    double A = 1. + 2./rr;              // 0PN + 1PN, always present
    if (pnorder_>=2) A += 1./r2;
    double val = A*delta;
    if (pnorder_>=2) val += xv[mu-1]*xv[nu-1]/r4;
    return val;
  }

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

// int Kerr2PN::christoffel(double dst[4][4][4], double const pos[4]) const
// {
//   int a, mu, nu;
//   for (a=0; a<4; ++a)
//     for(mu=0; mu<4; ++mu)
//       for(nu=0; nu<4; ++nu)
// 	dst[a][mu][nu]=0.;

//   double xx = pos[1], yy = pos[2], zz = pos[3];
//   double xv[3] = {xx, yy, zz};
//   double Jv[3] = {0., 0., spin_};
//   double ev[3] = {0., 0., 1.};
//   double r2 = xx*xx+yy*yy+zz*zz, rr = sqrt(r2);
//   double r3 = r2*rr, r4 = r2*r2, r5 = r4*rr, r6 = r4*r2;
//   double ne = zz/rr;
//   double chi2 = spin_*spin_;

//   // Gamma^0_00 = 0 (Eq. 103) -- already set above, nothing to do.

//   // Gamma^0_0i = Gamma^0_i0 = x^i / r^3  (Eq. 108)
//   for (int ii=0; ii<3; ++ii) {
//     double val = xv[ii]/r3;
//     dst[0][0][ii+1] = val;
//     dst[0][ii+1][0] = val;
//   }

//   // Gamma^0_ij = -1/2 (d_i g0j + d_j g0i)  (Eq. 111), 1.5PN spin term
//   for (int ii=0; ii<3; ++ii)
//     for (int jj=0; jj<3; ++jj)
//       dst[0][ii+1][jj+1] =
// 	-0.5*( dg0(ii,jj,xv,Jv,r3,r5) + dg0(jj,ii,xv,Jv,r3,r5) );

//   // Gamma^i_00  (Eq. 123)
//   for (int ii=0; ii<3; ++ii) {
//     double xi = xv[ii], ei = ev[ii];
//     double val = xi/r3*(1. - 4./rr + 11./r2)
//       + chi2/2.*( 6.*ne*ei/r4 - 3.*xi/r5*(5.*ne*ne - 1.) );
//     dst[ii+1][0][0] = val;
//   }

//   // Gamma^i_0j = Gamma^i_j0 = 1/2 (d_j g0i - d_i g0j)  (Eq. 128), spin-orbit
//   for (int ii=0; ii<3; ++ii) {
//     for (int jj=0; jj<3; ++jj) {
//       double val = 0.5*( dg0(jj,ii,xv,Jv,r3,r5) - dg0(ii,jj,xv,Jv,r3,r5) );
//       dst[ii+1][0][jj+1] = val;
//       dst[ii+1][jj+1][0] = val;
//     }
//   }

//   // Gamma^i_jk  (Eq. 139)
//   for (int ii=0; ii<3; ++ii) {
//     for (int jj=0; jj<3; ++jj) {
//       for (int kk=0; kk<3; ++kk) {
// 	double dij = (ii==jj) ? 1. : 0.;
// 	double dik = (ii==kk) ? 1. : 0.;
// 	double djk = (jj==kk) ? 1. : 0.;
// 	double val =
// 	  1./r3*( dij*xv[kk] + dik*xv[jj] - djk*xv[ii] )*(-1. + 2./rr)
// 	  + 1./r4*djk*xv[ii]
// 	  - 2./r6*xv[ii]*xv[jj]*xv[kk];
// 	dst[ii+1][jj+1][kk+1] = val;
//       }
//     }
//   }

//   return 0;
// }


int Kerr2PN::christoffel(double dst[4][4][4], double const pos[4]) const
{
  int a, mu, nu;
  for (a=0; a<4; ++a)
    for(mu=0; mu<4; ++mu)
      for(nu=0; nu<4; ++nu)
	dst[a][mu][nu]=0.;

  double xx = pos[1], yy = pos[2], zz = pos[3];
  double xv[3] = {xx, yy, zz};
  double Jv[3] = {0., 0., spin_};
  double ev[3] = {0., 0., 1.};
  double r2 = xx*xx+yy*yy+zz*zz, rr = sqrt(r2);
  double r3 = r2*rr, r4 = r2*r2, r5 = r4*rr, r6 = r4*r2;
  double ne = zz/rr;
  double chi2 = spin_*spin_;

  // Gamma^0_00 = 0 (Eq. 103) -- unchanged, always zero.

  // Gamma^0_0i = x^i / r^3  (Eq. 108) -- already the complete needed
  // expression at every order (its own 2PN correction cancels identically,
  // see the note in the Christoffel document) -- no gating needed.
  for (int ii=0; ii<3; ++ii) {
    double val = xv[ii]/r3;
    dst[0][0][ii+1] = val;
    dst[0][ii+1][0] = val;
  }

  // Gamma^0_ij -- 1.5PN frame-dragging, needs the full 2PN treatment
  if (pnorder_>=2) {
    for (int ii=0; ii<3; ++ii)
      for (int jj=0; jj<3; ++jj)
        dst[0][ii+1][jj+1] =
          -0.5*( dg0(ii,jj,xv,Jv,r3,r5) + dg0(jj,ii,xv,Jv,r3,r5) );
  }

  // Gamma^i_00  (Eq. 123)  <-- the "11." bug is fixed back to "3." here
  for (int ii=0; ii<3; ++ii) {
    double xi = xv[ii], ei = ev[ii];
    double bracket = 1. - 4./rr;         // 0PN + 1PN, always present
    if (pnorder_>=2) bracket += 9./r2;
    double val = xi/r3*bracket;
    if (pnorder_>=2)
      val += chi2/2.*( 6.*ne*ei/r4 - 3.*xi/r5*(5.*ne*ne - 1.) );
    dst[ii+1][0][0] = val;
  }

  // Gamma^i_0j -- spin-orbit, needs the full 2PN treatment
  if (pnorder_>=2) {
    for (int ii=0; ii<3; ++ii) {
      for (int jj=0; jj<3; ++jj) {
        double val = 0.5*( dg0(jj,ii,xv,Jv,r3,r5) - dg0(ii,jj,xv,Jv,r3,r5) );
        dst[ii+1][0][jj+1] = val;
        dst[ii+1][jj+1][0] = val;
      }
    }
  }

  // Gamma^i_jk  (Eq. 140, corrected)
  for (int ii=0; ii<3; ++ii) {
    for (int jj=0; jj<3; ++jj) {
      for (int kk=0; kk<3; ++kk) {
	double dij = (ii==jj) ? 1. : 0.;
	double dik = (ii==kk) ? 1. : 0.;
	double djk = (jj==kk) ? 1. : 0.;
	double val = -1./r3*( dij*xv[kk] + dik*xv[jj] - djk*xv[ii] );  // 1PN, always present
	if (pnorder_>=2)
	  val += 1./r4*( dij*xv[kk] + dik*xv[jj] ) - 2./r6*xv[ii]*xv[jj]*xv[kk];
	dst[ii+1][jj+1][kk+1] = val;
      }
    }
  }

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
