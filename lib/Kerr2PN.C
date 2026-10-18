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
#include <cstring>
#include <cmath>

using namespace Gyoto;
using namespace Gyoto::Metric;
using namespace std;

/*
  ---------------------------------------------------------------------
  Kerr metric in harmonic Cartesian coordinates, truncated at 1PN / 2PN.

  Conventions (see the document "Christoffel Symbols in the Kerr Metric,
  Harmonic Coordinates, Second-Order Post-Newtonian Approximation", 
  K. Abd El Dayem):

   - Geometrized units G = c = M = 1.
   - Cartesian harmonic coordinates x^mu = (t, x, y, z), pos[0..3].
   - Spin axis along z: e = (0,0,1), J = (0,0,spin_).
   - n = x/r, n.e = z/r.

  PERFORMANCE NOTES
  -----------------
  christoffel() is the hot path of the integrator (called several times
  per step), gmunu(g,pos) is called by every norm / scalar-product check.
  Both are therefore written as straight-line code:

   * one sqrt and one division per call (everything else is multiplication
     by powers of ir = 1/r),
   * no Levi-Civita symbol, no function calls, no branches inside loops:
     with J = spin*z_hat the frame-dragging terms have a closed form
     (derived below) and only the non-zero components are written,
   * the PN order and spin tests are hoisted out of all loops.

  Closed form of the gravitomagnetic part, h_i := g_{0i} = 2 (x X J)_i / r^3
  with J = (0,0,chi):
      h = ( 2 chi y / r^3 , -2 chi x / r^3 , 0 )
      d_x h_x = -6 chi x y / r^5          d_y h_x =  2chi/r^3 - 6 chi y^2/r^5
      d_z h_x = -6 chi y z / r^5          d_x h_y = -2chi/r^3 + 6 chi x^2/r^5
      d_y h_y =  6 chi x y / r^5          d_z h_y =  6 chi x z / r^5
      d_m h_z = 0
  and
      Gamma^0_ij = -1/2 (d_i h_j + d_j h_i)
      Gamma^i_0j =  1/2 (d_j h_i - d_i h_j)      (antisymmetric in i,j)
  ---------------------------------------------------------------------
*/

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



// accessors

Gyoto::Metric::Kerr2PN::Kerr2PN() :
Generic(GYOTO_COORDKIND_CARTESIAN, "Kerr2PN"),
  spin_(0.), Q_(0.), quadrupole_is_kerr_(1), pnorder_(2)   // default: full 2PN
{
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

/*
  Full 4x4 covariant metric in one shot.

  Without this override Generic::gmunu(g,pos) calls the scalar
  gmunu(pos,mu,nu) 16 times, each one redoing the sqrt, the divisions
  and the branch cascade.  This version is used by ScalarProd() and thus
  by every norm check of the integrator.
*/
void Kerr2PN::gmunu(double g[4][4], const double * pos) const {
  const double xx = pos[1], yy = pos[2], zz = pos[3];
  const double r2 = xx*xx + yy*yy + zz*zz;
  const double ir = 1./sqrt(r2);
  const double ir2 = ir*ir, ir3 = ir2*ir, ir4 = ir2*ir2;

  double g00 = -1. + 2.*ir - 2.*ir2;      // 0PN + 1PN
  double A   =  1. + 2.*ir;               // 0PN + 1PN
  double g01 = 0., g02 = 0., B = 0.;

  if (pnorder_>=2) {
    const double ne = zz*ir;
    g00 += ir3*( 2. - spin_*spin_*(3.*ne*ne - 1.) );
    A   += ir2;
    B    = ir4;                           // coefficient of x_i x_j
    const double w = 2.*spin_*ir3;        // g_0i = 2 (x X J)_i / r^3
    g01 =  w*yy;
    g02 = -w*xx;
  }

  const double x[3] = {xx, yy, zz};
  g[0][0] = g00;
  g[0][1] = g[1][0] = g01;
  g[0][2] = g[2][0] = g02;
  g[0][3] = g[3][0] = 0.;
  for (int i=0; i<3; ++i) {
    g[i+1][i+1] = A + B*x[i]*x[i];
    for (int j=i+1; j<3; ++j)
      g[i+1][j+1] = g[j+1][i+1] = B*x[i]*x[j];
  }
}

double Kerr2PN::gmunu(const double * pos, int mu, int nu) const {
  const double xx = pos[1], yy = pos[2], zz = pos[3];
  const double r2 = xx*xx + yy*yy + zz*zz;
  const double ir = 1./sqrt(r2);
  const double ir2 = ir*ir, ir3 = ir2*ir;
  const bool pn2 = (pnorder_>=2);

  if (mu==0 && nu==0) {
    double val = -1. + 2.*ir - 2.*ir2;
    if (pn2) {
      const double ne = zz*ir;
      val += ir3*( 2. - spin_*spin_*(3.*ne*ne - 1.) );
    }
    return val;
  }

  if (mu==0 || nu==0) {
    if (!pn2) return 0.;
    const int ii = (mu==0) ? nu : mu;
    if (ii==1) return  2.*spin_*yy*ir3;
    if (ii==2) return -2.*spin_*xx*ir3;
    return 0.;
  }

  // spatial block
  const double x[4] = {0., xx, yy, zz};
  double val = (mu==nu) ? 1. + 2.*ir + (pn2 ? ir2 : 0.) : 0.;
  if (pn2) val += x[mu]*x[nu]*ir2*ir2;
  return val;
}


int Kerr2PN::christoffel(double dst[4][4][4], double const pos[4]) const
{
  std::memset(dst, 0, 64*sizeof(double));

  const double x[3] = {pos[1], pos[2], pos[3]};
  const double r2  = x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
  const double ir  = 1./sqrt(r2);                 // the only sqrt + division
  const double ir2 = ir*ir,  ir3 = ir2*ir, ir4 = ir2*ir2;
  const double ir5 = ir4*ir, ir6 = ir4*ir2;
  const bool   pn2 = (pnorder_>=2);

  // ---- Gamma^0_0i = x^i/r^3 (complete at every order) --------------
  for (int i=0; i<3; ++i) {
    const double v = x[i]*ir3;
    dst[0][0][i+1] = v;
    dst[0][i+1][0] = v;
  }

  // ---- Gamma^i_00  (monopole part) ---------------------------------
  {
    double b = 1. - 4.*ir;                  // 0PN + 1PN
    if (pn2) b += 9.*ir2;
    const double f = ir3*b;
    for (int i=0; i<3; ++i) dst[i+1][0][0] = x[i]*f;
  }

  // ---- Gamma^i_jk = c1 (d_ij x_k + d_ik x_j) + c2 d_jk x_i + c3 x_i x_j x_k
  {
    const double c1 = pn2 ? (ir4 - ir3) : -ir3;
    const double c2 = ir3;
    const double c3 = pn2 ? -2.*ir6 : 0.;
    for (int i=0; i<3; ++i) {
      const int a = (i==2) ? 0 : i+1;           // the two other indices
      const int b = (a==2) ? 0 : a+1;
      const double xi = x[i], xa = x[a], xb = x[b];
      const double c3xi = c3*xi;
      const double v_ii = (2.*c1 + c2)*xi + c3xi*xi*xi;
      const double v_ia = (c1 + c3xi*xi)*xa;    // Gamma^i_{i a}
      const double v_ib = (c1 + c3xi*xi)*xb;    // Gamma^i_{i b}
      const double w_aa = xi*(c2 + c3*xa*xa);   // Gamma^i_{a a}
      const double w_bb = xi*(c2 + c3*xb*xb);   // Gamma^i_{b b}
      const double w_ab = c3xi*xa*xb;           // Gamma^i_{a b}
      dst[i+1][i+1][i+1] = v_ii;
      dst[i+1][i+1][a+1] = dst[i+1][a+1][i+1] = v_ia;
      dst[i+1][i+1][b+1] = dst[i+1][b+1][i+1] = v_ib;
      dst[i+1][a+1][a+1] = w_aa;
      dst[i+1][b+1][b+1] = w_bb;
      dst[i+1][a+1][b+1] = dst[i+1][b+1][a+1] = w_ab;
    }
  }

  // ---- 2PN-only terms ------------------------------------------------
  if (pn2 && spin_!=0.) {
    const double chi = spin_, chi2 = chi*chi;

    // quadrupole part of Gamma^i_00
    const double ne = x[2]*ir;
    const double h  = 0.5*chi2;
    const double q  = -3.*h*ir5*(5.*ne*ne - 1.);
    dst[1][0][0] += q*x[0];
    dst[2][0][0] += q*x[1];
    dst[3][0][0] += q*x[2] + 6.*h*ne*ir4;

    // frame dragging (1.5PN): Gamma^0_ij and Gamma^i_0j, J = chi z_hat
    const double s   = chi*ir5;
    const double xy  = x[0]*x[1], yz = x[1]*x[2], xz = x[0]*x[2];
    const double rho2 = x[0]*x[0] + x[1]*x[1];

    dst[0][1][1] =  6.*s*xy;
    dst[0][2][2] = -6.*s*xy;
    dst[0][1][2] = dst[0][2][1] = -3.*s*(x[0]*x[0] - x[1]*x[1]);
    dst[0][1][3] = dst[0][3][1] =  3.*s*yz;
    dst[0][2][3] = dst[0][3][2] = -3.*s*xz;

    const double a = chi*(2.*ir3 - 3.*rho2*ir5);
    dst[1][0][2] = dst[1][2][0] =  a;
    dst[2][0][1] = dst[2][1][0] = -a;
    dst[1][0][3] = dst[1][3][0] = -3.*s*yz;
    dst[3][0][1] = dst[3][1][0] =  3.*s*yz;
    dst[2][0][3] = dst[2][3][0] =  3.*s*xz;
    dst[3][0][2] = dst[3][2][0] = -3.*s*xz;
  }

  return 0;
}


int Kerr2PN::isStopCondition(double const * const coord) const {
  // Compare r^2 to rsink^2: no pow(), no sqrt on the position.
  double rsink2 = 0.;
  if (spin_*spin_<1.){ // Black hole solution with event horizon
    const double rsink = 1. + sqrt(1. - spin_*spin_) + GYOTO_KERR_HORIZON_SECURITY;
    rsink2 = rsink*rsink;
  }
  const double xx=coord[1], yy=coord[2], zz=coord[3];
  return (xx*xx + yy*yy + zz*zz) < rsink2;
}