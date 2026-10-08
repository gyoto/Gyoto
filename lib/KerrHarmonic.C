/*
  Copyright 2026 Frederic Vincent, Thibaut Paumard & Karim Abd El Dayem
  
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

#include "GyotoKerrHarmonic.h"

#include <iostream>
#include <cstdlib>
#include <cmath>

using namespace Gyoto;
using namespace Gyoto::Metric;
using namespace std;

/// Properties
#include "GyotoProperty.h"
GYOTO_PROPERTY_START(KerrHarmonic, "Kerr spacetime in Cartesian harmonic coordinates")
GYOTO_PROPERTY_DOUBLE(KerrHarmonic, Spin, spin,
		      "Spin parameter (dimensionless, 0).")
GYOTO_PROPERTY_END(KerrHarmonic, Generic::properties)

// accessors

Gyoto::Metric::KerrHarmonic::KerrHarmonic() :
Generic(GYOTO_COORDKIND_CARTESIAN, "KerrHarmonic"),
  spin_(0.)
{
  GYOTO_DEBUG << endl;
}

Gyoto::Metric::KerrHarmonic::KerrHarmonic(const KerrHarmonic & orig)
  : Generic(orig), spin_(orig.spin_)
{
  GYOTO_DEBUG << endl;
}

// default copy constructor should be fine 
KerrHarmonic * KerrHarmonic::clone() const {
  return new KerrHarmonic(*this); }

Gyoto::Metric::KerrHarmonic::~KerrHarmonic()
{
  GYOTO_DEBUG << endl;
}

// Mutators
void KerrHarmonic::spin(const double spin) {
  spin_=spin;
}

// Accessors
double KerrHarmonic::spin() const { return spin_ ; }

void KerrHarmonic::gmunu(double g[4][4], const double * pos) const {

  // DEFINE QUANTITIES
  /*
   The details of the computations of the metric coefficients
   are provided in a pdf note by FV (ask me if interested).
   */
  double aa = spin_, a2 = aa*aa;
  double mcal2 = 1-a2; // this is \mathcal{M}^2 = M^2 - a^2 in FV notes 
  // Harmonic Cartesian coordinates and harmonic radius
  double xx = pos[1], yy = pos[2], zz = pos[3],
    //double xx = 10., yy = 9., zz = 8., // test!!!!
    x2 = xx*xx, y2 = yy*yy, z2 = zz*zz,
    rr = pow(x2+y2+z2,0.5), r2 = rr*rr;
  double Rcal2 = pow((r2 - a2)*(r2 - a2) + 4.*a2*z2,0.5); // useful quantity,
                                            // this is \mathcal{R} in FV notes.

  // BL coordinates derived from harmonic
  double rBL = 1. + 1/pow(2.,0.5)*pow(r2 - a2 + Rcal2,0.5), // BL radius
    uu = rBL-1., // u = rBL-M is a useful quantity, see FV notes
    u2 = uu*uu,
    costhBL = zz/uu, costhBL2 = costhBL*costhBL, // BL cos(theta)
    sinthBL = pow(1 - costhBL2,0.5), sinthBL2 = sinthBL*sinthBL,
    rho2BL = rBL*rBL + a2 * costhBL2, // usual rho2 quantity in BL metric
    DeltaBL = rBL*rBL - 2.*rBL + a2; // usual Delta quantity in BL metric

  // Covariant BL metric coefficients expressed in
  // harmonic coordinates
  double gtt_BL = -(1. - 2*rBL/rho2BL),
    gtp_BL = -2*aa*rBL*sinthBL2/rho2BL,
    grr_BL = rho2BL/DeltaBL, 
    gthth_BL = rho2BL,
    gpp_BL = (rBL*rBL + a2 + 2.*a2*rBL*sinthBL2/rho2BL)*sinthBL2;

  // Derivative of BL coord wrt harmonic ones
  double
    drdx = uu*xx/Rcal2,
    drdy = uu*yy/Rcal2,
    drdz = (u2 + a2)*zz/(uu*Rcal2),
    dthdx = xx*zz/(Rcal2*pow(u2-z2,0.5)),
    dthdy = yy*zz/(Rcal2*pow(u2-z2,0.5)),
    dthdz = 1./pow(u2-z2,0.5) * (z2*(u2+a2)/(Rcal2*u2) - 1.),
    dphidx = aa/(a2+u2) * drdx - yy/(x2+y2) - drdx*aa/(u2 - mcal2), 
    dphidy = aa/(a2+u2) * drdy + xx/(x2+y2) - drdy*aa/(u2 - mcal2),
    dphidz = aa/(a2+u2) * drdz - drdz*aa/(u2 - mcal2);

  // FILL GMUNU

  g[0][0] = gtt_BL;

  g[0][1] = g[1][0] = gtp_BL * dphidx;
  g[0][2] = g[2][0] = gtp_BL * dphidy;
  g[0][3] = g[3][0] = gtp_BL * dphidz;

  g[1][1] = grr_BL * drdx * drdx + gthth_BL * dthdx * dthdx
      + gpp_BL * dphidx * dphidx;
  g[2][2] = grr_BL * drdy * drdy + gthth_BL * dthdy * dthdy
      + gpp_BL * dphidy * dphidy;
  g[3][3] = grr_BL * drdz * drdz + gthth_BL * dthdz * dthdz
      + gpp_BL * dphidz * dphidz;

  g[1][2] = g[2][1] = grr_BL * drdx * drdy + gthth_BL * dthdx * dthdy
      + gpp_BL * dphidx * dphidy;
  g[1][3] = g[3][1] = grr_BL * drdx * drdz + gthth_BL * dthdx * dthdz
      + gpp_BL * dphidx * dphidz;
  g[2][3] = g[3][2] = grr_BL * drdy * drdz + gthth_BL * dthdy * dthdz
      + gpp_BL * dphidy * dphidz;
  
}

double KerrHarmonic::gmunu(const double * pos, int mu, int nu) const {
  /*
   The details of the computations of the metric coefficients
   are provided in a pdf note by FV (ask me if interested).
   */
  double aa = spin_, a2 = aa*aa;
  double mcal2 = 1-a2; // this is \mathcal{M}^2 = M^2 - a^2 in FV notes 
  // Harmonic Cartesian coordinates and harmonic radius
  double xx = pos[1], yy = pos[2], zz = pos[3],
    //double xx = 10., yy = 9., zz = 8., // test!!!!
    x2 = xx*xx, y2 = yy*yy, z2 = zz*zz,
    rr = pow(x2+y2+z2,0.5), r2 = rr*rr;
  double Rcal2 = pow((r2 - a2)*(r2 - a2) + 4.*a2*z2,0.5); // useful quantity,
                                            // this is \mathcal{R} in FV notes.

  // BL coordinates derived from harmonic
  double rBL = 1. + 1/pow(2.,0.5)*pow(r2 - a2 + Rcal2,0.5), // BL radius
    uu = rBL-1., // u = rBL-M is a useful quantity, see FV notes
    u2 = uu*uu,
    costhBL = zz/uu, costhBL2 = costhBL*costhBL, // BL cos(theta)
    sinthBL = pow(1 - costhBL2,0.5), sinthBL2 = sinthBL*sinthBL,
    rho2BL = rBL*rBL + a2 * costhBL2, // usual rho2 quantity in BL metric
    DeltaBL = rBL*rBL - 2.*rBL + a2; // usual Delta quantity in BL metric

  // Covariant BL metric coefficients expressed in
  // harmonic coordinates
  double gtt_BL = -(1. - 2*rBL/rho2BL),
    gtp_BL = -2*aa*rBL*sinthBL2/rho2BL,
    grr_BL = rho2BL/DeltaBL, 
    gthth_BL = rho2BL,
    gpp_BL = (rBL*rBL + a2 + 2.*a2*rBL*sinthBL2/rho2BL)*sinthBL2;

  // Derivative of BL coord wrt harmonic ones
  double
    drdx = uu*xx/Rcal2,
    drdy = uu*yy/Rcal2,
    drdz = (u2 + a2)*zz/(uu*Rcal2),
    dthdx = xx*zz/(Rcal2*pow(u2-z2,0.5)),
    dthdy = yy*zz/(Rcal2*pow(u2-z2,0.5)),
    dthdz = 1./pow(u2-z2,0.5) * (z2*(u2+a2)/(Rcal2*u2) - 1.),
    dphidx = aa/(a2+u2) * drdx - yy/(x2+y2) - drdx*aa/(u2 - mcal2), 
    dphidy = aa/(a2+u2) * drdy + xx/(x2+y2) - drdy*aa/(u2 - mcal2),
    dphidz = aa/(a2+u2) * drdz - drdz*aa/(u2 - mcal2);

  if ((mu==0) && (nu==0)) return gtt_BL;
  
  if (((mu==0) && (nu==1)) || ((mu==1) && (nu==0)))
    return gtp_BL * dphidx;
  if (((mu==0) && (nu==2)) || ((mu==2) && (nu==0)))
    return gtp_BL * dphidy;
  if (((mu==0) && (nu==3)) || ((mu==3) && (nu==0)))
    return gtp_BL * dphidz;
  
  if ((mu==1) && (nu==1))
    return grr_BL * drdx * drdx + gthth_BL * dthdx * dthdx
      + gpp_BL * dphidx * dphidx;
  if ((mu==2) && (nu==2))
    return grr_BL * drdy * drdy + gthth_BL * dthdy * dthdy
      + gpp_BL * dphidy * dphidy;
  if ((mu==3) && (nu==3))
    return grr_BL * drdz * drdz + gthth_BL * dthdz * dthdz
      + gpp_BL * dphidz * dphidz;

  if (((mu==1) && (nu==2)) || ((mu==2) && (nu==1)))
    return grr_BL * drdx * drdy + gthth_BL * dthdx * dthdy
      + gpp_BL * dphidx * dphidy;
  if (((mu==1) && (nu==3)) || ((mu==3) && (nu==1)))
    return grr_BL * drdx * drdz + gthth_BL * dthdx * dthdz
      + gpp_BL * dphidx * dphidz;
  if (((mu==2) && (nu==3)) || ((mu==3) && (nu==2)))
    return grr_BL * drdy * drdz + gthth_BL * dthdy * dthdz
      + gpp_BL * dphidy * dphidz;

  GYOTO_ERROR("Should have returned before!");
  return 0.;
  
}

void KerrHarmonic::gmunu_up(double gup[4][4], const double * pos) const {

  //cout << "KH gmunuup" << endl;

  double aa = spin_, a2 = aa*aa;
  double mcal2 = 1-a2; // this is \mathcal{M}^2 = M^2 - a^2 in FV notes
  // Harmonic Cartesian coordinates and harmonic radius
  double xx = pos[1], yy = pos[2], zz = pos[3],
    //double xx = 10., yy = 9., zz = 8., // test!!!!
    x2 = xx*xx, y2 = yy*yy, z2 = zz*zz,
    rr = pow(x2+y2+z2,0.5), r2 = rr*rr;
  double Rcal2 = pow((r2 - a2)*(r2 - a2) + 4.*a2*z2,0.5); // useful quantity,
                                            // this is \mathcal{R} in FV notes.

  // BL coordinates derived from harmonic
  double rBL = 1. + 1/pow(2.,0.5)*pow(r2 - a2 + Rcal2,0.5), // BL radius
    uu = rBL-1., // u = rBL-M is a useful quantity, see FV notes
    u2 = uu*uu,
    costhBL = zz/uu, costhBL2 = costhBL*costhBL, // BL cos(theta)
    sinthBL = pow(1 - costhBL2,0.5), sinthBL2 = sinthBL*sinthBL,
    rho2BL = rBL*rBL + a2 * costhBL2, // usual rho2 quantity in BL metric
    DeltaBL = rBL*rBL - 2.*rBL + a2, // usual Delta quantity in BL metric
    cosphBL = (uu*xx + aa*yy)/(x2+y2) * sinthBL, // playing with the x ans y relations
    sinphBL = (uu*yy - aa*xx)/(x2+y2) * sinthBL;
  
  // BL inverse metric
  double gtt_BL = -1./DeltaBL * (rBL*rBL + a2 + 2.*a2*rBL*sinthBL2/rho2BL), // contravariant indices throughout
    grr_BL = DeltaBL/rho2BL,
    gthth_BL = 1./rho2BL,
    gpp_BL = 1./(DeltaBL*sinthBL2) * (1. - 2.*rBL/rho2BL),
    gtp_BL = -2.*aa*rBL/(DeltaBL*rho2BL);

  // Derivative of H coord wrt BL ones
  double
    dxdr = cosphBL * sinthBL - aa*yy/(u2-mcal2),
    dydr = sinphBL * sinthBL + aa*xx/(u2-mcal2),
    dzdr = costhBL,
    dxdth = xx*costhBL/sinthBL,
    dydth = yy*costhBL/sinthBL,
    dzdth = -uu*sinthBL,
    dxdp = -yy,
    dydp = xx,
    dzdp = 0.;

  gup[0][0] = gtt_BL; // g^tt

  gup[0][1] = gup[1][0] = gtp_BL*dxdp; // g^tx
  gup[0][2] = gup[2][0] = gtp_BL*dydp; // g^ty
  gup[0][3] = gup[3][0] = gtp_BL*dzdp; // g^tz

  gup[1][1] = grr_BL * dxdr*dxdr + gthth_BL * dxdth*dxdth + gpp_BL * dxdp*dxdp;
  gup[2][2] = grr_BL * dydr*dydr + gthth_BL * dydth*dydth + gpp_BL * dydp*dydp;
  gup[3][3] = grr_BL * dzdr*dzdr + gthth_BL * dzdth*dzdth + gpp_BL * dzdp*dzdp;
  
  gup[1][2] = gup[2][1] = grr_BL * dxdr*dydr + gthth_BL * dxdth*dydth + gpp_BL * dxdp*dydp;
  gup[1][3] = gup[3][1] = grr_BL * dxdr*dzdr + gthth_BL * dxdth*dzdth + gpp_BL * dxdp*dzdp;
  gup[2][3] = gup[3][2] = grr_BL * dzdr*dydr + gthth_BL * dzdth*dydth + gpp_BL * dzdp*dydp;

}

void KerrHarmonic::jacobian(double jac[4][4][4], const double * pos) const {

  double aa = spin_, a2 = aa*aa;
  double mcal2 = 1-a2; // this is \mathcal{M}^2 = M^2 - a^2 in FV notes 
  // Harmonic Cartesian coordinates and harmonic radius
  double xx = pos[1], yy = pos[2], zz = pos[3],
    //xx = 10., yy=9., zz=8., // TEST
    x2 = xx*xx, y2 = yy*yy, z2 = zz*zz,
    rr = pow(x2+y2+z2,0.5), r2 = rr*rr;
  double Rcal2 = pow((r2 - a2)*(r2 - a2) + 4.*a2*z2,0.5),
    Rcal = pow(Rcal2,0.5), Rcal3 = Rcal*Rcal2;

  // BL coordinates derived from harmonic
  double rBL = 1. + 1/pow(2.,0.5)*pow(r2 - a2 + Rcal2,0.5), // BL radius
    rBL2 = rBL*rBL,
    uu = rBL-1., // u = rBL-M is a useful quantity, see FV notes
    u2 = uu*uu,
    u2_mcal2 = u2 - mcal2,
    costhBL = zz/uu, costhBL2 = costhBL*costhBL, // BL cos(theta)
    sinthBL = pow(1 - costhBL2,0.5), sinthBL2 = sinthBL*sinthBL,
    rho2BL = rBL2 + a2 * costhBL2, // usual rho2 quantity in BL metric
    DeltaBL = rBL2 - 2.*rBL + a2; // usual Delta quantity in BL metric

  // Covariant BL metric coefficients expressed in
  // harmonic coordinates
  double //gtt_BL = -(1. - 2*rBL/rho2BL),
    gtp_BL = -2*aa*rBL*sinthBL2/rho2BL,
    grr_BL = rho2BL/DeltaBL, 
    gthth_BL = rho2BL,
    gpp_BL = (rBL2 + a2 + 2.*a2*rBL*sinthBL2/rho2BL)*sinthBL2;

  // Derivative of BL coord wrt harmonic ones
  double
    drdx = uu*xx/Rcal2,
    drdy = uu*yy/Rcal2,
    drdz = (u2 + a2)*zz/(uu*Rcal2),
    dthdx = xx*zz/(Rcal2*pow(u2-z2,0.5)),
    dthdy = yy*zz/(Rcal2*pow(u2-z2,0.5)),
    dthdz = 1./pow(u2-z2,0.5) * (z2*(u2+a2)/(Rcal2*u2) - 1.),
    dphidx = aa/(a2+u2) * drdx - yy/(x2+y2) - drdx*aa/(u2 - mcal2), 
    dphidy = aa/(a2+u2) * drdy + xx/(x2+y2) - drdy*aa/(u2 - mcal2),
    dphidz = aa/(a2+u2) * drdz - drdz*aa/(u2 - mcal2);

  // Derivatives of Rcal
  double dRcaldx = xx*(r2-a2) / Rcal3,
    dRcaldy = yy*(r2-a2) / Rcal3,
    dRcaldz = zz * (r2+a2) / Rcal3; // it is indeed a plus here unlike the 2 previous ones

  // Second derivative of BL coord wrt harmonic ones // all checked
  double
    d2rdx2 = xx/Rcal2*drdx + uu/Rcal2 - 2.*uu*xx/Rcal3 * dRcaldx,
    d2rdy2 = yy/Rcal2*drdy + uu/Rcal2 - 2.*uu*yy/Rcal3 * dRcaldy,
    d2rdz2 = drdz * (1. - a2/u2) * zz/Rcal2 + (u2+a2)/(uu*Rcal2) - 2.*dRcaldz*zz/Rcal3 * (uu + a2/uu),
    d2rdydx = xx/yy*d2rdy2 - xx/y2 * drdy,
    d2rdzdx = xx/Rcal2 * drdz - 2.*uu*xx/Rcal3 * dRcaldz,
    d2rdzdy = yy/Rcal2 * drdz - 2.*uu*yy/Rcal3 * dRcaldz,
    d2thdx2 = zz/(Rcal2 * pow(u2-z2,0.5)) * (1. - xx * (2./Rcal * dRcaldx + uu/(u2-z2) * drdx)),
    d2thdy2 = zz/(Rcal2 * pow(u2-z2,0.5)) * (1. - yy * (2./Rcal * dRcaldy + uu/(u2-z2) * drdy)),
    d2thdz2 = 1./pow(u2-z2,0.5) * ((zz - uu*drdz)/(u2-z2) * (z2*(u2+a2)/(Rcal2*u2) - 1.)
				   + 2.*zz*(u2+a2)/(Rcal2*u2) + 2.*z2/(Rcal2*uu)*drdz
				   - 2.*z2*(u2+a2)/(Rcal3*u2)*dRcaldz - 2.*z2*(u2+a2)/(Rcal2*u2*uu)*drdz), // this one is really pretty, isn't it?
    d2thdydx = -xx*zz/(Rcal2*pow(u2-z2,0.5)) * (2./Rcal * dRcaldy + uu/(u2-z2) * drdy),
    d2thdzdx = xx/(Rcal2*pow(u2-z2,0.5)) * (1. - 2.*zz/Rcal * dRcaldz - zz/(u2-z2) * (drdz*uu - zz)),
    d2thdzdy = yy/(Rcal2*pow(u2-z2,0.5)) * (1. - 2.*zz/Rcal * dRcaldz - zz/(u2-z2) * (drdz*uu - zz)),
    d2phidx2 = aa/(a2+u2) * d2rdx2 - 2.*aa * uu / ((a2+u2)*(a2+u2)) * drdx*drdx + 2.*xx*yy/((x2+y2)*(x2+y2)) - d2rdx2*aa/u2_mcal2 + 2.*drdx*drdx*aa*uu/(u2_mcal2 * u2_mcal2),
    d2phidy2 = aa/(a2+u2) * d2rdy2 - 2.*aa * uu / ((a2+u2)*(a2+u2)) * drdy*drdy - 2.*xx*yy/((x2+y2)*(x2+y2)) - d2rdy2*aa/u2_mcal2 + 2.*drdy*drdy*aa*uu/(u2_mcal2 * u2_mcal2),
    d2phidz2 = aa/(a2+u2) * d2rdz2 - 2.*aa * uu / ((a2+u2)*(a2+u2)) * drdz*drdz - d2rdz2*aa/u2_mcal2 + 2.*drdz*drdz*aa*uu/(u2_mcal2 * u2_mcal2),
    d2phidydx = aa/(a2+u2) * d2rdydx - 2.*aa * uu / ((a2+u2)*(a2+u2)) * drdx*drdy - 1./(x2+y2) + 2.*y2/((x2+y2)*(x2+y2)) - d2rdydx*aa/u2_mcal2 + 2.*drdx*drdy*aa*uu/(u2_mcal2 * u2_mcal2),
    d2phidzdx = aa/(a2+u2) * d2rdzdx - 2.*aa * uu / ((a2+u2)*(a2+u2)) * drdz*drdx - d2rdzdx*aa/u2_mcal2 + 2.*drdx*drdz*aa*uu/(u2_mcal2 * u2_mcal2),
    d2phidzdy = aa/(a2+u2) * d2rdzdy - 2.*aa * uu / ((a2+u2)*(a2+u2)) * drdz*drdy - d2rdzdy*aa/u2_mcal2 + 2.*drdy*drdz*aa*uu/(u2_mcal2 * u2_mcal2);

  // Derivatives of BL metric coefs
  double drho2dx = 2.*(rBL*drdx - a2*costhBL*sinthBL*dthdx),
    drho2dy = 2.*(rBL*drdy - a2*costhBL*sinthBL*dthdy),
    drho2dz = 2.*(rBL*drdz - a2*costhBL*sinthBL*dthdz),
    dDeltadx = 2.*drdx*(rBL-1.),
    dDeltady = 2.*drdy*(rBL-1.),
    dDeltadz = 2.*drdz*(rBL-1.),
    drho2_over_Delta_dx = (drho2dx * DeltaBL - rho2BL * dDeltadx)/(DeltaBL*DeltaBL),
    drho2_over_Delta_dy = (drho2dy * DeltaBL - rho2BL * dDeltady)/(DeltaBL*DeltaBL),
    drho2_over_Delta_dz = (drho2dz * DeltaBL - rho2BL * dDeltadz)/(DeltaBL*DeltaBL),
    dr_over_rho2_dx = (drdx*rho2BL - rBL*drho2dx)/(rho2BL*rho2BL),
    dr_over_rho2_dy = (drdy*rho2BL - rBL*drho2dy)/(rho2BL*rho2BL),
    dr_over_rho2_dz = (drdz*rho2BL - rBL*drho2dz)/(rho2BL*rho2BL),
    dsinth2dx = 2.*costhBL*sinthBL*dthdx,
    dsinth2dy = 2.*costhBL*sinthBL*dthdy,
    dsinth2dz = 2.*costhBL*sinthBL*dthdz,
    dr_over_rho2_sin2th_dx = dr_over_rho2_dx * sinthBL2 + rBL/rho2BL * dsinth2dx,
    dr_over_rho2_sin2th_dy = dr_over_rho2_dy * sinthBL2 + rBL/rho2BL * dsinth2dy,
    dr_over_rho2_sin2th_dz = dr_over_rho2_dz * sinthBL2 + rBL/rho2BL * dsinth2dz;

  //cout << "droverrho2dx= " << dr_over_rho2_dx << " " << drho2dx << endl;
    
  double dgtt_dx = 2. * dr_over_rho2_dx,
    dgtt_dy = 2. * dr_over_rho2_dy,
    dgtt_dz = 2. * dr_over_rho2_dz,
    dgtp_dx = -2.*aa*dr_over_rho2_sin2th_dx,
    dgtp_dy = -2.*aa*dr_over_rho2_sin2th_dy,
    dgtp_dz = -2.*aa*dr_over_rho2_sin2th_dz,
    dgrr_dx = drho2_over_Delta_dx,
    dgrr_dy = drho2_over_Delta_dy,
    dgrr_dz = drho2_over_Delta_dz,
    dgthth_dx = drho2dx,
    dgthth_dy = drho2dy,
    dgthth_dz = drho2dz,
    dgpp_dx = (2.*rBL*drdx + 2.*a2*dr_over_rho2_sin2th_dx)*sinthBL2 + (rBL2 + a2 + 2.*a2*rBL*sinthBL2/rho2BL) * dsinth2dx,
    dgpp_dy = (2.*rBL*drdy + 2.*a2*dr_over_rho2_sin2th_dy)*sinthBL2 + (rBL2 + a2 + 2.*a2*rBL*sinthBL2/rho2BL) * dsinth2dy,
    dgpp_dz = (2.*rBL*drdz + 2.*a2*dr_over_rho2_sin2th_dz)*sinthBL2 + (rBL2 + a2 + 2.*a2*rBL*sinthBL2/rho2BL) * dsinth2dz;

  for (int mu=0; mu<4; ++mu){
    for (int nu=0; nu<4; ++nu){
      jac[0][mu][nu] = 0.; // \partial_t = 0
    }
  }

  // Formulas below are checked
  jac[1][0][0] = dgtt_dx ; // \partial_x g_tt
  jac[2][0][0] = dgtt_dy; // \partial_y g_tt
  jac[3][0][0] = dgtt_dz; // \partial_z g_tt

  jac[1][0][1] = jac[1][1][0] = dgtp_dx * dphidx + gtp_BL * d2phidx2;
  jac[2][0][1] = jac[2][1][0] = dgtp_dy * dphidx + gtp_BL * d2phidydx; 
  jac[3][0][1] = jac[3][1][0] = dgtp_dz * dphidx + gtp_BL * d2phidzdx;

  jac[1][0][2] = jac[1][2][0] = dgtp_dx * dphidy + gtp_BL * d2phidydx;
  jac[2][0][2] = jac[2][2][0] = dgtp_dy * dphidy + gtp_BL * d2phidy2; 
  jac[3][0][2] = jac[3][2][0] = dgtp_dz * dphidy + gtp_BL * d2phidzdy;

  jac[1][0][3] = jac[1][3][0] = dgtp_dx * dphidz + gtp_BL * d2phidzdx;
  jac[2][0][3] = jac[2][3][0] = dgtp_dy * dphidz + gtp_BL * d2phidzdy; 
  jac[3][0][3] = jac[3][3][0] = dgtp_dz * dphidz + gtp_BL * d2phidz2; 

  jac[1][1][1] = dgrr_dx * drdx * drdx + 2. * grr_BL * d2rdx2 * drdx
    + dgthth_dx * dthdx * dthdx + 2. * gthth_BL * d2thdx2 * dthdx
    + dgpp_dx * dphidx * dphidx + 2. * gpp_BL * d2phidx2 * dphidx;
  jac[2][1][1] = dgrr_dy * drdx * drdx + 2. * grr_BL * d2rdydx * drdx
    + dgthth_dy * dthdx * dthdx + 2. * gthth_BL * d2thdydx * dthdx
    + dgpp_dy * dphidx * dphidx + 2. * gpp_BL * d2phidydx * dphidx;
  jac[3][1][1] = dgrr_dz * drdx * drdx + 2. * grr_BL * d2rdzdx * drdx
    + dgthth_dz * dthdx * dthdx + 2. * gthth_BL * d2thdzdx * dthdx
    + dgpp_dz * dphidx * dphidx + 2. * gpp_BL * d2phidzdx * dphidx;
  jac[1][2][2] = dgrr_dx * drdy * drdy + 2. * grr_BL * d2rdydx * drdy
    + dgthth_dx * dthdy * dthdy + 2. * gthth_BL * d2thdydx * dthdy
    + dgpp_dx * dphidy * dphidy + 2. * gpp_BL * d2phidydx * dphidy;
  jac[2][2][2] = dgrr_dy * drdy * drdy + 2. * grr_BL * d2rdy2 * drdy
    + dgthth_dy * dthdy * dthdy + 2. * gthth_BL * d2thdy2 * dthdy
    + dgpp_dy * dphidy * dphidy + 2. * gpp_BL * d2phidy2 * dphidy;
  jac[3][2][2] = dgrr_dz * drdy * drdy + 2. * grr_BL * d2rdzdy * drdy
    + dgthth_dz * dthdy * dthdy + 2. * gthth_BL * d2thdzdy * dthdy
    + dgpp_dz * dphidy * dphidy + 2. * gpp_BL * d2phidzdy * dphidy;
  jac[1][3][3] = dgrr_dx * drdz * drdz + 2. * grr_BL * d2rdzdx * drdz
    + dgthth_dx * dthdz * dthdz + 2. * gthth_BL * d2thdzdx * dthdz
    + dgpp_dx * dphidz * dphidz + 2. * gpp_BL * d2phidzdx * dphidz;
  jac[2][3][3] = dgrr_dy * drdz * drdz + 2. * grr_BL * d2rdzdy * drdz
    + dgthth_dy * dthdz * dthdz + 2. * gthth_BL * d2thdzdy * dthdz
    + dgpp_dy * dphidz * dphidz + 2. * gpp_BL * d2phidzdy * dphidz;
  jac[3][3][3] = dgrr_dz * drdz * drdz + 2. * grr_BL * d2rdz2 * drdz
    + dgthth_dz * dthdz * dthdz + 2. * gthth_BL * d2thdz2 * dthdz
    + dgpp_dz * dphidz * dphidz + 2. * gpp_BL * d2phidz2 * dphidz;

  jac[1][1][2] = jac[1][2][1] = dgrr_dx * drdx * drdy
    + grr_BL * d2rdx2 * drdy + grr_BL * drdx * d2rdydx
    + dgthth_dx * dthdx * dthdy
    + gthth_BL * d2thdx2 * dthdy + gthth_BL * dthdx * d2thdydx
    + dgpp_dx * dphidx * dphidy
    + gpp_BL * d2phidx2 * dphidy + gpp_BL * dphidx * d2phidydx;
  jac[2][1][2] = jac[2][2][1] = dgrr_dy * drdx * drdy
    + grr_BL * d2rdydx * drdy + grr_BL * drdx * d2rdy2
    + dgthth_dy * dthdx * dthdy
    + gthth_BL * d2thdydx * dthdy + gthth_BL * dthdx * d2thdy2
    + dgpp_dy * dphidx * dphidy
    + gpp_BL * d2phidydx * dphidy + gpp_BL * dphidx * d2phidy2;
  jac[3][1][2] = jac[3][2][1] = dgrr_dz * drdx * drdy
    + grr_BL * d2rdzdx * drdy + grr_BL * drdx * d2rdzdy
    + dgthth_dz * dthdx * dthdy
    + gthth_BL * d2thdzdx * dthdy + gthth_BL * dthdx * d2thdzdy
    + dgpp_dz * dphidx * dphidy
    + gpp_BL * d2phidzdx * dphidy + gpp_BL * dphidx * d2phidzdy;
  
  jac[1][1][3] = jac[1][3][1] = dgrr_dx * drdx * drdz 
    + grr_BL * d2rdx2 * drdz + grr_BL * drdx * d2rdzdx
    + dgthth_dx * dthdx * dthdz
    + gthth_BL * d2thdx2 * dthdz + gthth_BL * dthdx * d2thdzdx
    + dgpp_dx * dphidx * dphidz
    + gpp_BL * d2phidx2 * dphidz + gpp_BL * dphidx * d2phidzdx;
  jac[2][1][3] = jac[2][3][1] = dgrr_dy * drdx * drdz
    + grr_BL * d2rdydx * drdz + grr_BL * drdx * d2rdzdy
    + dgthth_dy * dthdx * dthdz
    + gthth_BL * d2thdydx * dthdz + gthth_BL * dthdx * d2thdzdy
    + dgpp_dy * dphidx * dphidz
    + gpp_BL * d2phidydx * dphidz + gpp_BL * dphidx * d2phidzdy;
  jac[3][1][3] = jac[3][3][1] = dgrr_dz * drdx * drdz 
    + grr_BL * d2rdzdx * drdz + grr_BL * drdx * d2rdz2
    + dgthth_dz * dthdx * dthdz
    + gthth_BL * d2thdzdx * dthdz + gthth_BL * dthdx * d2thdz2
    + dgpp_dz * dphidx * dphidz
    + gpp_BL * d2phidzdx * dphidz + gpp_BL * dphidx * d2phidz2;

  jac[1][2][3] = jac[1][3][2] = dgrr_dx * drdy * drdz
    + grr_BL * d2rdydx * drdz + grr_BL * drdy * d2rdzdx
    + dgthth_dx * dthdy * dthdz
    + gthth_BL * d2thdydx * dthdz + gthth_BL * dthdy * d2thdzdx
    + dgpp_dx * dphidy * dphidz
    + gpp_BL * d2phidydx * dphidz + gpp_BL * dphidy * d2phidzdx;
  jac[2][2][3] = jac[2][3][2] = dgrr_dy * drdy * drdz
    + grr_BL * d2rdy2 * drdz + grr_BL * drdy * d2rdzdy
    + dgthth_dy * dthdy * dthdz
    + gthth_BL * d2thdy2 * dthdz + gthth_BL * dthdy * d2thdzdy
    + dgpp_dy * dphidy * dphidz
    + gpp_BL * d2phidy2 * dphidz + gpp_BL * dphidy * d2phidzdy;
  jac[3][2][3] = jac[3][3][2] = dgrr_dz * drdy * drdz
    + grr_BL * d2rdzdy * drdz + grr_BL * drdy * d2rdz2
    + dgthth_dz * dthdy * dthdz
    + gthth_BL * d2thdzdy * dthdz + gthth_BL * dthdy * d2thdz2
    + dgpp_dz * dphidy * dphidz
    + gpp_BL * d2phidzdy * dphidz + gpp_BL * dphidy * d2phidz2;

}


void KerrHarmonic::gmunu_up_and_jacobian(double gup[4][4], double jac[4][4][4], const double * pos) const {

  // DEFINE ALL QUANTITIES
  double aa = spin_, a2 = aa*aa;
  double mcal2 = 1-a2; // this is \mathcal{M}^2 = M^2 - a^2 in FV notes 
  // Harmonic Cartesian coordinates and harmonic radius
  double xx = pos[1], yy = pos[2], zz = pos[3],
    //xx = 10., yy=9., zz=8., // TEST
    x2 = xx*xx, y2 = yy*yy, z2 = zz*zz,
    rr = pow(x2+y2+z2,0.5), r2 = rr*rr;
  double Rcal2 = pow((r2 - a2)*(r2 - a2) + 4.*a2*z2,0.5),
    Rcal = pow(Rcal2,0.5), Rcal3 = Rcal*Rcal2;

  // BL coordinates derived from harmonic
  double rBL = 1. + 1/pow(2.,0.5)*pow(r2 - a2 + Rcal2,0.5), // BL radius
    rBL2 = rBL*rBL,
    uu = rBL-1., // u = rBL-M is a useful quantity, see FV notes
    u2 = uu*uu,
    u2_mcal2 = u2 - mcal2,
    costhBL = zz/uu, costhBL2 = costhBL*costhBL, // BL cos(theta)
    sinthBL = pow(1 - costhBL2,0.5), sinthBL2 = sinthBL*sinthBL,
    rho2BL = rBL2 + a2 * costhBL2, // usual rho2 quantity in BL metric
    DeltaBL = rBL2 - 2.*rBL + a2, // usual Delta quantity in BL metric
    cosphBL = (uu*xx + aa*yy)/(x2+y2) * sinthBL, // playing with the x ans y relations
    sinphBL = (uu*yy - aa*xx)/(x2+y2) * sinthBL;

  // Covariant BL metric coefficients expressed in
  // harmonic coordinates
  double //gtt_BL = -(1. - 2*rBL/rho2BL),
    gtp_BL = -2*aa*rBL*sinthBL2/rho2BL,
    grr_BL = rho2BL/DeltaBL, 
    gthth_BL = rho2BL,
    gpp_BL = (rBL2 + a2 + 2.*a2*rBL*sinthBL2/rho2BL)*sinthBL2;

  // Contravariant (inverse) BL metric expressed in
  // harmonic coordinates
  double gtt_inv_BL = -1./DeltaBL * (rBL*rBL + a2 + 2.*a2*rBL*sinthBL2/rho2BL), // contravariant indices throughout
    grr_inv_BL = DeltaBL/rho2BL,
    gthth_inv_BL = 1./rho2BL,
    gpp_inv_BL = 1./(DeltaBL*sinthBL2) * (1. - 2.*rBL/rho2BL),
    gtp_inv_BL = -2.*aa*rBL/(DeltaBL*rho2BL);

  // Derivative of BL coord wrt harmonic ones
  double
    drdx = uu*xx/Rcal2,
    drdy = uu*yy/Rcal2,
    drdz = (u2 + a2)*zz/(uu*Rcal2),
    dthdx = xx*zz/(Rcal2*pow(u2-z2,0.5)),
    dthdy = yy*zz/(Rcal2*pow(u2-z2,0.5)),
    dthdz = 1./pow(u2-z2,0.5) * (z2*(u2+a2)/(Rcal2*u2) - 1.),
    dphidx = aa/(a2+u2) * drdx - yy/(x2+y2) - drdx*aa/(u2 - mcal2), 
    dphidy = aa/(a2+u2) * drdy + xx/(x2+y2) - drdy*aa/(u2 - mcal2),
    dphidz = aa/(a2+u2) * drdz - drdz*aa/(u2 - mcal2);

  // Derivative of harmonic coord wrt BL ones
  double
    dxdr = cosphBL * sinthBL - aa*yy/(u2-mcal2),
    dydr = sinphBL * sinthBL + aa*xx/(u2-mcal2),
    dzdr = costhBL,
    dxdth = xx*costhBL/sinthBL,
    dydth = yy*costhBL/sinthBL,
    dzdth = -uu*sinthBL,
    dxdp = -yy,
    dydp = xx,
    dzdp = 0.;

  // Derivatives of Rcal
  double dRcaldx = xx*(r2-a2) / Rcal3,
    dRcaldy = yy*(r2-a2) / Rcal3,
    dRcaldz = zz * (r2+a2) / Rcal3; // it is indeed a plus here unlike the 2 previous ones

  // Second derivative of BL coord wrt harmonic ones
  double
    d2rdx2 = xx/Rcal2*drdx + uu/Rcal2 - 2.*uu*xx/Rcal3 * dRcaldx,
    d2rdy2 = yy/Rcal2*drdy + uu/Rcal2 - 2.*uu*yy/Rcal3 * dRcaldy,
    d2rdz2 = drdz * (1. - a2/u2) * zz/Rcal2 + (u2+a2)/(uu*Rcal2) - 2.*dRcaldz*zz/Rcal3 * (uu + a2/uu),
    d2rdydx = xx/yy*d2rdy2 - xx/y2 * drdy,
    d2rdzdx = xx/Rcal2 * drdz - 2.*uu*xx/Rcal3 * dRcaldz,
    d2rdzdy = yy/Rcal2 * drdz - 2.*uu*yy/Rcal3 * dRcaldz,
    d2thdx2 = zz/(Rcal2 * pow(u2-z2,0.5)) * (1. - xx * (2./Rcal * dRcaldx + uu/(u2-z2) * drdx)),
    d2thdy2 = zz/(Rcal2 * pow(u2-z2,0.5)) * (1. - yy * (2./Rcal * dRcaldy + uu/(u2-z2) * drdy)),
    d2thdz2 = 1./pow(u2-z2,0.5) * ((zz - uu*drdz)/(u2-z2) * (z2*(u2+a2)/(Rcal2*u2) - 1.)
				   + 2.*zz*(u2+a2)/(Rcal2*u2) + 2.*z2/(Rcal2*uu)*drdz
				   - 2.*z2*(u2+a2)/(Rcal3*u2)*dRcaldz - 2.*z2*(u2+a2)/(Rcal2*u2*uu)*drdz), // this one is really pretty, isn't it?
    d2thdydx = -xx*zz/(Rcal2*pow(u2-z2,0.5)) * (2./Rcal * dRcaldy + uu/(u2-z2) * drdy),
    d2thdzdx = xx/(Rcal2*pow(u2-z2,0.5)) * (1. - 2.*zz/Rcal * dRcaldz - zz/(u2-z2) * (drdz*uu - zz)),
    d2thdzdy = yy/(Rcal2*pow(u2-z2,0.5)) * (1. - 2.*zz/Rcal * dRcaldz - zz/(u2-z2) * (drdz*uu - zz)),
    d2phidx2 = aa/(a2+u2) * d2rdx2 - 2.*aa * uu / ((a2+u2)*(a2+u2)) * drdx*drdx + 2.*xx*yy/((x2+y2)*(x2+y2)) - d2rdx2*aa/u2_mcal2 + 2.*drdx*drdx*aa*uu/(u2_mcal2 * u2_mcal2),
    d2phidy2 = aa/(a2+u2) * d2rdy2 - 2.*aa * uu / ((a2+u2)*(a2+u2)) * drdy*drdy - 2.*xx*yy/((x2+y2)*(x2+y2)) - d2rdy2*aa/u2_mcal2 + 2.*drdy*drdy*aa*uu/(u2_mcal2 * u2_mcal2),
    d2phidz2 = aa/(a2+u2) * d2rdz2 - 2.*aa * uu / ((a2+u2)*(a2+u2)) * drdz*drdz - d2rdz2*aa/u2_mcal2 + 2.*drdz*drdz*aa*uu/(u2_mcal2 * u2_mcal2),
    d2phidydx = aa/(a2+u2) * d2rdydx - 2.*aa * uu / ((a2+u2)*(a2+u2)) * drdx*drdy - 1./(x2+y2) + 2.*y2/((x2+y2)*(x2+y2)) - d2rdydx*aa/u2_mcal2 + 2.*drdx*drdy*aa*uu/(u2_mcal2 * u2_mcal2),
    d2phidzdx = aa/(a2+u2) * d2rdzdx - 2.*aa * uu / ((a2+u2)*(a2+u2)) * drdz*drdx - d2rdzdx*aa/u2_mcal2 + 2.*drdx*drdz*aa*uu/(u2_mcal2 * u2_mcal2),
    d2phidzdy = aa/(a2+u2) * d2rdzdy - 2.*aa * uu / ((a2+u2)*(a2+u2)) * drdz*drdy - d2rdzdy*aa/u2_mcal2 + 2.*drdy*drdz*aa*uu/(u2_mcal2 * u2_mcal2);

  // Derivatives of BL metric coefs
  double drho2dx = 2.*(rBL*drdx - a2*costhBL*sinthBL*dthdx),
    drho2dy = 2.*(rBL*drdy - a2*costhBL*sinthBL*dthdy),
    drho2dz = 2.*(rBL*drdz - a2*costhBL*sinthBL*dthdz),
    dDeltadx = 2.*drdx*(rBL-1.),
    dDeltady = 2.*drdy*(rBL-1.),
    dDeltadz = 2.*drdz*(rBL-1.),
    drho2_over_Delta_dx = (drho2dx * DeltaBL - rho2BL * dDeltadx)/(DeltaBL*DeltaBL),
    drho2_over_Delta_dy = (drho2dy * DeltaBL - rho2BL * dDeltady)/(DeltaBL*DeltaBL),
    drho2_over_Delta_dz = (drho2dz * DeltaBL - rho2BL * dDeltadz)/(DeltaBL*DeltaBL),
    dr_over_rho2_dx = (drdx*rho2BL - rBL*drho2dx)/(rho2BL*rho2BL),
    dr_over_rho2_dy = (drdy*rho2BL - rBL*drho2dy)/(rho2BL*rho2BL),
    dr_over_rho2_dz = (drdz*rho2BL - rBL*drho2dz)/(rho2BL*rho2BL),
    dsinth2dx = 2.*costhBL*sinthBL*dthdx,
    dsinth2dy = 2.*costhBL*sinthBL*dthdy,
    dsinth2dz = 2.*costhBL*sinthBL*dthdz,
    dr_over_rho2_sin2th_dx = dr_over_rho2_dx * sinthBL2 + rBL/rho2BL * dsinth2dx,
    dr_over_rho2_sin2th_dy = dr_over_rho2_dy * sinthBL2 + rBL/rho2BL * dsinth2dy,
    dr_over_rho2_sin2th_dz = dr_over_rho2_dz * sinthBL2 + rBL/rho2BL * dsinth2dz;
    
  double dgtt_dx = 2. * dr_over_rho2_dx,
    dgtt_dy = 2. * dr_over_rho2_dy,
    dgtt_dz = 2. * dr_over_rho2_dz,
    dgtp_dx = -2.*aa*dr_over_rho2_sin2th_dx,
    dgtp_dy = -2.*aa*dr_over_rho2_sin2th_dy,
    dgtp_dz = -2.*aa*dr_over_rho2_sin2th_dz,
    dgrr_dx = drho2_over_Delta_dx,
    dgrr_dy = drho2_over_Delta_dy,
    dgrr_dz = drho2_over_Delta_dz,
    dgthth_dx = drho2dx,
    dgthth_dy = drho2dy,
    dgthth_dz = drho2dz,
    dgpp_dx = (2.*rBL*drdx + 2.*a2*dr_over_rho2_sin2th_dx)*sinthBL2 + (rBL2 + a2 + 2.*a2*rBL*sinthBL2/rho2BL) * dsinth2dx,
    dgpp_dy = (2.*rBL*drdy + 2.*a2*dr_over_rho2_sin2th_dy)*sinthBL2 + (rBL2 + a2 + 2.*a2*rBL*sinthBL2/rho2BL) * dsinth2dy,
    dgpp_dz = (2.*rBL*drdz + 2.*a2*dr_over_rho2_sin2th_dz)*sinthBL2 + (rBL2 + a2 + 2.*a2*rBL*sinthBL2/rho2BL) * dsinth2dz;
  
  // FILL GMUNU_UP

  gup[0][0] = gtt_inv_BL; // g^tt
  
  gup[0][1] = gup[1][0] = gtp_inv_BL*dxdp; // g^tx
  gup[0][2] = gup[2][0] = gtp_inv_BL*dydp; // g^ty
  gup[0][3] = gup[3][0] = gtp_inv_BL*dzdp; // g^tz
  
  gup[1][1] = grr_inv_BL * dxdr*dxdr + gthth_inv_BL * dxdth*dxdth
    + gpp_inv_BL * dxdp*dxdp;
  gup[2][2] = grr_inv_BL * dydr*dydr + gthth_inv_BL * dydth*dydth
    + gpp_inv_BL * dydp*dydp;
  gup[3][3] = grr_inv_BL * dzdr*dzdr + gthth_inv_BL * dzdth*dzdth
    + gpp_inv_BL * dzdp*dzdp;
  
  gup[1][2] = gup[2][1] = grr_inv_BL * dxdr*dydr + gthth_inv_BL * dxdth*dydth
    + gpp_inv_BL * dxdp*dydp;
  gup[1][3] = gup[3][1] = grr_inv_BL * dxdr*dzdr + gthth_inv_BL * dxdth*dzdth
    + gpp_inv_BL * dxdp*dzdp;
  gup[2][3] = gup[3][2] = grr_inv_BL * dzdr*dydr + gthth_inv_BL * dzdth*dydth
    + gpp_inv_BL * dzdp*dydp;

  // FILL JACOBIAN

    for (int mu=0; mu<4; ++mu){
    for (int nu=0; nu<4; ++nu){
      jac[0][mu][nu] = 0.; // \partial_t = 0
    }
  }

  // Formulas below are checked
  jac[1][0][0] = dgtt_dx ; // \partial_x g_tt
  jac[2][0][0] = dgtt_dy; // \partial_y g_tt
  jac[3][0][0] = dgtt_dz; // \partial_z g_tt

  jac[1][0][1] = jac[1][1][0] = dgtp_dx * dphidx + gtp_BL * d2phidx2;
  jac[2][0][1] = jac[2][1][0] = dgtp_dy * dphidx + gtp_BL * d2phidydx; 
  jac[3][0][1] = jac[3][1][0] = dgtp_dz * dphidx + gtp_BL * d2phidzdx;

  jac[1][0][2] = jac[1][2][0] = dgtp_dx * dphidy + gtp_BL * d2phidydx;
  jac[2][0][2] = jac[2][2][0] = dgtp_dy * dphidy + gtp_BL * d2phidy2; 
  jac[3][0][2] = jac[3][2][0] = dgtp_dz * dphidy + gtp_BL * d2phidzdy;

  jac[1][0][3] = jac[1][3][0] = dgtp_dx * dphidz + gtp_BL * d2phidzdx;
  jac[2][0][3] = jac[2][3][0] = dgtp_dy * dphidz + gtp_BL * d2phidzdy; 
  jac[3][0][3] = jac[3][3][0] = dgtp_dz * dphidz + gtp_BL * d2phidz2; 

  jac[1][1][1] = dgrr_dx * drdx * drdx + 2. * grr_BL * d2rdx2 * drdx
    + dgthth_dx * dthdx * dthdx + 2. * gthth_BL * d2thdx2 * dthdx
    + dgpp_dx * dphidx * dphidx + 2. * gpp_BL * d2phidx2 * dphidx;
  jac[2][1][1] = dgrr_dy * drdx * drdx + 2. * grr_BL * d2rdydx * drdx
    + dgthth_dy * dthdx * dthdx + 2. * gthth_BL * d2thdydx * dthdx
    + dgpp_dy * dphidx * dphidx + 2. * gpp_BL * d2phidydx * dphidx;
  jac[3][1][1] = dgrr_dz * drdx * drdx + 2. * grr_BL * d2rdzdx * drdx
    + dgthth_dz * dthdx * dthdx + 2. * gthth_BL * d2thdzdx * dthdx
    + dgpp_dz * dphidx * dphidx + 2. * gpp_BL * d2phidzdx * dphidx;
  jac[1][2][2] = dgrr_dx * drdy * drdy + 2. * grr_BL * d2rdydx * drdy
    + dgthth_dx * dthdy * dthdy + 2. * gthth_BL * d2thdydx * dthdy
    + dgpp_dx * dphidy * dphidy + 2. * gpp_BL * d2phidydx * dphidy;
  jac[2][2][2] = dgrr_dy * drdy * drdy + 2. * grr_BL * d2rdy2 * drdy
    + dgthth_dy * dthdy * dthdy + 2. * gthth_BL * d2thdy2 * dthdy
    + dgpp_dy * dphidy * dphidy + 2. * gpp_BL * d2phidy2 * dphidy;
  jac[3][2][2] = dgrr_dz * drdy * drdy + 2. * grr_BL * d2rdzdy * drdy
    + dgthth_dz * dthdy * dthdy + 2. * gthth_BL * d2thdzdy * dthdy
    + dgpp_dz * dphidy * dphidy + 2. * gpp_BL * d2phidzdy * dphidy;
  jac[1][3][3] = dgrr_dx * drdz * drdz + 2. * grr_BL * d2rdzdx * drdz
    + dgthth_dx * dthdz * dthdz + 2. * gthth_BL * d2thdzdx * dthdz
    + dgpp_dx * dphidz * dphidz + 2. * gpp_BL * d2phidzdx * dphidz;
  jac[2][3][3] = dgrr_dy * drdz * drdz + 2. * grr_BL * d2rdzdy * drdz
    + dgthth_dy * dthdz * dthdz + 2. * gthth_BL * d2thdzdy * dthdz
    + dgpp_dy * dphidz * dphidz + 2. * gpp_BL * d2phidzdy * dphidz;
  jac[3][3][3] = dgrr_dz * drdz * drdz + 2. * grr_BL * d2rdz2 * drdz
    + dgthth_dz * dthdz * dthdz + 2. * gthth_BL * d2thdz2 * dthdz
    + dgpp_dz * dphidz * dphidz + 2. * gpp_BL * d2phidz2 * dphidz;

  jac[1][1][2] = jac[1][2][1] = dgrr_dx * drdx * drdy
    + grr_BL * d2rdx2 * drdy + grr_BL * drdx * d2rdydx
    + dgthth_dx * dthdx * dthdy
    + gthth_BL * d2thdx2 * dthdy + gthth_BL * dthdx * d2thdydx
    + dgpp_dx * dphidx * dphidy
    + gpp_BL * d2phidx2 * dphidy + gpp_BL * dphidx * d2phidydx;
  jac[2][1][2] = jac[2][2][1] = dgrr_dy * drdx * drdy
    + grr_BL * d2rdydx * drdy + grr_BL * drdx * d2rdy2
    + dgthth_dy * dthdx * dthdy
    + gthth_BL * d2thdydx * dthdy + gthth_BL * dthdx * d2thdy2
    + dgpp_dy * dphidx * dphidy
    + gpp_BL * d2phidydx * dphidy + gpp_BL * dphidx * d2phidy2;
  jac[3][1][2] = jac[3][2][1] = dgrr_dz * drdx * drdy
    + grr_BL * d2rdzdx * drdy + grr_BL * drdx * d2rdzdy
    + dgthth_dz * dthdx * dthdy
    + gthth_BL * d2thdzdx * dthdy + gthth_BL * dthdx * d2thdzdy
    + dgpp_dz * dphidx * dphidy
    + gpp_BL * d2phidzdx * dphidy + gpp_BL * dphidx * d2phidzdy;
  
  jac[1][1][3] = jac[1][3][1] = dgrr_dx * drdx * drdz 
    + grr_BL * d2rdx2 * drdz + grr_BL * drdx * d2rdzdx
    + dgthth_dx * dthdx * dthdz
    + gthth_BL * d2thdx2 * dthdz + gthth_BL * dthdx * d2thdzdx
    + dgpp_dx * dphidx * dphidz
    + gpp_BL * d2phidx2 * dphidz + gpp_BL * dphidx * d2phidzdx;
  jac[2][1][3] = jac[2][3][1] = dgrr_dy * drdx * drdz
    + grr_BL * d2rdydx * drdz + grr_BL * drdx * d2rdzdy
    + dgthth_dy * dthdx * dthdz
    + gthth_BL * d2thdydx * dthdz + gthth_BL * dthdx * d2thdzdy
    + dgpp_dy * dphidx * dphidz
    + gpp_BL * d2phidydx * dphidz + gpp_BL * dphidx * d2phidzdy;
  jac[3][1][3] = jac[3][3][1] = dgrr_dz * drdx * drdz 
    + grr_BL * d2rdzdx * drdz + grr_BL * drdx * d2rdz2
    + dgthth_dz * dthdx * dthdz
    + gthth_BL * d2thdzdx * dthdz + gthth_BL * dthdx * d2thdz2
    + dgpp_dz * dphidx * dphidz
    + gpp_BL * d2phidzdx * dphidz + gpp_BL * dphidx * d2phidz2;

  jac[1][2][3] = jac[1][3][2] = dgrr_dx * drdy * drdz
    + grr_BL * d2rdydx * drdz + grr_BL * drdy * d2rdzdx
    + dgthth_dx * dthdy * dthdz
    + gthth_BL * d2thdydx * dthdz + gthth_BL * dthdy * d2thdzdx
    + dgpp_dx * dphidy * dphidz
    + gpp_BL * d2phidydx * dphidz + gpp_BL * dphidy * d2phidzdx;
  jac[2][2][3] = jac[2][3][2] = dgrr_dy * drdy * drdz
    + grr_BL * d2rdy2 * drdz + grr_BL * drdy * d2rdzdy
    + dgthth_dy * dthdy * dthdz
    + gthth_BL * d2thdy2 * dthdz + gthth_BL * dthdy * d2thdzdy
    + dgpp_dy * dphidy * dphidz
    + gpp_BL * d2phidy2 * dphidz + gpp_BL * dphidy * d2phidzdy;
  jac[3][2][3] = jac[3][3][2] = dgrr_dz * drdy * drdz
    + grr_BL * d2rdzdy * drdz + grr_BL * drdy * d2rdz2
    + dgthth_dz * dthdy * dthdz
    + gthth_BL * d2thdzdy * dthdz + gthth_BL * dthdy * d2thdz2
    + dgpp_dz * dphidy * dphidz
    + gpp_BL * d2phidzdy * dphidz + gpp_BL * dphidy * d2phidz2;
}
  
int KerrHarmonic::isStopCondition(double const * const coord) const {
  double rsinkBL = 0.;
  if (spin_*spin_<1.){ // Black hole solution with event horizon
    double rhorBL = 1 + sqrt(1 - spin_*spin_); // BL radius EH
    rsinkBL = rhorBL + GYOTO_KERR_HORIZON_SECURITY;
  }

  double xx=coord[1], yy=coord[2], zz=coord[3],
    rr = pow(xx*xx+yy*yy+zz*zz,0.5), r2 = rr*rr,
    a2 = spin_*spin_,
    Rcal2 = pow((r2 - a2)*(r2 - a2) + 4.*a2*zz*zz,0.5),
    rBL = 1. + 1/pow(2.,0.5)*pow(r2 - a2 + Rcal2,0.5); // current BL radius

  return rBL < rsinkBL ;
}

void KerrHarmonic::circularVelocity(double const * coor, double* vel,
				    double dir) const {

  double xx=coor[1], yy=coor[2], zz=coor[3],
    rr = pow(xx*xx+yy*yy+zz*zz,0.5), r2 = rr*rr,
    a2 = spin_*spin_,
    Rcal2 = pow((r2 - a2)*(r2 - a2) + 4.*a2*zz*zz,0.5),
    rBL = 1. + 1/pow(2.,0.5)*pow(r2 - a2 + Rcal2,0.5);

  double Omega = dir/(pow(rBL,3./2.) + spin_); // BL dphi/dt

  vel[1] = -coor[2]*Omega;
  vel[2] =  coor[1]*Omega;
  vel[3] = 0.;
  vel[0] = SysPrimeToTdot(coor, vel+1);
  vel[1] *= vel[0];
  vel[2] *= vel[0];
}

