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

double KerrHarmonic::gmunu(const double * pos, int mu, int nu) const {
  /*
   The details of the computations of the metric coefficients
   are provided in a pdf note by FV (ask me if interested).
   */
  double aa = spin_, a2 = aa*aa;
  // Harmonic Cartesian coordinates and harmonic radius
  double xx = pos[1], yy = pos[2], zz = pos[3],
    //double xx = 10., yy = 9., zz = 8., // test!!!!
    x2 = xx*xx, y2 = yy*yy, z2 = zz*zz,
    rr = pow(x2+y2+z2,0.5), r2 = rr*rr;
  double Rcal2 = pow((r2 - a2)*(r2 - a2) + 4.*a2*z2,0.5); // useful quantity,
                                            // this is \mathcal{R} in FV notes.

  // KS coordinates derived from harmonic
  double rKS = 1. + 1/pow(2.,0.5)*pow(r2 - a2 + Rcal2,0.5), // KS radius
    uu = rKS-1., // u = rKS-M is a useful quantity, see FV notes
    u2 = uu*uu,
    costhKS = zz/uu, costhKS2 = costhKS*costhKS, // KS cos(theta)
    sinthKS = pow(1 - costhKS2,0.5), sinthKS2 = sinthKS*sinthKS,
    rho2KS = rKS*rKS + a2 * costhKS2; // usual rho2 quantity in KS metric

  //cout << "All stuff H coord in KerrH: " << xx << " " << yy << " " << zz << " " << rr << " " << Rcal2 << " " << rKS << " " << costhKS << " " << sinthKS << endl;

  // Covariant spherical KS metric coefficients expressed in
  // harmonic coordinates
  double gtt_KS = -(1. - 2*rKS/rho2KS),
    gtr_KS = 2*rKS/rho2KS, gtp_KS = -2*aa*rKS*sinthKS2/rho2KS,
    grr_KS = 1.+2*rKS/rho2KS, grp_KS = -aa*(1.+2*rKS/rho2KS)*sinthKS2,
    gthth_KS = rho2KS,
    gpp_KS = (rKS*rKS + a2 + 2.*a2*rKS*sinthKS2/rho2KS)*sinthKS2;

  // Derivative of KS coord wrt harmonic ones
  double drdx = uu*xx/Rcal2, drdy = uu*yy/Rcal2,
    drdz = (u2 + a2)*zz/(uu*Rcal2),//2.*(u2 + a2)*zz/(uu*Rcal2),
    dphidz = aa/(a2+u2) * drdz,
    dthdx = 0.,
    dthdy = 0.,
    dthdz = 0.,
    dphidx = 0., 
    dphidy = 0.;
  /*
    Special cases along the z axis.

    1. The general formulas for dthdx^i (see below)
    are undef along the axis. See FV notes for a simple demo
    that dthdx and dthdy are equal to 1/z along the axis,
    while dthdz=0 there. These dthdx^i multiply nonzero gmunu
    coef so they must be defined.
    
    2. dphidx and dphidy are undef along the z axis (dphidz is fine).
    However, they multiply gmunu coef that are \propro \sin \theta,
    and that thus cancel along the axis. So we can put them by
    default to zero, and update them to their value below *only*
    away from the z axis, where the denominator of the second term
    is well defined.
  */
  if (x2+y2>0.){
    dthdx = xx*zz/(Rcal2*pow(u2-z2,0.5)); // deno is zero along axis
    dthdy = yy*zz/(Rcal2*pow(u2-z2,0.5));
    dthdz = 1./pow(u2-z2,0.5) * (z2*(u2+a2)/(Rcal2*u2) - 1.),
      
    dphidx = aa/(a2+u2) * drdx - yy/(x2+y2); // second deno zero along axis
    dphidy = aa/(a2+u2) * drdy + xx/(x2+y2);
  }else{
    dthdx = 1./zz; // specific computation along the axis
    dthdy = 1./zz;
    // dthdz is left to its default value of zero.
  }

  //cout << "dthdx= " << xx*zz << " " << Rcal2 << " " << u2 << " " <<  pow(u2-z2,0.5) << endl;

  //cout << "All stuff gmunuKS in KerrH: " << gtt_KS << " " << gtr_KS << " " << gtp_KS << " " << grr_KS << " " << grp_KS << " " << gthth_KS << " " << gpp_KS << endl;
  //cout << "All stuff d/dx in KerrH: " << drdx << " " << drdy << " " << drdz << " " << dthdx << " " << dthdy << " " << dthdz << " " << dphidx << " " << dphidy << " " << dphidz << endl;

  /*cout << "all t gmunu= " << gtt_KS << " " << gtr_KS * drdx + gtp_KS * dphidx << " " << gtr_KS * drdy + gtp_KS * dphidy << " " << gtr_KS * drdz + gtp_KS * dphidz << " " << endl;

  cout << "all space gmunu= " << grr_KS * drdx * drdx + gthth_KS * dthdx * dthdx
    + gpp_KS * dphidx * dphidx + 2. * grp_KS * drdx * dphidx << " " << grr_KS * drdy * drdy + gthth_KS * dthdy * dthdy
    + gpp_KS * dphidy * dphidy + 2. * grp_KS * drdy * dphidy << " " <<
     grr_KS * drdz * drdz + gthth_KS * dthdz * dthdz
    + gpp_KS * dphidz * dphidz + 2. * grp_KS * drdz * dphidz << " " <<
    grr_KS * drdx * drdy + gthth_KS * dthdx * dthdy
      + gpp_KS * dphidx * dphidy + grp_KS * (drdx * dphidy + drdy * dphidx)
       << " " << grr_KS * drdx * drdz + gthth_KS * dthdx * dthdz
    + gpp_KS * dphidx * dphidz + grp_KS * (drdx * dphidz + drdz * dphidx) << " " << grr_KS * drdy * drdz + gthth_KS * dthdy * dthdz
    + gpp_KS * dphidy * dphidz + grp_KS * (drdy * dphidz + drdz * dphidy) << endl;

    GYOTO_ERROR("test KerrH gmunu");*/

  //cout << "gzz= " << grr_KS << " " << drdz << " " << gthth_KS << " " << dthdz << " " << gpp_KS << " " << dphidz << " " << grp_KS << endl;

  if ((mu==0) && (nu==0)) return gtt_KS;
  
  if (((mu==0) && (nu==1)) || ((mu==1) && (nu==0)))
    return gtr_KS * drdx + gtp_KS * dphidx;
  if (((mu==0) && (nu==2)) || ((mu==2) && (nu==0)))
    return gtr_KS * drdy + gtp_KS * dphidy;
  if (((mu==0) && (nu==3)) || ((mu==3) && (nu==0)))
    return gtr_KS * drdz + gtp_KS * dphidz;
  //cout << "gxx = " << grr_KS << " " << drdx << " " << gthth_KS << " " << dthdx << " " << gpp_KS << " " << dphidx << " " << grp_KS << endl;
  if ((mu==1) && (nu==1))
    return grr_KS * drdx * drdx + gthth_KS * dthdx * dthdx
      + gpp_KS * dphidx * dphidx + 2. * grp_KS * drdx * dphidx;
  if ((mu==2) && (nu==2))
    return grr_KS * drdy * drdy + gthth_KS * dthdy * dthdy
      + gpp_KS * dphidy * dphidy + 2. * grp_KS * drdy * dphidy;
  if ((mu==3) && (nu==3))
    return grr_KS * drdz * drdz + gthth_KS * dthdz * dthdz
      + gpp_KS * dphidz * dphidz + 2. * grp_KS * drdz * dphidz;

  if (((mu==1) && (nu==2)) || ((mu==2) && (nu==1)))
    return grr_KS * drdx * drdy + gthth_KS * dthdx * dthdy
      + gpp_KS * dphidx * dphidy + grp_KS * (drdx * dphidy + drdy * dphidx);
  if (((mu==1) && (nu==3)) || ((mu==3) && (nu==1)))
    return grr_KS * drdx * drdz + gthth_KS * dthdx * dthdz
      + gpp_KS * dphidx * dphidz + grp_KS * (drdx * dphidz + drdz * dphidx);
  if (((mu==2) && (nu==3)) || ((mu==3) && (nu==2)))
    return grr_KS * drdy * drdz + gthth_KS * dthdy * dthdz
      + gpp_KS * dphidy * dphidz + grp_KS * (drdy * dphidz + drdz * dphidy);

  // TEST: Implement KS spin zero gmunu

  /*cout << "all gmunu= " << -1.+2./rr << " " << 2./r2*xx << " " << 2./r2*yy << " " << 2./r2*zz << " " <<
    1.+2*x2/(r2*rr) << " " << 1.+2*y2/(r2*rr) << " " << 1.+2*z2/(r2*rr) <<" " << (1.+2./rr)*xx*yy/r2 + xx*yy*z2/(r2*(x2+y2)) - xx*yy/(x2+y2) << " " <<
    (1.+2./rr)*xx*zz/r2 - xx*zz/r2 << " " << (1.+2./rr)*yy*zz/r2 - yy*zz/r2 << endl;
    GYOTO_ERROR("test in gmunu");*/
  
  // if ((mu==0) && (nu==0)) return (-1.+2./rr);
  
  // if (((mu==0) && (nu==1)) || ((mu==1) && (nu==0)))
  //   return 2./r2*xx;
  // if (((mu==0) && (nu==2)) || ((mu==2) && (nu==0)))
  //   return 2./r2*yy;
  // if (((mu==0) && (nu==3)) || ((mu==3) && (nu==0)))
  //   return 2./r2*zz;

  // if ((mu==1) && (nu==1))
  //   return 1.+2*x2/(r2*rr);
  // if ((mu==2) && (nu==2))
  //   return 1.+2*y2/(r2*rr);
  // if ((mu==3) && (nu==3))
  //   return 1.+2*z2/(r2*rr);

  // if (((mu==1) && (nu==2)) || ((mu==2) && (nu==1)))
  //   return (1.+2./rr)*xx*yy/r2 + xx*yy*z2/(r2*(x2+y2)) - xx*yy/(x2+y2);
  // if (((mu==1) && (nu==3)) || ((mu==3) && (nu==1)))
  //   return (1.+2./rr)*xx*zz/r2 - xx*zz/r2;
  // if (((mu==2) && (nu==3)) || ((mu==3) && (nu==2)))
  //   return (1.+2./rr)*yy*zz/r2 - yy*zz/r2;
  
  GYOTO_ERROR("Should have returned before!");
  return 0.;
  
}

void KerrHarmonic::gmunu_up(double gup[4][4], const double * pos) const {

  //cout << "KH gmunuup" << endl;

  double aa = spin_, a2 = aa*aa;
  // Harmonic Cartesian coordinates and harmonic radius
  double xx = pos[1], yy = pos[2], zz = pos[3],
    //double xx = 10., yy = 9., zz = 8., // test!!!!
    x2 = xx*xx, y2 = yy*yy, z2 = zz*zz,
    rr = pow(x2+y2+z2,0.5), r2 = rr*rr;
  double Rcal2 = pow((r2 - a2)*(r2 - a2) + 4.*a2*z2,0.5); // useful quantity,
                                            // this is \mathcal{R} in FV notes.

  // KS coordinates derived from harmonic
  double rKS = 1. + 1/pow(2.,0.5)*pow(r2 - a2 + Rcal2,0.5), // KS radius
    uu = rKS-1., // u = rKS-M is a useful quantity, see FV notes
    //u2 = uu*uu,
    costhKS = zz/uu, costhKS2 = costhKS*costhKS, // KS cos(theta)
    sinthKS = pow(1 - costhKS2,0.5), sinthKS2 = sinthKS*sinthKS,
    rho2KS = rKS*rKS + a2 * costhKS2, // usual rho2 quantity in KS metric
    cosphKS = (uu*xx + aa*yy)/(x2+y2) * sinthKS, // playing with the x ans y relations
    sinphKS = (uu*yy - aa*xx)/(x2+y2) * sinthKS;
  
  // Spherical KS inverse metric
  double gtt_KS = -1.-2.*rKS/rho2KS, // contravariant indices throughout
    gtr_KS = 2.*rKS/rho2KS, grr_KS = (rKS*rKS - 2*rKS + a2)/rho2KS,
    gthth_KS = 1./rho2KS, gpp_KS = 1./(rho2KS*sinthKS2),
    grp_KS = aa/rho2KS;

  // Derivative of H coord wrt KS ones
  double dxdr = cosphKS*sinthKS, dydr = sinphKS*sinthKS, dzdr = costhKS,
    dxdth = xx*costhKS/sinthKS, dydth = yy*costhKS/sinthKS, dzdth = -uu*sinthKS,
    dxdp = -yy, dydp = xx, dzdp = 0.;

  
  gup[0][0] = gtt_KS; // g^tt

  gup[0][1] = gup[1][0] = gtr_KS*dxdr; // g^tx
  gup[0][2] = gup[2][0] = gtr_KS*dydr; // g^ty
  gup[0][3] = gup[3][0] = gtr_KS*dzdr; // g^tz

  gup[1][1] = grr_KS * dxdr*dxdr + gthth_KS * dxdth*dxdth + gpp_KS * dxdp*dxdp + 2.*grp_KS * dxdr*dxdp;
  gup[2][2] = grr_KS * dydr*dydr + gthth_KS * dydth*dydth + gpp_KS * dydp*dydp + 2.*grp_KS * dydr*dydp;
  gup[3][3] = grr_KS * dzdr*dzdr + gthth_KS * dzdth*dzdth + gpp_KS * dzdp*dzdp + 2.*grp_KS * dzdr*dzdp;
  gup[1][2] = gup[2][1] = grr_KS * dxdr*dydr + gthth_KS * dxdth*dydth + gpp_KS * dxdp*dydp + grp_KS * (dxdr*dydp+dydr*dxdp);
  gup[1][3] = gup[3][1] = grr_KS * dxdr*dzdr + gthth_KS * dxdth*dzdth + gpp_KS * dxdp*dzdp + grp_KS * (dxdr*dzdp+dzdr*dxdp);
  gup[2][3] = gup[3][2] = grr_KS * dzdr*dydr + gthth_KS * dzdth*dydth + gpp_KS * dzdp*dydp + grp_KS * (dzdr*dydp+dydr*dzdp);

  /* cout << "all gmunuup " << gup[0][0] << " " << gup[0][1] << " " << gup[0][2] << " " << gup[0][3] << " " <<
    gup[1][1] << " " << gup[2][2] << " " << gup[3][3] << " " << gup[1][2] << " " << gup[1][3] << " " << gup[2][3] << endl;
    GYOTO_ERROR("test gmunuup");*/
}

void KerrHarmonic::jacobian(double jac[4][4][4], const double * pos) const {

  double aa = spin_, a2 = aa*aa;
  // Harmonic Cartesian coordinates and harmonic radius
  double xx = pos[1], yy = pos[2], zz = pos[3],
    //xx = 10., yy=9., zz=8., // TEST
    x2 = xx*xx, y2 = yy*yy, z2 = zz*zz,
    rr = pow(x2+y2+z2,0.5), r2 = rr*rr;
  double Rcal2 = pow((r2 - a2)*(r2 - a2) + 4.*a2*z2,0.5),
    Rcal = pow(Rcal2,0.5), Rcal3 = Rcal*Rcal2;

  // KS coordinates derived from harmonic
  double rKS = 1. + 1/pow(2.,0.5)*pow(r2 - a2 + Rcal2,0.5), // KS radius
    rKS2 = rKS*rKS,
    uu = rKS-1., // u = rKS-M is a useful quantity, see FV notes
    u2 = uu*uu,
    costhKS = zz/uu, costhKS2 = costhKS*costhKS, // KS cos(theta)
    sinthKS = pow(1 - costhKS2,0.5), sinthKS2 = sinthKS*sinthKS,
    rho2KS = rKS*rKS + a2 * costhKS2; // usual rho2 quantity in KS metric

  // Covariant spherical KS metric coefficients expressed in
  // harmonic coordinates
  double //gtt_KS = -(1. - 2*rKS/rho2KS),
    gtr_KS = 2*rKS/rho2KS, gtp_KS = -2*aa*rKS*sinthKS2/rho2KS,
    grr_KS = 1.+2*rKS/rho2KS, grp_KS = -aa*(1.+2*rKS/rho2KS)*sinthKS2,
    gthth_KS = rho2KS,
    gpp_KS = (rKS*rKS + a2 + 2.*a2*rKS*sinthKS2/rho2KS)*sinthKS2;

  // Derivative of KS coord wrt harmonic ones
  double drdx = uu*xx/Rcal2, drdy = uu*yy/Rcal2,
    drdz = (u2 + a2)*zz/(uu*Rcal2),//2.*(u2 + a2)*zz/(uu*Rcal2),
    dthdx = xx*zz/(Rcal2*pow(u2-z2,0.5)),
    dthdy = yy*zz/(Rcal2*pow(u2-z2,0.5)),
    dthdz = 1./pow(u2-z2,0.5) * (z2*(u2+a2)/(Rcal2*u2) - 1.),
    dphidx = aa/(a2+u2) * drdx - yy/(x2+y2),
    dphidy = aa/(a2+u2) * drdy + xx/(x2+y2),
    dphidz = aa/(a2+u2) * drdz;

  // Derivatives of Rcal
  double dRcaldx = xx*(r2-a2) / Rcal3,
    dRcaldy = yy*(r2-a2) / Rcal3,
    dRcaldz = zz * (r2+a2) / Rcal3; // it is indeed a plus here unlike the 2 previous ones

  // Second derivative of KS coord wrt harmonic ones // all checked
  double d2rdx2 = xx/Rcal2*drdx + uu/Rcal2 - 2.*uu*xx/Rcal3 * dRcaldx,
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
    d2phidx2 = aa/(a2+u2) * d2rdx2 - 2.*aa * uu / ((a2+u2)*(a2+u2)) * drdx*drdx + 2.*xx*yy/((x2+y2)*(x2+y2)),
    d2phidy2 = aa/(a2+u2) * d2rdy2 - 2.*aa * uu / ((a2+u2)*(a2+u2)) * drdy*drdy - 2.*xx*yy/((x2+y2)*(x2+y2)),
    d2phidz2 = aa/(a2+u2) * d2rdz2 - 2.*aa * uu / ((a2+u2)*(a2+u2)) * drdz*drdz,
    d2phidydx = aa/(a2+u2) * d2rdydx - 2.*aa * uu / ((a2+u2)*(a2+u2)) * drdx*drdy - 1./(x2+y2) + 2.*y2/((x2+y2)*(x2+y2)),
    d2phidzdx = aa/(a2+u2) * d2rdzdx - 2.*aa * uu / ((a2+u2)*(a2+u2)) * drdz*drdx,
    d2phidzdy = aa/(a2+u2) * d2rdzdy - 2.*aa * uu / ((a2+u2)*(a2+u2)) * drdz*drdy;

  // Derivatives of spherical KS metric coefs
  double drho2dx = 2.*(rKS*drdx - a2*costhKS*sinthKS*dthdx),
    drho2dy = 2.*(rKS*drdy - a2*costhKS*sinthKS*dthdy),
    drho2dz = 2.*(rKS*drdz - a2*costhKS*sinthKS*dthdz),
    dr_over_rho2_dx = (drdx*rho2KS - rKS*drho2dx)/(rho2KS*rho2KS),
    dr_over_rho2_dy = (drdy*rho2KS - rKS*drho2dy)/(rho2KS*rho2KS),
    dr_over_rho2_dz = (drdz*rho2KS - rKS*drho2dz)/(rho2KS*rho2KS),
    dsinth2dx = 2.*costhKS*sinthKS*dthdx,
    dsinth2dy = 2.*costhKS*sinthKS*dthdy,
    dsinth2dz = 2.*costhKS*sinthKS*dthdz,
    dr_over_rho2_sin2th_dx = dr_over_rho2_dx * sinthKS2 + rKS/rho2KS * dsinth2dx,
    dr_over_rho2_sin2th_dy = dr_over_rho2_dy * sinthKS2 + rKS/rho2KS * dsinth2dy,
    dr_over_rho2_sin2th_dz = dr_over_rho2_dz * sinthKS2 + rKS/rho2KS * dsinth2dz;

  //cout << "droverrho2dx= " << dr_over_rho2_dx << " " << drho2dx << endl;
    
  double dgtt_dx = 2. * dr_over_rho2_dx,
    dgtt_dy = 2. * dr_over_rho2_dy,
    dgtt_dz = 2. * dr_over_rho2_dz,
    dgtr_dx = 2. * dr_over_rho2_dx,
    dgtr_dy = 2. * dr_over_rho2_dy,
    dgtr_dz = 2. * dr_over_rho2_dz,
    dgtp_dx = -2.*aa*dr_over_rho2_sin2th_dx,
    dgtp_dy = -2.*aa*dr_over_rho2_sin2th_dy,
    dgtp_dz = -2.*aa*dr_over_rho2_sin2th_dz,
    dgrr_dx = 2. * dr_over_rho2_dx,
    dgrr_dy = 2. * dr_over_rho2_dy,
    dgrr_dz = 2. * dr_over_rho2_dz,
    dgrp_dx = -aa * (2. * dr_over_rho2_dx) * sinthKS2 - aa*(1.+2.*rKS/rho2KS)*dsinth2dx,
    dgrp_dy = -aa * (2. * dr_over_rho2_dy) * sinthKS2 - aa*(1.+2.*rKS/rho2KS)*dsinth2dy,
    dgrp_dz = -aa * (2. * dr_over_rho2_dz) * sinthKS2 - aa*(1.+2.*rKS/rho2KS)*dsinth2dz,
    dgthth_dx = drho2dx,
    dgthth_dy = drho2dy,
    dgthth_dz = drho2dz,
    dgpp_dx = (2.*rKS*drdx + 2.*a2*dr_over_rho2_sin2th_dx)*sinthKS2 + (rKS2 + a2 + 2.*a2*rKS*sinthKS2/rho2KS) * dsinth2dx,
    dgpp_dy = (2.*rKS*drdy + 2.*a2*dr_over_rho2_sin2th_dy)*sinthKS2 + (rKS2 + a2 + 2.*a2*rKS*sinthKS2/rho2KS) * dsinth2dy,
    dgpp_dz = (2.*rKS*drdz + 2.*a2*dr_over_rho2_sin2th_dz)*sinthKS2 + (rKS2 + a2 + 2.*a2*rKS*sinthKS2/rho2KS) * dsinth2dz;

  for (int mu=0; mu<4; ++mu){
    for (int nu=0; nu<4; ++nu){
      jac[0][mu][nu] = 0.; // \partial_t = 0
    }
  }

  // Formulas below are checked
  jac[1][0][0] = dgtt_dx ; // \partial_x g_tt
  jac[2][0][0] = dgtt_dy; // \partial_y g_tt
  jac[3][0][0] = dgtt_dz; // \partial_z g_tt

  jac[1][0][1] = jac[1][1][0] = dgtr_dx * drdx + gtr_KS * d2rdx2
    + dgtp_dx * dphidx + gtp_KS * d2phidx2;
  jac[2][0][1] = jac[2][1][0] = dgtr_dy * drdx + gtr_KS * d2rdydx
    + dgtp_dy * dphidx + gtp_KS * d2phidydx; 
  jac[3][0][1] = jac[3][1][0] = dgtr_dz * drdx + gtr_KS * d2rdzdx
    + dgtp_dz * dphidx + gtp_KS * d2phidzdx;

  jac[1][0][2] = jac[1][2][0] = dgtr_dx * drdy + gtr_KS * d2rdydx
    + dgtp_dx * dphidy + gtp_KS * d2phidydx;
  jac[2][0][2] = jac[2][2][0] = dgtr_dy * drdy + gtr_KS * d2rdy2
    + dgtp_dy * dphidy + gtp_KS * d2phidy2; 
  jac[3][0][2] = jac[3][2][0] = dgtr_dz * drdy + gtr_KS * d2rdzdy
    + dgtp_dz * dphidy + gtp_KS * d2phidzdy;

  jac[1][0][3] = jac[1][3][0] = dgtr_dx * drdz + gtr_KS * d2rdzdx
    + dgtp_dx * dphidz + gtp_KS * d2phidzdx;
  jac[2][0][3] = jac[2][3][0] = dgtr_dy * drdz + gtr_KS * d2rdzdy
    + dgtp_dy * dphidz + gtp_KS * d2phidzdy; 
  jac[3][0][3] = jac[3][3][0] = dgtr_dz * drdz + gtr_KS * d2rdz2
    + dgtp_dz * dphidz + gtp_KS * d2phidz2; 

  jac[1][1][1] = dgrr_dx * drdx * drdx + 2. * grr_KS * d2rdx2 * drdx
    + dgthth_dx * dthdx * dthdx + 2. * gthth_KS * d2thdx2 * dthdx
    + dgpp_dx * dphidx * dphidx + 2. * gpp_KS * d2phidx2 * dphidx
    + 2. * (dgrp_dx * drdx * dphidx
	    + grp_KS * d2rdx2 * dphidx + grp_KS * drdx * d2phidx2);
  jac[2][1][1] = dgrr_dy * drdx * drdx + 2. * grr_KS * d2rdydx * drdx
    + dgthth_dy * dthdx * dthdx + 2. * gthth_KS * d2thdydx * dthdx
    + dgpp_dy * dphidx * dphidx + 2. * gpp_KS * d2phidydx * dphidx
    + 2. * (dgrp_dy * drdx * dphidx
	    + grp_KS * d2rdydx * dphidx + grp_KS * drdx * d2phidydx);
  jac[3][1][1] = dgrr_dz * drdx * drdx + 2. * grr_KS * d2rdzdx * drdx
    + dgthth_dz * dthdx * dthdx + 2. * gthth_KS * d2thdzdx * dthdx
    + dgpp_dz * dphidx * dphidx + 2. * gpp_KS * d2phidzdx * dphidx
    + 2. * (dgrp_dz * drdx * dphidx
	    + grp_KS * d2rdzdx * dphidx + grp_KS * drdx * d2phidzdx);
  jac[1][2][2] = dgrr_dx * drdy * drdy + 2. * grr_KS * d2rdydx * drdy
    + dgthth_dx * dthdy * dthdy + 2. * gthth_KS * d2thdydx * dthdy
    + dgpp_dx * dphidy * dphidy + 2. * gpp_KS * d2phidydx * dphidy
    + 2. * (dgrp_dx * drdy * dphidy
	  + grp_KS * d2rdydx * dphidy + grp_KS * drdy * d2phidydx);
  jac[2][2][2] = dgrr_dy * drdy * drdy + 2. * grr_KS * d2rdy2 * drdy
    + dgthth_dy * dthdy * dthdy + 2. * gthth_KS * d2thdy2 * dthdy
    + dgpp_dy * dphidy * dphidy + 2. * gpp_KS * d2phidy2 * dphidy
    + 2. * (dgrp_dy * drdy * dphidy
	  + grp_KS * d2rdy2 * dphidy + grp_KS * drdy * d2phidy2);
  jac[3][2][2] = dgrr_dz * drdy * drdy + 2. * grr_KS * d2rdzdy * drdy
    + dgthth_dz * dthdy * dthdy + 2. * gthth_KS * d2thdzdy * dthdy
    + dgpp_dz * dphidy * dphidy + 2. * gpp_KS * d2phidzdy * dphidy
    + 2.* (dgrp_dz * drdy * dphidy
	  + grp_KS * d2rdzdy * dphidy + grp_KS * drdy * d2phidzdy);
  jac[1][3][3] = dgrr_dx * drdz * drdz + 2. * grr_KS * d2rdzdx * drdz
    + dgthth_dx * dthdz * dthdz + 2. * gthth_KS * d2thdzdx * dthdz
    + dgpp_dx * dphidz * dphidz + 2. * gpp_KS * d2phidzdx * dphidz
    + 2.* (dgrp_dx * drdz * dphidz
	  + grp_KS * d2rdzdx * dphidz + grp_KS * drdz * d2phidzdx);
  jac[2][3][3] = dgrr_dy * drdz * drdz + 2. * grr_KS * d2rdzdy * drdz
    + dgthth_dy * dthdz * dthdz + 2. * gthth_KS * d2thdzdy * dthdz
    + dgpp_dy * dphidz * dphidz + 2. * gpp_KS * d2phidzdy * dphidz
    + 2.* (dgrp_dy * drdz * dphidz
	  + grp_KS * d2rdzdy * dphidz + grp_KS * drdz * d2phidzdy);
  jac[3][3][3] = dgrr_dz * drdz * drdz + 2. * grr_KS * d2rdz2 * drdz
    + dgthth_dz * dthdz * dthdz + 2. * gthth_KS * d2thdz2 * dthdz
    + dgpp_dz * dphidz * dphidz + 2. * gpp_KS * d2phidz2 * dphidz
    + 2. *(dgrp_dz * drdz * dphidz
	  + grp_KS * d2rdz2 * dphidz + grp_KS * drdz * d2phidz2);

  jac[1][1][2] = jac[1][2][1] = dgrr_dx * drdx * drdy
    + grr_KS * d2rdx2 * drdy + grr_KS * drdx * d2rdydx
    + dgthth_dx * dthdx * dthdy
    + gthth_KS * d2thdx2 * dthdy + gthth_KS * dthdx * d2thdydx
    + dgpp_dx * dphidx * dphidy
    + gpp_KS * d2phidx2 * dphidy + gpp_KS * dphidx * d2phidydx
    + dgrp_dx * (drdx * dphidy + drdy * dphidx)
    + grp_KS * (d2rdx2 * dphidy + drdx * d2phidydx + d2rdydx * dphidx + drdy * d2phidx2);
  jac[2][1][2] = jac[2][2][1] = dgrr_dy * drdx * drdy
    + grr_KS * d2rdydx * drdy + grr_KS * drdx * d2rdy2
    + dgthth_dy * dthdx * dthdy
    + gthth_KS * d2thdydx * dthdy + gthth_KS * dthdx * d2thdy2
    + dgpp_dy * dphidx * dphidy
    + gpp_KS * d2phidydx * dphidy + gpp_KS * dphidx * d2phidy2
    + dgrp_dy * (drdx * dphidy + drdy * dphidx)
    + grp_KS * (d2rdydx * dphidy + drdx * d2phidy2 + d2rdy2 * dphidx + drdy * d2phidydx);
  jac[3][1][2] = jac[3][2][1] = dgrr_dz * drdx * drdy
    + grr_KS * d2rdzdx * drdy + grr_KS * drdx * d2rdzdy
    + dgthth_dz * dthdx * dthdy
    + gthth_KS * d2thdzdx * dthdy + gthth_KS * dthdx * d2thdzdy
    + dgpp_dz * dphidx * dphidy
    + gpp_KS * d2phidzdx * dphidy + gpp_KS * dphidx * d2phidzdy
    + dgrp_dz * (drdx * dphidy + drdy * dphidx)
    + grp_KS * (d2rdzdx * dphidy + drdx * d2phidzdy + d2rdzdy * dphidx + drdy * d2phidzdx);

  jac[1][1][3] = jac[1][3][1] = dgrr_dx * drdx * drdz 
    + grr_KS * d2rdx2 * drdz + grr_KS * drdx * d2rdzdx
    + dgthth_dx * dthdx * dthdz
    + gthth_KS * d2thdx2 * dthdz + gthth_KS * dthdx * d2thdzdx
    + dgpp_dx * dphidx * dphidz
    + gpp_KS * d2phidx2 * dphidz + gpp_KS * dphidx * d2phidzdx
    + dgrp_dx * (drdx * dphidz + drdz * dphidx)
    + grp_KS * (d2rdx2 * dphidz + drdx * d2phidzdx + d2rdzdx * dphidx + drdz * d2phidx2);
  jac[2][1][3] = jac[2][3][1] = dgrr_dy * drdx * drdz
    + grr_KS * d2rdydx * drdz + grr_KS * drdx * d2rdzdy
    + dgthth_dy * dthdx * dthdz
    + gthth_KS * d2thdydx * dthdz + gthth_KS * dthdx * d2thdzdy
    + dgpp_dy * dphidx * dphidz
    + gpp_KS * d2phidydx * dphidz + gpp_KS * dphidx * d2phidzdy
    + dgrp_dy * (drdx * dphidz + drdz * dphidx)
    + grp_KS * (d2rdydx * dphidz + drdx * d2phidzdy + d2rdzdy * dphidx + drdz * d2phidydx);
  jac[3][1][3] = jac[3][3][1] = dgrr_dz * drdx * drdz 
    + grr_KS * d2rdzdx * drdz + grr_KS * drdx * d2rdz2
    + dgthth_dz * dthdx * dthdz
    + gthth_KS * d2thdzdx * dthdz + gthth_KS * dthdx * d2thdz2
    + dgpp_dz * dphidx * dphidz
    + gpp_KS * d2phidzdx * dphidz + gpp_KS * dphidx * d2phidz2
    + dgrp_dz * (drdx * dphidz + drdz * dphidx)
    + grp_KS * (d2rdzdx * dphidz + drdx * d2phidz2 + d2rdz2 * dphidx + drdz * d2phidzdx);

  jac[1][2][3] = jac[1][3][2] = dgrr_dx * drdy * drdz
    + grr_KS * d2rdydx * drdz + grr_KS * drdy * d2rdzdx
    + dgthth_dx * dthdy * dthdz
    + gthth_KS * d2thdydx * dthdz + gthth_KS * dthdy * d2thdzdx
    + dgpp_dx * dphidy * dphidz
    + gpp_KS * d2phidydx * dphidz + gpp_KS * dphidy * d2phidzdx
    + dgrp_dx * (drdy * dphidz + drdz * dphidy)
    + grp_KS * (d2rdydx * dphidz + drdy * d2phidzdx + d2rdzdx * dphidy + drdz * d2phidydx);
  jac[2][2][3] = jac[2][3][2] = dgrr_dy * drdy * drdz
    + grr_KS * d2rdy2 * drdz + grr_KS * drdy * d2rdzdy
    + dgthth_dy * dthdy * dthdz
    + gthth_KS * d2thdy2 * dthdz + gthth_KS * dthdy * d2thdzdy
    + dgpp_dy * dphidy * dphidz
    + gpp_KS * d2phidy2 * dphidz + gpp_KS * dphidy * d2phidzdy
    + dgrp_dy * (drdy * dphidz + drdz * dphidy)
    + grp_KS * (d2rdy2 * dphidz + drdy * d2phidzdy + d2rdzdy * dphidy + drdz * d2phidy2);
  jac[3][2][3] = jac[3][3][2] = dgrr_dz * drdy * drdz
    + grr_KS * d2rdzdy * drdz + grr_KS * drdy * d2rdz2
    + dgthth_dz * dthdy * dthdz
    + gthth_KS * d2thdzdy * dthdz + gthth_KS * dthdy * d2thdz2
    + dgpp_dz * dphidy * dphidz
    + gpp_KS * d2phidzdy * dphidz + gpp_KS * dphidy * d2phidz2
    + dgrp_dz * (drdy * dphidz + drdz * dphidy)
    + grp_KS * (d2rdzdy * dphidz + drdy * d2phidz2 + d2rdz2 * dphidy + drdz * d2phidzdy);

  /*cout << "all jac" ;
  for (int a=0; a<4;++a){
    for (int mu=0; mu<4;++mu){
      for (int nu=0; nu<4;++nu){
	cout << "a mu nu=" << a << " " << mu << " " << nu << " " << jac[a][mu][nu] << endl ;
      }
    }
  }
  //cout << "all gmunuup " << gup[0][0] << " " << gup[0][1] << " " << gup[0][2] << " " << gup[0][3] << " " <<
  //gup[1][1] << " " << gup[2][2] << " " << gup[3][3] << " " << gup[1][2] << " " << gup[1][3] << " " << gup[2][3] << endl;
  GYOTO_ERROR("test jac"); */
  
}
  
int KerrHarmonic::isStopCondition(double const * const coord) const {
  double rsinkKS = 0.;
  if (spin_*spin_<1.){ // Black hole solution with event horizon
    double rhorKS = 1 + sqrt(1 - spin_*spin_); // KS radius EH
    rsinkKS = rhorKS + GYOTO_KERR_HORIZON_SECURITY;
  }

  double xx=coord[1], yy=coord[2], zz=coord[3],
    rr = pow(xx*xx+yy*yy+zz*zz,0.5), r2 = rr*rr,
    a2 = spin_*spin_,
    Rcal2 = pow((r2 - a2)*(r2 - a2) + 4.*a2*zz*zz,0.5),
    rKS = 1. + 1/pow(2.,0.5)*pow(r2 - a2 + Rcal2,0.5); // current KS radius

  return rKS < rsinkKS ;
}

void KerrHarmonic::circularVelocity(double const * coor, double* vel,
				    double dir) const {

  double xx=coor[1], yy=coor[2], zz=coor[3],
    rr = pow(xx*xx+yy*yy+zz*zz,0.5), r2 = rr*rr,
    a2 = spin_*spin_,
    Rcal2 = pow((r2 - a2)*(r2 - a2) + 4.*a2*zz*zz,0.5),
    rKS = 1. + 1/pow(2.,0.5)*pow(r2 - a2 + Rcal2,0.5);

  double Omega = dir/(pow(rKS,3./2.) + spin_); // same value as in KerrBL

  vel[1] = -coor[2]*Omega;
  vel[2] =  coor[1]*Omega;
  vel[3] = 0.;
  vel[0] = SysPrimeToTdot(coor, vel+1);
  vel[1] *= vel[0];
  vel[2] *= vel[0];
}

// // PolishDoughnut specific functions
// double KerrHarmonic::getPotential(double const pos[4], double l_cst) const {
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
// double KerrHarmonic::getSpecificAngularMomentum(double rr) const {
//   // this is l = -u_phi/u_t for a circular equatorial 4-velocity
//   double qq=charge_, q2=qq*qq, r2=rr*rr;
//   if (rr<q2 or 1. - 2./rr + q2/r2 == 0.){
//     cout << "r, q2, term= " << rr << " " << q2 << " " << 1. - 2./rr + q2/r2 << endl;
//     GYOTO_ERROR("bad values in l !");
//   }
//   return sqrt(rr-q2)/(1. - 2./rr + q2/r2);
// }
