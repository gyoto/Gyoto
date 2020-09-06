/*
    Copyright 2011-2016 Frederic Vincent, Thibaut Paumard, Odele Straub

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
#include "GyotoXmasSKS.h"
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


//// Property list:
////////////////////////////////////////////////////////
GYOTO_PROPERTY_START(XmasSKS,
         "Metric around a Schwarzschild-like black hole, in (cartesian) Kerr-Schild coordinates.")
// INPUT metric properties 
GYOTO_PROPERTY_DOUBLE(XmasSKS, Mdm_over_Mbh, Mdm_over_Mbh,
          "Ratio between dark matter and BH mass (q:= M_dm/M_bh)")
GYOTO_PROPERTY_DOUBLE(XmasSKS, Rdm, Rdm,
          "Boundary radius of the DM distribution.")
GYOTO_PROPERTY_DOUBLE(XmasSKS, gamma, gamma,
          "Slope of the DM distribution.")
// Integration properties
GYOTO_PROPERTY_DOUBLE(XmasSKS, HorizonSecurity, horizonSecurity,
          "Thickness of sink layer around horizon (geometrical units, 0.01).")
GYOTO_PROPERTY_BOOL(XmasSKS, GenericIntegrator, SpecificIntegrator,
        genericIntegrator,
        "Which version of the Legacy integrator should be used (specific).")
GYOTO_PROPERTY_DOUBLE(XmasSKS, DiffTol, difftol,
          "Tuning parameter for the specific Legacy integrator (0.01).")
//
GYOTO_PROPERTY_END(XmasSKS, Generic::properties)


//#define GYOTO_MIN_THETA 1e-5 //seems too big
#define GYOTO_MIN_THETA 1e-10



// Construct
// The minimal constructor just sets the coordinate kind and
// the metric kind name.     
// //////////////////////            
XmasSKS::XmasSKS() :
  Generic(GYOTO_COORDKIND_CARTESIAN, "XmasSKS"),
  Mdm_over_Mbh_(0.), 
  Rdm_(0.),
  gamma_(0.),
  difftol_(GYOTO_XmasSKS_DEFAULT_DIFFTOL),
  rsink_(2.+GYOTO_KERR_HORIZON_SECURITY),
  drhor_(GYOTO_KERR_HORIZON_SECURITY), 
  generic_integrator_(false)
{}



// Clone 
// If the metric class is not trivial (e.g. contains arrays), 
// it may be necessary to implement the copy constructor as well. 
// //////////////////////////////////////////////////////////////
XmasSKS * XmasSKS::clone () const {
 return new XmasSKS(*this); }







// Accessors
// /////////

GYOTO_PROPERTY_ACCESSORS(XmasSKS, double, Mdm_over_Mbh_, Mdm_over_Mbh)
GYOTO_PROPERTY_ACCESSORS(XmasSKS, double, Rdm_, Rdm)
GYOTO_PROPERTY_ACCESSORS(XmasSKS, double, gamma_, gamma)

double XmasSKS::difftol() const { return difftol_;}
void XmasSKS::difftol(double t) {difftol_=t;}

void XmasSKS::horizonSecurity(const double drhor) {
  drhor_=drhor;
  rsink_=2.+drhor_;
  tellListeners();
}
double XmasSKS::horizonSecurity() const {return drhor_; }

void XmasSKS::genericIntegrator(bool t)
{
  generic_integrator_=t;
  tellListeners();
}
bool XmasSKS::genericIntegrator() const {return generic_integrator_;}









// METRIC 
// (cartesian KS coord)
// ------------------------------------------------
// This Schwarzschild-like metric takes into account 
// a smooth dark matter distribution of
// mass=Mdm, radius=Rdm, and slope=beta.
// REFERENCE : Lacroix + Silk (2013)
// //////////////////////////////////////////////////
void XmasSKS::gmunu(double g[4][4], const double * pos) const {
  double x = pos[1], y=pos[2], z=pos[3];
  double x2=x*x, y2=y*y, z2=z*z;
  double r2= x2 + y2 + z2;
  double r = sqrt(r2);
  double r_2 = r-2., r_2_2 = r_2*r_2;

  double beta = 3. - gamma_; //gamma = 7./3., 
  double B = 1.- 2./r * (1. + Mdm_over_Mbh_ * pow(r/Rdm_,beta));
  double tmp1 = 1./B - 4.*B/r_2_2;

  int mu, nu;
  for (mu=0; mu<4; ++mu)
    for (nu=0; nu<4; ++nu)
      g[mu][nu]=0.;

  
  //diagonal metric coefficients
  g[0][0] = - B;                                         
  g[1][1] = (y2+z2)/r2 + x2/r2*tmp1;                                
  g[2][2] = (x2+z2)/r2 + y2/r2*tmp1;                                      
  g[3][3] = (x2+y2)/r2 + z2/r2*tmp1;   
  // off diagonal metric coefficients
  g[0][1] = g[1][0] = 2.*x/(r*r_2)*B;
  g[0][2] = g[2][0] = 2.*y/(r*r_2)*B;
  g[0][3] = g[3][0] = 2.*z/(r*r_2)*B;
  g[1][2] = g[2][1] = -x*y/r2 + x*y/r2*tmp1;
  g[1][3] = g[3][1] = -x*z/r2 + x*z/r2*tmp1;
  g[2][3] = g[3][2] = -y*z/r2 + y*z/r2*tmp1;



  GYOTO_DEBUG<<"g_tt = "<< g[0][0] <<endl;
  GYOTO_DEBUG<<"g_xx = "<< g[1][1] <<endl;
  GYOTO_DEBUG<<"g_yy = "<< g[2][2] <<endl;
  GYOTO_DEBUG<<"g_zz = "<< g[3][3] <<endl;
  GYOTO_DEBUG<<"g_tx = "<< g[0][1] <<endl;
  GYOTO_DEBUG<<"g_ty = "<< g[0][2] <<endl;
  GYOTO_DEBUG<<"g_tz = "<< g[0][3] <<endl;
  GYOTO_DEBUG<<"g_xy = "<< g[1][2] <<endl;
  GYOTO_DEBUG<<"g_xz = "<< g[1][3] <<endl;
  GYOTO_DEBUG<<"g_yz = "<< g[2][3] <<endl;
  
  GYOTO_DEBUG<<"g_munu done"<<endl;
}

double XmasSKS::gmunu(const double * pos, int mu, int nu) const {
  double x = pos[1], y=pos[2], z=pos[3];
  double x2=x*x, y2=y*y, z2=z*z;
  double r2= x2 + y2 + z2;
  double r = sqrt(r2);
  double r_2 = r-2., r_2_2 = r_2*r_2;

  double beta = 3. - gamma_; //gamma = 7./3., 
  double B = 1.- 2./r * (1. + Mdm_over_Mbh_ * pow(r/Rdm_,beta));
  double tmp1 = 1./B - 4.*B/r_2_2;
  
  //diagonal metric coefficients
  if ((mu==0) && (nu==0)) return - B;
  if ((mu==1) && (nu==1)) return (y2+z2)/r2 + x2/r2*tmp1;
  if ((mu==2) && (nu==2)) return (x2+z2)/r2 + y2/r2*tmp1;
  if ((mu==3) && (nu==3)) return (x2+y2)/r2 + z2/r2*tmp1; 
  // off diagonal metric coefficients
  if (((mu==0) && (nu==1)) || ((mu==1) && (nu==0))) return 2.*x/(r*r_2)*B;
  if (((mu==0) && (nu==2)) || ((mu==2) && (nu==0))) return 2.*y/(r*r_2)*B;
  if (((mu==0) && (nu==3)) || ((mu==3) && (nu==0))) return 2.*z/(r*r_2)*B;
  if (((mu==1) && (nu==2)) || ((mu==2) && (nu==1))) return -x*y/r2 + x*y/r2*tmp1;
  if (((mu==1) && (nu==3)) || ((mu==3) && (nu==1))) return -x*z/r2 + x*z/r2*tmp1;
  if (((mu==2) && (nu==3)) || ((mu==3) && (nu==2))) return -y*z/r2 + y*z/r2*tmp1;
  
  return 0.;
} 






// INVERSE METRIC
// //////////////
void XmasSKS::gmunu_up(double gup[4][4], const double * pos) const {
  double x = pos[1], y=pos[2], z=pos[3];
  double x2=x*x, y2=y*y, z2=z*z;
  double r2= x2 + y2 + z2;
  double r = sqrt(r2);
  double r_2 = r-2., r_2_2 = r_2*r_2;

  double beta = 3. - gamma_; //gamma = 7./3., 
  double B = 1.- 2./r * (1. + Mdm_over_Mbh_ * pow(r/Rdm_,beta));
  double tmp2 = 1./B + 4.*B/r_2_2;
  

  int mu, nu;
  for (mu=0; mu<4; ++mu) for (nu=0; nu<4; ++nu) gup[mu][nu]=0.;

  //diagonal metric coefficients
  gup[0][0] = - tmp2; 
  gup[1][1] = (B*x2 + y2 + z2)/r2;
  gup[2][2] = (x2 + B*y2 + z2)/r2;
  gup[3][3] = (x2 + y2 + B*z2)/r2;
  //off diagonal metric coefficients
  gup[0][1] = gup[1][0] = 2.*B*x/(r*r_2);
  gup[0][2] = gup[2][0] = 2.*B*y/(r*r_2);
  gup[0][3] = gup[3][0] = 2.*B*z/(r*r_2);
  gup[1][2] = gup[2][1] = x*y/r2*(B - 1.);
  gup[1][3] = gup[3][1] = x*z/r2*(B - 1.);
  gup[2][3] = gup[3][2] = y*z/r2*(B - 1.);

}

double XmasSKS::gmunu_up(const double * pos, int mu, int nu) const {
  double x = pos[1], y=pos[2], z=pos[3];
  double x2=x*x, y2=y*y, z2=z*z;
  double r2= x2 + y2 + z2;
  double r = sqrt(r2);
  double r_2 = r-2., r_2_2 = r_2*r_2;

  double beta = 3. - gamma_; //gamma = 7./3., 
  double B = 1.- 2./r * (1. + Mdm_over_Mbh_ * pow(r/Rdm_,beta));
  double tmp2 = 1./B + 4.*B/r_2_2;

  //diagonal metric coefficients
  if ((mu==0) && (nu==0)) return - tmp2; 
  if ((mu==1) && (nu==1)) return (B*x2 + y2 + z2)/r2;
  if ((mu==2) && (nu==2)) return (x2 + B*y2 + z2)/r2;
  if ((mu==3) && (nu==3)) return (x2 + y2 + B*z2)/r2;
  // off diagonal metric coefficients
  if (((mu==0) && (nu==1)) || ((mu==1) && (nu==0))) return 2.*B*x/(r*r_2);
  if (((mu==0) && (nu==2)) || ((mu==2) && (nu==0))) return 2.*B*y/(r*r_2);
  if (((mu==0) && (nu==3)) || ((mu==3) && (nu==0))) return 2.*B*z/(r*r_2);
  if (((mu==1) && (nu==2)) || ((mu==2) && (nu==1))) return x*y/r2*(B - 1.);
  if (((mu==1) && (nu==3)) || ((mu==3) && (nu==1))) return x*z/r2*(B - 1.);
  if (((mu==2) && (nu==3)) || ((mu==3) && (nu==2))) return y*z/r2*(B - 1.);

  return 0.;
} 





// CHRISTOFFELS
// /////////////
int XmasSKS::christoffel(double dst[4][4][4], double const pos[4]) const
{

 // derivatives d/dx of the metric
 // d/dt = 0, no t-dependence 
  int a, mu, nu;
  for (a=0; a<4; ++a)
    for(mu=0; mu<4; ++mu)
      for(nu=0; nu<4; ++nu)
  dst[a][mu][nu]=0.;

  double x = pos[1], y=pos[2], z=pos[3];
  double x2=x*x, y2=y*y, z2=z*z;
  double x3=x*x2, y3=y*y2, z3=z*z2, x4=x2*x2, y4=y2*y2, z4=z2*z2;
  double x2y2 = x2 + y2, x2z2 = x2 + z2, y2z2 = y2 + z2;
  double r2= x2 + y2 + z2, r4=r2*r2, r6=r4*r2;
  double r = sqrt(r2), r3=r*r2;
  double r_2 = r - 2., r_2_2 = r_2*r_2, r_2_3 = r_2*r_2_2;

  double beta = 3. - gamma_; //gamma = 7./3., 
  double BB     = 2./r * beta * Mdm_over_Mbh_ * pow(r/Rdm_,beta);
  double B      = 1.- 2./r - BB/beta, B2 = B*B;
  double dB1    = (-1. + B + BB);
  double dB2    = (-2. + 2.*B + BB);
  double tmp1   = 1./B - 4.*B/r_2_2;
  double tmp2   = 1./B + 4.*B/r_2_2;
  double Btmp1  = tmp1 - 1.;
  double Btmp2  = 8.*B/r_2_3 + dB1/(r*B)*tmp2;
  double Bterm  = -2.*Btmp1 + r*Btmp2; 

  double Term01=1./(B*r_2_2*r_2_2*r*r4) * (8.*(-1. + BB 
    + B*(-1. + B*(-3. + 3.*B + BB)))*r2 - 4.*(3.*(-1. + BB) 
    + B*(-3. + B*(-7. + 5.*B + BB)))*r3 + 4.*(-3. + 8.*B + 6.*B2 + 3.*BB)*r4 
    - 2.*(-1. + B*(3. + 2.*B) + BB)*r*r4 + r3*(24.*B*B2 + B*(-56. + r2) 
    - (-1. + BB)*(24. + r2) + 2.*B2*(-28. + 4.*BB + r2)) + 2.*r2*(-16.*B*B2 
    - B*(-16. + 3.*r2) + (-1. + BB)*(8. + 3.*r2) - 2.*B2*(4.*BB + 3.*(-4. + r2))));
 

  // 0
  dst[0][0][0] = - B*dB1/(r_2*r);

  dst[0][0][1] = dst[0][1][0] = x * dB1*(4.*B2 - r_2_2)/(2.*B*r_2_2*r2);
  
  dst[0][0][2] = dst[0][2][0] = y * dB1*(4.*B2 - r_2_2)/(2.*B*r_2_2*r2);
  
  dst[0][0][3] = dst[0][3][0] = z * dB1*(4.*B2 - r_2_2)/(2.*B*r_2_2*r2);
  
  dst[0][1][1] = 1./(r_2_2*r*r4)*(-2.*r*tmp1*(-(-1. + BB)*r_2*r*x2 
    + B*r*(-2.*r2 + r3 + 4.*x2 - 3.*r*x2)) 
    + B*r_2*y2z2*(2.*r2*(-2. + tmp1) + r*tmp2*x2 + 2.*(-(-2. + tmp1)*x2 + y2z2)) 
    + B*r_2*x2*(2.*r2*tmp1 + dB1*x2/B2 + 8.*B*r*x2/r_2_3 + 4.*dB1*r4*r2*x2/r_2_2 - 
     2.*(tmp1*x2 + y2z2)));
  
  dst[0][1][2] = dst[0][2][1] = - x*y*Term01; 
 
  dst[0][1][3] = dst[0][3][1] = - x*z*Term01;

  dst[0][2][2] = 1./(r_2_2*r*r4)*(-2.*r*tmp1*(-(-1. + BB)*r_2*r*y2 
    + B*r*(-2.*r2 + r3 + 4.*y2 - 3.*r*y2)) 
    + B*r_2*x2z2*(2.*r2*(-2. + tmp1) + r*tmp2*y2 + 2.*(-(-2. + tmp1)*y2 + x2z2)) 
    + B*r_2*y2*(2.*r2*tmp1 + dB1*y2/B2 + 8.*B*r*y2/r_2_3 + 4.*dB1*r4*r2*y2/r_2_2 - 
     2.*(tmp1*y2 + x2z2)));

  dst[0][2][3] = dst[0][3][2] =  - y*z*Term01;

  dst[0][3][3] = 1./(r_2_2*r*r4)*(-2.*r*tmp1*(-(-1. + BB)*r_2*r*z2 
    + B*r*(-2.*r2 + r3 + 4.*z2 - 3.*r*z2)) 
    + B*r_2*x2y2*(2.*r2*(-2. + tmp1) + r*tmp2*z2 + 2.*(-(-2. + tmp1)*z2 + x2y2)) 
    + B*r_2*z2*(2.*r2*tmp1 + dB1*z2/B2 + 8.*B*r*z2/r_2_3 + 4.*dB1*r4*r2*z2/r_2_2 - 
     2.*(tmp1*z2 + x2y2)));
 

 // 1
  dst[1][0][0] = -B*dB1*x/(2.*r2);

  dst[1][0][1] = dst[1][1][0] = B*dB1*x2/(r_2*r3);
  
  dst[1][0][2] = dst[1][2][0] = B*dB1*x*y/(r_2*r3);
  
  dst[1][0][3] = dst[1][3][0] = B*dB1*x*z/(r_2*r3);

  dst[1][1][1] = -1./(2.*B*r_2_2*r2*r4) * x * (r_2_2*r*(r - BB*r)*x2 - 
    4.*B*B2*r2*(2.*r2 - 3.*x2 - 2.*y2z2) + 2.*B2*r2*(2.*(-1. + BB)*x2 + r_2_2*y2z2) 
    - B*r_2_2*(-x4 - 5.*x2*y2z2 - 4.*y2z2*y2z2 + 2.*r2*(x2 + 3.*y2z2)));
  
  dst[1][1][2] = dst[1][2][1] = -1./(2.*B*r_2_2*r2*r4) * y * (r_2_2*r*(r - BB*r)*x2 
    + 4.*B*B2*r2*x2 - 2.*B2*r2*(6. - 2.*BB - 4.*r + r2)*x2 
    + B*r_2_2*(-x4 + 2.*r2*(x2 - y2z2) + x2*y2z2 + 2.*y2z2*y2z2));
  
  dst[1][1][3] = dst[1][3][1] =  -1./(2.*B*r_2_2*r2*r4) * z * (r_2_2*r2*(1. - BB)*x2 
    + 4.*B*B2*r2*x2 - 2.*B2*r2*(6. - 2.*BB - 4.*r + r2)*x2 
    + B*r_2_2*(-x4 + 2.*r2*(x2 - y2z2) + x2*y2z2 + 2.*y2z2*y2z2));

  dst[1][2][2] = -1./(2.*B*r_2_2*r4) * x * (4.*B*B2*y2 - (-1. + BB)*r_2_2*y2 
    - B*r_2_2*(r2 + x2z2) + 2.*B2*(r_2_2*x2 + 2.*(-1. + BB)*y2 + r_2_2*z2));
  
  dst[1][2][3] = dst[1][3][2] = -1./(2.*B*r_2_2*r4) * x*y*z * (4.*B*B2 + B*r_2_2 
    - (-1. + BB)*r_2_2 + 2.*B2*(-6. + 2.*BB - (-4. + r)*r));
  
  dst[1][3][3] = 1./(2.*B*r_2_2*r4) * x * ((-1. + BB)*r_2_2*z2 + 4.*B*B2*(2.*r2 
    - 2.*x2y2 - 3.*z2) + B*r_2_2*(2.*x2y2 + z2) + 2.* B2*(-r_2_2*x2y2 - 2.* (-1. + BB)*z2));


  // 2
  dst[2][0][0] = -B*dB1*y/(2.*r2);

  dst[2][0][1] = dst[2][1][0] = B*dB1*x*y/(r_2*r3);
  
  dst[2][0][2] = dst[2][2][0] = B*dB1*y2/(r_2*r3);
  
  dst[2][0][3] = dst[2][3][0] = B*dB1*y*z/(r_2*r3);
    
  dst[2][1][1] = 1./(2.*B*r_2_2*r4) * y * (-dB1*(2. + 2.*B - r)*(-2. + 2.*B + r)*x2 
    - 2.*(-1. + B)*B*r_2_2*x2z2);
  
  dst[2][1][2] = dst[2][2][1] = - 1./(2.*B*r_2_2*r4) * x*y2 * (4.*B*B2 
    + B*r_2_2 - (-1. + BB)*r_2_2 + 2.*B2*(-6. + 2.*BB + 4.*r - r2)); 
  
  dst[2][1][3] = dst[2][3][1] = - 1./(2.*B*r_2_2*r4) * x*y*z * (4.*B*B2  
    + B*r_2_2 - (-1. + BB)*r_2_2 + 2.*B2*(-6. + 2.*BB + 4.*r - r2));  

  dst[2][2][2] = 1./(2.*B*r_2_2*r4) * y * (4.*B*B2*y2 + (-1. + BB)*r_2_2*y2 
    + B*r_2_2*(r2 + x2z2) - 2.*B2*(r_2_2*x2 + 2.*(-1. + BB)*y2 + r_2_2*z2));
  
  dst[2][2][3] = dst[2][3][2] = - 1./(2.*B*r_2_2*r4) * y2*z * (4.*B*B2  
    + B*r_2_2 - (-1. + BB)*r_2_2 + 2.*B2*(-6. + 2.*BB + 4.*r - r2));
  
  dst[2][3][3] = 1./(2.*B*r_2_2*r4) * y * (B*r_2_2*(r2 - (-1. + 2.*B)*x2y2) 
    - 4.*z2 + (-4.*B2*dB1 - 4.*BB*r + (-1. + BB)*r2 + 4.*(BB + r))*z2);



  // 3
  dst[3][0][0] = - B*dB1*z/(2.*r2);

  dst[3][0][1] = dst[3][1][0] = B*dB1*x*z/(r_2*r3);
  
  dst[3][0][2] = dst[3][2][0] = B*dB1*y*z/(r_2*r3);

  dst[3][0][3] = dst[3][3][0] = B*dB1*z2/(r_2*r3);
    
  dst[3][1][1] = 1./(2.*B*r_2_2*r4) * z * ((-1. + BB)*r_2_2*x2 - 4.*B*B2*x2 
    + B*r_2_2*(r2 + y2z2) - 6.*B2*(2.*(-1. + BB)*x2 + r_2_2*y2z2));
  
  dst[3][1][2] = dst[3][2][1] = - 1./(2.*B*r_2_2*r4) * x*y*z * (4.*B*B2 
    + B*r_2_2 - (-1. + BB)*r_2_2 + 6.* B2*(-6. + 2.*BB - (-4. + r)*r));
  
  dst[3][1][3] = dst[3][3][1] = - 1./(2.*B*r_2_2*r4) * x*z2 * (4.*B*B2 
    + B*r_2_2 - (-1. + BB)*r_2_2 + 2.*B2*(-6. + 2.*BB - (-4. + r)*r)); 

  dst[3][2][2] = - 1./(2.*B*r_2_2*r4) * z * (4.*B*B2*y2 - (-1. + BB)*r_2_2*y2 
    - B*r_2_2*(r2 + x2z2) + 2.*B2*(r_2_2*x2 + 2.*(-1. + BB)*y2 + r_2_2*z2));
  
  dst[3][2][3] = dst[3][3][2] = - 1./(2.*B*r_2_2*r4) * y*z2 * (4.*B*B2 
    + B*r_2_2 - (-1. + BB)*r_2_2 + 2.*B2*(-6. + 2.*BB - (-4. + r)*r));

  dst[3][3][3] = - 1./(2.*B*r_2_2*r4) * z * (4.*B*B2*z2 - B*r_2_2*(r2 + x2y2) 
    - (-1. + BB)*r_2_2*z2 + 2.*B2*(r_2_2*x2y2 + 2.*(-1. + BB)*z2));

  return 0;

} 







void XmasSKS::circularVelocity(double const coor[4], double vel[4],
			      double dir) const {

  if (keplerian_) {
    // If keplerian_ is true, let the generic implementation return
    // the Keplerian velocity instead of the true circular velocity
    Generic::circularVelocity(coor, vel, dir);
    return;
  }

  double rcross=sqrt(coor[1]*coor[1] + coor[2]*coor[2]);
  double Omega=dir*pow(rcross*rcross*rcross, -0.5); //angular Keplerian velocity
  
  vel[1] = -coor[2]*Omega;
  vel[2] =  coor[1]*Omega;
  vel[3] = 0.;
  vel[0] = SysPrimeToTdot(coor, vel+1);
  vel[1] *= vel[0];
  vel[2] *= vel[0];
  
}




int XmasSKS::isStopCondition(double const * const coord) const {
  double
    x=coord[1], y=coord[2], z=coord[3],
    Tdot=coord[4], xdot=coord[5], ydot=coord[6], zdot=coord[7],
    x2=x*x, y2=y*y, z2=z*z,
    r2=x2+y2+z2,
    r=sqrt(r2);
  double rdot=(x*xdot+y*ydot+z*zdot)/r;

  return (r<rsink_);
}





