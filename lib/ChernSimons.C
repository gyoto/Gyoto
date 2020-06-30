/*
    Copyright 2013, 2016 Frederic Vincent & Thibaut Paumard

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

#include "GyotoChernSimons.h"

#include <iostream>
#include <cmath>
#include <cstdlib>


using namespace Gyoto;
using namespace Gyoto::Metric;
using namespace std;

/// Properties
#include "GyotoProperty.h"
GYOTO_PROPERTY_START(ChernSimons, "Chern-Simons 1st order perturbation to KerrBL metric")
GYOTO_PROPERTY_DOUBLE(ChernSimons, DzetaCS, dzetaCS, "Chern-Simons coupling constant")
GYOTO_PROPERTY_END(ChernSimons, KerrBL::properties)

// accessor
void ChernSimons::dzetaCS(double d) {dzetaCS_=d;}
double ChernSimons::dzetaCS() const {return dzetaCS_;}

///

#define drhor 2e-1 // 1e-1 leads to problem close to horizon

Gyoto::Metric::ChernSimons::ChernSimons()
  : KerrBL(), dzetaCS_(0.)
{
  // The constructor should only initialize to default values, the
  // subcontractor will later call setParameter() to set each
  // parameter.
  kind("ChernSimons");
  GYOTO_DEBUG << "Building ChernSimons" << endl;
}

Gyoto::Metric::ChernSimons::ChernSimons(const ChernSimons &o)
  : KerrBL(o), dzetaCS_(o.dzetaCS_)
{
  // don't forget to copy each member here
  kind("ChernSimons");
  GYOTO_DEBUG << "Copying ChernSimons" << endl;
}
ChernSimons * ChernSimons::clone() const { return new ChernSimons(*this); }

Gyoto::Metric::ChernSimons::~ChernSimons()
{
  // delete each member pointers or decremente smartpointers by
  // setting them to NULL
  GYOTO_DEBUG << "Destroying ChernSimons";
}

double ChernSimons::gmunu(const double * pos, int mu, int nu) const {
  double r = pos[1];
  double sth2, cth2;
  sincos(pos[2], &sth2, &cth2);
  sth2*=sth2; cth2*=cth2;
  double r2=r*r, r3=r2*r;
  double a2=spin_*spin_;
  double sigma=r2+a2*cth2;
  double delta=r2-2.*r+a2;

  double ff=1.-2./r;

  if ((mu==0) && (nu==0)) {
    return -ff-2.*a2/r3*cth2;
  }
  if ((mu==1) && (nu==1)) return 1./ff+a2/(ff*r2)*(cth2-1./ff);
  if ((mu==2) && (nu==2)) return sigma;
  if ((mu==3) && (nu==3))
    return r2*sth2+a2*sth2*(1.+2./r*sth2);
  if (((mu==0) && (nu==3)) || ((mu==3) && (nu==0))){
    return -2*spin_/r*sth2
      +5./8.*dzetaCS_*spin_/(r2*r2)*(1.+12./7.*1./r+27./10.*1./r2)*sth2;
  }

  return 0.;
}

void ChernSimons::gmunu(double ARGOUT_ARRAY2[4][4], double const IN_ARRAY1[4]) const {
  // Let's make sure the Generic version is called, not the KerrBL one
  Generic::gmunu(ARGOUT_ARRAY2, IN_ARRAY1);
}

int ChernSimons::christoffel(double ARGOUT_ARRAY3[4][4][4], double const IN_ARRAY1[4]) const {
  // Let's make sure the Generic version is called, not the KerrBL one
  return Generic::christoffel(ARGOUT_ARRAY3, IN_ARRAY1);
}

double ChernSimons::christoffel(const double coord[4], const int alpha,
				const int mu, const int nu) const {
  return Generic::christoffel(coord, alpha, mu, nu);
}

void ChernSimons::gmunu_up(double ARGOUT_ARRAY2[4][4], double const IN_ARRAY1[4]) const {
  double g[4][4];
  gmunu(g, IN_ARRAY1);
  matrix4CircularInvert(ARGOUT_ARRAY2, g);
  return;
}

double ChernSimons::gmunu_up(const double * pos, int mu, int nu) const {
  return Generic::gmunu_up(pos, mu, nu);
}

int ChernSimons::diff(const double* coordGen, const double* cst,
		      double* res) const{
  double a2=spin_*spin_;

  //int width=25;//15;
  //int prec=15;//8;

  double rsink=1.+sqrt(1.-a2)+drhor;

  double r = coordGen[1];

  if (r < 0.) {
    cerr << "r= " << r << endl;
    GYOTO_ERROR( "ChernSimons.C: r negative!!!!! the horizon"
		" may have been crossed..." );
  }

  if (r < rsink) {
#   if GYOTO_DEBUG_ENABLED
    GYOTO_DEBUG << "Too close to horizon in ChernSimons::diff at r= " 
		<< r << endl;
#   endif
    return 1;
  }

  double r2 = r*r;
  double r3 = r2*r;
  double r4 = r2*r2;
  double r5=r4*r;
  double ff = 1.-2./r;
  // compute and store efficiently sin, cos, cotan
  double theta=coordGen[2];
  double costheta, sintheta;
  sincos(theta, &sintheta, &costheta);
  double costheta2=costheta*costheta;
  if (sintheta==0.) GYOTO_ERROR("sintheta==0");
  double cotantheta=costheta/sintheta;
  double cotantheta2=cotantheta*cotantheta;
  double cotantheta3=cotantheta2*cotantheta;
  double sin2theta=2.*sintheta*costheta;
  double cos2theta=2.*costheta2-1.;

  double pr=coordGen[5];
  double ptheta=coordGen[6];

  double a3=a2*spin_;

  double Sigma=r2+a2*costheta2;
  if (Sigma==0) GYOTO_ERROR("In ChernSimons::diff(): Sigma==0");
  double Sigmam1=1./Sigma;
  double Sigmam2=Sigmam1*Sigmam1;

  double Delta=r2-2*r+a2;

  double E=cst[1];
  double E2=E*E;
  double L=cst[2];
  double L2=L*L;

  double tmp1=(2.*Delta*Sigma);
  if (tmp1==0)  GYOTO_ERROR("In ChernSimons::diff(): 2.*Delta*Sigma==0");
  double tmp1m1=1./tmp1;

  if (Delta==0) GYOTO_ERROR("In ChernSimons::diff(): Delta==0");

  //NB: equations of motion are independent of Carter constant in this
  //form. However, the dependency of the dynamic on this constant
  //appears when transforming from principal momenta to coordinate
  //derivatives (e.g. p_theta -> thetadot)


  /*
    ---> Standard Kerr equations of geodesics, slow rotation approx.
  */
  res[0] =
    1./(2.*ff*r4)
    *(2.*(
	  r*(-2.*spin_*L+E*r3+a2*E*(2.+r))+a2*E*(a2+r*(-2.+r))*costheta2
	  )
      )
    -
    1./(2.*ff*r4)*a2/r2*(costheta2-1./ff)
    *(2.*(
	  E*r4 // to order zero in a2
	  )
      );// tdot

  res[1] = (ff+a2/r2*(1.-ff*costheta2))*pr; //rdot

  res[2] = 1./r2*(1.-a2/r2*costheta2)*ptheta; //thetadot

  res[3] = -1./(2.*ff*r4)
    *(-2.*(
	   r*(2.*spin_*E+L*(-2.+r))+L*(a2+r*(-2.+r))*cotantheta2
	   )
      )
    +1./(2.*ff*r4)*a2/r2*(costheta2-1./ff)
    *(-2.*(r*(2.*spin_*E+L*(-2.+r))+L*(r*(-2.+r))
	   *cotantheta2)
      )
    ; //phidot

  res[4] = 0.;// ptdot: pt = cst = -E

  double tmp2=r2+a2*costheta2;
  if (tmp2==0) GYOTO_ERROR("r2+a2*costheta2==0");
  double tmp2m2=1./(tmp2*tmp2);

  double tmp3=a2+r*(-2.+r);
  double tmp3_2=tmp3*tmp3;

  res[5] =
    (-1./r4*(r*(r-a2)-a2*(1.-r)*costheta2)+2.*a2*costheta2/r4)*pr*pr
    +1./r3*(1.-2.*a2/r2*costheta2)*ptheta*ptheta
    +(1./(r4*r4*ff*ff)
      *(
	costheta2*a2*E2*r3*(r-4.)
	-2.*r3*spin_*E*L*(4.-3.*r)
	-r2*a2*(L2+2*E2*r*(r-2.))
	-r3*(E2*r3-L2*(r-2.)*(r-2.))
	+L2*cotantheta2*r4*r*ff*ff*(1.+2*a2/(ff*r2))
	)
      -2.*a2/(r4*r4*r2*ff*ff)*(costheta2+1./ff)
      *(
	-r3*(E2*r3-L2*(r-2.)*(r-2.))
	+L2*cotantheta2*r5*ff*ff*(1.+2.*a2/(ff*r2))
	)
      );// prdot

  res[6]=
    -0.5*(a2*sin2theta*ff/r2)*pr*pr
    -0.5*(a2*sin2theta*1./r4)*ptheta*ptheta
    +(
      1./r4
      *(
	L2*r2*cotantheta
	+0.5*L2*(a2+2.*r2+a2*cos2theta)*cotantheta3
	+a2/(ff*r)*(L2*(2.-r)+2.*E2*r2)*costheta*sintheta
	)
      -2.*a2*costheta2/(r4*r2)
      *(
	L2*r2*cotantheta+L2*r2*cotantheta3
	)
      ); // pthetadot

  res[7] = 0.;//pphidot: pphi = cst = L

  /*
    ---> Chern-Simons modifications to 1st order
    Checked Jul 10: OK
  */
  res[0]+=
    1./(2.*ff*r4)*spin_*L*(189.+120.*r+70.*r2)*dzetaCS_/(56.*r2*r2); //tdot

  res[3]+=
    -1./(2.*ff*r4)*spin_*E*(189.+120.*r+70.*r2)*dzetaCS_/(56.*r2*r2); //phidot

  res[5]+=-spin_*E*L*dzetaCS_/(56.*r4*r2*(r-2.))*
    (
     -1323.+36.*r+70.*r2+210.*r2*r
     )
    /
    (
     r2*(r-2.)+2.*a2*costheta2*(r-2.)+2.*a2*r
     );

  res[6]+=a3*E*L*dzetaCS_*
    (189.+120.*r+70.*r2)*costheta*sintheta
    /
    (56.*r4*r4*(a2+r*(r-2.)*(1.+2.*a2/r2*costheta2))); //pthdot

  return 0;
}

void ChernSimons::circularVelocity(double const coor[4], double vel[4],
			      double dir) const {

  if (keplerian_) {
    // If keplerian_ is true, let the generic implementation return
    // the Keplerian velocity instead of the true circular velocity
    Generic::circularVelocity(coor, vel, dir);
    return;
  }

# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG<<"coor=["<<coor[0]<<", "<<coor[1]<<", "<<coor[2]<<", "<<coor[3]
	     <<"], dir="<<dir<<endl;
# endif
  double sinth = sin(coor[2]);
  double coord[4] = {coor[0], coor[1]*sinth, M_PI*0.5, coor[3]};

  vel[1] = vel[2] = 0.;
  double rr=coord[1], r2=rr*rr, r3=r2*rr, r4=r3*rr, r5=r4*rr, aa=spin_,
    a2=aa*aa;
  double fact=-112.*r5+567.*dzetaCS_+300.*rr*dzetaCS_+140.*r2*dzetaCS_;
  // This is Omega_CS:
  vel[3] =
    (
     aa*fact+56.*r5*r2*sqrt(4.*(r3-a2)/r4+a2*fact*fact/(3136.*r5*r5*r4))
     )
    /
    (
     112.*r5*(r3-a2)
     );

  vel[0] = SysPrimeToTdot(coor, vel+1);
  vel[3] *= vel[0];
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG_ARRAY(vel,4);
# endif
}
