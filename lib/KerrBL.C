/*
    Copyright 2011 Frederic Vincent, Thibaut Paumard

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
#include "GyotoKerrBL.h"
#include "GyotoWorldline.h"
#include "GyotoError.h"

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

/*
  NB: points delicats de KerrBL:
  - z-axis problems (d'ou interet de KS? a creuser)
  - probleme pres de l'horizon (d'ou interet de KS? a creuser)
  - changement de thetadot dans CheckCons (mais: test dans RK4ada pour
    eviter de "trop" changer)
  - changement de rdot dans Normalize4v (mais: test dans RK4ada pour
    eviter de "trop" changer)
*/

/*
Comment on z-axis problem:
z-axis pb is never launched by the condition theta<Gyoto_min_theta, it is always launched by the derlim, thetaaol tests in diff() (see comments there) 


To prevent any infinite loop, a test is done (see: if (countbis > 50
&& zaxis)....)  which should never be read : the integration is
stopped brutally if it's read.  A message is printed on the screen and
the case must be investigated.  The parameter 50 (=countbis_max) can
also be changed (put higher) in order to try to avoid the manu
militari stop.

[*] see also the same .C just changing delta to pi/4N. With same value of i.
*/
//#define GYOTO_MIN_THETA 1e-5 //seems too big
#define GYOTO_MIN_THETA 1e-10

#define drhor 1e-1
  //stop integration at r_horizon+drhor ; it speeds up a lot if high
  //enough (prevents geodesic from "accumulating" next to horizon) ;
  //but must not be too big not to erase multiple images photons who
  //can travel close to horizon In any case, for debug calculation,
  //chose drhor high enough (1e-1 eg, not 1e-3) If ray-tracing ISCO,
  //be at least sure that r_ISCO > rhor+drhor (eg for drhor=1e-1, it's
  //OK for a<0.999)
					       
KerrBL::KerrBL() :
  Generic(GYOTO_COORDKIND_SPHERICAL), spin_(0.)
{
  setKind("KerrBL");
}

KerrBL::KerrBL(double a, double m) :
  Generic(m, GYOTO_COORDKIND_SPHERICAL), spin_(a) 
{
  //DEBUG!!!
  //spin_=0.;
  //setMass(0.);
  //******************
  setKind("KerrBL");

}

// default copy constructor should be fine 
KerrBL::KerrBL(const KerrBL& gg) : Metric::Generic(gg), spin_(gg.spin_)
{setKind("KerrBL");}
KerrBL * KerrBL::clone () const { return new KerrBL(*this); }



KerrBL::~KerrBL() {
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << endl;
# endif
}

// Output
/*
std::ostream& Gyoto::operator<<( std::ostream& o, const KerrBL& met ) {
  return  met.print(o);
}

std::ostream& KerrBL::print( std::ostream& o) const {
  o << "spin=" << spin_ << ", " ;
  Metric::print(o);
  return o;
}
*/

// Mutators
void KerrBL::setSpin(const double spin) {
  spin_=spin;
  tellListeners();
}

// Accessors
double KerrBL::getSpin() const { return spin_ ; }

//Prograde marginally stable orbit
double KerrBL::getRms() const {
  double aa=spin_, aa2=aa*aa;
  double  z1 = 1. + pow((1. - aa2),1./3.)*(pow((1. + aa),1./3.) + pow((1. - aa),1./3.)); 
  double  z2 = pow(3.*aa2 + z1*z1,1./2.);

  return (3. + z2 - pow((3. - z1)*(3. + z1 + 2.*z2),1./2.));
}

//Prograde marginally bound orbit
double KerrBL::getRmb() const {
  return 2.-spin_+2.*sqrt(1.-spin_);
}

//Computation of metric coefficients
double KerrBL::gmunu(const double * pos, int mu, int nu) const {
  double r = pos[1];
  double sth2, cth2;
  sincos(pos[2], &sth2, &cth2);
  sth2*=sth2; cth2*=cth2;
  double r2=r*r;
  double a2=spin_*spin_;
  double rho2=r2+a2*cth2;
//   if (debug()) {
//     cerr << "DEBUG: KerrBL::gmunu: SPIN="<< spin_ << ", R=" << r << ", THETA=" << theta<< endl;
//   }

  if ((mu==0) && (nu==0)) return -(1.-2.*r/rho2); // 2*r*mass
  if ((mu==1) && (nu==1)) return rho2/(r2-2.*r+a2);
  if ((mu==2) && (nu==2)) return rho2;
  if ((mu==3) && (nu==3))
    return (r2+a2+2.*r*a2*sth2/rho2)*sth2; // 2*r*mass
  if (((mu==0) && (nu==3)) || ((mu==3) && (nu==0)))
    return -2*spin_*r*sth2/rho2; //2*r*mass

  return 0.;
} 

//Computation of metric coefficients in contravariant form
double KerrBL::gmunu_up(const double * pos, int mu, int nu) const {
  double r = pos[1];
  double sth2, cth2, cscth; //sin, cos, cosec theta
  sincos(pos[2], &sth2, &cth2);
  cscth=1/sth2; sth2*=sth2; cth2*=cth2;
  double c2th=2.*cth2-1.; //cos 2*theta
  double r2=r*r;
  double a2=spin_*spin_;
  double rho2=r2+a2*cth2;
  double delta=r2-2.*r+a2;

  if ((mu==0) && (nu==0)) return 
			    -(a2*a2+2.*r2*r2+a2*r*(2.+3.*r)
			      +a2*delta*c2th)
			    / (delta*(a2+2.*r2+a2*c2th));
  if ((mu==1) && (nu==1)) return (r2-2.*r+a2)/rho2;
  if ((mu==2) && (nu==2)) return 1./rho2;
  if ((mu==3) && (nu==3))
    return 
      2.*(r2-2.*r+a2*cth2)*cscth*cscth/(delta*(a2+2.*r2+a2*c2th));
  if (((mu==0) && (nu==3)) || ((mu==3) && (nu==0)))
    return -4*spin_*r/(delta*(a2+2.*r2+a2*c2th)); //2*r*mass

  return 0.;
} 

double KerrBL::christoffel(const double[8],
		   const int, const int, const int) const{
  throwError( "KerrBL.C : should never come here to find christoffel!!" );
  return 0.; // avoid warning
}

// Optimized version
double KerrBL::ScalarProd(const double* pos,
			const double* u1, const double* u2) const {
# if GYOTO_DEBUG_ENABLED
  GYOTO_IF_DEBUG
    GYOTO_DEBUG_ARRAY(pos, 4);
    GYOTO_DEBUG_ARRAY(u1, 4);
    GYOTO_DEBUG_ARRAY(u2, 4);
    GYOTO_DEBUG   << "ScalarProd(pos, u1, u2)="
		  <<  (gmunu(pos,0,0)*u1[0]*u2[0]
		      +gmunu(pos,1,1)*u1[1]*u2[1]
		      +gmunu(pos,2,2)*u1[2]*u2[2]
		      +gmunu(pos,3,3)*u1[3]*u2[3]
		      +gmunu(pos,0,3)*u1[0]*u2[3]
		      +gmunu(pos,3,0)*u1[3]*u2[0])
	 << endl;
  GYOTO_ENDIF_DEBUG
# endif
    
  return gmunu(pos,0,0)*u1[0]*u2[0]
    +gmunu(pos,1,1)*u1[1]*u2[1]
    +gmunu(pos,2,2)*u1[2]*u2[2]
    +gmunu(pos,3,3)*u1[3]*u2[3]
    +gmunu(pos,0,3)*u1[0]*u2[3]
    +gmunu(pos,3,0)*u1[3]*u2[0];
  
}

/*For integration of KerrBL geodesics.
diff is such that : y_dot=diff(y,cst) where cst are constants of motion (mu,E,L,Q in KerrBL)
and y contains [r,theta,phi,t,pr,ptheta] (pr and ptheta are canonical momentum)
and y_dot is [rdot,thetadot,phidot,tdot,prdot,pthetadot]
*/
int KerrBL::diff(const double* coordGen, const double* cst, double* res) const{
  double a=spin_;
  double a2=a*a;

  //int width=25;//15;
  //int prec=15;//8;

  double rsink=1.+sqrt(1.-a2)+drhor;

  double r = coordGen[1] ; 

  if (r < 0.) {
    cerr << "r= " << r << endl;
    throwError( "KerrBL.C : r negative!!!!! the horizon may have been crossed..." );
    
  }

  if (r < rsink) {
#   if GYOTO_DEBUG_ENABLED
    GYOTO_DEBUG << "Too close to horizon in KerrBL::diff at r= " << r << endl;
#   endif
    return 1;
  }

  double r2 = r*r ; 
  double r3 = r2*r ;  

  // compute and store efficiently sin, cos, cotan
  double theta=coordGen[2];
  double costheta, sintheta;
  sincos(theta, &sintheta, &costheta);
  double costheta2=costheta*costheta;
  if (sintheta==0.) throwError("sintheta==0");
  double cotantheta=costheta/sintheta;
  double cotantheta2=cotantheta*cotantheta;
  double cotantheta3=cotantheta2*cotantheta;
  double sin2theta=2.*sintheta*costheta;
  double cos2theta=2.*costheta2-1.;

  double pr=coordGen[5];
  double ptheta=coordGen[6];

  double a3=a2*a;

  double Sigma=r2+a2*costheta2;
  if (Sigma==0) throwError("In KerrBL::diff(): Sigma==0");
  double Sigmam1=1./Sigma;
  double Sigmam2=Sigmam1*Sigmam1;

  double Delta=r2-2*r+a2;

  double E=cst[1];
  double E2=E*E;
  double L=cst[2];
  double L2=L*L;

  double tmp1=(2.*Delta*Sigma);
  if (tmp1==0)  throwError("In KerrBL::diff(): 2.*Delta*Sigma==0");
  double tmp1m1=1./tmp1;

  if (Delta==0) throwError("In KerrBL::diff(): Delta==0");

  //NB: equations of motion are independent of Carter constant in this
  //form. However, the dependency of the dynamic on this constant
  //appears when transforming from principal momenta to coordinate
  //derivatives (e.g. p_theta -> thetadot)

  res[0] = tmp1m1*(2.*(r*(-2.*a*L+E*r3+a2*E*(2.+r))+a2*E*(a2+r*(-2.+r))*costheta2));// tdot

  res[1] = Delta*Sigmam1*pr; //rdot

  res[2] = Sigmam1*ptheta; //thetadot

  res[3] = -tmp1m1*(-2.*(r*(2.*a*E+L*(-2.+r))+L*(a2+r*(-2.+r))*cotantheta2)); //phidot
  
  res[4] = 0.;// ptdot : pt = cst = -E

  double tmp2=r2+a2*costheta2;
  if (tmp2==0) throwError("r2+a2*costheta2==0");
  double tmp2m2=1./(tmp2*tmp2);

  double tmp3=a2+r*(-2.+r);
  double tmp3_2=tmp3*tmp3;

  res[5] =
    -0.5*(2.*(r*(r-a2)-a2*(1.-r)*costheta2)*tmp2m2)*pr*pr
    -0.5*(-2.*r*tmp2m2)*ptheta*ptheta
    +(tmp2m2/tmp3_2
      *(a2*(a2*a2*E2-2.*a3*E*L+2.*a*E*L*r2+E2*r3*(-4.+r)
	    +a2*(L2*(1.-r)+2*E2*r2))*costheta2
	+r*(-r*(a2*a2*E2-2.*a3*E*L+2.*a*E*L*(4.-3.*r)*r
		+a2*(L2+2.*E2*r*(-2.+r))+r*(E2*r3-L2*(-2.+r)*(-2.+r)))
	    +L2*tmp3_2*cotantheta2)));// prdot

  res[6]=
    -0.5*(a2*Delta*sin2theta*Sigmam2)*pr*pr
    -0.5*(a2*sin2theta*Sigmam2)*ptheta*ptheta
    +(
      Sigmam2
      *(
	L2*r2*cotantheta
	+0.5*L2*(a2+2.*r2+a2*cos2theta)*cotantheta3
	+a2*r*(2.*a2*E2-4.*a*E*L+L2*(2.-r)+2.*E2*r2)*costheta*sintheta/Delta
	)
      ); // pthetadot

  res[7] = 0.;//pphi = cst = L

  return 0;
}

void KerrBL::circularVelocity(double const coor[4], double vel[4],
			      double dir) const {
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG<<"coor=["<<coor[0]<<", "<<coor[1]<<", "<<coor[2]<<", "<<coor[3]
	     <<"], dir="<<dir<<endl;
# endif
  double sinth = sin(coor[2]);
  double coord[4] = {coor[0], coor[1]*sinth, M_PI*0.5, coor[3]};

  vel[1] = vel[2] = 0.;
  vel[3] = 1./((dir*pow(coord[1], 1.5) + spin_)*sinth);
  vel[0] = SysPrimeToTdot(coor, vel+1);
  vel[3] *= vel[0];
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG_ARRAY(vel,4);
# endif
}

//Runge Kutta to order 4
int KerrBL::myrk4(Worldline * line, const double coordin[8],
		  double h, double res[8]) const
{
  
  /*
    For external use only (i.e. from WlIntegState::nextstep) ; coor
    must be [t,r,th,ph,tdot,rdot,thdot,phdot]
   */
  
  /*
    Returns 0 if everything is fine, 1 is theta is too close to 0, 2
    if a value of r < horizon is used.
   */
  
  /*Switch BL -> principal momenta*/
  double coor[8], res_mom[8] ;
  double const * const cst = line -> getCst();
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG_ARRAY(cst,3);
#endif
  MakeMomentum(coordin,cst,coor);

  double k1[8] ; 
  double k2[8] ; 
  double k3[8] ; 
  double k4[8] ; 
  double coor_plus_halfk1[8] ; 
  double sixth_k1[8] ; 
  double coor_plus_halfk2[8] ; 
  double third_k2[8] ; 
  double coor_plus_k3[8] ; 
  double third_k3[8] ; 
  double sixth_k4[8] ; 
	  
  if (fabs(fmod(coor[2]+M_PI/2, M_PI)-M_PI/2) < GYOTO_MIN_THETA){
    return 1;
  }

  if (diff(coor,cst, k1)) //if diff returns 1, the geodesic is too
			  //close to horizon --> stop integration with
			  //error code 2
    return 2; 

  for (int i=0;i<8;i++) {
    k1[i]=h*k1[i];
    coor_plus_halfk1[i]=coor[i]+0.5*k1[i];
    sixth_k1[i]=1./6.*k1[i];
  }
  
  if (fabs(fmod(coor_plus_halfk1[2]+M_PI/2, M_PI)-M_PI/2) < GYOTO_MIN_THETA){
    return 1;
  }

  if (diff(coor_plus_halfk1,cst,k2))
    return 2; 

  for (int i=0;i<8;i++) {
    k2[i]=h*k2[i];
    coor_plus_halfk2[i]=coor[i]+0.5*k2[i];
    third_k2[i]=1./3.*k2[i];
  }
	
  if (fabs(fmod(coor_plus_halfk2[2]+M_PI/2, M_PI)-M_PI/2) < GYOTO_MIN_THETA){
    return 1;
  }

  if (diff(coor_plus_halfk2,cst,k3))
    return 2 ;

  for (int i=0;i<8;i++) {
    k3[i]=h*k3[i];
    coor_plus_k3[i]=coor[i]+k3[i];
    third_k3[i]=1./3.*k3[i];
  }

  if (fabs(fmod(coor_plus_k3[2]+M_PI/2, M_PI)-M_PI/2) < GYOTO_MIN_THETA){
    return 1;
  }

  if (diff(coor_plus_k3,cst,k4))
    return 2 ;

  for (int i=0;i<8;i++) {
    k4[i]=h*k4[i];
    sixth_k4[i]=1./6.*k4[i];
  }

  for (int i=0;i<8;i++)
    res_mom[i]=coor[i]+sixth_k1[i]+third_k2[i]+third_k3[i]+sixth_k4[i];
  
  /*Switch principal momenta -> BL  */
  MakeCoord(res_mom,cst,res);
 
  return 0;
}

int KerrBL::myrk4(const double coor[8], const double cst[5],
		  double h, double res[8]) const{
  //  cout << "In KerrBL::myrk4 - 2" << endl;
  /*
    For internal use only (i.e. from adaptive rk4) ; coor must be
    [t,r,th,ph,pt,pr,pth,pph]
   */
    
  /*
    Returns 0 if everything is fine, 1 is theta is too close to 0, 2
    if a value of r < horizon is used.
   */
  
  double k1[8] ; 
  double k2[8] ; 
  double k3[8] ; 
  double k4[8] ; 
  double coor_plus_halfk1[8] ; 
  double sixth_k1[8] ; 
  double coor_plus_halfk2[8] ; 
  double third_k2[8] ; 
  double coor_plus_k3[8] ; 
  double third_k3[8] ; 
  double sixth_k4[8] ; 
  double derlim_hor=1e5, derlim_gen=1e6, derlim;
     // throw a "z-axis problem" (return 1) if prdot or pthdot becomes higher
  double aa=spin_;
  double rhor=1.+sqrt(1.-aa*aa), factrtol=5.;
  double thetatol_hor=1e-1, thetatol_gen=1e-3, thetatol;
  /*
    If theta is closer to 0 than thetatol and if the derivatives
    returned by diff are higher than derlim, a z-axis problem is
    launched (return 1) which leads to increasing the integration step
    in adaptive RK4 and stride across the z axis.

    What "close" means depends on whether we are close or not to the
    horizon. When close to the horizon, thetatol should be bigger to
    prevent from infinite integration (the worldline "accumulating"
    towards theta=0[pi]) and derlim should be smaller to prevent error
    in diff (typically r<0 error).

    These factrtol, thetatol_* and derlim_* parameters are to be
    played with to smooth the integration.
   */

  if (coor[1] < factrtol*rhor) {
    thetatol=thetatol_hor;derlim=derlim_hor;
  }else{
    thetatol=thetatol_gen;derlim=derlim_gen;
  }

  double thetacompare = fabs(fmod(coor[2]+M_PI/2, M_PI)-M_PI/2);

  /*  if (thetacompare < GYOTO_MIN_THETA) {
    return 1; // will lead to z-axis problem in RK4 adaptive
    }*/
  //I think these tests are not useful, the derivatives tests on
  //k1,...,k4 are stronger (and necessary!)
 
  if (diff(coor,cst, k1)) // Too close to horizon
    return 2; 

  // NB: if derivatives are high close to z-axis -> throw z-axis pb by
  // returning 1;
  if (   (thetacompare < thetatol)
	 && (fabs(k1[5]) > derlim || fabs(k1[6]) > derlim) ){
    return 1;
  }

  for (int i=0;i<8;i++) {
    k1[i]=h*k1[i];
    coor_plus_halfk1[i]=coor[i]+0.5*k1[i];
    sixth_k1[i]=1./6.*k1[i];
  }

  if (diff(coor_plus_halfk1,cst,k2))
    return 2;

  if (   (thetacompare < thetatol)
	 && (fabs(k2[5]) > derlim || fabs(k2[6]) > derlim) ){
    return 1;
  }

  for (int i=0;i<8;i++) {
    k2[i]=h*k2[i];
    coor_plus_halfk2[i]=coor[i]+0.5*k2[i];
    third_k2[i]=1./3.*k2[i];
  }
	
  if (diff(coor_plus_halfk2,cst,k3))
    return 2 ;

  if (   (thetacompare < thetatol) 
	 && (fabs(k3[5]) > derlim || fabs(k3[6]) > derlim) ){
    return 1;
  }


  for (int i=0;i<8;i++) {
    k3[i]=h*k3[i];
    coor_plus_k3[i]=coor[i]+k3[i];
    third_k3[i]=1./3.*k3[i];
  }

  if (diff(coor_plus_k3,cst,k4))
    return 2 ;

  if (   (thetacompare < thetatol)
	 && (fabs(k4[5]) > derlim || fabs(k4[6]) > derlim) ){
    return 1;
  }

  for (int i=0;i<8;i++) {
    k4[i]=h*k4[i];
    sixth_k4[i]=1./6.*k4[i];
  }

  for (int i=0;i<8;i++)
    res[i]=coor[i]+sixth_k1[i]+third_k2[i]+third_k3[i]+sixth_k4[i];
  
  return 0;
}

int KerrBL::myrk4_adaptive(Worldline * line, const double coordin[8],
			   double , double , double coordout1[8],
			   double h0, double& h1) const
{
  
  /*Switch BL -> principal momenta*/
  double coor[8], coor1[8], cstest[5], coorhalf[8], coor2[8],
    coor1bis[8], mycoor[8], delta1[8];
  double const * const cst = line -> getCst();
  MakeMomentum(coordin,cst,coor);

  double delta0[8], dcoor[8];
  double delta0min=1e-15, eps=0.0001, S=0.9, errmin=1e-6, hbis=0.5*h0,
    err, h1min=0.01, h1max=coor[1]*0.5, diffr, diffth, difftol=0.01, normtemp,
    cstol_gen=1e-3, cstol_hor=1e-2, cstol, div, QCarter;
  int countbis=0, countbislim=50, zaxis=0; // for z-axis problem in myrk4
  //int norm1=0, normhalf=0, norm2=0, rk1=0, rkhalf=0, rk2=0, update, makerr=0.;
  int norm1=0, rk1=0, rkhalf=0, rk2=0, update, makerr=0;
  double a=spin_, factrtol=3.,
    rtol=factrtol*(1.+sqrt(1.-a*a)), rlimitol=10.;

  if (coor[1] < rtol) cstol = cstol_hor;
        // for tests of cst of motion conservation; don't ask too much
        // if near horizon...
  else cstol = cstol_gen;

  if (fabs(fmod(coor[2]+M_PI/2, M_PI)-M_PI/2) < GYOTO_MIN_THETA) {
    GYOTO_WARNING << "Too close to Z axis: stopping integration"<< endl;
    return 1;
  }

  if (diff(coor,cst,dcoor)) return 1 ;

  for (int i = 0;i<8;i++) delta0[i]=delta0min+eps*(fabs(h0*dcoor[i]));

  while (1){
   
    err=0.;
    countbis++;

    //*** z-axis problem ***

    /*
      This while loop tests whether theta is too close to zero
      (coordinate singularity in BL) and whether we are outside
      horizon.
    */
    while
      ( (( (rk1   =myrk4(coor,cst,h0,coor1)      ) == 1) || (rk1    == 2)) ||
        (( (rkhalf=myrk4(coor,cst,hbis,coorhalf) ) == 1) || (rkhalf == 2)) ||
        (( (rk2   =myrk4(coorhalf,cst,hbis,coor2)) == 1) || (rk2    == 2))   ) 
      { 
	
	if (rk1==2 || rkhalf==2 || rk2==2){
	  //if any of rks is 2, stop integration, we're inside horizon
	  return 1;
	}
	
	zaxis=1;
	h0*=1.1; hbis*=1.1;

	GYOTO_INFO << "NOTE: Passing close to z-axis at theta= "
		   << coor[2] << " and r= " << coor[1]
		   << ", jumping ahead with h0= " << h0 << endl;
	
	//throwError("stop");
        
      }

    
    if (countbis > countbislim && zaxis) {
      // If the "jumping ahead" trick above does not work and leads to
      // infinite jumping ahead behavior (should never happen of
      // course)
      
      GYOTO_INFO << "WARNING: " << endl
		 << "in KerrBL.C couldn't solve z-axis problem ; stopping..."
		 << endl;
      return 1;
    }

    //*** Error determination: ***
    
    for (int i = 0;i<8;i++){
      delta1[i]=coor2[i]-coor1[i];
      if ((err<fabs(delta1[i]/delta0[i]))) {
	err=fabs(delta1[i]/delta0[i]);
      }
    }
    
    //*** What next, as a function of error value: ***
    
    if (err>1) {
      h0=S*h0*pow(err,-0.25);
      hbis=0.5*h0;
    }else{
      h1=(err > errmin ? S*h0*pow(err,-0.2) : 4.*h0);//pour éviter les explosions
      if (fabs(h1)<h1min) {
	h1= (h1>0)?h1min:-h1min;
      }
      if (fabs(h1)>h1max) h1=(h1>0.)?h1max:-h1max;

      //*** Normalizing ***
      update=1;
      norm1=CheckCons(coor1,cst,coor1bis); 

      // pr and ptheta relative difference (NB: only pr and ptheta are modified in CheckCons)
      if (coor1[5]) diffr = fabs(coor1[5]-coor1bis[5])/fabs(coor1[5]);
      else if (coor1bis[5]) diffr = fabs(coor1[5]-coor1bis[5])/fabs(coor1[5]);
      else diffr = 0.;
      if (coor1[6]) diffth = fabs(coor1[6]-coor1bis[6])/fabs(coor1[6]);
      else if (coor1bis[6]) diffth = fabs(coor1[6]-coor1bis[6])/fabs(coor1[6]);
      else diffth = 0.;

      if ((diffr > difftol || diffth > difftol)) {
	// don't update coordinates if the relative differences are
	//	more than 1% --> consequence = norm and Carter cst
	//	won't be exactly the same for this step --> below,
	//	test to be sure they're not "too much violated"; NB:
	//	if this test is not performed, the "corrected"
	//	worldline can diverge from the "true" one after a long
	//	time of integration.
	update=0;
	MakeCoord(coor1,cst,mycoor);
	normtemp = ScalarProd(mycoor, mycoor+4, mycoor+4);
	computeCst(mycoor,cstest);
	if (cst[3]>cstol) {
	  QCarter=cst[3];
	  div = cst[4]; // cst[4] = cst[3]==0? 1. : 1./cst[3]
	} else {
	  QCarter=0.;
	  div=1.;
	}
	if ( fabs(normtemp+cst[0])>cstol ) makerr=1; // cst[0] == -real_norm
	if ( makerr && (fabs(cstest[3]-QCarter)*div>cstol) ) makerr=3;
	if ( fabs(cstest[3]-QCarter)*div>cstol ) makerr=2;
	
	if (makerr) {
	  if (verbose() >= GYOTO_SEVERE_VERBOSITY) {
	    cerr << "WARNING:" << endl;
	    if (makerr==1)
	      cerr << "Real norm, current norm= " << (cst[0]?-1.:0.) << " " << normtemp << endl;
	    else if (makerr==2){
	      cerr << "Carter cst error= (" 
		   << QCarter << "-" << cstest[3] << ")*" << div << "*100.= "
		   << fabs(QCarter-cstest[3])*div*100. << " %, cstol=" << cstol
		   << endl;
	    }else{
	      cerr << "Real norm, current norm= " << (cst[0]?-1.:0.) << " " 
		   << normtemp << endl
		   << "Carter cst error= (" 
		   << QCarter << "-" << cstest[3] << ")*" << div << "*100.= "
		   << fabs(QCarter-cstest[3])*div*100. << " %"
		   << endl;
	    }
	  }
	  if (coor1[1]<rlimitol) {
#           if GYOTO_DEBUG_ENABLED
	    GYOTO_DEBUG << "Probable cause of warning:"
			<< "z-axis problem badly treated in "
			<< "KerrBL::myrk4_adaptive" << endl;
#           endif
	    // some rare cases can end up with bad cst conservation
	    // even at r = a few rhor...
	  }else{
	    GYOTO_SEVERE << "This warning occured at r= " << coor1[1] << endl
			 << "i.e. far from horizon --> to be investigated"
			 << ", or maybe increase parameter cstol" 
			 << "in KerrBL.C" << endl;
	  }
	}

      }
	
      //Update coord
      if (update && !norm1){ // norm1=1 if impossible to normalize in CheckCons due to z-axis pb
	for (int i=0;i<8;i++) coor1[i]=coor1bis[i];
      }
      //      cout << "KerrBL Used h0= " << h0 << endl;
      //*** Switch principal momenta -> BL: ***
      
      MakeCoord(coor1,cst,coordout1);
      break;
      
    } //err>1 if-loop end
    
    
  } // while loop end
  
  return 0;
}

int KerrBL::CheckCons(const double coor_init[8], const double cst[5], double coor_fin[8]) const {
  /*
    Ensures that the cst of motion are conserved.
    E and L always are (see diff). But Q and norm may not be.
  */

  double mycoor[8];
  //int thetalert=0;

  /* *** From momenta to normal BL coordinates *** */
  
  MakeCoord(coor_init,cst,mycoor); //Computes mycoor=[t,r,theta,phi,tdot,rdot,thetadot,phidot] from coor_init=[t,r,theta,phi,pt=-E,pr,ptheta,pphi=L] and cst


  /*
    *** Carter constant's check ***
    As the equations of motion (cf diff) are independent of Q, it is necessary to check whether this constant is conserved.
   */

  double argsqrt, limarg=1e-5, limargbis=0.1;
  double costh, sinth, a2=spin_*spin_;
  sincos(mycoor[2], &sinth, &costh);
  double sinthm2=1./(sinth*sinth), costh2=costh*costh;
  double mu=cst[0], EE=cst[1], LL=cst[2], QQ=cst[3], QQm1=cst[4];
  double mu2=mu*mu, EE2=EE*EE, LL2=LL*LL;
  double Sigma=mycoor[1]*mycoor[1]+a2*costh2;
  double Sigma2=Sigma*Sigma;
  double Qtest=
    Sigma2*mycoor[6]*mycoor[6]+costh2*(a2*(mu2-EE2)+LL2*sinthm2);
  //this should be equal to constant QQ
  int thdotpos=1;
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG
    << "mu="<<mu<<", EE="<<EE<<", LL="<<LL<<", QQ="<<QQ<<", QQm1="<<QQm1
    << ", Qtest="<<Qtest<< ", fabs(Qtest-QQ)/QQm1="<< fabs(Qtest-QQ)/QQm1
    << endl;
# endif
  if (fabs(Qtest-QQ)*QQm1 > 1e-6){//Then change thetadot to allow Qtest=QQ
                                  // Note: if Q==0, QQm1==1.

    if (mycoor[6]<0.) thdotpos=0; //to preserve the sign of thetadot
    argsqrt=QQ-costh2*(a2*(mu2-EE2)+LL2*sinthm2);//thetadot should be the sqrt of this quantity

    if (argsqrt<0 && fabs(argsqrt)<=limarg) {//if argsqrt is <0, but "small enough", put it to zero
      argsqrt=0.;
    }else if (argsqrt<0 && fabs(argsqrt)>limarg) {//if it is <0 and "big", two cases: either we are too close to z-axis which causes theta coordinate to behave badly (then return 1: "z-axis problem" in RK4 adaptive); or Error
      if (fabs(fmod(coor_init[2]+M_PI/2, M_PI)-M_PI/2) < M_PI/50.) {
	return 1; // z-axis problem 
      }else{
	if (fabs(argsqrt)>limargbis)
	  throwError("In KerrBL::CheckCons Impossible to determine thetadot; "
		     "maybe try to increase parameter limarg");
	GYOTO_INFO << "KerrBL::CheckCons argsqrt= " << argsqrt << " at theta= " << coor_init[2] << ". Putting it to 0..." << endl;
	argsqrt=0.;
      }
    }
    
    mycoor[6]=sqrt(argsqrt)/Sigma;//update thetadot to impose QQ constant
    
    if (!thdotpos) mycoor[6]*=-1.;
  }
  
  //*** Normalizing 4v: ***

  Normalize4v(mycoor,mu);

  //*** Back to generalized momenta: ***

  MakeMomentum(mycoor,cst,coor_fin);

  return 0;
}

void KerrBL::Normalize4v(double coord[8], const double part_mass) const {
  /*
    Changes rdot to allow norm conservation.
  */
  
  double aa=spin_;
  double rhor=1.+sqrt(1.-aa*aa);
  
  double gtt=gmunu(coord, 0, 0), gtph=gmunu(coord, 0, 3), grr=gmunu(coord, 1, 1), gthth=gmunu(coord, 2, 2), gphph=gmunu(coord, 3, 3);
  double rr=coord[1], tdot=coord[4], rdot=coord[5], thdot=coord[6], phdot=coord[7];
  
  int valuepos;//to preserve rdot sign
  if (rdot>0.) valuepos=1;
  else valuepos=0;
  
  if (part_mass==0.){// *** ZERO - NORM CASE ***
    
    double argrac=-(gtt*tdot*tdot+2.*gtph*phdot*tdot+gthth*thdot*thdot+gphph*phdot*phdot)/grr;//rdot should be equal to the sqrt of this quantity
    double arglim=1e-4;//generates lots of warning with argrac /approx 1e-5...
    if (argrac<0 && fabs(argrac)<arglim) {//if the quantity is <0 but "small enough" --> put it to 0
      argrac=0.;
    }
    if (argrac<0 && fabs(argrac)>arglim) {
      if (rr/rhor < 2.) {//A AFFINER?? //if the quantity is <0 and "big", but we are close to horizon --> put it to zero with Warning message
	if (verbose() >= GYOTO_WARNING_VERBOSITY) {
	  GYOTO_WARNING << "0-NORM CLOSE TO HORIZON : "
               << "in KerrBL::Normalize4v impossible to normalize 0-mass "
	       << "particule next to horizon. Putting argrac to 0. "
	       << "Effective value of argrac= " << argrac << endl
	       << "with coord= ";
	  for (int myind=0;myind<8;myind++) cerr << coord[myind] << " ";
	  cerr << endl;
	}
	argrac=0.;
      }else{//if the quantity is <0 and "big", and we are not close to horizon --> error
	throwError( "In KerrBL::Normalize4v impossible to normalize 0-mass particle outside horizon!" );
      }	  
    }
    
    coord[5]=sqrt(argrac);//update rdot
    if (!valuepos) coord[5]*=-1.;
    
  }else if (part_mass>0.){ // *** -1 - NORM CASE ***
 
    double argrac=-(1.+gtt*tdot*tdot+2.*gtph*phdot*tdot+gthth*thdot*thdot+gphph*phdot*phdot)/grr;
    double arglim=1e-7;
    if (argrac<0 && fabs(argrac)<arglim) argrac=0.;
    if (argrac<0 && fabs(argrac)>arglim) {
      if (rr/rhor < 2.) {//A AFFINER??
	if (verbose()>=GYOTO_WARNING_VERBOSITY) {
	  cerr << "WARNING -1 - NORM CLOSE TO HORIZON : "
	       << "in KerrBL::Normalize4v impossible to normalize massive "
	       << "particle next to horizon. Putting argrac to 0. "
	       << "Effective value of argrac= " << argrac << endl
	       << "with coord= ";
	  for (int myind=0;myind<8;myind++) cerr << coord[myind] << " ";
	  cerr << endl;
	}
	argrac=0.;
      }else{
	throwError( "In KerrBL::Normalize4v impossible to normalize massive particule outside horizon!" );
      }
    }
    
    coord[5]=sqrt(argrac);//update rdot
    if (!valuepos) coord[5]*=-1.;
    
  }else{
    throwError("In KerrBL::Normalize4v: negative mass!");
  }

}


void KerrBL::MakeCoord(const double coordin[8], const double cst[5], double coord[8]) const {

  double tt=coordin[0], rr = coordin[1], theta=coordin[2], phi=coordin[3], 
    pr=coordin[5], ptheta=coordin[6];
 
  double r2 = rr*rr ;
 
  double sinth, costh;
  sincos(theta, &sinth, &costh);
  double costheta2=costh*costh, sintheta2=sinth*sinth, 
    aa=spin_, a2=aa*aa;

  double Sigma=r2+a2*costheta2, Delta=r2-2*rr+a2, lambda=1.-2*rr/Sigma, 
    xi=2.*aa*rr*sintheta2/Sigma, xi2=xi*xi, 
    gamma = sintheta2*(r2+a2+2.*a2*rr*sintheta2/Sigma), 
    fact=1./(gamma*lambda+xi2);

  double EE=cst[1], LL=cst[2];
    
  double rdot=Delta/Sigma*pr, thetadot=1./Sigma*ptheta, 
    phidot=lambda*fact*LL+xi*fact*EE, tdot=-xi*fact*LL+gamma*fact*EE;

  coord[0]=tt;coord[1]=rr;coord[2]=theta;coord[3]=phi;coord[4]=tdot;coord[5]=rdot;coord[6]=thetadot;coord[7]=phidot;
 
}

void KerrBL::MakeMomentum(const double coord[8], const double cst[5], double coordout[8]) const{

  double EE=cst[1], LL=cst[2];
  
  // double tt=coord[0], rr = coord[1], theta=coord[2], phi=coord[3], tdot=coord[4],
  //  rdot=coord[5], thetadot=coord[6], phidot=coord[7];
  double tt=coord[0], rr = coord[1], theta=coord[2], phi=coord[3],
    rdot=coord[5], thetadot=coord[6];

  double r2 = rr*rr, costheta2=cos(theta); costheta2*=costheta2;
    //    sintheta2=sin(theta)*sin(theta);
  
  double aa=spin_, a2=aa*aa;

  double Sigma=r2+a2*costheta2, Delta=r2-2*rr+a2;

  double pr=Sigma/Delta*rdot, ptheta=Sigma*thetadot;

  coordout[0]=tt;coordout[1]=rr;coordout[2]=theta;coordout[3]=phi;coordout[4]=-EE;coordout[5]=pr;coordout[6]=ptheta;coordout[7]=LL;
}

void KerrBL::nullifyCoord(double coord[8]) const {
  double tdot2;
  nullifyCoord(coord,tdot2);
}

void KerrBL::nullifyCoord(double coord[4], double & tdot2) const {

  int i;
  double a, b=0., c=0.;
  
  a=gmunu(coord, 0, 0);
  b=gmunu(coord, 0, 3)*coord[7];

  for (i=1;i<=3;++i){
    c+=gmunu(coord, i, i)*coord[4+i]*coord[4+i];
  }
  
  double sDelta=sqrt(b*b-a*c), am1=1./a;
  tdot2=(-b+sDelta)*am1;
  coord[4]=(-b-sDelta)*am1;
}

void KerrBL::computeCst(const double coord[8], double cst[5]) const{
  
  //double tt=coord[0], rr = coord[1], theta=coord[2], phi=coord[3], tdot=coord[4],
  //  rdot=coord[5], thetadot=coord[6], phidot=coord[7];
  double rr = coord[1], theta=coord[2], tdot=coord[4],
    thetadot=coord[6], phidot=coord[7];
  
  double norm=ScalarProd(coord, coord+4, coord+4);
  
  double sinth, costh;
  sincos(theta, &sinth, &costh);
  double r2 = rr*rr, costheta2=costh*costh, 
    sintheta2=sinth*sinth;
  
  double a2=spin_*spin_;
  
  double Sigma=r2+a2*costheta2, fact=2.*spin_*rr*sintheta2/Sigma;
  
  double mu;//Particule mass: 0 or 1
  if (fabs(norm)<fabs(norm+1.)){
    mu=0.;
  }else{
    mu=1.;
  }
  //See e.g. Levin&PerezGiz 07 and MTW
  double EE=(1-2*rr/Sigma)*tdot+fact*phidot, 
                   //OK for particule mass = 0 or 1 
    LL=sintheta2*(r2+a2+spin_*fact)*phidot-fact*tdot, 
                   //OK for particule mass = 0 or 1
    QQ=Sigma*thetadot*Sigma*thetadot+costheta2
                           *(a2*(mu*mu-EE*EE)+LL*LL/sintheta2); 
                   //different for a 0-mass and a 1-mass particule
  
  cst[0]=mu;cst[1]=EE;cst[2]=LL;cst[3]=QQ;
  cst[4]= cst[3]==0. ? 1. : 1./cst[3]; // cache 1/Q or 1. if Q==0
}

void KerrBL::setParticleProperties(Worldline * line, const double* coord) const
{
  double cst[5];
  computeCst(coord,cst);
  line -> setCst(cst,5);
}

#ifdef GYOTO_USE_XERCES
void KerrBL::fillElement(Gyoto::FactoryMessenger *fmp) {
  fmp -> setParameter("Spin", spin_);
  Metric::Generic::fillElement(fmp);
}

void KerrBL::setParameter(string name, string content, string unit) {
  if(name=="Spin") setSpin(atof(content.c_str()));
  else Generic::setParameter(name, content, unit);
}

#endif
