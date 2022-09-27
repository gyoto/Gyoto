/*
    Copyright 2011-2016, 2018-2020 Frederic Vincent, Thibaut Paumard

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
#include <iostream>
#include <cstdlib>
#include "GyotoMetric.h"
#include "GyotoFactoryMessenger.h"
#include "GyotoConverters.h"
#include "GyotoProperty.h"
#include "GyotoWorldline.h"
#include <cmath>
#include <string>
#include <sstream>

using namespace std ; 
using namespace Gyoto;

Register::Entry* Metric::Register_ = NULL;

//// Gyoto::Object API
GYOTO_PROPERTY_START(Metric::Generic,
		     "The geometry of space-time at this end of the Universe.")
GYOTO_PROPERTY_DOUBLE_UNIT(Metric::Generic, Mass, mass,
			   "Mass for scaling geometrical units to meters (kg, 1.).")
GYOTO_PROPERTY_BOOL(Metric::Generic, Keplerian, NonKeplerian, keplerian,
		    "Whether to use the Keplerian approximation in circularVelocity().")
GYOTO_PROPERTY_DOUBLE(Metric::Generic, DeltaMin, deltaMin,
		      "Minimum step for Legacy integrator (geometrical units, DBL_MIN).")
GYOTO_PROPERTY_DOUBLE(Metric::Generic, DeltaMax, deltaMax,
		      "Maximum step for Legacy integrator (geometrical units, DBL_MAX).")
GYOTO_PROPERTY_DOUBLE(Metric::Generic, DeltaMaxOverR, deltaMaxOverR,
		      "Max of step/r coordinate for Legacy integrator (geometrical units, 1)")
GYOTO_PROPERTY_END(Metric::Generic, Object::properties)

///
/*
  User code is free to provide any or none of the various versions of
  gmunu_up(). The default implementations call one another to try and
  find user-provided code, but a default implementation needs to be
  triggered if none are provided.  In order to avoid infinite
  recursion as well as for efficiency, several of those methose set a
  flag in __defaultfeatures if they are called to inform the other
  methods. This is what each method will try:

    - coefficient gmunu_up:
      + matrix gmunu_up;
    - matrix gmunu_up:
      + coefficient gmunu_up
      + fall back to inverting gmunu matrix.
 */

#define __default_gmunu_up_coef 1
#define __default_gmunu_up_matrix 2
#define __default_jacobian 4
#define __default_gmunu_up_and_jacobian 8
#define __default_christoffel_coef 16
#define __default_christoffel_matrix 32

Metric::Generic::Generic() :
  SmartPointee(), Object("anonymous metric"), mass_(1.), coordkind_(GYOTO_COORDKIND_UNSPECIFIED),
  __defaultfeatures(0),
  delta_min_(GYOTO_DEFAULT_DELTA_MIN),
  delta_max_(GYOTO_DEFAULT_DELTA_MAX),
  delta_max_over_r_(GYOTO_DEFAULT_DELTA_MAX_OVER_R),
  keplerian_(false)
{
}

Metric::Generic::Generic(const int coordkind, const std::string &name) :
  SmartPointee(), Object(name), mass_(1.), coordkind_(coordkind),
  __defaultfeatures(0),
  delta_min_(GYOTO_DEFAULT_DELTA_MIN),
  delta_max_(GYOTO_DEFAULT_DELTA_MAX),
  delta_max_over_r_(GYOTO_DEFAULT_DELTA_MAX_OVER_R),
  keplerian_(false)
{
# if GYOTO_DEBUG_ENABLED
  GYOTO_IF_DEBUG
  GYOTO_DEBUG_EXPR(coordkind_);
  GYOTO_DEBUG_EXPR(kind_);
  GYOTO_ENDIF_DEBUG
# endif
}

Metric::Generic::Generic(Generic const &o):
  SmartPointee(o), Object(o), mass_(o.mass_), coordkind_(o.coordkind_),
  __defaultfeatures(o.__defaultfeatures),
  delta_min_(o.delta_min_), delta_max_(o.delta_max_),
  delta_max_over_r_(o.delta_max_over_r_), keplerian_(o.keplerian_)
{}

Metric::Generic * Metric::Generic::clone() const {
  string msg = "Metric::Generic::clone() called: cloning not supported for metric kind ";
  msg += kind();
  GYOTO_ERROR (msg);
  return const_cast<Metric::Generic*>(this); // to avoid warning
}

Metric::Generic::~Generic(){
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << endl;
# endif
}

// Output

/***************Definition of the physical scene**************/

void Metric::Generic::mass(const double mas)        {
  mass_=mas;
  tellListeners();
}
void Metric::Generic::mass(const double mas, const string &unit) {
# ifdef GYOTO_DEBUG_ENABLED
  GYOTO_IF_DEBUG
  GYOTO_DEBUG_EXPR(mas);
  GYOTO_DEBUG_EXPR(unit);
  GYOTO_ENDIF_DEBUG
# endif
  mass(Units::ToKilograms(mas, unit));
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << "(mass="<<mas<<", unit=\"" << unit << "\") : mass_="
	      << mass_ <<" kg"<< endl;
# endif
}

int Metric::Generic::coordKind()               const { return coordkind_; }
void Metric::Generic::coordKind(int coordkind)
{
  coordkind_=coordkind;
  tellListeners();
}

double Metric::Generic::mass()                 const { return mass_; }
double Metric::Generic::mass(const string &unit) const {
  return Units::FromKilograms(mass(), unit);
}

double Metric::Generic::deltaMin() const {return delta_min_;}
double Metric::Generic::deltaMax() const {return delta_max_;}
void  Metric::Generic::deltaMin(double h1) {delta_min_=h1;}
void  Metric::Generic::deltaMax(double h1) {delta_max_=h1;}
double Metric::Generic::deltaMaxOverR() const { return delta_max_over_r_;}
void Metric::Generic::deltaMaxOverR(double t) {delta_max_over_r_=t;}

double Metric::Generic::deltaMax(double const pos[8], double h1max) const
{
  double h1max_at_r=abs(pos[1]);
  if (coordkind_==GYOTO_COORDKIND_CARTESIAN) {
    if (abs(pos[2])>h1max_at_r) h1max_at_r=abs(pos[2]);
    if (abs(pos[3])>h1max_at_r) h1max_at_r=abs(pos[3]);
  }
  h1max_at_r *= delta_max_over_r_;
  if (h1max > h1max_at_r) h1max = h1max_at_r;
  if (h1max>delta_max_) h1max=delta_max_;
  if (h1max<delta_min_) h1max=delta_min_;
  return h1max;
}

bool Metric::Generic::keplerian() const {return keplerian_;}
void Metric::Generic::keplerian(bool t) {keplerian_=t;}

double Metric::Generic::SysPrimeToTdot(const double pos[4], const double v[3]) const {
  double sum=0.,xpr[4];
  double g[4][4];
  int i,j;
# if GYOTO_DEBUG_ENABLED
  GYOTO_IF_DEBUG
    GYOTO_DEBUG_ARRAY(pos,4);
    GYOTO_DEBUG_ARRAY(v,3);
  GYOTO_ENDIF_DEBUG
# endif


  xpr[0]=1.; // dt/dt=1;
  for (i=0;i<3;++i) xpr[i+1]=v[i];
  gmunu(g, pos);
  for (i=0;i<4;++i) {
    for (j=0;j<4;++j) {
      sum+=g[i][j]*xpr[i]*xpr[j];
    }
  }
  if (sum>=0) {
    GYOTO_WARNING << "v>c\n";
    return 0.;
  }
  return pow(-sum, -0.5);
}


void Metric::Generic::nullifyCoord(double coord[8]) const {
  double tdot2;
  nullifyCoord(coord, tdot2);
}

void Metric::Generic::nullifyCoord(double coord[8], double& tdot2) const {
  int i, j;
  double a, b=0., c=0.;
  double g[4][4];
  gmunu(g, coord);
  a=g[0][0];
  for (i=1;i<=3;++i) {
    b+=g[0][i]*coord[4+i];
    for (j=1;j<=3;++j) {
      c+=g[i][j]*coord[4+i]*coord[4+j];
    }
  }
  double sDelta=sqrt(b*b-a*c), am1=1./a;
  tdot2=(-b+sDelta)*am1;
  coord[4]=(-b-sDelta)*am1;
}

void Metric::Generic::normalizeFourVel(double coord[8]) const {
  normalizeFourVel(coord, coord+4);
}

void Metric::Generic::normalizeFourVel(double const pos[4],
				       double fourvel[4]) const {
  double tdotin=fourvel[0];
  double threevel[3]={fourvel[1]/tdotin,
		      fourvel[2]/tdotin,
		      fourvel[3]/tdotin};
  fourvel[0]=SysPrimeToTdot(pos, threevel);
  for (int k = 0; k<3; ++k)
    fourvel[k+1]=threevel[k]*fourvel[0];
}

double Metric::Generic::ScalarProd(const double pos[4],
			  const double u1[4], const double u2[4]) const {
  double res=0.;
  double g[4][4];
  gmunu(g, pos);
  for (int i=0;i<4;i++) {
    for (int j=0;j<4;j++) {
      res+=g[i][j]*u1[i]*u2[j];
    }
  }
  return res;
}

double Metric::Generic::norm(const double pos[4],
			     const double u1[4]) const {
  double norm2=ScalarProd(pos, u1, u1);
  return ((norm2>0.)-(norm2<0.))*sqrt(abs(norm2));
}

void Metric::Generic::multiplyFourVect(double vect[4],
				       double a) const {
  for (int k=0; k<4; ++k) vect[k] *= a;
}

void Metric::Generic::addFourVect(double u1[4],
				  double const u2[4]) const {
  for (int k=0; k<4; ++k) u1[k] += u2[k];
}

void Metric::Generic::projectFourVect(double const pos[4],
				      double u1[4],
				      double const u2[4]) const {
  double n=norm(pos, u2);
  double s=(n>0)-(n<0);
  double absinvn=1./abs(n);
  double u[4]={u2[0]*absinvn, u2[1]*absinvn,
	       u2[2]*absinvn, u2[3]*absinvn};
  multiplyFourVect(u, -s*ScalarProd(pos, u, u1));
  addFourVect(u1, u);
}

void Metric::Generic::dualOneForm(double const IN_ARRAY1[4],
				  double const IN_ARRAY2[4],
				  double ARGOUT_ARRAY3[4]) const {
  double g[4][4];
  gmunu(g, IN_ARRAY1);
  for (int nu=0; nu<4; nu++) {
    ARGOUT_ARRAY3[nu]=0.;
    for (int mu=0; mu<4; mu++) {
      ARGOUT_ARRAY3[nu] += g[mu][nu]*IN_ARRAY2[mu];
    }
  }
}


/***************Geodesic most general integration**************/


double Metric::Generic::gmunu(const double x[4], int mu, int nu) const {
  double g[4][4];
  gmunu(g, x);
  return g[mu][nu];
}

void Metric::Generic::gmunu(double g[4][4], const double x[4]) const {
  int mu, nu;
  for (mu=0; mu<4; ++mu) {
    g[mu][mu]=gmunu(x, mu, mu);
    for (nu=mu+1; nu<4; ++nu)
      g[mu][nu]=g[nu][mu]=gmunu(x, mu, nu);
  }
}

double Metric::Generic::gmunu_up(const double x[4], int mu, int nu) const {
  const_cast<Generic*>(this)->__defaultfeatures |= __default_gmunu_up_coef;
  double g[4][4];
  gmunu_up(g, x);
  return g[mu][nu];
}

void Metric::Generic::gmunu_up(double gup[4][4], const double x[4]) const {
  const_cast<Generic*>(this)->__defaultfeatures |= __default_gmunu_up_matrix;
  if (!(__defaultfeatures & __default_gmunu_up_coef)) {
    // if gmunu_up(x, mu, nu) is not the default (or we don't know), use it
    int mu, nu;
    for (mu=0; mu<4; ++mu) {
      gup[mu][mu]=gmunu_up(x, mu, mu);
      for (nu=mu+1; nu<4; ++nu)
	gup[mu][nu]=gup[nu][mu]=gmunu_up(x, mu, nu);
    }
  } else if (!(__defaultfeatures & __default_gmunu_up_and_jacobian)){
    // if gmunu_up_and_jacobian is not the default (or we don't know), use it
    double jac[4][4][4];
    gmunu_up_and_jacobian(gup, jac, x);
  } else {
    // else call g and invert it
    double g[4][4];
    gmunu(g, x);
    Gyoto::matrix4Invert(gup, g);
  }
}

void Metric::Generic::jacobian(double jac[4][4][4], const double x0[4]) const {
  const_cast<Generic*>(this)->__defaultfeatures |= __default_jacobian;

  if (__defaultfeatures & __default_gmunu_up_and_jacobian) {
    // If gmunu_up_and_jacobian is the default, do the same
    // thing. This saves computing gmunu_up needlessly.
    double g0[4][4], gx[4][4], h=1e-7, x[4]={x0[0], x0[1], x0[2], x0[3]};

    gmunu(g0, x0);
    for (int alpha=0; alpha<4; ++alpha) {
      x[alpha]=x0[alpha]+h;
      gmunu(gx, x);
      for (int mu=0; mu<4;mu++) {
	jac[alpha][mu][mu] = (gx[mu][mu]-g0[mu][mu])/h;
	for (int nu=mu+1; nu<4;nu++) {
	  jac[alpha][nu][mu] = jac[alpha][mu][nu] = (gx[mu][nu]-g0[mu][nu])/h;
	}
      }
      x[alpha]=x0[alpha];
    }

  } else {
    double gup[4][4];
    gmunu_up_and_jacobian(gup, jac, x0) ;
  }
}

void Metric::Generic::gmunu_up_and_jacobian(double gup[4][4], double jac[4][4][4], const double x0[4]) const {
  const_cast<Generic*>(this)->__defaultfeatures |= __default_gmunu_up_and_jacobian;

  double g0[4][4];
  if ( ((__defaultfeatures & __default_gmunu_up_coef)
	&&(__defaultfeatures & __default_gmunu_up_matrix))
       || (__defaultfeatures & __default_jacobian) ) {
    // If (all gmunu_up flavors) or (jacobian) is the default
    // implementation , we will need the covariant metric coefficient
    // matrix at this point.
    gmunu(g0, x0);
  }

  if ((__defaultfeatures & __default_gmunu_up_coef)
	&&(__defaultfeatures & __default_gmunu_up_matrix)) {
    // If all flavors of gmunu_up are default implementations,
    // don't call them: invert g0.
    Gyoto::matrix4Invert(gup, g0);
  } else {
    gmunu_up(gup, x0);
  }

  if (__defaultfeatures & __default_jacobian) {
    // If jacobian is the default implementation, derive g numerically.
    double gx[4][4], h=1e-7, x[4]={x0[0], x0[1], x0[2], x0[3]};

    for (int alpha=0; alpha<4; ++alpha) {
      x[alpha]=x0[alpha]+h;
      gmunu(gx, x);
      for (int mu=0; mu<4;mu++) {
	jac[alpha][mu][mu] = (gx[mu][mu]-g0[mu][mu])/h;
	for (int nu=mu+1; nu<4;nu++) {
	  jac[alpha][nu][mu] = jac[alpha][mu][nu] = (gx[mu][nu]-g0[mu][nu])/h;
	}
      }
      x[alpha]=x0[alpha];
    }

    return;
    
  } else {
    jacobian(jac, x0);
  }

 }

void Metric::Generic::computeNBeta(const double coord[4],
					 double &NN,double beta[3]) const
{
  throwError("In Metric::computeNBeta not implemented");
}

double Metric::Generic::christoffel(const double * x, int alpha, int mu, int nu) const {
  const_cast<Generic*>(this)->__defaultfeatures |= __default_christoffel_coef;
  double dst[4][4][4];
  christoffel(dst, x);
  return dst[alpha][mu][nu];

}

int Metric::Generic::christoffel(double dst[4][4][4], const double * x) const {
  const_cast<Generic*>(this)->__defaultfeatures |= __default_christoffel_matrix;
  int a, mu, nu, i;
  if (__defaultfeatures & __default_christoffel_coef) {
    // if christoffel(x, a, mu, nu) is the default implementation,
    // rely on gmunu_up_and_jacobian
    double gup[4][4], jac[4][4][4];
    gmunu_up_and_jacobian(gup, jac, x);
    // computing Gamma^a_mu_nu
    for (a=0; a<4; ++a) {
      for (mu=0; mu<4; ++mu) {
	for (nu=mu; nu<4; ++nu) {
	  dst[a][mu][nu]=0.;
	  for (i=0; i<4; ++i) {
	    dst[a][mu][nu]+=0.5*gup[i][a]*
	      (jac[mu][i][nu]+jac[nu][mu][i]-jac[i][mu][nu]);
	  }
	  if (mu!=nu) dst[a][nu][mu] = dst[a][mu][nu];
	}
      }
    }
  } else {
    // else get the coefficients one by one
    for (a=0; a<4; ++a) {
      for (mu=0; mu<4; ++mu) {
	dst[a][mu][mu]=christoffel(x, a, mu, mu);
	for (nu=mu+1; nu<4; ++nu)
	  dst[a][mu][nu]=dst[a][nu][mu]=christoffel(x, a, mu, nu);
      }
    }
  }
  return 0;
}

/*
Let : Y=[x0,x1,x2,x3,x0_dot,x1_dot,x2_dot,x3_dot] (dot=d/dtau, tau=proper time)
diff is such as : Y_dot=diff(Y)
The general equation of geodesics is used.
 */
int Metric::Generic::diff(const state_t &x,
			  state_t &dxdt,
			  double /* mass */) const {
  if (x.size()<8) GYOTO_ERROR("x should have at least 8 elements");
  if (x.size() != dxdt.size()) GYOTO_ERROR("x.size() should be the same as dxdt.size()");
  if (x[4]<1e-6) return 1;
  int nvec = (x.size()-4)/4;
  dxdt[0]=x[4];
  dxdt[1]=x[5];
  dxdt[2]=x[6];
  dxdt[3]=x[7];
  double dst[4][4][4];
  int retval=christoffel(dst, x.data());
  if (retval) return retval;
  for(int alpha=0; alpha<4; ++alpha) {
    for (int v=1; v<=nvec; ++v)
      dxdt[alpha+4*v]=0.;
    for (int i=0;i<4;i++)
      for (int j=0;j<4;j++)
	for (int v=1; v<=nvec; ++v) 
	  dxdt[alpha+v*4] -= dst[alpha][i][j]*x[4+i]*x[v*4+j];
  }
  return 0;
}

int Metric::Generic::diff31(const state_t &x,
			  state_t &dxdt,
			  double /* mass */) const {
  throwError("In Metric::diff31 not implemented");
}

/*Runge Kutta to order 4

 */

//int Metric::Generic::myrk4(const double y[6], const double* cst , double h, double* res) const{
int Metric::Generic::myrk4(Worldline * line, const state_t &coord, double h, state_t &res) const{
  //cout << "In Metric::Generic::myrk4" << endl;
  size_t sz = coord.size();
  state_t k1(sz) ;
  state_t k2(sz) ;
  state_t k3(sz) ;
  state_t k4(sz) ;
  state_t coord_plus_halfk1(sz) ;
  state_t sixth_k1(sz) ;
  state_t coord_plus_halfk2(sz) ;
  state_t third_k2(sz) ;
  state_t coord_plus_k3(sz) ;
  state_t third_k3(sz) ;
  state_t sixth_k4(sz) ;
  double mass=line->getMass();

  if (diff(coord, k1, mass)) return 1 ;
  
  for (int i=0;i<8;i++) {
    k1[i]=h*k1[i];
    coord_plus_halfk1[i]=coord[i]+0.5*k1[i];
    sixth_k1[i]=1./6.*k1[i];
  }
  
  if (diff(coord_plus_halfk1, k2, mass)) return 1 ;
  for (int i=0;i<8;i++) {
    k2[i]=h*k2[i];
    coord_plus_halfk2[i]=coord[i]+0.5*k2[i];
    third_k2[i]=1./3.*k2[i];
  }
	
  if (diff(coord_plus_halfk2, k3, mass)) return 1 ;
  for (int i=0;i<8;i++) {
    k3[i]=h*k3[i];
    coord_plus_k3[i]=coord[i]+k3[i];
    third_k3[i]=1./3.*k3[i];
  }

  if (diff(coord_plus_k3, k4, mass)) return 1 ;
  for (int i=0;i<8;i++) {
    k4[i]=h*k4[i];
    sixth_k4[i]=1./6.*k4[i];
  }


  for (int i=0;i<8;i++) {
    res[i]=coord[i]+sixth_k1[i]+third_k2[i]+third_k3[i]+sixth_k4[i]; 
  }

  return 0;
}

void Metric::Generic::circularVelocity(double const * coor, double* vel,
				       double dir) const {
  if (!keplerian_) {
    stringstream ss;
    ss << kind_
       << "::circularVelocity() is not implemented. "
       <<"Use \"<Keplerian/>\" for the Keplerian approximation.";
    GYOTO_ERROR(ss.str());
  }

  if (coordkind_==GYOTO_COORDKIND_SPHERICAL) {
    double sinth = sin(coor[2]);
    double coord[4] = {coor[0], coor[1]*sinth, M_PI*0.5, coor[3]};

    vel[1] = vel[2] = 0.;
    vel[3] = 1./(dir*pow(coord[1], 1.5));

    vel[0] = SysPrimeToTdot(coor, vel+1);
    vel[3] *= vel[0];
  } else if (coordkind_==GYOTO_COORDKIND_CARTESIAN) {
    double rcross=sqrt ( coor[1]*coor[1] + coor[2]*coor[2] );
    double Omega=dir*pow(rcross*rcross*rcross, -0.5);
    //angular Keplerian velocity
  
    vel[1] = -coor[2]*Omega;
    vel[2] =  coor[1]*Omega;
    vel[3] = 0.;
    vel[0] = SysPrimeToTdot(coor, vel+1);
    vel[1] *= vel[0];
    vel[2] *= vel[0];
  } else GYOTO_ERROR("Unknownn COORDKIND");
}

void Metric::Generic::zamoVelocity(double const * pos, double* vel) const {
  double ephi[4] = {0., 0., 0., 1.};
  vel [1] = vel [2] = vel [3] = 0.;
  vel [0] = 1.;

  if (coordkind_==GYOTO_COORDKIND_CARTESIAN) {
    // ephi in Cartesian
    double phi=atan2(pos[2], pos[1]);
    double cp, sp;
    sincos(phi, &sp, &cp);
    ephi[0]=0.;
    ephi[1]=-sp;
    ephi[2]=cp;
    ephi[3]=0.;
  }

  projectFourVect(pos, vel, ephi);
  multiplyFourVect(vel, 1./fabs(norm(pos, vel)));
}

void Metric::Generic::cartesianVelocity(double const coord[8], double vel[3]) {
  double tauprime;
  switch(coordkind_) {
  case GYOTO_COORDKIND_SPHERICAL:
    {
      double r = coord[1];
      double costheta = cos(coord[2]), sintheta = sin(coord[2]); 
      double cosphi   = cos(coord[3]), sinphi   = sin(coord[3]);

      tauprime   = 1./coord[4];
      double rprime     = coord[5]*tauprime;
      double thetaprime = coord[6]*tauprime;
      double phiprime   = coord[7]*tauprime;
      vel[0] = rprime * sintheta * cosphi
	+ r * thetaprime * costheta * cosphi
	- r * phiprime * sintheta * sinphi;
      vel[1] = rprime * sintheta * sinphi
	+ r * thetaprime * costheta * sinphi
	+ r * phiprime * sintheta * cosphi;
      vel[2] = rprime * costheta
	- r * thetaprime * sintheta
	;
    }
    break;
  case GYOTO_COORDKIND_CARTESIAN:
	tauprime = 1./coord[4];
	vel[0] = coord[5]*tauprime;
	vel[1] = coord[6]*tauprime;
	vel[2] = coord[7]*tauprime;
	break;
  default:
    GYOTO_ERROR
      ("Metric::Generic::cartesianVelocity: unknown coordinate kind");
  }
}

int Metric::Generic::myrk4_adaptive(Worldline* line, state_t const &coord, double lastnorm , double normref, state_t &coordnew, double h0, double& h1, double h1max) const{ 
  
  double delta0[8];
  double delta0min=1e-15;
  state_t dcoord(coord.size());
  double eps=0.0001;
  double S=0.9;
  double errmin=1e-6;
  double factnorm=2.;
  double mass=line->getMass();

  h1max=deltaMax(coord.data(), h1max);
 
  //cout << "1st diff" << endl;
  diff(coord, dcoord, mass) ;

  for (int i = 0;i<8;i++) delta0[i]=delta0min+eps*(fabs(h0*dcoord[i]));

  double hbis=0.5*h0;
  state_t coordhalf(coord.size());
  state_t coord2(coord.size());
  double delta1[8];
  
  double err;
  int count=0;
  double newnorm;

  /*cout << "coord= ";
  for (int jj=0;jj<8;jj++) {
    cout << coord[jj] << " " ;
  }
  cout << endl;*/
  
  while (1){
    count++;
    //cout << "count in rk Met= " << count << endl;
    err=0.;
    //cout << "then diff" << endl;
    if (
	myrk4(line,coord,h0,coordnew) |
	myrk4(line,coord,hbis,coordhalf)| 
	myrk4(line,coordhalf,hbis,coord2)
	)
      return 1;
    
    //cout << "end then diff" << endl;

    /* cout << "coordnew= ";
    for (int jj=0;jj<8;jj++) {
      cout << coordnew[jj] << " " ;
    }
    cout << endl;*/
    
    for (int i = 0;i<8;i++){
      delta1[i]=coord2[i]-coordnew[i];
      double err_i=fabs(delta1[i]/delta0[i]);
      if (err<err_i) err=err_i;
    }

    if (err>1) {
      h0=S*h0*pow(err,-0.25);
      hbis=0.5*h0;
    }else{
      h1=(err > errmin ? S*h0*pow(err,-0.2) : 4.*h0);//pour éviter les explosions
      if (fabs(h1)<delta_min_) h1=(h0>0.)?delta_min_:-delta_min_;
      if (fabs(h1)>h1max) h1=(h0>0.)?h1max:-h1max;
      
      //Testing tangent vector norm stays next to 0 :
      
      newnorm=ScalarProd(coordnew.data(), coordnew.data()+4, coordnew.data()+4);

      if ( fabs(newnorm-normref) > factnorm*fabs(lastnorm-normref) ) {
	//cout << "norm big!" << endl;
	//cout << "newnorm= " << newnorm << endl;
      	//myrk4(coord,h0/10.,coordnew);
      	//h1/=10.;
      }
      GYOTO_DEBUG << "step used= " << h0 << endl;
      break;
    

    }
  }

  return 0;

}

double Metric::Generic::unitLength() const { 
  return mass_ * GYOTO_G_OVER_C_SQUARE; 
}

double Metric::Generic::unitLength(const string &unit) const { 
  return Units::FromMeters(unitLength(), unit);
}

int Metric::Generic::isStopCondition(double const * const ) const {
  return 0;
}

void Metric::Generic::setParticleProperties(Worldline*, const double*) const {
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << endl;
# endif
}

void Metric::Generic::observerTetrad(obskind_t obskind,
				     double const coord[4], double fourvel[4],
				     double screen1[4], double screen2[4],
				     double screen3[4]) const{
  if (obskind == GYOTO_OBSKIND_ZAMO) {
    zamoVelocity(coord, fourvel);
  } else if (obskind== GYOTO_OBSKIND_KEPLERIAN) {
    circularVelocity(coord, fourvel);
  } else if (obskind != GYOTO_OBSKIND_FULLYSPECIFIED) {
    normalizeFourVel(coord, fourvel);
  }
  if (obskind != GYOTO_OBSKIND_FULLYSPECIFIED)
    observerTetrad(coord, fourvel, screen1, screen2, screen3);

  // No general way to define the tetrad, should be defined
  // in specific metrics. Test below will obviously fail for
  // a machine-initialized tetrad.
  double normtol=1e-10;
  if (fabs(ScalarProd(coord,fourvel,fourvel)+1.)>normtol ||
      fabs(ScalarProd(coord,screen1,screen1)-1.)>normtol ||
      fabs(ScalarProd(coord,screen2,screen2)-1.)>normtol ||
      fabs(ScalarProd(coord,screen3,screen3)-1.)>normtol){
    GYOTO_SEVERE << "In Metric:observerTetrad: observer's local "
		  << "basis is not properly normalized "
		  << "norm-(-1, 1, 1, 1)= "
		  << ScalarProd(coord,fourvel,fourvel)+1. << " "
		  << ScalarProd(coord,screen1,screen1)-1. << " "
		  << ScalarProd(coord,screen2,screen2)-1. << " "
		  << ScalarProd(coord,screen3,screen3)-1. << endl;
  }
  
  if (fabs(ScalarProd(coord,fourvel,screen1))>normtol ||
      fabs(ScalarProd(coord,fourvel,screen2))>normtol ||
      fabs(ScalarProd(coord,fourvel,screen3))>normtol ||
      fabs(ScalarProd(coord,screen1,screen2))>normtol ||
      fabs(ScalarProd(coord,screen1,screen3))>normtol ||
      fabs(ScalarProd(coord,screen2,screen3))>normtol){
    GYOTO_SEVERE << "In Metric:observerTetrad: observer's local "
		 << "basis is not orthogonal" << endl;
  }
}

void Metric::Generic::GramSchmidt(double const pos[4], double u0[4],
				  double u1[4], double u2[4], double u3[4]) const {
  // This is the Gram-Schmidt process according to
  // https://en.wikipedia.org/wiki/Gram%E2%80%93Schmidt_process

  // normalize u0
  multiplyFourVect(u0, 1./abs(norm(pos, u0)));

  // project u1 onto hyperplane othogonal to u0
  // then normalize u1
  projectFourVect(pos, u1, u0);
  multiplyFourVect(u1, 1./abs(norm(pos, u1)));

  // project u2 along u0, then along u1, then normalize it
  projectFourVect(pos, u2, u0);
  projectFourVect(pos, u2, u1);
  multiplyFourVect(u2, 1./abs(norm(pos, u2)));


  // project u3 along u0, u1 and u2 then normalize it
  projectFourVect(pos, u3, u0);
  projectFourVect(pos, u3, u1);
  projectFourVect(pos, u3, u2);
  multiplyFourVect(u3, 1./abs(norm(pos, u3)));

}

void Metric::Generic::observerTetrad(double const pos[4], double fourvel[4],
				     double screen1[4], double screen2[4],
				     double screen3[4]) const{
  // following Krolik & Hawley 2004
  // https://iopscience.iop.org/article/10.1086/427932/fulltext/
  // Start with U, ephi, er, etheta and us Gram-Schmidt orthonormalization
  // Warning, this is not really what Krolik & Hawley did.

  switch(coordkind_) {
  case GYOTO_COORDKIND_SPHERICAL:
    screen1[0]=0.;
    screen1[1]=0.;
    screen1[2]=0.;
    screen1[3]=-1.;

    screen2[0]=0.;
    screen2[1]=0.;
    screen2[2]=-1.;
    screen2[3]=0.;

    screen3[0]=0.;
    screen3[1]=-1.;
    screen3[2]=0.;
    screen3[3]=0.;
    break;

  case GYOTO_COORDKIND_CARTESIAN:
    {
      double rp=sqrt(pos[1]*pos[1]+pos[2]*pos[2]);
      double theta=atan2(rp, pos[3]);
      double phi=atan2(pos[2], pos[1]);
      double ct, st, cp, sp;
      sincos(phi, &sp, &cp);
      sincos(theta, &st, &ct);
      screen1[0]=0.;
      screen1[1]=sp;
      screen1[2]=-cp;
      screen1[3]=0.;

      screen2[0]=0.;
      screen2[1]=-ct*cp;
      screen2[2]=-ct*sp;
      screen2[3]=st;

      screen3[0]=0.;
      screen3[1]=-pos[1];
      screen3[2]=-pos[2];
      screen3[3]=-pos[3];
    }
    break;

  default:
    GYOTO_ERROR
      ("Metric::Generic::observerTetrad: unknown coordinate kind");
  }

  GramSchmidt(pos, fourvel, screen2, screen3, screen1);

}

double Metric::Generic::getRmb() const{
  GYOTO_ERROR("In Metric::getRmb: should be implemented "
	     "in the derived metric");
  return 0.; // silence warning
}

double Metric::Generic::getRms() const{
  GYOTO_ERROR("In Metric::getRms: should be implemented "
	     "in the derived metric");
  return 0.; // silence warning
}

double Metric::Generic::getSpecificAngularMomentum(double rr) const{
  GYOTO_ERROR("In Metric::getSpecificAngularMomentum: should be implemented "
	     "in the derived metric");
  return 0.; // silence warning
}

double Metric::Generic::getPotential(double const pos[4], double l_cst) const{
  GYOTO_ERROR("In Metric::getPotential: should be implemented "
	     "in the derived metric");
  return 0.; // silence warning
}

/***************For SmartPointers**************/

int Metric::Generic::getRefCount() { return SmartPointee::getRefCount(); }

void Metric::initRegister() {
  if (Gyoto::Metric::Register_) delete Gyoto::Metric::Register_;
  Gyoto::Metric::Register_ = NULL;
}

void Gyoto::Metric::Register(std::string name, Metric::Subcontractor_t* scp) {
  Register::Entry* ne =
    new Register::Entry(name,
			(Gyoto::SmartPointee::Subcontractor_t*) scp,
			Gyoto::Metric::Register_);
  Gyoto::Metric::Register_ = ne;
}

GYOTO_GETSUBCONTRACTOR(Metric)
