/*
    Copyright 2011-2016, 2018-2019 Frederic Vincent, Thibaut Paumard

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

Metric::Generic::Generic(const int coordkind, const std::string &name) :
  SmartPointee(), Object(name), mass_(1.), coordkind_(coordkind),
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

const string Metric::Generic::kind() const {return kind_;}
void Metric::Generic::kind(const string src) { kind_ = src;}

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


double Metric::Generic::christoffel(const double * x, int alpha, int mu, int nu) const {
  double dst[4][4][4];
  christoffel(dst, x);
  return dst[alpha][mu][nu];
}

int Metric::Generic::christoffel(double dst[4][4][4], const double * x) const {
  int alpha, mu, nu;
  for (alpha=0; alpha<4; ++alpha) {
    for (mu=0; mu<4; ++mu) {
      dst[alpha][mu][mu]=christoffel(x, alpha, mu, mu);
      for (nu=mu+1; nu<4; ++nu)
	dst[alpha][mu][nu]=dst[alpha][nu][mu]=christoffel(x, alpha, mu, nu);
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
			  state_t &dxdt) const {
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

/*Runge Kutta to order 4

 */

//int Metric::Generic::myrk4(const double y[6], const double* cst , double h, double* res) const{
int Metric::Generic::myrk4(Worldline * /* line */ , const state_t &coord, double h, state_t &res) const{
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
  
  if (diff(coord, k1)) return 1 ; 
  
  for (int i=0;i<8;i++) {
    k1[i]=h*k1[i];
    coord_plus_halfk1[i]=coord[i]+0.5*k1[i];
    sixth_k1[i]=1./6.*k1[i];
  }
  
  if (diff(coord_plus_halfk1,k2)) return 1 ; 
  for (int i=0;i<8;i++) {
    k2[i]=h*k2[i];
    coord_plus_halfk2[i]=coord[i]+0.5*k2[i];
    third_k2[i]=1./3.*k2[i];
  }
	
  if (diff(coord_plus_halfk2,k3)) return 1 ;
  for (int i=0;i<8;i++) {
    k3[i]=h*k3[i];
    coord_plus_k3[i]=coord[i]+k3[i];
    third_k3[i]=1./3.*k3[i];
  }

  if (diff(coord_plus_k3,k4)) return 1 ; 
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

void Metric::Generic::zamoVelocity(double const * coor, double* vel) const {
  static bool needwarning = true;
  if (needwarning) {
    GYOTO_WARNING <<
      "Trivial implementation of Metric::Generic::zamoVelocity() may be wrong for you Metric\n";
    needwarning = false;
  }
  vel [1] = vel [2] = vel [3] = 0.;
  vel [0] = SysPrimeToTdot(coor, vel+1) ;
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
	+ r * phiprime * cosphi;
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


  h1max=deltaMax(coord.data(), h1max);
 
  //cout << "1st diff" << endl;
  diff(coord,dcoord) ;

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

void Metric::Generic::observerTetrad(string const obskind,
				     double const coord[4], double fourvel[4],
				     double screen1[4], double screen2[4],
				     double screen3[4]) const{
  if (obskind == "ZAMO") {
    zamoVelocity(coord, fourvel);
  } else if (obskind=="KeplerianObserver") {
    circularVelocity(coord, fourvel);
  } else {
    // Let's renormalize fourvel
    double tdotin=fourvel[0];
    double threevel[3]={fourvel[1]/tdotin,
			fourvel[2]/tdotin,
			fourvel[3]/tdotin};
    fourvel[0]=SysPrimeToTdot(coord, threevel);
  }
  if (obskind != "FullySpecified")
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

void Metric::Generic::observerTetrad(double const coord[4], double const fourvel[4],
				     double screen1[4], double screen2[4],
				     double screen3[4]) const{
  GYOTO_ERROR("Metric::<ThisKind>::observerTetrad not implemented.");
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
