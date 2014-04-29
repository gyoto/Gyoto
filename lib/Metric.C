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
#include <iostream>
#include <cstdlib>
#include "GyotoMetric.h"
#include "GyotoFactoryMessenger.h"
#include "GyotoConverters.h"
#include <cmath>
#include <string>
#include <sstream>

using namespace std ; 
using namespace Gyoto;

Register::Entry* Metric::Register_ = NULL;

// Default constructor
Metric::Generic::Generic() :
  mass_(1.), coordkind_(GYOTO_COORDKIND_UNSPECIFIED),
  delta_min_(GYOTO_DEFAULT_DELTA_MIN),
  delta_max_(GYOTO_DEFAULT_DELTA_MAX),
  delta_max_over_r_(GYOTO_DEFAULT_DELTA_MAX_OVER_R)
{
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << endl;
#endif
  kind("Unspecified");
}

Metric::Generic::Generic(const double mas, const int coordkind) :
  mass_(mas), coordkind_(coordkind),
  delta_min_(GYOTO_DEFAULT_DELTA_MIN),
  delta_max_(GYOTO_DEFAULT_DELTA_MAX),
  delta_max_over_r_(GYOTO_DEFAULT_DELTA_MAX_OVER_R)
{
# if GYOTO_DEBUG_ENABLED
  GYOTO_IF_DEBUG;
  GYOTO_DEBUG_EXPR(mass_);
  GYOTO_DEBUG_EXPR(coordkind_);
  GYOTO_ENDIF_DEBUG;
# endif
  kind("Unspecified");
}

Metric::Generic::Generic(const int coordkind) :
  mass_(1.), coordkind_(coordkind),
  delta_min_(GYOTO_DEFAULT_DELTA_MIN),
  delta_max_(GYOTO_DEFAULT_DELTA_MAX),
  delta_max_over_r_(GYOTO_DEFAULT_DELTA_MAX_OVER_R)
{
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG_EXPR(coordkind_);
# endif
  kind("Unspecified");
}

// No copy constructor needed, default is fine
Metric::Generic * Metric::Generic::clone() const {
  string msg = "Metric::Generic::clone() called: cloning not supported for metric kind ";
  msg += kind();
  throwError (msg);
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
void Metric::Generic::coordKind(int coordkind)       { coordkind_=coordkind; }

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
  double h1max_at_r=pos[1];
  if (coordkind_==GYOTO_COORDKIND_CARTESIAN) {
    if (pos[2]>h1max_at_r) h1max_at_r=pos[2];
    if (pos[3]>h1max_at_r) h1max_at_r=pos[3];
  }
  h1max_at_r *= delta_max_over_r_;
  if (h1max > h1max_at_r) h1max = h1max_at_r;
  if (h1max>delta_max_) h1max=delta_max_;
  if (h1max<delta_min_) h1max=delta_min_;
  return h1max;
}

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

double Metric::Generic::Norm3D(double* pos) const {
  throwError("Check Norm3D");
  double res=0.;
  double g[4][4];
  gmunu(g, pos);
  for (int i=0;i<3;i++) {
    for (int j=0;j<3;j++) {
      res+=g[i+1][j+1]*pos[i]*pos[j];
    }
  }
  return sqrt(res);
}


/***************Geodesic most general integration**************/


double Metric::Generic::gmunu(const double * x, int mu, int nu) const {
  double g[4][4];
  gmunu(g, x);
  return g[mu][nu];
}

void Metric::Generic::gmunu(double g[4][4], const double * x) const {
  size_t mu, nu;
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

void Metric::Generic::christoffel(double dst[4][4][4], const double * x) const {
  size_t alpha, mu, nu;
  for (alpha=0; alpha<4; ++alpha) {
    for (mu=0; mu<4; ++mu) {
      dst[alpha][mu][mu]=christoffel(x, alpha, mu, mu);
      for (nu=mu+1; nu<4; ++nu)
	dst[alpha][mu][nu]=dst[alpha][nu][mu]=christoffel(x, alpha, mu, nu);
    }
  }
}

/*
Let : Y=[x0,x1,x2,x3,x0_dot,x1_dot,x2_dot,x3_dot] (dot=d/dtau, tau=proper time)
diff is such as : Y_dot=diff(Y)
The general equation of geodesics is used.
 */
int Metric::Generic::diff(const double coord[8], double res[8]) const{
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << endl;
# endif
  res[0]=coord[4];
  res[1]=coord[5];
  res[2]=coord[6];
  res[3]=coord[7];
  double dst[4][4][4];
  christoffel(dst, coord);
  for(int alpha=0; alpha<4; ++alpha) {
    res[alpha+4]=0.;
    for (int i=0;i<4;i++)
      for (int j=0;j<4;j++)
	res[alpha+4] -= dst[alpha][i][j]*coord[4+i]*coord[4+j];
  }
  return 0;
}

/*Runge Kutta to order 4

 */

//int Metric::Generic::myrk4(const double y[6], const double* cst , double h, double* res) const{
int Metric::Generic::myrk4(Worldline *, const double coord[8], double h, double res[8]) const{
  //cout << "In Metric::Generic::myrk4" << endl;
  double k1[8] ; 
  double k2[8] ; 
  double k3[8] ; 
  double k4[8] ; 
  double coord_plus_halfk1[8] ; 
  double sixth_k1[8] ; 
  double coord_plus_halfk2[8] ; 
  double third_k2[8] ; 
  double coord_plus_k3[8] ; 
  double third_k3[8] ; 
  double sixth_k4[8] ; 
	  
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

void Metric::Generic::circularVelocity(double const * , double*, double) const {
  stringstream ss;
  ss << kind_ << "::circularVelocity() is not implemented";
  throwError(ss.str());
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
    Gyoto::throwError
      ("Metric::Generic::cartesianVelocity: unknown coordinate kind");
  }
}

int Metric::Generic::myrk4_adaptive(Worldline* line, const double * coord, double lastnorm , double normref, double* coordnew , double h0, double& h1, double h1max) const{ 
  
  double delta0[8];
  double delta0min=1e-15;
  double dcoord[8];
  double eps=0.0001;
  double S=0.9;
  double errmin=1e-6;
  double factnorm=2.;

  if (h1max>delta_max_) h1max=delta_max_;
  if (h1max<delta_min_) h1max=delta_min_;
 
  //cout << "1st diff" << endl;
  diff(coord,dcoord) ;

  for (int i = 0;i<8;i++) delta0[i]=delta0min+eps*(fabs(h0*dcoord[i]));

  double hbis=0.5*h0;
  double coordhalf[8];
  double coord2[8];
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
    myrk4(line,coord,h0,coordnew);
    myrk4(line,coord,hbis,coordhalf);
    myrk4(line,coordhalf,hbis,coord2);
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
      
      newnorm=ScalarProd(coordnew, coordnew+4, coordnew+4);

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



/***************For SmartPointers**************/

int Metric::Generic::getRefCount() { return SmartPointee::getRefCount(); }

#ifdef GYOTO_USE_XERCES
void Metric::Generic::fillElement(Gyoto::FactoryMessenger *fmp) {
  fmp -> setSelfAttribute("kind", kind_);
  fmp -> setParameter("Mass", mass());
  if (delta_min_!=GYOTO_DEFAULT_DELTA_MIN)
    fmp -> setParameter("DeltaMin", delta_min_);
  if (delta_max_!=GYOTO_DEFAULT_DELTA_MAX)
    fmp -> setParameter("DeltaMax", delta_max_);
  if (delta_max_over_r_ != GYOTO_DEFAULT_DELTA_MAX_OVER_R)
    fmp -> setParameter("DeltaMaxOverR", delta_max_over_r_);
}

void Metric::Generic::setParameter(string name, string content, string unit) {
  if      (name=="Mass")     mass(atof(content.c_str()), unit);
  else if (name=="DeltaMin") deltaMin(atof(content.c_str()));
  else if (name=="DeltaMax") deltaMax(atof(content.c_str()));
  else if (name=="DeltaMaxOverR") deltaMaxOverR (atof(content.c_str()));
}

void Metric::Generic::setParameters(Gyoto::FactoryMessenger *fmp)  {
  string name="", content="", unit="";
  if (fmp)
    while (fmp->getNextParameter(&name, &content, &unit))
      setParameter(name, content, unit);
}

void Metric::Generic::processGenericParameters(Gyoto::FactoryMessenger *fmp)  {
  if (!fmp) return;
  string name="", content="";
  fmp -> reset();
  while (fmp->getNextParameter(&name, &content)) {
    if(name=="Mass")
      mass(atof(content.c_str()), fmp -> getAttribute("unit"));
  }
}

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

Metric::Subcontractor_t*
Metric::getSubcontractor(std::string name, int errmode) {
  if (!Gyoto::Metric::Register_) throwError("No Metric kind registered!");
  return (Metric::Subcontractor_t*)Gyoto::Metric::Register_
    -> getSubcontractor(name, errmode);
}

#endif
