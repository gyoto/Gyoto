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
#include <cmath>
#include <string>
#include <sstream>

#ifdef HAVE_UDUNITS
# include <udunits2.h>
#endif

using namespace std ; 
using namespace Gyoto;

Register::Entry* Metric::Register_ = NULL;

// Default constructor
Metric::Generic::Generic() :
  mass_(1.), coordkind_(GYOTO_COORDKIND_UNSPECIFIED)
{
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << endl;
#endif
  setKind("Unspecified");
}

Metric::Generic::Generic(const double mass, const int coordkind) :
  mass_(mass), coordkind_(coordkind)
{
# if GYOTO_DEBUG_ENABLED
  GYOTO_IF_DEBUG;
  GYOTO_DEBUG_EXPR(mass_);
  GYOTO_DEBUG_EXPR(coordkind_);
  GYOTO_ENDIF_DEBUG;
# endif
  setKind("Unspecified");
}

Metric::Generic::Generic(const int coordkind) :
  mass_(1.), coordkind_(coordkind)
{
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG_EXPR(coordkind_);
# endif
  setKind("Unspecified");
}

// No copy constructor needed, default is fine
Metric::Generic * Metric::Generic::clone() const {
  string msg = "Metric::Generic::clone() called: cloning not supported for metric kind ";
  msg += getKind();
  throwError (msg);
  return const_cast<Metric::Generic*>(this); // to avoid warning
}

Metric::Generic::~Generic(){
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << endl;
# endif
}

// Output

const string Metric::Generic::getKind() const {return kind_;}
void Metric::Generic::setKind(const string src) { kind_ = src;}

/***************Definition of the physical scene**************/

void Metric::Generic::setMass(const double mass)        { mass_=mass;}
void Metric::Generic::setMass(const double mass, string unit) {
# ifdef GYOTO_DEBUG_ENABLED
  GYOTO_IF_DEBUG
  GYOTO_DEBUG_EXPR(mass);
  GYOTO_DEBUG_EXPR(unit);
  GYOTO_ENDIF_DEBUG
# endif
# ifdef HAVE_UDUNITS
  if (unit=="") {
    mass_=mass;
    return;
  }
  ut_system * utsys = ut_read_xml(NULL);
  ut_unit * kg = ut_get_unit_by_symbol(utsys, "kg");
  if (!kg) throwError("BUG");
  ut_unit * from = NULL;
  if (unit=="sunmass") {
    from = ut_scale(GYOTO_SUN_MASS, kg);
  }
  if (!from) from = ut_get_unit_by_symbol(utsys, unit.c_str());
  if (!from) from = ut_get_unit_by_name(utsys, unit.c_str());
  if (!ut_are_convertible(from, kg)) {
    stringstream ss;
    ss << "Unsupported mass unit: \"" << unit;
    throwError(ss.str());
  }
  cv_converter * converter = ut_get_converter(from, kg);
  mass_ = cv_convert_double(converter, mass);
  ut_free(kg);
  ut_free(from);
  cv_free(converter);
  ut_free_system(utsys);
# else
  mass_ = mass;
  if (unit=="" || unit=="kg") ; // do nothing !
  else if (unit=="sunmass") mass_ *= GYOTO_SUN_MASS;
  else if (unit=="g") mass_ *= 1e-3;
  else {
    stringstream ss;
    ss << "Unsupported mass unit: \"" << unit
       << "\". Supported units: [kg] g sunmass";
    throwError(ss.str());
  }
# endif
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << "(mass="<<mass<<", unit=\"" << unit << "\") : mass_="
	      << mass_ <<" kg"<< endl;
# endif
}

int Metric::Generic::getCoordKind()               const { return coordkind_; }
void Metric::Generic::setCoordKind(int coordkind)       { coordkind_=coordkind; }

double Metric::Generic::getMass()                 const { return mass_; }


double Metric::Generic::SysPrimeToTdot(const double pos[4], const double v[3]) const {
  double sum=0.,xpr[4];
  int i,j;
# if GYOTO_DEBUG_ENABLED
  GYOTO_IF_DEBUG
    GYOTO_DEBUG_ARRAY(pos,4);
    GYOTO_DEBUG_ARRAY(v,3);
  GYOTO_ENDIF_DEBUG
# endif


  xpr[0]=1.; // dt/dt=1;
  for (i=0;i<3;++i) xpr[i+1]=v[i];
  for (i=0;i<4;++i) {
    for (j=0;j<4;++j) {
      sum+=gmunu(pos,i, j)*xpr[i]*xpr[j];
    }
  }
  if (sum>=0) {
    stringstream ss;
    ss << "In Metric.C: Impossible condition (v>c): tdot^2="<< -1/sum
       << " at pos=["<<pos[0] << ", " << pos[1] <<", "<<pos[2] <<", "<< pos[3]
       << "], vel=[" << v[0] << ", " << v[1] <<", "<<v[2]<<"]";
    throwError(ss.str());

  }
  return sqrt(-1./sum);
}


void Metric::Generic::nullifyCoord(double coord[8]) const {
  double tdot2;
  nullifyCoord(coord, tdot2);
}

void Metric::Generic::nullifyCoord(double coord[8], double& tdot2) const {
  //int width=25;//15;
  //int prec=15;//8;
  //cout << "DNAS NUllIFY" << endl;

  // if (debug()) cout << "dans Metric::Generic::nullifyCoord coord au debut:  " << coord[0] << " " << coord[1] << " " << coord[2] << " " << coord[3] << " " << coord[4] << " " << coord[5] << " " << coord[6] << " " << coord[7] << " " << endl;

  int i, j;
  double a, b=0., c=0.;
  a=gmunu(coord, 0, 0);
  //cout << "gtt = " << gmunu(coord, 0, 0) << endl;
  //cout << "coord' dans nullify:" << endl;
  //for (int myi=4;myi<8;myi++) cout << coord[myi] << " ";
  //cout << endl;
  for (i=1;i<=3;++i) {
    b+=gmunu(coord, 0, i)*coord[4+i];
    //cout << "i g0i= " << i << " " << gmunu(coord, 0, i) << endl;
    for (j=1;j<=3;++j) {
      //cout << "add to c: " << gmunu(coord, i, j)*coord[4+i]*coord[4+j] << endl;
      //cout << "i j gij= " << i << " " << j << " " << gmunu(coord, i, j) << endl;
      c+=gmunu(coord, i, j)*coord[4+i]*coord[4+j];
    }
  }
  //cout << "a b c in Metric::Generic::nullify= " << a << " " << b << " " << c << endl;
  // a*tdot*tdot + 2*b*tdot + c ==0.
  double Delta=b*b-a*c;
  tdot2=(-b+sqrt(Delta))/a;
  //  cout << "Delta a b c= " << setprecision(prec) << setw(width) << Delta << " " << setprecision(prec) << setw(width) << a << " " << setprecision(prec) << setw(width) << b << " " << setprecision(prec) << setw(width) << c << endl;
  coord[4]=(-b-sqrt(Delta))/a;
  // cout << "Tdot= " << coord[4] << endl;
  /*
    A creuser : en un point x0,x1,x2,x3 je spécifie la 3vitesse. Il existe deux géodésiques vérifiant ces conditions : celle qui pointe vers le futur (x0dot > 0) et celle qui pointe vers le passé (x0dot < 0).
    La solution (-b-sqrt(Delta))/a semble être tjrs positive. C'est celle qu'on veut.

   */
}

double Metric::Generic::ScalarProd(const double pos[4],
			  const double u1[4], const double u2[4]) const {
  // cout << "a comparer a = " << gmunu(pos,0,0) << " " << gmunu(pos,2,2) << endl;
  double res=0.;
  for (int i=0;i<4;i++) {
    for (int j=0;j<4;j++) {
      
      res+=gmunu(pos, i, j)*u1[i]*u2[j];
      //cout << "dans ScalarP: " << i << " " << j << " " << gmunu(pos, i, j) << " " << u1[i] <<  " " << u2[j] << " " << gmunu(pos, i, j)*u1[i]*u2[j] << endl;

    }
  }

  /*for (int kk=0;kk<4;kk++) cout << "Metric: In Scalar: " << setprecision(10) << gmunu(pos, kk, kk) << " " << u1[kk] << " " << u2[kk] << endl;
    cout << "** pour comp!" << gmunu(pos, 0, 0)*u1[0]*u2[0]+gmunu(pos, 1, 1)*u1[1]*u2[1]+gmunu(pos, 2, 2)*u1[2]*u2[2]+gmunu(pos, 3, 3)*u1[3]*u2[3] << endl;*/
  
  return res;
}

double Metric::Generic::Norm3D(double* pos) const {
  throwError("Check Norm3D");
  double res=0.;
  for (int i=0;i<3;i++) {
    for (int j=0;j<3;j++) {
      
      res+=gmunu(pos,i+1, j+1)*pos[i]*pos[j];
      // if (debug()) cout << "i,j,res= " << i << " " << j << " " << res << endl;
    }
  }
  return sqrt(res);
}


/***************Geodesic most general integration**************/

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
  res[4]=res[5]=res[6]=res[7]=0.;
  for (int i=0;i<4;i++) {
    for (int j=0;j<4;j++) {
      res[4]-=christoffel(coord,0,i,j)*coord[4+i]*coord[4+j];
      res[5]-=christoffel(coord,1,i,j)*coord[4+i]*coord[4+j];
      res[6]-=christoffel(coord,2,i,j)*coord[4+i]*coord[4+j];
      res[7]-=christoffel(coord,3,i,j)*coord[4+i]*coord[4+j];
    }
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

int Metric::Generic::myrk4_adaptive(Worldline* line, const double * coord, double lastnorm , double normref, double* coordnew , double h0, double& h1) const{ 
  
  double delta0[8];
  double delta0min=1e-15;
  double dcoord[8];
  double eps=0.0001;
  double S=0.9;
  double errmin=1e-6;
  double h1min=0.001;
  double h1max=1e6;
  double factnorm=2.;
 
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
      if (err<fabs(delta1[i]/delta0[i])) err=fabs(delta1[i]/delta0[i]);
    }

    if (err>1) {
      h0=S*h0*pow(err,-0.25);
      hbis=0.5*h0;
    }else{
      h1=(err > errmin ? S*h0*pow(err,-0.2) : 4.*h0);//pour éviter les explosions
      if (fabs(h1)<h1min) h1=(h0>0.)?h1min:-h1min;
      if (fabs(h1)>h1max) h1=(h0>0.)?h1max:-h1max;
      
      //Testing tangent vector norm stays next to 0 :
      
      newnorm=ScalarProd(coordnew, coordnew+4, coordnew+4);

      if ( fabs(newnorm-normref) > factnorm*fabs(lastnorm-normref) ) {
	//cout << "norm big!" << endl;
	//cout << "newnorm= " << newnorm << endl;
      	//myrk4(coord,h0/10.,coordnew);
      	//h1/=10.;
      }
      //cout << "Metric.C: step used= " << h0 << endl;
      break;
    

    }
  }

  return 0;

}

double Metric::Generic::unitLength() const { 
  return mass_ * GYOTO_G_OVER_C_SQUARE; 
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
  fmp -> setParameter("Mass", getMass());
}

void Metric::Generic::processGenericParameters(Gyoto::FactoryMessenger *fmp)  {
  string name="", content="";
  fmp -> reset();
  while (fmp->getNextParameter(&name, &content)) {
    if(name=="Mass")
      setMass(atof(content.c_str()), fmp -> getAttribute("unit"));
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

Metric::Subcontractor_t* Metric::getSubcontractor(std::string name) {
  if (!Gyoto::Metric::Register_) throwError("No Metric kind registered!");
  return (Metric::Subcontractor_t*)Gyoto::Metric::Register_
    -> getSubcontractor(name);
}

#endif
