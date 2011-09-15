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
#include <GyotoKerrKS.h>
#include <GyotoWorldline.h>
#include <GyotoError.h>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <string>
#include <cstring>

using namespace std ; 
using namespace Gyoto ; 

#define drhor 0.5
//stop integration at r_horizon+drhor for outgoing geodesics

/*
NOTA BENE: to improve KerrKS
So far (March 2011) KerrKS is just a stub ; to improve it, it's necessary to imporve myrk4_adaptive, which is really obsolete (and longer!) as compared to KerrBL ; a check of cst of motion conservation should be implemented: so far, the norm behaves badly when approaching (not-so-)close to horizon, that's why drhor is chosen so big.
In particular, don't trust too much the result with spin>0
 */

KerrKS::KerrKS() : Metric(GYOTO_COORDKIND_CARTESIAN), spin_(0.) {setKind("KerrKS");}

KerrKS::KerrKS(double a, double m) : Metric(m, GYOTO_COORDKIND_CARTESIAN), spin_(a) {setKind("KerrKS");}

// default copy constructor should be fine
// KerrKS::KerrKS(const KerrKS& gg) : 
//          Metric(gg), spin_(gg.spin_) {setKind("KerrKS");}
// "virtual" copy constructor:
KerrKS * KerrKS::clone () const { return new KerrKS(*this); }


KerrKS::~KerrKS() {}

// Output
/*
std::ostream& Gyoto::operator<<( std::ostream& o, const KerrKS& met ) {
  return  met.print(o);
}

std::ostream& KerrKS::print( std::ostream& o) const {
  o << "spin=" << spin_ << ", " ;
  Metric::print(o);
  return o;
}
*/

// Mutators
void KerrKS::setSpin(const double spin) { spin_=spin; }

// Accessors
double KerrKS::getSpin() const { return spin_ ; }

double KerrKS::christoffel(const double[8],
		   const int, const int, const int) const{
  throwError( "KerrKS.C : should never come here to find christoffel!!" );
  return 0.;
}

double KerrKS::gmunu(const double * pos, int mu, int nu) const {
  if (mu<0 || nu<0 || mu>3 || nu>3) throwError ("KerrKS::gmunu: incorrect value for mu or nu");
  //double x=pos[0], y=pos[1], z=pos[2];
  double x=pos[1], y=pos[2], z=pos[3];
  double x2=x*x;
  double y2=y*y;
  double z2=z*z;
  double a2=spin_*spin_;
  double temp=x2+y2+z2-a2;
  double rr=sqrt(0.5*(temp+sqrt(temp*temp+4*a2*z2)));
  double r2=rr*rr;
  double r3=rr*r2;
  double r4=rr*r3;
  double fact=2.*r3/(r4+a2*z2);

  double res=0.;
  if (mu==nu) {
    if ((mu==0) && (nu==0)) res=fact-1.;
    if ((mu==1) && (nu==1)) res= 1.+fact*pow((rr*x+spin_*y)/(r2+a2),2);
    if ((mu==2) && (nu==2)) res= 1.+fact*pow((rr*y-spin_*x)/(r2+a2),2);
    if ((mu==3) && (nu==3)) res= 1.+fact*z2/r2;
  }
  if (nu<mu) {int vu=nu; nu=mu; mu=vu;}
  if (mu==0) {
    if (nu==1) res= fact/(r2+a2)*(rr*x+spin_*y);
    if (nu==2) res= fact/(r2+a2)*(rr*y-spin_*x);
    if (nu==3) res= fact*z/rr;
  }
  if (mu==1) {
    if (nu==2) res= fact/pow(r2+a2,2)*(x*y*(r2-a2)+spin_*rr*(y2-x2));
    if (nu==3) res= fact/(r2+a2)*(rr*x+spin_*y)*z/rr;
  }
  if ((mu==2) && (nu==3)) res= fact/(r2+a2)*(rr*y-spin_*x)*z/rr;

  return res;
  
} 

int KerrKS::diff(const double* coord, const double* cst, double* res) const{

  //Important note: coord and res are double[7], not double[8] because
  //there's no Tdotdot equation Here coord MUST BE:
  //(T,x,y,z,xdot,ydot,zdot)

  //So far the following equations are only true for a 0-mass
  //  particule if (cst[0]!=0. && debug()) cout << " ****WARNING**** :
  //  KS equations used for a non 0-mass particule!" << endl;
  if (cst[0]!=0. && debug()) throwError("Kerr-Schild equations used for a non 0-mass particle!");
  //  int width=15;
  //  int prec=20;
  
  double xx=coord[1]; 
  double yy=coord[2];
  double zz=coord[3];
  double xdot=coord[4]; //4, not 5, see above comment
  double ydot=coord[5];
  double zdot=coord[6];

  //Using Eqs. in Hameury, Marck, Pelat 94
  double spin2=spin_*spin_;
  //to determine rr and rdot, use x^2+y^2+z^2=r^2+a^2*(1-z^2/r^2)
  double temp=xx*xx+yy*yy+zz*zz-spin2;
  double rr=sqrt(0.5*(temp+sqrt(temp*temp+4*spin2*zz*zz)));
  double r2=rr*rr;

  //  double theta=acos(zz/rr);
  double rdot=(xx*xdot+yy*ydot+zz*zdot+spin2*zz*zdot/r2)/(rr+spin2*zz*zz/(rr*r2));

  double sico=(xx*rr+spin_*yy)/(r2+spin2);// = sin(theta)*cos(psi)
  double sisi=(rr*yy-spin_*xx)/(r2+spin2);// = sin(theta)*sin(psi)

  double EE=cst[1];
  double LL=cst[2];
  double KK=cst[3]+pow(LL-spin_*EE,2);

  double Sigma=r2+spin2*zz*zz/r2; // NB: cos(theta)=z/r
  double BigE=EE*(r2+spin2)-spin_*LL;
  double BigQ=Sigma*rdot+BigE;
  double Delta=r2-2.*rr+spin2;
  double Sigma3=Sigma*Sigma*Sigma;
  double QoverD=BigQ/Delta;

  if (BigE == Sigma*rdot) {
    if (debug()) cout << "WARNING: Outgoing geodesic can't cross the horizon! Stopping..." << endl;
    return 1;
  }
  
  double Tdot = 2*KK*rr/(Sigma*(BigE-Sigma*rdot))+EE;

  //Horizon test:
  double rsink=1.+sqrt(1.-spin2)+drhor;

  if (rr<rsink && rdot >0 && Tdot>0) {// even if KS is okay at the horizon, a geodesic can't escape the horizon, so the eq of geodesic fail there --> must check we're far enough
    if (debug()) 
      cerr << "Too close to horizon in KerrKS::diff at r= " << rr << endl;
    return 1;
  }

  res[0]=Tdot;
  res[1]=xdot;
  res[2]=ydot;
  res[3]=zdot;
  
  res[4]=1./Sigma3*(-4.*spin_*rr*QoverD*Sigma*ydot
		    +(Sigma-4.*r2)*sico*(KK-spin2*QoverD*QoverD)
		    -spin_*rr*sisi*QoverD*(4*(EE*Sigma-BigQ)+(4.*spin2-Sigma)*QoverD));//xdotdot

  res[5]=1./Sigma3*(4.*spin_*rr*QoverD*Sigma*ydot
  		    +(Sigma-4.*r2)*sisi*(KK-spin2*QoverD*QoverD)
  		    +spin_*rr*sico*QoverD*(4*(EE*Sigma-BigQ)+(4.*spin2-Sigma)*QoverD));//ydotdot
  
  res[6]=-1./Sigma3*KK*zz/rr*(3.*r2-spin2*zz*zz/r2);//zdotdot ; NB: cos(theta)=z/r
 
  return 0;
}

int KerrKS::diff(const double*, double*) const{
  throwError("In KerrKS::diff should never get here!");
  return 0;
}
int KerrKS::myrk4_adaptive(Worldline* line, const double * coord, double , double , double* coord1, double h0, double& h1) const{

  double const * const cst = line -> getCst();

  /*  cout << "Cst motion in KerrKS= ";
  for (int ii=0;ii<4;ii++) cout << cst[ii] << " " ;
  cout << endl;*/

  //  int width=15;
  //  int prec=8;

  double delta0[7]; //caution, diff's result is a 7-sized vector, there's no Tdotdot equation
  double delta0min=1e-15;
  double dcoord[7];
  double eps=1e-4;//0.0001;
  double S=0.9;
  double errmin=1e-6;
  double h1min=0.001;
  double h1max=1.; //NB : .1 is necessary to obtain nice images! but long...
  //1e6; //low h1max necessary to integrate properly the thin disk
  //double factnorm=2.;
  //  double newnorm1, newnorm2;
 
  double coordtemp[7]={coord[0],coord[1],coord[2],coord[3],coord[5],coord[6],coord[7]};
  //Caution!! diff must be fed with 7-sized vectors! see comment in diff
  if (diff(coordtemp,cst,dcoord)) return 1;

  for (int i = 0;i<7;i++){
    delta0[i]=eps*(fabs(coord[i])+fabs(h0*dcoord[i]));
    delta0[i]=delta0min+eps*(fabs(h0*dcoord[i]));
  }

  double hbis=0.5*h0;
  
  double coordhalf[8];
  double coord2[8];
  double delta1[7];
  
  double err;
  int count=0;

  double newnorm;
  
  while (1){
    count++;
    //cout << "count= " << count << endl;
    err=0.;
    if (myrk4(coord,cst,h0,coord1)) return 1; 
                   //NB: myrk4 is fed with 8-sized vector, no pb here
    if (myrk4(coord,cst,hbis,coordhalf)) return 1;
    if (myrk4(coordhalf,cst,hbis,coord2)) return 1;

    for (int i = 0;i<4;i++){
      delta1[i]=coord2[i]-coord1[i];
      //cout << "del= " << delta1[i] << endl;
      if (err<fabs(delta1[i]/delta0[i])) err=fabs(delta1[i]/delta0[i]);
    }
    for (int i = 5;i<8;i++){//NB: delta0 is a 7-sized vector, see comment in diff
      delta1[i]=coord2[i]-coord1[i];
      //cout << "del= " << delta1[i] << endl;
      if (err<fabs(delta1[i]/delta0[i-1])) err=fabs(delta1[i]/delta0[i-1]);
    }

    /*cout << "coord in RK= " ;
    for (int ii=0;ii<8;ii++) cout << setprecision(GYOTO_PREC) << setw(GYOTO_WIDTH) << coord[ii] << " ";
    cout << endl;
    cout << "coord1 in RK= " ;
    for (int ii=0;ii<8;ii++) cout << setprecision(GYOTO_PREC) << setw(GYOTO_WIDTH) << coord1[ii] << " ";
    cout << endl;
    cout << "coord2 in RK= " ;
    for (int ii=0;ii<8;ii++) cout << setprecision(GYOTO_PREC) << setw(GYOTO_WIDTH) << coord2[ii] << " ";
    cout << endl;
    
    cout << "err= " << err << endl;*/
    
    if (err>1) {
      h0=S*h0*pow(err,-0.25);
      hbis=0.5*h0;
            
    }else{
      h1=(err > errmin ? S*h0*pow(err,-0.2) : 4.*h0);//pour éviter les explosions
      //      cout << "h0,h1= " << h0 << " " << h1 << endl;
      if (fabs(h1)<h1min) h1=h0>0?h1min:-h1min;
      if (fabs(h1)>h1max) h1=h0>0?h1max:-h1max;
 
      newnorm=ScalarProd(coord1, coord1+4, coord1+4);

      //This "norm check" is trash, it's of no use, except to slow down drastically the computation (even with it, the norm can be very bad)
      /*if (fabs(newnorm-normref) > factnorm*fabs(lastnorm-normref)) {      
	
      	if (myrk4(coord,cst,h0/10.,coord1)) return 1;
      	h1/=10.;
	
	}*/
      
      break;
    }
    
    
  }
  /*cout << "coordnew in RK= " ;
  for (int ii=0;ii<8;ii++) cout << setprecision(GYOTO_PREC) << setw(GYOTO_WIDTH) << coord1[ii] << " ";
  cout << endl;*/
  //cout << "KerrKS: norm, h1= " << setprecision(GYOTO_PREC) << setw(GYOTO_WIDTH) << newnorm << " " << ScalarProd(coord1, coord1+4, coord1+4) << endl;
 
  /*double xx=coord1[1], yy=coord1[2], zz=coord1[3];
  double aa=spin_, a2=aa*aa;
  double temp=xx*xx+yy*yy+zz*zz-a2;
  double rr=sqrt(0.5*(temp+sqrt(temp*temp+4*a2*zz*zz)));
  cout << "In KerrKS::ada r= " << rr << endl;*/
  return 0;
    
}

int KerrKS::myrk4(Worldline *line, const double coord[8], double hh, double res[8]) const {
  double const * const cst = line -> getCst();
  return myrk4(coord,cst,hh,res);
}

int KerrKS::myrk4(const double* coord, const double* cst , double h, double* res) const{
  
  //cout << "In KerrKS rk4 h= " << h << endl;
  double coordtemp[7]={coord[0],coord[1],coord[2],coord[3],coord[5],coord[6],coord[7]};
  double restemp[7];
  double k1[7] ; 
  double k2[7] ; 
  double k3[7] ; 
  double k4[7] ; 
  double coord_plus_halfk1[7] ; 
  double sixth_k1[7] ; 
  double coord_plus_halfk2[7] ; 
  double third_k2[7] ; 
  double coord_plus_k3[7] ; 
  double third_k3[7] ; 
  double sixth_k4[7] ; 
	  
  if (diff(coordtemp,cst, k1)) return 1; 
  //cout << "TEST dans rk= " << coordtemp[4] << endl;
  
  for (int i=0;i<7;i++) {
    k1[i]=h*k1[i];
    coord_plus_halfk1[i]=coordtemp[i]+0.5*k1[i];
    sixth_k1[i]=1./6.*k1[i];
  }
  
  if (diff(coord_plus_halfk1,cst,k2)) return 1; 
  for (int i=0;i<7;i++) {
    k2[i]=h*k2[i];
    coord_plus_halfk2[i]=coordtemp[i]+0.5*k2[i];
    third_k2[i]=1./3.*k2[i];
  }
	
  if (diff(coord_plus_halfk2,cst,k3)) return 1;
  for (int i=0;i<7;i++) {
    k3[i]=h*k3[i];
    coord_plus_k3[i]=coordtemp[i]+k3[i];
    third_k3[i]=1./3.*k3[i];
  }

  if (diff(coord_plus_k3,cst,k4)) return 1; 
  for (int i=0;i<7;i++) {
    k4[i]=h*k4[i];
    sixth_k4[i]=1./6.*k4[i];
  }


  for (int i=0;i<7;i++) {
    restemp[i]=coordtemp[i]+sixth_k1[i]+third_k2[i]+third_k3[i]+sixth_k4[i]; 
  }

  for (int i=0;i<4;i++) {
    res[i]=restemp[i];
  }
  
  for (int i=5;i<8;i++) {
    res[i]=restemp[i-1];
  }

  //computing Tdot 
  double ktemp[7];
  if (diff(restemp,cst,ktemp)) return 1;
  res[4]=ktemp[0];

  return 0;

}


void KerrKS::MakeCst(const double* coord, double cst[4]) const{
  /*
    To keep in mind : check the cst of integration by comparing to BL case, it must be the same for the same initial condition.
    Pay attention to the fact that it is NORMAL if the integration fails near the horizon when rdot<0 and Tdot<0 (see Eq 7 in Marck 96, if E=rdot - sign of E is related to sign of Tdot ; and Marck 96 reminds that at the horizon, E=abs(rdot)): a geodesic can't escape the horizon! However, if rdot<0 and Tdot >0 (which was the case in the earliest version of the code), the horizon can be passed. But it is an ingoing geodesic that is being integrated (here E=rdot never happens, E=-rdot at the horizon).
    Don't forget that near the horizon, the behavior of a photon can't be reversed.

   */
  //for the transformation between KS and BL see Hameury, Marck, Pelat 94
  // int width=15, prec=10;

  double xx=coord[1], yy=coord[2], zz=coord[3], Tdot=coord[4], xdot=coord[5], 
    ydot=coord[6], zdot=coord[7];

  double norm=ScalarProd(coord, coord+4, coord+4);
  // cout << "norm in MakeCst= " << norm << endl;

  double aa=spin_, a2=aa*aa;
  double temp=xx*xx+yy*yy+zz*zz-a2;
  double rr=sqrt(0.5*(temp+sqrt(temp*temp+4*a2*zz*zz))), r2=rr*rr;
  //double z2=zz*zz;

  
  double costh=zz/rr, costheta2=costh*costh, sintheta2=1.-costheta2, theta=acos(costh), sinth=sin(theta); //NB: acos gives a value between 0 and pi, which is exactly what we want for theta 

  if (sinth==0) {
    throwError("KerrKS::MakeCst : Initial condition on z-axis : BL coordinates singular at this point with theta=0");
  }
  
  double rdot=(xx*xdot+yy*ydot+zz*zdot+a2*zz*zdot/r2)/(rr+a2*zz*zz/(rr*r2)), thetadot=(rdot*costh-zdot)/(rr*sinth);
 
  //See Hameury, Mark & Pelat 94
  double Sigma=r2+a2*costheta2, Delta=r2-2*rr+a2, fact=2.*aa*rr*sintheta2/Sigma;
  double BLtdot=Tdot-2.*rr/Delta*rdot;

  double cosp=(rr*xx+aa*yy)/(sinth*(a2+r2)), sinp=(rr*yy-aa*xx)/(sinth*(a2+r2));//cos(psi) and sin(psi)
  double psidot;
  
  if (aa != 0.) {
    psidot=(xdot*cosp+ydot*sinp-rdot*sinth-rr*thetadot*costh)/(-aa*sinth);
  }else{
    if (cosp!=0.) {
      psidot=(ydot - sinp*(rdot*sinth+rr*thetadot*costh))/(rr*sinth*cosp);
    }else{
      psidot=(xdot - cosp*(rdot*sinth+rr*thetadot*costh))/(-rr*sinth*sinp);
    }
  }

  double phidot=psidot-aa/Delta*rdot;
  //cout << "Make Cst psidot= " << setprecision(GYOTO_PREC) << setw(GYOTO_WIDTH) << psidot << endl;
  
  double mu;//Particule mass: 0 or 1
  if (fabs(norm)<fabs(norm+1.)){
    mu=0.;
  }else{
    mu=1.;
  }

  double E=(1-2*rr/Sigma)*BLtdot+fact*phidot, //OK for particule mass = 0 or 1 
    L=sintheta2*(r2+a2+aa*fact)*phidot-fact*BLtdot, //OK for particule mass = 0 or 1
    Q=Sigma*thetadot*Sigma*thetadot+costheta2*(a2*(mu*mu-E*E)+L*L/sintheta2); //different for a 0-mass and a 1-mass particule

  cst[0]=mu;cst[1]=E;cst[2]=L;cst[3]=Q;
}

void KerrKS::nullifyCoord(double coord[8]) const {
  double tdot2;
  nullifyCoord(coord, tdot2);
}

void KerrKS::nullifyCoord(double coord[4], double& tdot2) const {
  Metric::nullifyCoord(coord, tdot2);
}

/*
  To compute the cst of motion for a given particle
 */
void KerrKS::setParticleProperties(Worldline * line, const double* coord) const {
  double cst[4];
  MakeCst(coord,cst);
  line -> setCst(cst,4);
}

int KerrKS::isStopCondition(double const * const coord) const {
  if (fabs(coord[5]) > 5000.) return 1;
  return 0;
}

#ifdef GYOTO_USE_XERCES
void KerrKS::fillElement(Gyoto::factoryMessenger *fmp) {
  fmp -> setParameter("Spin", spin_);
  Metric::fillElement(fmp);
}

SmartPointer<Metric> Gyoto::KerrKS::Subcontractor(factoryMessenger* fmp) {

  double spin=0., mass=1.; //default values
  string name="", content="";

  while (fmp->getNextParameter(&name, &content)) {
    if(name=="Spin") spin=atof(content.c_str());
  }
  SmartPointer<KerrKS> gg = new KerrKS(spin, mass);
  gg -> processGenericParameters(fmp);
  return gg;
}

void Gyoto::KerrKS::Init() {
  Gyoto::Metric::Register("KerrKS", &Subcontractor);
}
#endif
