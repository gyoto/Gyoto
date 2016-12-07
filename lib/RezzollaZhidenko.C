/*
    Copyright 2014 Frederic Vincent

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

#include "GyotoRezzollaZhidenko.h"

#include <iostream>
#include <cstdlib>
#include <cmath>

#define STRINGIFY(foo) #foo

using namespace Gyoto;
using namespace Gyoto::Metric;
using namespace std;

#define GYOTO_DRHOR 0.1
#define GYOTO_NBPARAM_MAX 4 // only this number of parameters is allowed, e.g. a0, a1, a2, a3 and not more for the time being; assumed the same for a and b.

/// Properties
#include "GyotoProperty.h"
GYOTO_PROPERTY_START(RezzollaZhidenko, "Spherically-symmetric parametrized metric of Rezzolla & Zhidenko 2014")
GYOTO_PROPERTY_DOUBLE(RezzollaZhidenko, Epsilon, epsilon)
GYOTO_PROPERTY_DOUBLE(RezzollaZhidenko, Rms, rms)
GYOTO_PROPERTY_DOUBLE(RezzollaZhidenko, Rmb, rmb)
GYOTO_PROPERTY_VECTOR_DOUBLE(RezzollaZhidenko, AParam, aparam,
			     "At most 4 elements")
GYOTO_PROPERTY_VECTOR_DOUBLE(RezzollaZhidenko, BParam, bparam,
			     "At most 4 elements")
GYOTO_PROPERTY_END(RezzollaZhidenko, Generic::properties)

// accessors
GYOTO_PROPERTY_ACCESSORS(RezzollaZhidenko, double, epsilon_, epsilon)
GYOTO_PROPERTY_ACCESSORS(RezzollaZhidenko, double, rms_, rms)
GYOTO_PROPERTY_ACCESSORS(RezzollaZhidenko, double, rmb_, rmb)

void RezzollaZhidenko::aparam(std::vector<double> const &v) {
  size_t n = v.size();
  if (n>GYOTO_NBPARAM_MAX)
    throwError("In RezzollaZhidenko: choose at most "
	       STRINGIFY(GYOTO_NBPARAM_MAX) " parameters");
  for (size_t i=0; i<n; ++i) {
    aparam_[i]=v[i];
    if (aparam_[i]<0.) throwError("In RezzollaZhidenko: param < 0!");
  }
  for (size_t i=n; i<GYOTO_NBPARAM_MAX; ++i) aparam_[i]=0.;
}
std::vector<double> RezzollaZhidenko::aparam() const {
  std::vector<double> v(GYOTO_NBPARAM_MAX, 0.);
  for (size_t i=0; i<GYOTO_NBPARAM_MAX; ++i) v[i]=aparam_[i];
  return v;
}
void RezzollaZhidenko::bparam(std::vector<double> const &v) {
  size_t n = v.size();
  if (n>GYOTO_NBPARAM_MAX)
    throwError("In RezzollaZhidenko: choose at most "
	       STRINGIFY(GYOTO_NBPARAM_MAX) " parameters");
  for (size_t i=0; i<n; ++i) {
    bparam_[i]=v[i];
    if (bparam_[i]<0.) throwError("In RezzollaZhidenko: param < 0!");
  }
  for (size_t i=n; i<GYOTO_NBPARAM_MAX; ++i) bparam_[i]=0.;
}
std::vector<double> RezzollaZhidenko::bparam() const {
  std::vector<double> v(GYOTO_NBPARAM_MAX, 0.);
  for (size_t i=0; i<GYOTO_NBPARAM_MAX; ++i) v[i]=bparam_[i];
  return v;
}

///

Gyoto::Metric::RezzollaZhidenko::RezzollaZhidenko()
  : Generic(GYOTO_COORDKIND_SPHERICAL, "RezzollaZhidenko"),
    epsilon_(0.), rms_(0.), rmb_(0.)
{
  GYOTO_DEBUG << endl;
  aparam_ = new double[GYOTO_NBPARAM_MAX];
  bparam_ = new double[GYOTO_NBPARAM_MAX];
  for (int ii=0;ii<GYOTO_NBPARAM_MAX;ii++){
    aparam_[ii]=0.;
    bparam_[ii]=0.;
  }
}

// default copy constructor should be fine 
RezzollaZhidenko * RezzollaZhidenko::clone() const { return new RezzollaZhidenko(*this); }

Gyoto::Metric::RezzollaZhidenko::~RezzollaZhidenko()
{
  GYOTO_DEBUG << endl;
  delete [] aparam_;
  delete [] bparam_;
}

double RezzollaZhidenko::getRms() const {
  return rms_;
}

double RezzollaZhidenko::getRmb() const {
  return rmb_;
}

double RezzollaZhidenko::getSpecificAngularMomentum(double rr) const {
  double NN=sqrt(N2(rr)), N3=NN*NN*NN;
  return sqrt(rr*rr*rr*Nprime(rr)/N3);
}

double RezzollaZhidenko::getPotential(double pos[4], double l_cst) const {
  double gtt = gmunu(pos,0,0);
  double gpp = gmunu(pos,3,3);
  if (gpp==0.) throwError("In RezzollaZhidenko: bad gpp");
  double rr = pos[1];
  double NN2=N2(rr), NN = sqrt(NN2);
  double  Omega = -l_cst * gtt/gpp ;
  return -2.*log(fabs(NN))+0.5*log(fabs(-NN2+gpp*Omega*Omega));
}

double RezzollaZhidenko::N2(const double rr) const{
  /*cout << "eps= " << epsilon_ << endl;
  cout << "aa= " ;
  for (int ii=0;ii<GYOTO_NBPARAM_MAX;ii++) cout << aparam_[ii] << " ";
  cout << endl;
  cout << "bb= " ;
  for (int ii=0;ii<GYOTO_NBPARAM_MAX;ii++) cout << bparam_[ii] << " ";
  cout << endl;*/

  /*
    RZ metric defines r0 as the event horizon location.
    It is related to epsilon through: epsilon = 2M/r0 - 1,
    where M is ADM mass of spacetime.
    Here it is assumed that units are such that M=1.
  */
  double r0 = 2./(1.+epsilon_);
  double xx = 1. - r0/rr, onemx = 1. - xx,
    onemx2 = onemx*onemx, onemx3 = onemx2*onemx;
  double Atilde = aparam_[1]/(1.+aparam_[2]*xx/(1.+aparam_[3]*xx));
  double N2 = xx*(1.-epsilon_*onemx+(aparam_[0]-epsilon_)*onemx2
		  *Atilde*onemx3);

  return N2;
}

double RezzollaZhidenko::B2(const double rr) const{
  double r0 = 2./(1.+epsilon_);
  double xx = 1. - r0/rr, onemx = 1. - xx,
    onemx2 = onemx*onemx, onemx3 = onemx2*onemx;
  double Btilde = bparam_[1]/(1.+bparam_[2]*xx/(1.+bparam_[3]*xx));
  double BB = 1+bparam_[0]*onemx+Btilde*onemx2,
    B2 = BB*BB;

  return B2;
}

double RezzollaZhidenko::Nprime(const double rr) const{
  double r0 = 2./(1.+epsilon_);
  double r2 = rr*rr;
  double xx = 1. - r0/rr, onemx = 1. - xx,
    onemx2 = onemx*onemx, onemx3 = onemx2*onemx;
  double Atilde = aparam_[1]/(1.+aparam_[2]*xx/(1.+aparam_[3]*xx));
  double AA = 1.-epsilon_*onemx+(aparam_[0]-epsilon_)*onemx2
    *Atilde*onemx3;
  double Atilde_der = -aparam_[1]*aparam_[2]/((1.+(aparam_[2]+aparam_[3])*xx)*(1.+(aparam_[2]+aparam_[3])*xx));
  double A_der = epsilon_ - 2.*(aparam_[0]-epsilon_)*onemx-3.*Atilde*onemx2
    +Atilde_der*onemx3;
  double NN = sqrt(N2(rr));
  
  return 1./(r2*NN)*(AA+xx*A_der);
  
}

double RezzollaZhidenko::Bprime(const double rr) const{
  double r0 = 2./(1.+epsilon_);
  double r2 = rr*rr;
  double xx = 1. - r0/rr, onemx = 1. - xx,
    onemx2 = onemx*onemx, onemx3 = onemx2*onemx;
  double Btilde = bparam_[1]/(1.+bparam_[2]*xx/(1.+bparam_[3]*xx));
  double Btilde_der = -bparam_[1]*bparam_[2]/((1.+(bparam_[2]+bparam_[3])*xx)*(1.+(bparam_[2]+bparam_[3])*xx));
  double B_der = -bparam_[0] - 2.*Btilde*onemx + Btilde_der*onemx2;
  
  return 1./r2*B_der;
  
}

double RezzollaZhidenko::gmunu(const double * pos, int mu, int nu) const {
  double rr = pos[1];
  if (rr<=0.) throwError("In RezzollaZhidenko::gmunu: r<0!");

  double sth2, cth2;
  sincos(pos[2], &sth2, &cth2);
  sth2*=sth2; cth2*=cth2;
  double r2=rr*rr;

  if ((mu==0) && (nu==0)) return -N2(rr);
  if ((mu==1) && (nu==1)) return B2(rr)/N2(rr);
  if ((mu==2) && (nu==2)) return r2;
  if ((mu==3) && (nu==3)) return r2*sth2;
  return 0.;
} 

int RezzollaZhidenko::christoffel(double dst[4][4][4], double const pos[4]) const
{
  int a, mu, nu;
  for (a=0; a<4; ++a)
    for(mu=0; mu<4; ++mu)
      for(nu=0; nu<4; ++nu)
	dst[a][mu][nu]=0.;

  double rr = pos[1];
  double sth, cth;
  sincos(pos[2], &sth, &cth);
  if (rr==0. || sth==0.) throwError("In RezzollaZhidenko::christoffel: "
				    "bad coord");
  double NN2=N2(rr), NN = sqrt(NN2), BB2=B2(rr), BB=sqrt(BB2);
  double NNprime=Nprime(rr), BBprime=Bprime(rr);

  dst[0][0][1]=dst[0][1][0]=NNprime/NN;
  dst[1][0][0]=NN*NN2/BB2*NNprime;
  dst[2][1][2]=dst[2][2][1]=1./rr;
  dst[3][1][3]=dst[3][3][1]=1./rr;
  dst[1][1][1]=BBprime/BB-NNprime/NN;
  dst[2][3][3]=-cth*sth;
  dst[3][2][3]=dst[3][3][2]=cth/sth;
  dst[1][2][2]=-rr*NN2/BB2;
  dst[1][3][3]=-rr*sth*sth*NN2/BB2;

  int compareSch=0;
  if (compareSch){
    // Comparison to Schwarzschild Christoffels for checking
    double dstS[4][4][4];
    for (a=0; a<4; ++a)
      for(mu=0; mu<4; ++mu)
	for(nu=0; nu<4; ++nu)
	  dstS[a][mu][nu]=0.;
    dstS[0][0][1]=dstS[0][1][0]=1./(rr*rr-2.*rr);
    dstS[1][0][0]=(rr-2.)/(rr*rr*rr);
    dstS[2][1][2]=dstS[2][2][1]=1./rr;
    dstS[3][1][3]=dstS[3][3][1]=1./rr;
    dstS[1][1][1]=1./(2*rr-rr*rr);
    dstS[2][3][3]=-cth*sth;
    dstS[3][2][3]=dstS[3][3][2]=cth/sth;
    dstS[1][2][2]=2.-rr;
    dstS[1][3][3]=(2.-rr)*sth*sth;
    for (a=0; a<4; ++a){
      for(mu=0; mu<4; ++mu){
	for(nu=0; nu<4; ++nu){
	  if(fabs(dst[a][mu][nu]-dstS[a][mu][nu])>1e-3) {
	    cout << "at r,th= " << rr <<  " " << pos[2] << endl;
	    cout << setprecision(20)<< a << " " << mu << " " << nu << " " << dst[a][mu][nu] << " " << dstS[a][mu][nu]<< endl;
	    throwError("test RZ");
	  }
	}
      }
    }
  }

  return 0;
}

int RezzollaZhidenko::isStopCondition(double const * const coord) const {
  double r0 = 2./(1.+epsilon_);
  double rsink = r0 + GYOTO_KERR_HORIZON_SECURITY;
  return coord[1] < rsink ;
}

int RezzollaZhidenko::diff(const double* coordGen, const double* cst, 
		      double* res) const{
  // OUTDATED TO BE CHECKED
  double rsink=2.+GYOTO_DRHOR;
  double r = coordGen[1] ; 

  if (r < 0.) {
    cerr << "r= " << r << endl;
    throwError( "RezzollaZhidenko.C : r negative!!!!! the horizon"
		" may have been crossed..." );
  }

  if (r < rsink) {
#   if GYOTO_DEBUG_ENABLED
    GYOTO_DEBUG << "Too close to horizon in RezzollaZhidenko::diff at r= " 
		<< r << endl;
#   endif
    return 1;
  }

  double r2 = r*r ; 
  double r3 = r2*r ;  

  double ff       = 1.+r2*1.*(1.-sqrt(1.+4./(1.*r3)));
  if (ff==0.) throwError("In RezzollaZhidenko::gmunu: ff is zero");
  double fprime   = 6./(r2*sqrt(1.+4./(1.*r3)))
    +2.*r*1.*(1.-sqrt(1.+4./(1.*r3))),
    fprimef2      = fprime/(ff*ff);

  double pr=coordGen[5];
  
  double EE=cst[1];
  double E2=EE*EE;
  double LL=cst[2];
  double L2=LL*LL;

  /*
    ---> Spherically symmetric EOM
  */
  res[0] = EE/ff; // tdot
  res[1] = ff*pr; // rdot
  res[2] = 0.;    // thdot (planar motion)
  res[3] = LL/r2; // phidot
  
  res[4] = 0.; // ptdot: pt = cst = -E
  res[5] = 0.5*fprime*pr*pr-L2/r3+0.5*E2*fprimef2;// prdot
  res[6]=  0.; // pthetadot
  res[7] = 0.; // pphidot: pphi = cst = L
  
  return 0;
}

void RezzollaZhidenko::circularVelocity(double const coor[4], double vel[4],
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
  double rr=coord[1], r2=rr*rr, r3=r2*rr;

  vel[3] = 1.; // to change
  
  vel[0] = SysPrimeToTdot(coor, vel+1);
  vel[3] *= vel[0];
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG_ARRAY(vel,4);
# endif
}
