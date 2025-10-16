/*
    Copyright 2025 Irene Urso

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

#include "GyotoKonoplyaRezzollaZhidenko.h"

#include <iostream>
#include <cstdlib>
#include <cmath>

#define STRINGIFY(foo) #foo

using namespace Gyoto;
using namespace Gyoto::Metric;
using namespace std;

#define GYOTO_NBPARAM_MAX 6 
// only this number of parameters is allowed, only the delta deformation parameters of Ni et al. 2016 (JCAP09(2016)014) for the time being
// Kerr metric retrieved for the functions given in Table 1 of Cárdenas-Avendaño & Held 2024 (PRD 109, 064052)

/// Properties
#include "GyotoProperty.h"
GYOTO_PROPERTY_START(KonoplyaRezzollaZhidenko, 
                     "Axisymmetric parametrized metric of Konoplya & Rezzolla & Zhidenko 2016")
GYOTO_PROPERTY_DOUBLE(KonoplyaRezzollaZhidenko, Spin, spin,
		      "Spin parameter (adimensioned, 0).")
GYOTO_PROPERTY_DOUBLE(KonoplyaRezzollaZhidenko, Rms, rms,
                      "Radius of the ISCO (geometrical units, 6).")
GYOTO_PROPERTY_VECTOR_DOUBLE(KonoplyaRezzollaZhidenko, DeltaParam, deltaparam,
			     "Delta deformation parameters of Ni et al. 2016 (at most 6 elements)")
GYOTO_PROPERTY_END(KonoplyaRezzollaZhidenko, Generic::properties)

// accessors
GYOTO_PROPERTY_ACCESSORS(KonoplyaRezzollaZhidenko, double, rms_, rms)

void KonoplyaRezzollaZhidenko::deltaparam(std::vector<double> const &v) {
  size_t n = v.size();
  if (n>GYOTO_NBPARAM_MAX)
    GYOTO_ERROR("In KonoplyaRezzollaZhidenko: choose at most "
	       STRINGIFY(GYOTO_NBPARAM_MAX) " parameters");
  double r0 = getr0(), r02 = r0*r0;
  for (size_t i=0; i<n; ++i) {
    deltaparam_[i]=v[i];
    //if (deltaparam_[i]<0.) GYOTO_ERROR("In KonoplyaRezzollaZhidenko: param < 0!");
    switch(i) {
      // See equations 9-12 of Nampalliwar et al. 2021 (PRD, 102, 124071) for the theoretical allowed ranges
      case 0: if (deltaparam_[i] < (4.*r0-3.*r02-spin2_)/r02) 
              GYOTO_ERROR("In KonoplyaRezzollaZhidenko: Restriction range not respected for δ1");
              break;
      case 1: 
      case 2: if (spin_ > 0.){
                if (deltaparam_[i] < -4./spin3_*(1.-sqrt(1.-spin2_))) 
                GYOTO_ERROR("In KonoplyaRezzollaZhidenko: Restriction range not respected for δ2 or δ3");
                break;
              }
              else if (spin_ < 0.){
                if (deltaparam_[i] > -4./spin3_*(1.-sqrt(1.-spin2_))) 
                GYOTO_ERROR("In KonoplyaRezzollaZhidenko: Restriction range not respected for δ2 or δ3");
                break;
              }
       case 3:
       case 4: if (deltaparam_[i] < -1.) 
              GYOTO_ERROR("In KonoplyaRezzollaZhidenko: Restriction range not respected for δ4 or δ5");
              break;
       case 5: if (deltaparam_[i] > r02/(4.-spin2_)) 
              GYOTO_ERROR("In KonoplyaRezzollaZhidenko: Restriction range not respected for δ6");
              // Irene: Given the disagreement in the literature, I don't know if this restriction range still applies...!!!
              break;
    }
  }
  for (size_t i=n; i<GYOTO_NBPARAM_MAX; ++i) deltaparam_[i]=0.;
}
std::vector<double> KonoplyaRezzollaZhidenko::deltaparam() const {
  std::vector<double> v(GYOTO_NBPARAM_MAX, 0.);
  for (size_t i=0; i<GYOTO_NBPARAM_MAX; ++i) v[i]=deltaparam_[i];
  return v;
}

///

Gyoto::Metric::KonoplyaRezzollaZhidenko::KonoplyaRezzollaZhidenko()
  : Generic(GYOTO_COORDKIND_SPHERICAL, "KonoplyaRezzollaZhidenko"),
    spin_(0.), spin2_(0.), spin3_(0.), spin4_(0.),
    rms_(6.), deltaparam_(NULL)
{
  GYOTO_DEBUG << endl;
  deltaparam_ = new double[GYOTO_NBPARAM_MAX];
  for (int ii=0;ii<GYOTO_NBPARAM_MAX;ii++){
    deltaparam_[ii]=0.;
  }
}

Gyoto::Metric::KonoplyaRezzollaZhidenko::KonoplyaRezzollaZhidenko(const KonoplyaRezzollaZhidenko & orig)
  : Generic(orig),
    spin_(orig.spin_), spin2_(orig.spin2_), spin3_(orig.spin4_), spin4_(orig.spin4_),
    rms_(orig.rms_), deltaparam_(NULL)
{
  GYOTO_DEBUG << endl;
  deltaparam_ = new double[GYOTO_NBPARAM_MAX];
  for (int ii=0;ii<GYOTO_NBPARAM_MAX;ii++){
    deltaparam_[ii]=orig.deltaparam_[ii];
  }
}

// default copy constructor should be fine 
KonoplyaRezzollaZhidenko * KonoplyaRezzollaZhidenko::clone() const {
  return new KonoplyaRezzollaZhidenko(*this); }
  
void KonoplyaRezzollaZhidenko::spin(const double val) {
  spin_  = val;
  spin2_ = spin_ * spin_;
  spin3_ = spin2_ * spin_;
  spin4_ = spin2_ * spin2_;
  tellListeners();
}
double KonoplyaRezzollaZhidenko::spin() const { return spin_ ; }

Gyoto::Metric::KonoplyaRezzollaZhidenko::~KonoplyaRezzollaZhidenko()
{
  GYOTO_DEBUG << endl;
  delete [] deltaparam_;
}

double KonoplyaRezzollaZhidenko::getRms() const {
  return rms_;
}

double KonoplyaRezzollaZhidenko::epsilon0Function() const{
  /*
    KRZ metric defines r0 as the exterior event horizon location.
    It is related to epsilon0 through: epsilon0 = 2M/r0 - 1 so that r0 = 2M/(1+espilon0),
    and in Kerr, r0 = r0Kerr = M + sqrt(M^2-a^2)
    where M is ADM mass of spacetime and a=J/M is the specific angular momentum.
    Here it is assumed that units are such that M=1.
  */
  double r0Kerr = 1.+sqrt(1.-spin2_);
  // Kerr case: 2./r0Kerr-1. or a^2/r0Kerr^2
  double epsilon0 = 2./r0Kerr-1.; 
  return epsilon0;
}

double KonoplyaRezzollaZhidenko::getr0() const{
  // Kerr case: 2./(1.+epsilon0Kerr) or 1.+sqrt(1.-a^2)
  double epsilon0 = epsilon0Function();
  double r0 = 2./(1.+epsilon0); 
  return r0;
}

double KonoplyaRezzollaZhidenko::k00Function(const double r0) const{
  // Kerr case: k00 = a^2/r0Kerr^2
  double r02 = r0*r0;
  double k00 = spin2_/r02;
  return k00;
}

double KonoplyaRezzollaZhidenko::k22Function(const double r0) const{
  // Kerr case: k22 = -a^2/r0Kerr^2
  // Irene: Disagreement in the literature!!!!! k22 = a^2/r0Kerr^2
  double r02 = r0*r0;
  double k22 = -spin2_/r02; //spin2_/r02; !!!
  return k22;
}

double KonoplyaRezzollaZhidenko::k23Function(const double r0) const{
  // Kerr case: k23 = -a^2/r0Kerr^2
  double r02 = r0*r0;
  double k23 = spin2_/r02;
  return k23;
}

double KonoplyaRezzollaZhidenko::k21Function(const double r0) const{
  // Kerr case: k21 = -a^2/r0Kerr^2
  // Irene: Disagreement in the literature!!!!! k21 = a^4/r0Kerr^4-2a^2/r0Kerr^3-ẟ5
  double r02 = r0*r0, r03 = r02*r0, r04 = r03*r0;
  double k21 = -spin2_/r02-deltaparam_[5]; //spin4_/r04-2.*spin2_/r03-deltaparam_[5]; !!!
  return k21;
}

double KonoplyaRezzollaZhidenko::a20Function(const double r0) const{
  // Kerr case: a20 = (1+a^2/r0Kerr^2)spin2/r0Kerr^2
  // Irene: Disagreement in the literature!!!!! 2a^2/r0Kerr^3
  double r02 = r0*r0, r03 = r02*r0;
  double a20 = (1.+spin2_/r02)*spin2_/r02; //2.*spin2_/r03; !!!
  return a20;
}

double KonoplyaRezzollaZhidenko::a21Function(const double r0) const{
  // Kerr case: a21 = -a^4/r0Kerr^4
  double r04 = r0*r0*r0*r0;
  double a21 = -spin4_/r04+deltaparam_[5];
  return a21;
}

double KonoplyaRezzollaZhidenko::w00Function(const double r0) const{
  // Kerr case: w00 = (1+a^2/r0Kerr^2)a/r0Kerr
  // Irene: Disagreement in the literature!!!!! 2a/r0Kerr^2
  double r02 = r0*r0;
  double w00 = (1.+spin2_/r02)*spin_/r0; //2.*spin_/r02; !!!
  return w00;
}

double KonoplyaRezzollaZhidenko::N2(const double rr, const double th) const{
  double r0 = getr0(), epsilon0 = epsilon0Function();
  double rr2 = rr*rr, rr3 = rr2*rr, rr4 = rr3*rr,
    r02 = r0*r0, r03 = r02*r0, r04 = r03*r0,
    costh2 = cos(th)*cos(th);
  double k00 = k00Function(r0),
         k22 = k22Function(r0),
         k23 = k23Function(r0),
         k21 = k21Function(r0),
         a20 = a20Function(r0),
         a21 = a21Function(r0);
  double N2 = (1.-r0/rr)*(1.-epsilon0*r0/rr+(k00-epsilon0)*r02/rr2+deltaparam_[0]*r03/rr3)
              +((k21/(1.+k22*(1.-r0/rr)/(1.+k23*(1.-r0/rr)))+a20)*r03/rr3+a21*r04/rr4)*costh2; 

  return N2;
}

double KonoplyaRezzollaZhidenko::B(const double rr, const double th) const{
  double r0 = getr0();
  double rr2 = rr*rr, r02 = r0*r0, costh2 = cos(th)*cos(th);
  double B = 1.+deltaparam_[3]*r02/rr2+deltaparam_[4]*r02/rr2*costh2; 

  return B;
}

double KonoplyaRezzollaZhidenko::Sigma(const double rr, const double th) const{
  double rr2 = rr*rr, costh2 = cos(th)*cos(th);
  double Sigma = 1.+spin2_*costh2/rr2;

  return Sigma;
}

double KonoplyaRezzollaZhidenko::W(const double rr, const double th) const{
  double r0 = getr0();
  double rr2 = rr*rr, rr3 = rr2*rr, 
         r02 = r0*r0, r03 = r02*r0,
         costh2 = cos(th)*cos(th);
  double w00 = w00Function(r0);
  double W = (w00*r02/rr2+deltaparam_[1]*r03/rr3+deltaparam_[2]*r03/rr3*costh2)/Sigma(rr,th);

  return W;
}

double KonoplyaRezzollaZhidenko::K2(const double rr, const double th) const{
  double r0 = getr0();
  double rr2 = rr*rr, rr3 = rr2*rr, 
         r02 = r0*r0, r03 = r02*r0,
         costh2 = cos(th)*cos(th);
  double k00 = k00Function(r0),
         k21 = k21Function(r0),
         k22 = k22Function(r0),
         k23 = k23Function(r0);
  double K2 = 1.+spin_*W(rr,th)/rr+(k00*r02/rr2+k21/(1.+k22*(1.-r0/rr)/(1.+k23*(1.-r0/rr)))*r03/rr3*costh2)/Sigma(rr,th);

  return K2;
}

double KonoplyaRezzollaZhidenko::gmunu(const double * pos, int mu, int nu) const {
  double rr = pos[1], th = pos[2];
  if (rr<=0.) GYOTO_ERROR("In KonoplyaRezzollaZhidenko::gmunu: r<0!");

  double rr2 = rr*rr, sinth2 = sin(th)*sin(th);

  if ((mu==0) && (nu==0)) return -(N2(rr,th)-W(rr,th)*W(rr,th)*sinth2)/K2(rr,th);
  if ((mu==1) && (nu==1)) return Sigma(rr,th)*B(rr,th)*B(rr,th)/N2(rr,th);
  if ((mu==2) && (nu==2)) return Sigma(rr,th)*rr2;
  if ((mu==3) && (nu==3)) return K2(rr,th)*rr2*sinth2;
  if ((mu==0) && (nu==3)) return -W(rr,th)*rr*sinth2;
  if ((mu==3) && (nu==0)) return -W(rr,th)*rr*sinth2;
  return 0.;
} 

double KonoplyaRezzollaZhidenko::gmunu_up(const double * pos, int mu, int nu) const {
  double rr = pos[1], th = pos[2];
  if (rr<=0.) GYOTO_ERROR("In KonoplyaRezzollaZhidenko::gmunu: r<0!");

  double rr2 = rr*rr, sinth2 = sin(th)*sin(th);

  if ((mu==0) && (nu==0)) return -K2(rr,th)/N2(rr,th);
  if ((mu==1) && (nu==1)) return N2(rr,th)/(B(rr,th)*B(rr,th)*Sigma(rr,th));
  if ((mu==2) && (nu==2)) return 1./(rr2*Sigma(rr,th));
  if ((mu==3) && (nu==3)) return -(W(rr,th)*W(rr,th)*sinth2-N2(rr,th))/(rr2*K2(rr,th)*N2(rr,th)*sinth2);
  if ((mu==0) && (nu==3)) return -W(rr,th)/(rr*N2(rr,th));
  if ((mu==3) && (nu==0)) return -W(rr,th)/(rr*N2(rr,th));
  return 0.;
} 

double KonoplyaRezzollaZhidenko::drN2(const double rr, const double th) const{
  double r0 = getr0(), epsilon0 = epsilon0Function();
  double rr2 = rr*rr, rr3 = rr2*rr, rr4 = rr3*rr, rr5 = rr4*rr,
    r02 = r0*r0, r03 = r02*r0, r04 = r03*r0,
    costh2 = cos(th)*cos(th);
  double k00 = k00Function(r0),
         k22 = k22Function(r0),
         k23 = k23Function(r0),
         k21 = k21Function(r0),
         a20 = a20Function(r0),
         a21 = a21Function(r0);
  double drN2 = r0/rr2*(1.-epsilon0*r0/rr+(k00-epsilon0)*r02/rr2+deltaparam_[0]*r03/rr3)
              +(1.-r0/rr)*(epsilon0*r0/rr2-2.*(k00-epsilon0)*r02/rr3-3.*deltaparam_[0]*r03/rr4)
              -(k21*k22*r04/rr5/((1.+(k23+k22)*(1.-r0/rr))*(1.+(k23+k22)*(1.-r0/rr)))
              +3.*(k21/(1.+k22*(1.-r0/rr)/(1.+k23*(1.-r0/rr)))+a20)*r03/rr4
              +4.*a21*r04/rr5)*costh2; 
  return drN2;
}

double KonoplyaRezzollaZhidenko::dthN2(const double rr, const double th) const{
  double r0 = getr0(), epsilon0 = epsilon0Function();
  double rr2 = rr*rr, rr3 = rr2*rr, rr4 = rr3*rr,
    r02 = r0*r0, r03 = r02*r0, r04 = r03*r0,
    costh = cos(th), sinth = sin(th);
  double k00 = k00Function(r0),
         k22 = k22Function(r0),
         k23 = k23Function(r0),
         k21 = k21Function(r0),
         a20 = a20Function(r0),
         a21 = a21Function(r0);
  double dthN2 = -2.*((k21/(1.+k22*(1.-r0/rr)/(1.+k23*(1.-r0/rr)))+a20)*r03/rr3+a21*r04/rr4)*costh*sinth; 
  return dthN2;
}

double KonoplyaRezzollaZhidenko::drB(const double rr, const double th) const{
  double r0 = getr0();
  double rr2 = rr*rr, rr3 = rr2*rr, r02 = r0*r0, costh2 = cos(th)*cos(th);
  double drB = -2.*deltaparam_[3]*r02/rr3-2.*deltaparam_[4]*r02/rr3*costh2; 

  return drB;
}

double KonoplyaRezzollaZhidenko::dthB(const double rr, const double th) const{
  double r0 = getr0();
  double rr2 = rr*rr, r02 = r0*r0, costh = cos(th), sinth = sin(th);
  double dthB = -2.*deltaparam_[4]*r02/rr2*costh*sinth; 

  return dthB;
}

double KonoplyaRezzollaZhidenko::drSigma(const double rr, const double th) const{
  double rr2 = rr*rr, rr3 = rr2*rr, costh2 = cos(th)*cos(th);
  double drSigma = -2.*spin2_*costh2/rr3;

  return drSigma;
}

double KonoplyaRezzollaZhidenko::dthSigma(const double rr, const double th) const{
  double rr2 = rr*rr, costh = cos(th), sinth = sin(th);
  double dthSigma = -2.*spin2_*costh*sinth/rr2;

  return dthSigma;
}

double KonoplyaRezzollaZhidenko::drW(const double rr, const double th) const{
  double r0 = getr0();
  double rr2 = rr*rr, rr3 = rr2*rr, rr4 = rr3*rr,
         r02 = r0*r0, r03 = r02*r0,
         costh2 = cos(th)*cos(th);
  double w00 = w00Function(r0);
  double drW = -drSigma(rr,th)*(w00*r02/rr2+deltaparam_[1]*r03/rr3+deltaparam_[2]*r03/rr3*costh2)/(Sigma(rr,th)*Sigma(rr,th))
               -(2.*w00*r02/rr3+3.*deltaparam_[1]*r03/rr4+3.*deltaparam_[2]*r03/rr4*costh2)/Sigma(rr,th);

  return drW;
}

double KonoplyaRezzollaZhidenko::dthW(const double rr, const double th) const{
  double r0 = getr0();
  double rr2 = rr*rr, rr3 = rr2*rr, 
         r02 = r0*r0, r03 = r02*r0,
         costh = cos(th), costh2 = costh*costh, sinth = sin(th);
  double w00 = w00Function(r0);
  double dthW = -dthSigma(rr,th)*(w00*r02/rr2+deltaparam_[1]*r03/rr3+deltaparam_[2]*r03/rr3*costh2)/(Sigma(rr,th)*Sigma(rr,th))
               -(2.*deltaparam_[2]*r03/rr3*costh*sinth)/Sigma(rr,th);

  return dthW;
}

double KonoplyaRezzollaZhidenko::drK2(const double rr, const double th) const{
  double r0 = getr0();
  double rr2 = rr*rr, rr3 = rr2*rr, rr4 = rr3*rr, rr5 = rr4*rr,
         r02 = r0*r0, r03 = r02*r0, r04 = r03*r0,
         costh2 = cos(th)*cos(th);
  double k00 = k00Function(r0),
         k21 = k21Function(r0),
         k22 = k22Function(r0),
         k23 = k23Function(r0);
  double drK2 = spin_*drW(rr,th)/rr-spin_*W(rr,th)/rr2
               -drSigma(rr,th)*(k00*r02/rr2+k21/(1.+k22*(1.-r0/rr)/(1.+k23*(1.-r0/rr)))*r03/rr3*costh2)/(Sigma(rr,th)*Sigma(rr,th))
               -(2.*k00*r02/rr3
               +k21*k22*costh2*r04/rr5/((1.+(k23+k22)*(1.-r0/rr))*(1.+(k23+k22)*(1.-r0/rr)))
               +3.*(k21/(1.+k22*(1.-r0/rr)/(1.+k23*(1.-r0/rr))))*costh2*r03/rr4)/Sigma(rr,th);

  return drK2;
}

double KonoplyaRezzollaZhidenko::dthK2(const double rr, const double th) const{
  double r0 = getr0();
  double rr2 = rr*rr, rr3 = rr2*rr, 
         r02 = r0*r0, r03 = r02*r0,
         costh = cos(th), costh2 = costh*costh, sinth = sin(th);
  double k00 = k00Function(r0),
         k21 = k21Function(r0),
         k22 = k22Function(r0),
         k23 = k23Function(r0);
  double dthK2 = spin_*dthW(rr,th)/rr
                 -dthSigma(rr,th)*(k00*r02/rr2+k21/(1.+k22*(1.-r0/rr)/(1.+k23*(1.-r0/rr)))*r03/rr3*costh2)/(Sigma(rr,th)*Sigma(rr,th))
                 -2.*k21/(1.+k22*(1.-r0/rr)/(1.+k23*(1.-r0/rr)))*sinth*costh*r03/rr3/Sigma(rr,th);

  return dthK2;
}

int KonoplyaRezzollaZhidenko::christoffel(double dst[4][4][4], double const pos[4]) const
{
  int a, mu, nu;
  for (a=0; a<4; ++a)
    for(mu=0; mu<4; ++mu)
      for(nu=0; nu<4; ++nu)
	dst[a][mu][nu]=0.;

  double rr = pos[1], rr2 = rr*rr,
  th = pos[2], sinth = sin(th), sinth2 = sinth*sinth, sinth3 = sinth2*sinth, costh = cos(th), costh2 = costh*costh;
  if (rr==0. || sinth==0.) GYOTO_ERROR("In KonoplyaRezzollaZhidenko::christoffel: bad coord");
				    
  double NN2=N2(rr,th), NN22=NN2*NN2, drNN2 = drN2(rr,th), dthNN2 = dthN2(rr,th),
  BB=B(rr,th), BB2=BB*BB, drBB = drB(rr,th), dthBB = dthB(rr,th),
  SSigma=Sigma(rr,th), drSSigma = drSigma(rr,th), dthSSigma = dthSigma(rr,th),
  WW=W(rr,th), WW2 = WW*WW, WW3 = WW2*WW, drWW = drW(rr,th), dthWW = dthW(rr,th),
  KK2=K2(rr,th), KK22=KK2*KK2, drKK2 = drK2(rr,th), dthKK2 = dthK2(rr,th);

  dst[0][0][1] = dst[0][1][0] = ((rr*WW2*drKK2-rr*KK2*WW*drWW+KK2*WW2)*sinth2-rr*NN2*drKK2+rr*KK2*drNN2)/(2.*rr*KK2*NN2); 
  dst[0][0][2] = dst[0][2][0] = ((WW2*dthKK2-KK2*WW*dthWW)*sinth2-NN2*dthKK2+KK2*dthNN2)/(2.*KK2*NN2);
  dst[0][1][3] = dst[0][3][1] = -(rr*WW*drKK2-rr*KK2*drWW+KK2*WW)*sinth2/(2.*NN2);
  dst[0][2][3] = dst[0][3][2] = -(rr*WW*dthKK2-rr*KK2*dthWW)*sinth2/(2.*NN2);
  dst[1][0][0] = ((NN2*WW2*drKK2-2.*KK2*NN2*WW*drWW)*sinth2-NN22*drKK2+KK2*NN2*drNN2)/(2.*BB2*KK22*SSigma); 
  dst[1][0][3] = dst[1][3][0] = (rr*NN2*drWW+NN2*WW)*sinth2/(2.*BB2*SSigma);
  dst[1][1][1] = (2.*NN2*SSigma*drBB-BB*SSigma*drNN2+BB*NN2*drSSigma)/(2.*BB*NN2*SSigma);
  dst[1][1][2] = dst[1][2][1] = (2.*NN2*SSigma*dthBB-BB*SSigma*dthNN2+BB*NN2*dthSSigma)/(2.*BB*NN2*SSigma);
  dst[1][2][2] = -(rr2*NN2*drSSigma+2.*rr*NN2*SSigma)/(2.*BB2*SSigma); 
  dst[1][3][3] = -(rr2*NN2*drKK2+2.*rr*KK2*NN2)*sinth2/(2.*BB2*SSigma); 
  dst[2][0][0] = -(2.*KK2*WW2*costh*sinth-(WW2*dthKK2-2.*KK2*WW*dthWW)*sinth2+NN2*dthKK2-KK2*dthNN2)/(2.*rr2*KK22*SSigma);
  dst[2][0][3] = dst[2][3][0] = (2.*WW*costh*sinth+sinth2*dthWW)/(2.*rr*SSigma);
  dst[2][1][1] = -(2.*BB*NN2*SSigma*dthBB-BB2*SSigma*dthNN2+BB2*NN2*dthSSigma)/(2.*rr2*NN22*SSigma);
  dst[2][1][2] = dst[2][2][1] = (rr*drSSigma+2.*SSigma)/(2.*rr*SSigma); 
  dst[2][2][2] = dthSSigma/(2.*SSigma);
  dst[2][3][3] = -(2.*KK2*costh*sinth+sinth2*dthKK2)/(2.*SSigma); 
  dst[3][0][1] = dst[3][1][0] = -(rr*NN2*WW*drKK2-rr*KK2*WW*drNN2+rr*KK2*NN2*drWW+KK2*NN2*WW-(rr*WW3*drKK2-rr*KK2*WW2*drWW+KK2*WW3)*sinth2)/(2.*rr2*KK22*NN2);
  dst[3][0][2] = dst[3][2][0] = -(2.*KK2*NN2*WW*costh-(WW3*dthKK2-KK2*WW2*dthWW)*sinth3+(NN2*WW*dthKK2-KK2*WW*dthNN2+KK2*NN2*dthWW)*sinth)/(2.*rr*KK22*NN2*sinth);
  dst[3][1][3] = dst[3][3][1] = -((rr*WW2*drKK2-rr*KK2*WW*drWW+KK2*WW2)*sinth2-rr*NN2*drKK2-2.*KK2*NN2)/(2.*rr*KK2*NN2); 
  dst[3][2][3] = dst[3][3][2] = -((WW2*dthKK2-KK2*WW*dthWW)*sinth3-2.*KK2*NN2*costh-NN2*sinth*dthKK2)/(2.*KK2*NN2*sinth);
  
  int compareKerr=0;
  if (compareKerr){
    double r2plusa2 = rr2+spin2_;
    double sin2th = 2.*sinth*costh, cotanth=costh/sinth, a2costhsinth=spin2_*costh*sinth;
    double SigmaK=rr2+spin2_*costh2, SigmaK2=SigmaK*SigmaK;
    double Delta=rr2-2.*rr+spin2_;
    double Deltam1=1./Delta,
      SigmaKm1=1./SigmaK, SigmaKm2=SigmaKm1*SigmaKm1, SigmaKm3=SigmaKm2*SigmaKm1;
    double rSigmaKm1=rr*SigmaKm1,
      Deltam1SigmaKm2=Deltam1*SigmaKm2;
    // Comparison to Kerr Christoffels for checking
    double dstK[4][4][4];
    for (a=0; a<4; ++a)
      for(mu=0; mu<4; ++mu)
	for(nu=0; nu<4; ++nu)
	  dstK[a][mu][nu]=0.;
    dstK[1][1][1]=(1.-rr)*Deltam1+rSigmaKm1;
    dstK[1][2][1]=dstK[1][1][2]=-a2costhsinth*SigmaKm1;
    dstK[1][2][2]=-Delta*rSigmaKm1;
    dstK[1][3][3]=-Delta*sinth2*(rr+(spin2_*(-2.*rr2+SigmaK)*sinth2)/SigmaK2)/SigmaK; 
    dstK[1][3][0]=dstK[1][0][3]=spin_*Delta*(-2*rr2+SigmaK)*sinth2*SigmaKm3;
    dstK[1][0][0]=-Delta*(-2.*rr2+SigmaK)*SigmaKm3;
    dstK[2][1][1]=a2costhsinth*Deltam1*SigmaKm1;
    dstK[2][2][1]=dstK[2][1][2]=rSigmaKm1;
    dstK[2][2][2]=-a2costhsinth*SigmaKm1;
    dstK[2][3][3]=-sinth*costh*SigmaKm3 * (Delta*SigmaK2 + 2.*rr*r2plusa2*r2plusa2);
    dstK[2][0][3]=dstK[2][3][0]=spin_*rr*r2plusa2*sin2th*SigmaKm3;
    dstK[2][0][0]=-2.*a2costhsinth*rr*SigmaKm3;
    dstK[3][3][1]=dstK[3][1][3]=Deltam1*SigmaKm2 * (rr*SigmaK*(SigmaK-2.*rr) + spin2_*(SigmaK-2.*rr2)*sinth2);
    dstK[3][3][2]=dstK[3][2][3]=SigmaKm2*cotanth * (-(SigmaK+Delta)*spin2_*sinth2 + r2plusa2*r2plusa2);
    dstK[3][0][1]=dstK[3][1][0]=spin_*(2.*rr2-SigmaK)*Deltam1SigmaKm2;
    dstK[3][0][2]=dstK[3][2][0]=-2.*spin_*rr*cotanth*SigmaKm2;
    dstK[0][3][1]=dstK[0][1][3]=-spin_*sinth2*Deltam1SigmaKm2 * (2.*rr2*r2plusa2 + SigmaK*(rr2-spin2_));
    dstK[0][3][2]=dstK[0][2][3]=SigmaKm2*spin_*spin2_*rr*sinth2*sin2th;
    dstK[0][0][1]=dstK[0][1][0]=(spin2_+rr2)*(2.*rr2-SigmaK)*Deltam1SigmaKm2;
    dstK[0][0][2]=dstK[0][2][0]=-spin2_*rr*sin2th*SigmaKm2;
    for (a=0; a<4; ++a){
      for(mu=0; mu<4; ++mu){
	for(nu=0; nu<4; ++nu){
	  if(fabs(dst[a][mu][nu]-dstK[a][mu][nu])>1e-4) {
	    cout << "at r,th= " << rr <<  " " << pos[2] << endl;
	    cout << setprecision(20)<< a << " " << mu << " " << nu << " " << dst[a][mu][nu] << " " << dstK[a][mu][nu]<< " " << fabs(dst[a][mu][nu]-dstK[a][mu][nu]) <<endl;
	    GYOTO_ERROR("test Kerr KRZ");
	  }
	}
      }
    }
  } 
  return 0;
}

int KonoplyaRezzollaZhidenko::isStopCondition(double const * const coord) const {
  double r0 = getr0();
  double rsink = r0 + GYOTO_KERR_HORIZON_SECURITY;
  return coord[1] < rsink ;
}

double KonoplyaRezzollaZhidenko::KeplerianSpecificAngularMomentum(const double rr, const double th) const {
  double rr2 = rr*rr, rr3 = rr2*rr, 
         sinth = sin(th), sinth2 = sinth*sinth;
  double NN2=N2(rr,th), NN22 = NN2*NN2, drNN2 = drN2(rr,th),
         WW=W(rr,th), WW2 = WW*WW, drWW = drW(rr,th),
         KK2=K2(rr,th), KK22=KK2*KK2, drKK2 = drK2(rr,th);
  double drguptt = -(drKK2*NN2-KK2*drNN2)/NN22,
         drguptp = -(drWW*rr*NN2-WW*(NN2+rr*drNN2))/(rr2*NN22),
         drguppp = -(2.*NN22*KK2+rr*NN22*drKK2-(2.*WW2*KK2*NN2+rr*WW2*(drKK2*NN2+KK2*drNN2)-2.*rr*KK2*NN2*WW*drWW)*sinth2)/(rr3*KK22*NN22*sinth2);
  double ell1 = (drguptp+sqrt(drguptp*drguptp-drguptt*drguppp))/drguppp,
         ell2 = (drguptp-sqrt(drguptp*drguptp-drguptt*drguppp))/drguppp; 
  double gtt = -(NN2-WW2*sinth2)/KK2,
         gtp = -WW*rr*sinth2,
         gpp = KK2*rr2*sinth2;
  // Positive Omega for for a prograde motion
  double Omega1 = -(gtp+ell1*gtt)/(gpp+ell1*gtp),
         Omega2 = -(gtp+ell2*gtt)/(gpp+ell2*gtp);
  double ell = (Omega1 > 0) ? ell1 : ell2;       
  return ell;
} 

double KonoplyaRezzollaZhidenko::KeplerianAngularVelocity(double const rr, double const th, const double ell) const {
  double rr2 = rr*rr, sinth = sin(th), sinth2 = sinth*sinth;
  double NN2=N2(rr,th), WW=W(rr,th), WW2 = WW*WW, KK2=K2(rr,th);
  double gtt = -(NN2-WW2*sinth2)/KK2,
         gtp = -WW*rr*sinth2,
         gpp = KK2*rr2*sinth2;
  double Omega = -(gtp+ell*gtt)/(gpp+ell*gtp);
  
  
  /*double rr3 = rr2*rr;
  double NN22 = NN2*NN2, drNN2 = drN2(rr,th),
         WW2 = WW*WW, drWW = drW(rr,th),
         KK22=KK2*KK2, drKK2 = drK2(rr,th);
  double drgtt = -((drNN2-2.*WW*drWW*sinth2)*KK2-(NN2-WW2*sinth2)*drKK2)/KK22,
         drgtp = -(drWW*rr+WW)*sinth2,
         drgpp = (drKK2*rr2+2.*KK2*rr)*sinth2;
  double Omega_test = (-drgtp+sqrt(drgtp*drgtp-drgtt*drgpp))/drgpp;
  cout << "Omega_test, Omega: " << Omega_test << ", " << Omega << endl;*/
  
  return Omega;
} 

double KonoplyaRezzollaZhidenko::KeplerianEnergy(double const rr, double const th, const double Omega) const { 
  double rr2 = rr*rr, sinth = sin(th), sinth2 = sinth*sinth;
  double NN2=N2(rr,th), WW=W(rr,th), WW2 = WW*WW, KK2=K2(rr,th);
  double gtt = -(NN2-WW2*sinth2)/KK2,
         gtp = -WW*rr*sinth2,
         gpp = KK2*rr2*sinth2;
  double E = -(gtt+Omega*gtp)/sqrt(-gtt-2.*gtp*Omega-gpp*Omega*Omega);
  return E;
} 

double KonoplyaRezzollaZhidenko::KeplerianAngularMomentum(double const rr, double const th, const double Omega) const {
  double rr2 = rr*rr, sinth = sin(th), sinth2 = sinth*sinth;
  double NN2=N2(rr,th), WW=W(rr,th), WW2 = WW*WW, KK2=K2(rr,th);
  double gtt = -(NN2-WW2*sinth2)/KK2,
         gtp = -WW*rr*sinth2,
         gpp = KK2*rr2*sinth2;
  double L = (gtp+Omega*gpp)/sqrt(-gtt-2.*gtp*Omega-gpp*Omega*Omega);
  return L;
} 

void KonoplyaRezzollaZhidenko::circularVelocity(double const coor[4], double vel[4],
					double dir) const {
  if (keplerian_) {
    // If keplerian_ is true, let the generic implementation return
    // the Newtonian Keplerian (i.e. r^-3/2, no spin) velocity 
    // instead of the true circular velocity
    Generic::circularVelocity(coor, vel, dir);
    return;
  }  
  
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG<<"coor=["<<coor[0]<<", "<<coor[1]<<", "<<coor[2]<<", "<<coor[3]
	     <<"], dir="<<dir<<endl;
# endif
  double th = coor[2], sinth = sin(th), rr = coor[1]*sinth;// rr projected on equat plane // Irene: why?? 
  double coord[4] = {coor[0],rr,M_PI*0.5,coor[3]};
  double ell = KeplerianSpecificAngularMomentum(rr,th),
         Omega = KeplerianAngularVelocity(rr,th,ell);

  vel[1] = vel[2] = 0.;
  // In Kerr: vel[3] = 1./((dir*pow(coord[1], 1.5) + spin_)); // Irene: dir?? 
  vel[3] = Omega; //sqrt(2./rr)*sqrt((drNN2-2.*WW*drWW-(NN2-WW2)*drKK2)/KK22); // this is Omega=dphi/dt

  vel[0] = SysPrimeToTdot(coord, vel+1); // dt/dtau
  vel[3] *= vel[0]; // dphi/dtau
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG_ARRAY(vel,4);
# endif
}
