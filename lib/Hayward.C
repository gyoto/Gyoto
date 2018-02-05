/*
 Copyright 2011-2016 Frederic Vincent, Thibaut Paumard

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
#include "GyotoHayward.h"
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

GYOTO_PROPERTY_START(Hayward,
                     "Metric around a rotating black-hole, in spherical Boyer-Lindquist coordinates.")
GYOTO_PROPERTY_DOUBLE(Hayward, Spin, spin,
                      "Spin parameter (adimensioned, 0).")
GYOTO_PROPERTY_DOUBLE(Hayward, Charge, charge,
                      "Charge parameter (adimensioned, 0).")
GYOTO_PROPERTY_DOUBLE(Hayward, HorizonSecurity, horizonSecurity,
                      "Thickness of sink layer around horizon (geometrical units, 0.01).")
GYOTO_PROPERTY_BOOL(Hayward, GenericIntegrator, SpecificIntegrator,
                    genericIntegrator,
                    "Which version of the Legacy integrator should be used (specific).")
GYOTO_PROPERTY_DOUBLE(Hayward, DiffTol, difftol,
                      "Tuning parameter for the specific Legacy integrator (0.01).")
GYOTO_PROPERTY_END(Hayward, Generic::properties)

/*
 NB: points delicats de Hayward:
 - z-axis problems (d'ou interet de KS? a creuser)
 - probleme pres de l'horizon (d'ou interet de KS? a creuser)
 - changement de thetadot dans CheckCons (mais: test dans RK4ada pour
 eviter de "trop" changer)
 - changement de rdot dans Normalize4v (mais: test dans RK4ada pour
 eviter de "trop" changer)
 */

/*
 Comment on z-axis problem:
 z-axis pb is never launched by the condition theta<Gyoto_min_theta,
 it is always launched by the derlim, thetaaol tests in diff()
 (see comments there)

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


Hayward::Hayward() :
Generic(GYOTO_COORDKIND_SPHERICAL, "Hayward"),
spin_(0.), a2_(0.), a3_(0.), a4_(0.), charge_(0.), b2_(0.),
difftol_(GYOTO_Hayward_DEFAULT_DIFFTOL),
rsink_(2.+GYOTO_KERR_HORIZON_SECURITY),
drhor_(GYOTO_KERR_HORIZON_SECURITY), generic_integrator_(false)
{}

// default copy constructor should be fine
Hayward * Hayward::clone () const { return new Hayward(*this); }

// Mutators
void Hayward::spin(const double a) {
    spin_=a;
    a2_=spin_*spin_;
    a3_=a2_*spin_;
    a4_=a2_*a2_;
//    rsink_=-100.+drhor_;    // rsink_=1.+sqrt(1.-a2_)+drhor_;
    tellListeners();
}

void Hayward::charge(const double b) {
    charge_=b;
    b2_=charge_*charge_;
    tellListeners();
}

// Accessors
double Hayward::spin() const { return spin_ ; }
double Hayward::charge() const { return charge_ ; }
double Hayward::difftol() const { return difftol_;}
void Hayward::difftol(double t) {difftol_=t;}

void Hayward::horizonSecurity(const double drhor) {
    drhor_=drhor;
//    rsink_=-100.+drhor_; //rsink_=1.+sqrt(1.-a2_)+drhor_;

    tellListeners();
}
double Hayward::horizonSecurity() const {return drhor_; }

void Hayward::genericIntegrator(bool t)
{
    generic_integrator_=t;
    tellListeners();
}
bool Hayward::genericIntegrator() const {return generic_integrator_;}

int Hayward::isStopCondition(double const * const coord) const {
    return 0.; // coord[1] < rsink_ ;
}

//Prograde marginally stable orbit
double Hayward::getRms() const {
    double aa=spin_;
    double  z1 = 1. + pow((1. - a2_),1./3.)*(pow((1. + aa),1./3.) + pow((1. - aa),1./3.));
    double  z2 = pow(3.*a2_ + z1*z1,1./2.);

    return (3. + z2 - pow((3. - z1)*(3. + z1 + 2.*z2),1./2.));
}

//Prograde marginally bound orbit
double Hayward::getRmb() const {
    return 2.-spin_+2.*sqrt(1.-spin_);
}

double Hayward::getSpecificAngularMomentum(double rr) const {
    // this is l = -u_phi/u_t for a circular equatorial 4-velocity
    double aa=spin_, sqrtr=sqrt(rr);
    return (rr*rr-2.*aa*sqrtr+aa*aa)/(pow(rr,1.5)-2.*sqrtr+aa);
}

double Hayward::getPotential(double const pos[4], double l_cst) const {
    // this is W = -ln(|u_t|) for a circular equatorial 4-velocity
    double  gtt = gmunu(pos,0,0);
    double  gtp = gmunu(pos,0,3);
    double  gpp = gmunu(pos,3,3);
    double  Omega = -(gtp + l_cst * gtt)/(gpp + l_cst * gtp) ;

    double  W = 0.5 * log(abs(gtt + 2. * Omega * gtp + Omega*Omega * gpp))
    - log(abs(gtt + Omega * gtp)) ;
    return  W ;
}

//Computation of metric coefficients in covariant form

void Hayward::gmunu(double g[4][4], const double * pos) const
{
    double r = pos[1];
    double sth2, cth2;
    sincos(pos[2], &sth2, &cth2);
    sth2*=sth2; cth2*=cth2;

    int mu, nu;
    for (mu=0; mu<4; ++mu)
        for (nu=0; nu<4; ++nu)
            g[mu][nu]=0.;


    if (r>=1.)
            {
            double u=1./r, u2=u*u, u3=u2*u, u4=u3*u, u5=u4*u, u6=u5*u, u7=u6*u;
            double a2b2=a2_*b2_;

            g[0][0] =-(2.*a2b2*u5*cth2+a2_*u2*cth2+2.*b2_*u3-2.*u+1.)/(a2_*u2*cth2+1.)/(2.*b2_*u3+1.);
            g[1][1] =(a2_*u2*cth2+1.)*(2.*b2_*u3+1.)/(2.*a2b2*u5+2.*b2_*u3+a2_*u2-2.*u+1.);
            g[2][2] =(a2_*u2*cth2+1.)/(u2); //pb de 1/u2
            g[3][3] =(2.*a4_*b2_*u7*cth2+2.*a2b2*u5*cth2+a4_*u4*cth2+2.*a2b2*u5+2.*a2_*u3*sth2+a2_*u2*cth2+2.*b2_*u3+a2_*u2+1.)*sth2/(a2_*u2*cth2+1.)/(2.*b2_*u3+1.)/(u2);
            g[0][3] = g[3][0] =-2.*spin_*u*sth2/(a2_*u2*cth2+1.)/(2.*b2_*u3+1.);
            }

    if (r>=0. && r<1.)
            {
            double r2=r*r, r3=r2*r;
            double m=r3/(r3+2.*b2_);
            double sigma=r2+a2_*cth2;
            double delta=r2-2.*m*r+a2_;

            g[0][0] = -1.+2.*m*r/sigma;
            g[1][1] = sigma/delta;
            g[2][2] = sigma;
            g[3][3] = (r2+a2_+2.*m*r*a2_*sth2/sigma)*sth2;
            g[0][3] = g[3][0] = -2.*spin_*m*r*sth2/sigma;
            }

    if (r<0.)
            {
            double r2=r*r, r3=r2*r;
            double m=-r3/(-r3+2.*b2_);
            double sigma=r2+a2_*cth2;
            double delta=r2-2.*m*r+a2_;

            g[0][0] = -1.+2.*m*r/sigma;
            g[1][1] = sigma/delta;
            g[2][2] = sigma;
            g[3][3] = (r2+a2_+2.*m*r*a2_*sth2/sigma)*sth2;
            g[0][3] = g[3][0] = -2.*spin_*m*r*sth2/sigma;
            }
}




double Hayward::gmunu(const double * pos, int mu, int nu) const
{
    double r = pos[1];
    double sth2, cth2;
           sincos(pos[2], &sth2, &cth2);
    sth2*=sth2; cth2*=cth2;

    if (r>=1.)
            {
            double u=1./r, u2=u*u, u3=u2*u, u4=u3*u, u5=u4*u, u6=u5*u, u7=u6*u;
            double a2b2=a2_*b2_;


            if ((mu==0) && (nu==0)) return -(2.*a2b2*u5*cth2+a2_*u2*cth2+2.*b2_*u3-2.*u+1.)/(a2_*u2*cth2+1.)/(2.*b2_*u3+1.);
            if ((mu==1) && (nu==1)) return (a2_*u2*cth2+1.)*(2.*b2_*u3+1.)/(2.*a2b2*u5+2.*b2_*u3+a2_*u2-2.*u+1.);
            if ((mu==2) && (nu==2)) return (a2_*u2*cth2+1.)/(u2);
            if ((mu==3) && (nu==3)) return (2.*a4_*b2_*u7*cth2+2.*a2b2*u5*cth2+a4_*u4*cth2+2.*a2b2*u5+2.*a2_*u3*sth2+a2_*u2*cth2+2.*b2_*u3+a2_*u2+1.)*sth2/(a2_*u2*cth2+1.)/(2.*b2_*u3+1.)/(u2);
            if (((mu==0) && (nu==3)) || ((mu==3) && (nu==0))) return -2.*spin_*u*sth2/(a2_*u2*cth2+1.)/(2.*b2_*u3+1.);
            }

    if (r>=0 && r<1.)
    {
            double r2=r*r, r3=r2*r, r4=r2*r2;
            double d=r3+2.*b2_;
            double m=r3/d;
            double sigma=r2+a2_*cth2;
            double delta=r2-2.*m*r+a2_;

            if ((mu==0) && (nu==0)) return -1.+2.*m*r/sigma;
            if ((mu==1) && (nu==1)) return sigma/delta;
            if ((mu==2) && (nu==2)) return sigma;
            if ((mu==3) && (nu==3)) return (r2+a2_+2.*m*r*a2_*sth2/sigma)*sth2;
            if (((mu==0) && (nu==3)) || ((mu==3) && (nu==0))) return  -2.*spin_*m*r*sth2/sigma;
    }

    if (r<0.)
            {
            double r2=r*r, r3=r2*r;
            double m=-r3/(-r3+2.*b2_);
            double sigma=r2+a2_*cth2;
            double delta=r2-2.*m*r+a2_;

            if ((mu==0) && (nu==0)) return -1.+2.*m*r/sigma;
            if ((mu==1) && (nu==1)) return sigma/delta;
            if ((mu==2) && (nu==2)) return sigma;
            if ((mu==3) && (nu==3)) return (r2+a2_+2.*m*r*a2_*sth2/sigma)*sth2;
            if (((mu==0) && (nu==3)) || ((mu==3) && (nu==0))) return  -2.*spin_*m*r*sth2/sigma;
            }


    return 0.;

}


//Computation of metric coefficients in contravariant form

void Hayward::gmunu_up(double gup[4][4], const double * pos) const
{
            double r = pos[1];
            double sth2, cth2;
                   sincos(pos[2], &sth2, &cth2);
                   sth2*=sth2; cth2*=cth2;
            double a2b2=a2_*b2_, a4b2=a2_*a2b2;

    int mu, nu;
    for (mu=0; mu<4; ++mu) for (nu=0; nu<4; ++nu) gup[mu][nu]=0.;

    if (r>=1.)
            {
            double u=1./r, u2=u*u, u3=u2*u, u4=u3*u, u5=u4*u, u6=u5*u, u7=u6*u;

            gup[0][0]=-(2.*a4_*b2_*u7*cth2+2.*a2b2*u5*cth2+a4_*u4*cth2+2.*a2b2*u5+2.*a2_*u3*sth2+a2_*u2*cth2+2.*b2_*u3+a2_*u2+1.)/(2.*a2b2*u5+2.*b2_*u3+a2_*u2-2.*u+1.)/(a2_*u2*cth2+1.);
            gup[1][1]=(2.*a2b2*u5+2.*b2_*u3+a2_*u2-2.*u+1.)/(a2_*u2*cth2+1.)/(2.*b2_*u3+1.);
            gup[2][2]=u2/(a2_*u2*cth2+1.);
            gup[3][3]=(2.*a2b2*u5*cth2+2.*b2_*u3+a2_*u2*cth2-2.*u+1.)*u2/(2.*a4_*b2_*u7*cth2+2.*a2b2*u5*cth2+a4_*u4*cth2+2.*a2b2*u5+2.*a2_*u3*sth2+a2_*u2*cth2-2.*a2_*u3+2.*b2_*u3+a2_*u2-2.*u+1.)/(sth2);
            gup[0][3] = gup[3][0] = -2.*spin_*u3/(2.*a2b2*u5+2.*b2_*u3+a2_*u2-2.*u+1.)/(a2_*u2*cth2+1.);
            }

    if (r>=0. && r<1.)
            {
            double r2=r*r, r3=r2*r, r4=r2*r2, r5=r4*r, r6=r5*r, r7=r6*r;
            double d=r3+2.*b2_;
            double m=r3/d;
            double sigma=r2+a2_*cth2;
            double delta=r2-2.*m*r+a2_;

            gup[0][0] = -(a4_*r3*cth2+a2_*r5*cth2+2.*a4b2*cth2+2.*a2b2*r2*cth2+2.*a2_*r4*sth2+a2_*r5+r7+2.*a2b2*r2+2.*b2_*r4)/sigma/(a2_*r3+r5+2.*a2b2+2.*b2_*r2-2.*r4);
            gup[1][1] = (a2_*r3+r5+2.*a2b2+2.*b2_*r2-2.*r4)/sigma/d;
            gup[2][2] = 1./sigma;
            gup[3][3] = (a2_*r3*cth2+2.*a2b2*cth2+r5+2.*b2_*r2-2.*r4)/(a4_*r3*cth2+a2_*r5*cth2+2.*a4b2*cth2+2.*a2b2*r2*cth2+2.*a2_*r4*sth2+a2_*r5+r7+2.*a2b2*r2-2.*a2_*r4+2.*b2_*r4-2.*r6)/sth2;
            gup[0][3] = gup[3][0] = -2.*spin_*r4/sigma/(a2_*r3+r5+2.*a2b2+2.*b2_*r2-2.*r4);
            }

    if (r<0.)
            {
            double r2=r*r, r3=r2*r, r4=r2*r2, r5=r4*r, r6=r5*r, r7=r6*r;
            double m=-r3/(-r3+2.*b2_);
            double sigma=r2+a2_*cth2;
            double delta=r2-2.*m*r+a2_;
            double dm=r3-2.*b2_;

            gup[0][0] = -(a4_*r3*cth2+a2_*r5*cth2-2.*a4b2*cth2-2.*a2b2*r2*cth2+2.*a2_*r4*sth2+a2_*r5+r7-2.*a2b2*r2-2.*b2_*r4)/sigma/(a2_*r3+r5-2.*a2b2-2.*b2_*r2-2.*r4);
            gup[1][1] = (a2_*r3+r5-2.*a2b2-2.*b2_*r2-2.*r4)/sigma/dm;
            gup[2][2] = 1./sigma;
            gup[3][3] = (a2_*r3*cth2-2.*a2b2*cth2+r5-2.*b2_*r2-2.*r4)/(a4_*r3*cth2+a2_*r5*cth2-2.*a4b2*cth2-2.*a2b2*r2*cth2+2.*a2_*r4*sth2+a2_*r5+r7-2.*a2b2*r2-2.*a2_*r4-2.*b2_*r4-2.*r6)/sth2;
            gup[0][3] = gup[3][0] = -2.*spin_*r4/sigma/(a2_*r3+r5-2.*a2b2-2.*b2_*r2-2.*r4);
            }


//    cout << "g^00=" << gup[0][0] << endl ;
//    cout << "g^11=" << gup[1][1] << endl ;
//    cout << "g^22=" << gup[2][2] << endl ;
//    cout << "g^33=" << gup[3][3] << endl ;
//    cout << "g^03=" << gup[0][3] << endl ;

}


double Hayward::gmunu_up(const double * pos, int mu, int nu) const
{
            double r = pos[1];
            double sth2, cth2;
                   sincos(pos[2], &sth2, &cth2);
                   sth2*=sth2; cth2*=cth2;
            double a2b2=a2_*b2_, a4b2=a2_*a2b2;

    if (r>=1.)
            {
            double u=1./r, u2=u*u, u3=u2*u, u4=u3*u, u5=u4*u, u6=u5*u, u7=u6*u;

            if ((mu==0) && (nu==0)) return -(2.*a4_*b2_*u7*cth2+2.*a2b2*u5*cth2+a4_*u4*cth2+2.*a2b2*u5+2.*a2_*u3*sth2+a2_*u2*cth2+2.*b2_*u3+a2_*u2+1.)/(2.*a2b2*u5+2.*b2_*u3+a2_*u2+2.*u+1.)/(a2_*u2*cth2+1.);
            if ((mu==1) && (nu==1)) return (2.*a2b2*u5+2.*b2_*u3+a2_*u2-2.*u+1.)/(a2_*u2*cth2+1.)/(2.*b2_*u3+1.);
            if ((mu==2) && (nu==2)) return u2/(a2_*u2*cth2+1.);
            if ((mu==3) && (nu==3)) return (2.*a2b2*u5*cth2+2.*b2_*u3+a2_*u2*cth2-2.*u+1.)*u2/(2.*a4_*b2_*u7*cth2+2.*a2b2*u5*cth2+a4_*u4*cth2+2.*a2b2*u5+2.*a2_*u3*sth2+a2_*u2*cth2-2.*a2_*u3+2.*b2_*u3+a2_*u2-2.*u+1.)/(sth2);
            if (((mu==0) && (nu==3)) || ((mu==3) && (nu==0)))return -2.*spin_*u3/(2.*a2b2*u5+2.*b2_*u3+a2_*u2-2.*u+1.)/(a2_*u2*cth2+1.);
            }

    if (r>=0. && r<1.)
            {
            double r2=r*r, r3=r2*r, r4=r2*r2, r5=r4*r, r6=r5*r, r7=r6*r;
            double d=r3+2.*b2_;
            double m=r3/d;
            double sigma=r2+a2_*cth2;
            double delta=r2-2.*m*r+a2_;

            if ((mu==0) && (nu==0)) return -(a4_*r3*cth2+a2_*r5*cth2+2.*a4b2*cth2+2.*a2b2*r2*cth2+2.*a2_*r4*sth2+a2_*r5+r7+2.*a2b2*r2+2.*b2_*r4)/sigma/(a2_*r3+r5+2.*a2b2+2.*b2_*r2-2.*r4);
            if ((mu==1) && (nu==1)) return (a2_*r3+r5+2.*a2b2+2.*b2_*r2-2.*r4)/sigma/d;
            if ((mu==2) && (nu==2)) return 1./sigma;
            if ((mu==3) && (nu==3)) return (a2_*r3*cth2+2.*a2b2*cth2+r5+2.*b2_*r2-2.*r4)/(a4_*r3*cth2+a2_*r5*cth2+2.*a4b2*cth2+2.*a2b2*r2*cth2+2.*a2_*r4*sth2+a2_*r5+r7+2.*a2b2*r2-2.*a2_*r4+2.*b2_*r4-2.*r6)/sth2;
            if (((mu==0) && (nu==3)) || ((mu==3) && (nu==0)))return -2*spin_*r4/sigma/(a2_*r3+r5+2.*a2b2+2.*b2_*r2-2.*r4);
            }

    if (r<0.)
            {
            double r2=r*r, r3=r2*r, r4=r2*r2, r5=r4*r, r6=r5*r, r7=r6*r;
            double m=-r3/(-r3+2.*b2_);
            double sigma=r2+a2_*cth2;
            double delta=r2-2.*m*r+a2_;
            double dm=r3-2.*b2_;

            if ((mu==0) && (nu==0)) return -(a4_*r3*cth2+a2_*r5*cth2-2.*a4b2*cth2-2.*a2b2*r2*cth2+2.*a2_*r4*sth2+a2_*r5+r7-2.*a2b2*r2-2.*b2_*r4)/sigma/(a2_*r3+r5-2.*a2b2-2.*b2_*r2-2.*r4);
            if ((mu==1) && (nu==1)) return (a2_*r3+r5-2.*a2b2-2.*b2_*r2-2.*r4)/sigma/dm;
            if ((mu==2) && (nu==2)) return 1./sigma;
            if ((mu==3) && (nu==3)) return (a2_*r3*cth2-2.*a2b2*cth2+r5-2.*b2_*r2-2.*r4)/(a4_*r3*cth2+a2_*r5*cth2-2.*a4b2*cth2-2.*a2b2*r2*cth2+2.*a2_*r4*sth2+a2_*r5+r7-2.*a2b2*r2-2.*a2_*r4-2.*b2_*r4-2.*r6)/sth2;
            if (((mu==0) && (nu==3)) || ((mu==3) && (nu==0)))return  -2.*spin_*r4/sigma/(a2_*r3+r5+2.*a2b2+2.*b2_*r2-2.*r4);
            }

    return 0.;
}



int Hayward::christoffel(double dst[4][4][4], double const pos[4]) const
{
    int a, mu, nu;
    for (a=0; a<4; ++a)
        for(mu=0; mu<4; ++mu)
            for(nu=0; nu<4; ++nu)
                dst[a][mu][nu]=0.;

                double r = pos[1];
                double sth, cth;
                       sincos(pos[2], &sth, &cth);
                double sth2 = sth*sth,
                       cth2 = cth*cth,
                       cth4=cth2*cth2,
                       sth4=sth2*sth2,
                       a2cthsth=a2_*cth*sth,
                       a2b2=a2_*b2_,
                       a4b2=a4_*b2_,
                       b4_=b2_*b2_,
                       a4b4=a4_*b4_,
                       a2b4=a2_*b4_;
                double a3=a2_*spin_,
                       a6=a4_*a2_;

    if (r>=1.)
            {

            double u=1./r, u2=u*u, u3=u2*u, u4=u3*u, u5=u4*u, u6=u5*u, u7=u6*u, u8=u7*u, u9=u8*u, u10=u9*u;

            dst[1][1][1]=-(4.*a2b2*b2_*u7*cth2-4.*a2b2*b2_*u7-8.*a2b2*u5*cth2+4.*a2b2*u4*cth2-4.*a2b2*u4-a2_*u2*cth2-4.*b2_*u3+a2_*u*cth2-a2_*u+1.)*u2/(2.*a2b2*u5+2.*b2_*u3+a2_*u2-2.*u+1.)/(a2_*u2*cth2+1.)/(2.*b2_*u3+1.);
            dst[1][2][1]=dst[1][1][2]=-a2cthsth*u2/(a2_*u2*cth2+1.);
            dst[1][2][2]=-(2.*a2b2*u5+2.*b2_*u3+a2_*u2-2.*u+1)/(a2_*u2*cth2+1.)/(2.*b2_*u3+1.)/u;
            dst[1][3][3]=-(4.*a4b4*u10*cth4+8.*a4b2*u8*cth2*sth2+4.*a4b2*u7*cth4+8.*a2b4*u8*cth2+a4_*u5*cth2*sth2+a4_*u4*cth4+4.*a2b2*u6*sth2+8.*a2b2*u5*cth2+4.*b4_*u6-a2_*u3*sth2+2.*a2_*u2*cth2+4.*b2_*u3+1.)*(2.*a2b2*u5+2.*b2_*u3+a2_*u2-2.*u+1)*sth2/(a2_*u2*cth2+1.)/(a2_*u2*cth2+1.)/(a2_*u2*cth2+1.)/(2.*b2_*u3+1.)/(2.*b2_*u3+1.)/(2.*b2_*u3+1.)/u;
            dst[1][3][0]=dst[1][0][3]=(8.*a2b2*u5*cth2+a2_*u2*cth2+4.*b2_*u3-1.)*(2.*a2b2*u5+2.*b2_*u3+a2_*u2-2.*u+1.)*spin_*u2*sth2/(a2_*u2*cth2+1.)/(a2_*u2*cth2+1.)/(a2_*u2*cth2+1.)/(2.*b2_*u3+1.)/(2.*b2_*u3+1.)/(2.*b2_*u3+1.);
            dst[1][0][0]=(8.*a2b2*u5*sth2-8.*a2b2*u5+a2_*u2*sth2-4.*b2_*u3-a2_*u2+1.)*(2.*a2b2*u5+2.*b2_*u3+a2_*u2-2.*u+1.)*u2/(a2_*u2*cth2+1.)/(a2_*u2*cth2+1.)/(a2_*u2*cth2+1.)/(2.*b2_*u3+1.)/(2.*b2_*u3+1.)/(2.*b2_*u3+1.);
            dst[2][1][1]=(2.*b2_*u3+1.)*a2cthsth*u4/(2.*a2b2*u5+2.*b2_*u3+a2_*u2-2.*u+1.)/(a2_*u2*cth2+1.);
            dst[2][2][1]=dst[2][1][2]=u/(a2_*u2*cth2+1.);
            dst[2][2][2]=-a2cthsth*u2/(a2_*u2*cth2+1.);
            dst[2][3][3]=-(2.*a4_*a2b2*u9*cth4+2.*a4b2*u7*cth4+a6*u6*cth4+4.*a4b2*u7*cth2-2.*a4_*u5*cth4+a4_*u4*cth4+4.*a2b2*u5*cth2+2.*a4_*u4*cth2+2.*a4_*u5+2.*a2b2*u5-4.*a2_*u3*cth2+2.*a2_*u2*cth2+4.*a2_*u3+2.*b2_*u3+a2_*u2+1.)*cth*sth/(a2_*u2*cth2+1.)/(a2_*u2*cth2+1.)/(a2_*u2*cth2+1.)/(2.*b2_*u3+1.);
            dst[2][0][3]=dst[2][3][0]=2.*(a2_*u2+1.)*spin_*u3*cth*sth/(a2_*u2*cth2+1.)/(a2_*u2*cth2+1.)/(a2_*u2*cth2+1.)/(2.*b2_*u3+1.);
            dst[2][0][0]=-2.*a2_*u5*cth*sth/(a2_*u2*cth2+1.)/(a2_*u2*cth2+1.)/(a2_*u2*cth2+1.)/(2.*b2_*u3+1.);
            dst[3][3][1]=dst[3][1][3]=(4.*a4b4*u10*cth4+8.*a4b2*u8*cth2*sth2+4.*a4b2*u7*cth4+8.*a2b4*u8*cth2+a4_*u5*cth2*sth2-4.*a2b2*u6*cth2+a4_*u4*cth4+4.*a2b2*u6*sth2+8.*a2b2*u5*cth2+4.*b4_*u6-2.*a2_*u3*cth2-a2_*u3*sth2-4.*b2_*u4+2.*a2_*u2*cth2+4.*b2_*u3-2.*u+1.)*u/(2.*a2b2*u5+2.*b2_*u3+a2_*u2-2.*u+1.)/(a2_*u2*cth2+1.)/(a2_*u2*cth2+1.)/(2.*b2_*u3+1.);
            dst[3][3][2]=dst[3][2][3]=(2.*a4b2*u7*cth4+a4_*u4*cth4+4.*a2b2*u5*cth2-2.*a2_*u3*cth2+2.*a2_*u2*cth2+2.*a2_*u3+2.*b2_*u3+1.)*cth/(a2_*u2*cth2+1.)/(a2_*u2*cth2+1.)/(2.*b2_*u3+1.)/sth;
            dst[3][0][1]=dst[3][1][0]=-(8.*a2b2*u5*cth2+a2_*u2*cth2+4.*b2_*u3-1.)*spin_*u4/(2.*a2b2*u5+2.*b2_*u3+a2_*u2-2.*u+1.)/(a2_*u2*cth2+1.)/(a2_*u2*cth2+1.)/(2.*b2_*u3+1.);
            dst[3][0][2]=dst[3][2][0]=-2.*spin_*u3*cth/(a2_*u2*cth2+1.)/(a2_*u2*cth2+1.)/(2.*b2_*u3+1.)/sth;
            dst[0][3][1]=dst[0][1][3]=(8.*a4b2*u7*cth2+4.*a2b2*u5*cth2+a4_*u4*cth2+4.*a2b2*u5-a2_*u2*cth2-a2_*u2-3.)*spin_*u2*sth2/(2.*a2b2*u5+2.*b2_*u3+a2_*u2-2.*u+1.)/(a2_*u2*cth2+1.)/(a2_*u2*cth2+1.)/(2.*b2_*u3+1.);
            dst[0][3][2]=dst[0][2][3]=-2.*(a2_*u2*sth2-a2_*u2-1.)*a3_*u3*cth*sth*sth2/(a2_*u2*cth2+1.)/(a2_*u2*cth2+1.)/(a2_*u2*cth2+1.)/(2.*b2_*u3+1.);
            dst[0][0][1]=dst[0][1][0]=(8.*a2b2*u5*sth2-8.*a2b2*u5+a2_*u2*sth2-4.*b2_*u3-a2_*u2+1.)*(a2_*u2+1.)*u2/(2.*a2b2*u5+2.*b2_*u3+a2_*u2-2.*u+1.)/(a2_*u2*cth2+1.)/(a2_*u2*cth2+1.)/(2.*b2_*u3+1.);
            dst[0][0][2]=dst[0][2][0]=-2.*a2_*u3*cth*sth/(a2_*u2*cth2+1.)/(a2_*u2*cth2+1.)/(2.*b2_*u3+1.);


                //cout << "t=" << pos[0] << " r=" << pos[1] << " theta=" << pos[2] << " dt=" << pos[4] << " dr=" << pos[5] << " dth=" << pos[6] << " dph=" << pos[7] << endl ;
            }

    if (r>=0. && r<1.)
            {
            double r2=r*r, r3=r2*r, r4=r3*r, r5=r4*r, r6=r5*r, r7=r6*r, r8=r7*r, r9=r8*r, r10=r9*r;
            double d=r3+2.*b2_;
            double Sigma=r2+a2_*cth2, Sigma2=Sigma*Sigma, Sigma3=Sigma2*Sigma;

            dst[1][1][1]=(a2_*r6*sth2+4.*a2b2*r3*sth2-a2_*r5*sth2+4.*a2b2*b2_*sth2-8.*a2b2*r2*sth2+a2_*r5-r7+8.*a2b2*r2+4.*b2_*r4)*r/(a2_*r3+r5+2.*a2b2+2.*b2_*r2-2.*r4)/Sigma/d;
            dst[1][2][1]=dst[1][1][2]=-a2cthsth/Sigma;
            dst[1][2][2]=-(a2_*r3+r5+2.*a2b2+2.*b2_*r2-2.*r4)*r/Sigma/d;
            dst[1][3][3]=-(a4_*r6*cth4+4.*a4_*b2_*r3*cth4+a4_*r5*cth2*sth2+2.*a2_*r8*cth2+4.*a4b4*cth4+8.*a4b2*r2*cth2*sth2+8.*a2b2*r5*cth2-a2_*r7*sth2+r10+8.*a2b4*r2*cth2+4.*a2b2*r4*sth2+4.*b2_*r7+4.*b4_*r4)*(a2_*r3+r5+2.*a2b2+2.*b2_*r2-2.*r4)*r*sth2/Sigma3/d/d/d;
            dst[1][3][0]=dst[1][0][3]=(a2_*r3*cth2+8.*a2b2*cth2-r5+4.*b2_*r2)*(a2_*r3+r5+2.*a2b2+2.*b2_*r2-2.*r4)*spin_*r3*sth2/Sigma3/d/d/d;
            dst[1][0][0]=-(a2_*r3*cth2+8.*a2b2*cth2-r5+4.*b2_*r2)*(a2_*r3+r5+2.*a2b2+2.*b2_*r2-2.*r4)*r3/Sigma3/d/d/d;
            dst[2][1][1]=d*a2cthsth/(a2_*r3+r5+2.*a2b2+2.*b2_*r2-2.*r4)/Sigma;
            dst[2][2][1]=dst[2][1][2]=r/Sigma;
            dst[2][2][2]=-a2cthsth/Sigma;
            dst[2][3][3]=-(a4_*a2_*r3*cth4+a4_*r5*cth4+2.*a4_*a2b2*cth4+2.*a4b2*r2*cth4-2*a4_*r4*cth4+2.*a4_*r5*cth2+2.*a2_*r7*cth2+4.*a2b2*a2_*r2*cth2+4.*a2b2*r4*cth2-4.*a2_*r6*cth2+a2_*r7+r9+2.*a4_*r4+2.*a2b2*r4+4.*a2_*r6+2.*b2_*r6)*sth*cth/Sigma3/d;
            dst[2][0][3]=dst[2][3][0]=2.*(a2_+r2)*spin_*r4*cth*sth/Sigma3/d;
            dst[2][0][0]=-2.*a2_*r4*cth*sth/Sigma3/d;
            dst[3][3][1]=dst[3][1][3]=(a4_*r6*cth4+4.*a4_*b2_*r3*cth4+a4_*r5*cth2*sth2+2.*a2_*r8*cth2+4.*a4b4*cth4+8.*a4_*b2_*r2*cth2*sth2+8.*a2b2*r5*cth2-2.*a2_*r7*cth2-a2_*r7*sth2+r10+8.*a2b4*r2*cth2-4.*a2b2*r4*cth2+4.*a2b2*r4*sth2+4.*b2_*r7-2.*r9+4.*b4_*r4-4.*b2_*r6)*r/Sigma2/(a2_*r3+r5+2.*a2b2+2.*b2_*r2-2.*r4)/d;
            dst[3][3][2]=dst[3][2][3]=(a4_*r3*cth4+2.*a4_*b2_*cth4+2.*a2_*r5*cth2+4.*a2b2*r2*cth2-2.*a2_*r4*cth2+r7+2.*a2_*r4+2.*b2_*r4)*cth/Sigma2/d/sth;
            dst[3][0][1]=dst[3][1][0]=-(a2_*r3*cth2+8.*a2b2*cth2-r5+4.*b2_*r2)*spin_*r3/Sigma2/(a2_*r3+r5+2.*a2b2+2.*b2_*r2-2.*r4)/d;
            dst[3][0][2]=dst[3][2][0]=-2.*spin_*r4*cth/Sigma2/d/sth;
            dst[0][3][1]=dst[0][1][3]=(a4_*r3*cth2-a2_*r5*cth2+8.*a4b2*cth2+4.*a2b2*r2*cth2-a2_*r5-3.*r7+4.*a2b2*r2)*spin_*r3*sth2/Sigma2/(a2_*r3+r5+2.*a2b2+2.*b2_*r2-2.*r4)/d;
            dst[0][3][2]=dst[0][2][3]=-2.*(a2_*sth2-a2_-r2)*a3_*r4*cth*sth2*sth/Sigma3/d;
            dst[0][0][1]=dst[0][1][0]=(a2_*r3*sth2+8.*a2b2*sth2-a2_*r3+r5-8.*a2b2-4.*b2_*r2)*(a2_+r2)*r3/(a2_*r3+r5+2.*a2b2+2.*b2_*r2-2.*r4)/Sigma2/d;
            dst[0][0][2]=dst[0][2][0]=-2.*a2_*r4*cth*sth/Sigma2/d;

                // cout.precision(16);
                // cout << "t=" << pos[0] << " r=" << pos[1] << " theta=" << pos[2] << " dt=" << pos[4] << " dr=" << pos[5] << " dth=" << pos[6] << " dph=" << pos[7] << endl ;
                
            }

    if (r<0.)
            {
            double r2=r*r, r3=r2*r, r4=r3*r, r5=r4*r, r6=r5*r, r7=r6*r, r8=r7*r, r9=r8*r, r10=r9*r;
            double sth3=sth2*sth;
            double Sigma=r2+a2_*cth2, Sigma2=Sigma*Sigma, Sigma3=Sigma2*Sigma;
            double dm=r3-2.*b2_, dm3=dm*dm*dm;


            dst[1][1][1]=(a2_*r6*sth2 - 4.*a2b2*r3*sth2 - a2_*r5*sth2 + 4.*a2b4*sth2 + 8.*a2b2*r2*sth2 + a2_*r5 - r7 - 8.*a2b2*r2 - 4.*b2_*r4)*r/(a2_*r3 + r5 - 2.*a2b2 - 2.*b2_*r2 - 2.*r4)/Sigma/dm;
            dst[1][2][1]=dst[1][1][2]=-a2cthsth/Sigma;
            dst[1][2][2]=-(a2_*r3 + r5 - 2.*a2b2 - 2.*b2_*r2 - 2.*r4)*r/Sigma/dm;
            dst[1][3][3]=-(a4_*r6*cth4 - 4.*a4_*b2_*r3*cth4 + a4_*r5*cth2*sth2 + 2.*a2_*r8*cth2 + 4.*a4_*b4_*cth4 - 8.*a4_*b2_*r2*cth2*sth2 - 8.*a2b2*r5*cth2 - a2_*r7*sth2 + r10 + 8.*a2b4*r2*cth2 - 4.*a2b2*r4*sth2 - 4.*b2_*r7 + 4.*b4_*r4)*(a2_*r3 + r5 - 2.*a2b2 - 2.*b2_*r2 - 2.*r4)*r*sth2/Sigma3/dm3;
            dst[1][3][0]=dst[1][0][3]=(a2_*r3*cth2 - 8.*a2b2*cth2 - r5 - 4.*b2_*r2)*(a2_*r3 + r5 - 2.*a2b2 - 2.*b2_*r2 - 2.*r4)*spin_*r3*sth2/Sigma3/dm3;
            dst[1][0][0]=-(a2_*r3*cth2 - 8.*a2b2*cth2 - r5 - 4.*b2_*r2)*(a2_*r3 + r5 - 2.*a2b2 - 2.*b2_*r2 - 2.*r4)*r3/Sigma3/dm3;
            dst[2][1][1]=dm*a2cthsth/(a2_*r3 + r5 - 2.*a2b2 - 2.*b2_*r2 - 2.*r4)/Sigma;
            dst[2][2][1]=dst[2][1][2]=r/Sigma;
            dst[2][2][2]=-a2cthsth/Sigma;
            dst[2][3][3]=-(a6*r3*cth4 + a4_*r5*cth4 - 2*a6*b2_*cth4 - 2*a4_*b2_*r2*cth4 - 2.*a4_*r4*cth4 + 2.*a4_*r5*cth2 + 2.*a2_*r7*cth2 - 4.*a4_*b2_*r2*cth2 - 4.*a2b2*r4*cth2 - 4.*a2_*r6*cth2 + a2_*r7 + r9 + 2.*a4_*r4 - 2.*a2b2*r4 + 4.*a2_*r6 - 2.*b2_*r6)*cth*sth/Sigma3/dm;
            dst[2][0][3]=dst[2][3][0]=2*(a2_ + r2)*spin_*r4*cth*sth/Sigma3/dm;
            dst[2][0][0]=-2*a2_*r4*cth*sth/Sigma3/dm;
            dst[3][3][1]=dst[3][1][3]=(a4_*r6*cth4 - 4.*a4_*b2_*r3*cth4 + a4_*r5*cth2*sth2 + 2.*a2_*r8*cth2 + 4.*a4_*b4_*cth4 - 8.*a4_*b2_*r2*cth2*sth2 - 8.*a2b2*r5*cth2 - 2.*a2_*r7*cth2 - a2_*r7*sth2 + r10 + 8.*a2b4*r2*cth2 + 4.*a2b2*r4*cth2 - 4.*a2b2*r4*sth2 - 4.*b2_*r7 - 2.*r9 + 4.*b4_*r4 + 4.*b2_*r6)*r/(a2_*r3 + r5 - 2.*a2b2 - 2.*b2_*r2 - 2.*r4)/Sigma2/dm;
            dst[3][3][2]=dst[3][2][3]=(a4_*r3*cth4 - 2.*a4_*b2_*cth4 + 2.*a2_*r5*cth2 - 4.*a2b2* r2*cth2 - 2.*a2_*r4*cth2 + r7 + 2.*a2_*r4 - 2.*b2_*r4)*cth/Sigma2/dm/sth;
            dst[3][0][1]=dst[3][1][0]=-(a2_*r3*cth2 - 8.*a2b2*cth2 - r5 - 4.*b2_*r2)*spin_*r3/(a2_*r3 + r5 - 2.*a2b2 - 2.*b2_*r2 - 2.*r4)/Sigma2/dm;
            dst[3][0][2]=dst[3][2][0]=-2.*spin_*r4*cth/Sigma2/dm/sth;
            dst[0][3][1]=dst[0][1][3]=(a4_*r3*cth2 - a2_*r5*cth2 - 8.*a4_*b2_*cth2 - 4.*a2b2*r2*cth2 - a2_*r5 - 3*r7 - 4.*a2b2*r2)*spin_*r3*sth2/(a2_*r3 + r5 - 2.*a2b2 - 2.*b2_*r2 - 2.*r4)/Sigma2/dm;
            dst[0][3][2]=dst[0][2][3]=-2.*(a2_*sth2 - a2_ - r2)*a3*r4*cth*sth3/Sigma3/dm;
            dst[0][0][1]=dst[0][1][0]=-(a2_*r3*cth2 - 8.*a2b2*cth2 - r5 - 4.*b2_*r2)*(a2_ + r2)*r3/(a2_*r3 + r5 - 2.*a2b2 - 2.*b2_*r2 - 2.*r4)/Sigma2/dm;
            dst[0][0][2]=dst[0][2][0]=-2.*a2_*r4*cth*sth/Sigma2/dm;

                //cout.precision(16);
                //cout << "t=" << pos[0] << " r=" << pos[1] << " theta=" << pos[2] << " dt=" << pos[4] << " dr=" << pos[5] << " dth=" << pos[6] << " dph=" << pos[7] << endl ;
            }

    return 0;
}



// Optimized version
double Hayward::ScalarProd(const double* pos,
                           const double* u1, const double* u2) const {
    double g[4][4];
    gmunu(g, pos);
    double res=
    g[0][0]*u1[0]*u2[0]+
    g[1][1]*u1[1]*u2[1]+
    g[2][2]*u1[2]*u2[2]+
    g[3][3]*u1[3]*u2[3]+
    g[0][3]*u1[0]*u2[3]+
    g[3][0]*u1[3]*u2[0];

# if GYOTO_DEBUG_ENABLED
    GYOTO_IF_DEBUG
    GYOTO_DEBUG_ARRAY(pos, 4);
    GYOTO_DEBUG_ARRAY(u1, 4);
    GYOTO_DEBUG_ARRAY(u2, 4);
    GYOTO_DEBUG   << "ScalarProd(pos, u1, u2)="
          <<  res
          << endl;
    GYOTO_ENDIF_DEBUG
# endif

    return res;

}

/*For integration of Hayward geodesics.
 diff is such that : y_dot=diff(y,cst) where cst are constants of motion (mu,E,L,Q in Hayward)
 and y contains [r,theta,phi,t,pr,ptheta] (pr and ptheta are canonical momentum)
 and y_dot is [rdot,thetadot,phidot,tdot,prdot,pthetadot]
 */
int Hayward::diff(const double* coordGen, const double* cst, double* res) const{
    double a=spin_;

    //int width=25;//15;
    //int prec=15;//8;

    double r = coordGen[1] ;

    if (r < 0.) {
        cerr << "r= " << r << endl;
        throwError( "Hayward.C : r negative!!!!! the horizon may have been crossed..." );

    }

    if (r < rsink_) {
#   if GYOTO_DEBUG_ENABLED
        GYOTO_DEBUG << "Too close to horizon in Hayward::diff at r= " << r << endl;
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

        double Sigma=r2+a2_*costheta2;
        if (Sigma==0) throwError("In Hayward::diff(): Sigma==0");
            double Sigmam1=1./Sigma;
            double Sigmam2=Sigmam1*Sigmam1;

            double Delta=r2-2*r+a2_;

            double E=cst[1];
            double E2=E*E;
            double L=cst[2];
            double L2=L*L;

            double tmp1=(2.*Delta*Sigma);
            if (tmp1==0)  throwError("In Hayward::diff(): 2.*Delta*Sigma==0");
                double tmp1m1=1./tmp1;

                if (Delta==0) throwError("In Hayward::diff(): Delta==0");

                    //NB: equations of motion are independent of Carter constant in this
                    //form. However, the dependency of the dynamic on this constant
                    //appears when transforming from principal momenta to coordinate
                    //derivatives (e.g. p_theta -> thetadot)


                /*
                 ---> Standard Kerr equations of geodesics
                 */
                    res[0] = tmp1m1*(2.*(r*(-2.*a*L+E*r3+a2_*E*(2.+r))+a2_*E*(a2_+r*(-2.+r))
                                         *costheta2));// tdot

                    res[1] = Delta*Sigmam1*pr; //rdot

                    res[2] = Sigmam1*ptheta; //thetadot

                    res[3] = -tmp1m1*(-2.*(r*(2.*a*E+L*(-2.+r))+L*(a2_+r*(-2.+r))
                                           *cotantheta2)); //phidot

                    res[4] = 0.;// ptdot: pt = cst = -E

                    double tmp2=r2+a2_*costheta2;
                    if (tmp2==0) throwError("r2+a2*costheta2==0");
                        double tmp2m2=1./(tmp2*tmp2);

                        double tmp3=a2_+r*(-2.+r);
                        double tmp3_2=tmp3*tmp3;

                        res[5] =
                        -0.5*(2.*(r*(r-a2_)-a2_*(1.-r)*costheta2)*tmp2m2)*pr*pr
                        -0.5*(-2.*r*tmp2m2)*ptheta*ptheta
                        +(tmp2m2/tmp3_2
                          *(a2_*(a4_*E2-2.*a3_*E*L+2.*a*E*L*r2+E2*r3*(-4.+r)
                                 +a2_*(L2*(1.-r)+2*E2*r2))*costheta2
                            +r*(-r*(a4_*E2-2.*a3_*E*L+2.*a*E*L*(4.-3.*r)*r
                                    +a2_*(L2+2.*E2*r*(-2.+r))+r*(E2*r3-L2*(-2.+r)*(-2.+r)))
                                +L2*tmp3_2*cotantheta2)));// prdot

                        res[6]=
                        -0.5*(a2_*Delta*sin2theta*Sigmam2)*pr*pr
                        -0.5*(a2_*sin2theta*Sigmam2)*ptheta*ptheta
                        +(
                          Sigmam2
                          *(
                            L2*r2*cotantheta
                            +0.5*L2*(a2_+2.*r2+a2_*cos2theta)*cotantheta3
                            +a2_*r*(2.*a2_*E2-4.*a*E*L+L2*(2.-r)+2.*E2*r2)*costheta*sintheta/Delta
                            )
                          ); // pthetadot

                        res[7] = 0.;//pphidot: pphi = cst = L

                        return 0;
}

void Hayward::circularVelocity(double const coor[4], double vel[4],
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
    double sinth = sin(coor[2]);
    double coord[4] = {coor[0], coor[1]*sinth, M_PI*0.5, coor[3]};

    vel[1] = vel[2] = 0.;
    //vel[3] = 1./((dir*pow(coord[1], 1.5) + spin_)*sinth);
    vel[3] = 1./((dir*pow(coord[1], 1.5) + spin_));

    vel[0] = SysPrimeToTdot(coor, vel+1);
    vel[3] *= vel[0];
# if GYOTO_DEBUG_ENABLED
    GYOTO_DEBUG_ARRAY(vel,4);
# endif
}

//Runge Kutta to order 4
int Hayward::myrk4(Worldline * line, const double coordin[8],
                   double h, double res[8]) const
{
    if (generic_integrator_) return Generic::myrk4(line, coordin, h, res);

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

        //  line -> checkPhiTheta(res);

        return 0;
}

int Hayward::myrk4(const double coor[8], const double cst[5],
                   double h, double res[8]) const{
    //  cout << "In Hayward::myrk4 - 2" << endl;
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
    double derlim_hor=1e4, derlim_gen=1e6, derlim;
    // throw a "z-axis problem" (return 1) if prdot or pthdot becomes higher
    double aa=spin_;
    double rhor=1.+sqrt(1.-aa*aa), factrtol=5.;
    double thetatol_hor=1e-1, thetatol_gen=2e-3, thetatol;
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
     played with to smooth the integration. In particular going to
     higher spin may impose using higher thetatol_gen
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

int Hayward::myrk4_adaptive(Worldline * line, const double coordin[8],
                            double lastnorm, double normref, double coordout1[8],
                            double h0, double& h1, double h1max) const
{
    if (generic_integrator_)
        return Generic::myrk4_adaptive(line, coordin, lastnorm, normref,
                                       coordout1, h0, h1, h1max);

    /*Switch BL -> principal momenta*/

        double coor[8], coor1[8], cstest[5], coorhalf[8], coor2[8],
        coor1bis[8], mycoor[8], delta1[8];
    double const * const cst = line -> getCst();
    MakeMomentum(coordin,cst,coor);
    double delta0[8], dcoor[8];
    double delta0min=1e-15, eps=0.0001, S=0.9, errmin=1e-6, hbis=0.5*h0,
    err, diffr, diffth, normtemp,
    cstol_gen=1e-3, cstol_hor=1e-2, cstol, div, QCarter;
    int countbis=0, countbislim=50, zaxis=0; // for z-axis problem in myrk4
    //int norm1=0, normhalf=0, norm2=0, rk1=0, rkhalf=0, rk2=0, update, makerr=0.;
    int norm1=0, rk1=0, rkhalf=0, rk2=0, update, makerr=0;
    double a=spin_, factrtol=3.,
    rtol=factrtol*(1.+sqrt(1.-a*a)), rlimitol=10.;

    h1max=deltaMax(coor, h1max);

    if (coor[1] < rtol) cstol = cstol_hor;
        // for tests of cst of motion conservation; don't ask too much
        // if near horizon...
        else cstol = cstol_gen;

            //NB: following test ok only if theta is 0-pi
            //    coz -pi/2[pi]=-pi/2 ...
        /*if (fabs(fmod(coor[2]+M_PI/2, M_PI)-M_PI/2) < GYOTO_MIN_THETA) {
         cout << "in kerr th= " << coor[2] << endl;
         GYOTO_WARNING << "Too close to Z axis: stopping integration"<< endl;
         return 1;
         }*/
            //Fred: Test removed Feb 2013, seems not necessary

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

            double hinit=h0, steptol=1e3;
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

                if (h0/hinit>steptol){
                    GYOTO_SEVERE << "In Hayward::myrk4_adaptive: "
                    "stop integration because unsolved z-axis problem" << endl;
                    return 1;
                }


            }


            if (countbis > countbislim && zaxis) {
                // If the "jumping ahead" trick above does not work and leads to
                // infinite jumping ahead behavior (should never happen of
                // course)

                GYOTO_INFO << "WARNING: " << endl
                << "in Hayward.C couldn't solve z-axis problem ; stopping..."
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
                if (fabs(h1)<delta_min_) {
                    GYOTO_SEVERE << "DeltaMin is too small" << endl;
                    h1= (h1>0)?delta_min_:-delta_min_;
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

                if ((diffr > difftol_ || diffth > difftol_)) {
                    // don't update coordinates if the relative differences are
                    //  more than 1% --> consequence = norm and Carter cst
                    //  won't be exactly the same for this step --> below,
                    //  test to be sure they're not "too much violated"; NB:
                    //  if this test is not performed, the "corrected"
                    //  worldline can diverge from the "true" one after a long
                    //  time of integration.
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
                            cerr << "WARNING: at r=" << coordin[1] << endl;
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
                            << "Hayward::myrk4_adaptive" << endl;
#           endif
                            // some rare cases can end up with bad cst conservation
                            // even at r = a few rhor...
                        }else{
                            GYOTO_SEVERE << "This warning occurred at r= " << coor1[1] << endl
                            << "i.e. far from horizon --> to be investigated"
                            << ", or maybe increase parameter cstol"
                            << "in Hayward.C" << endl;
                        }
                    }

                }

                //Update coord
                if (update && !norm1){ // norm1=1 if impossible to normalize in CheckCons due to z-axis pb
                    for (int i=0;i<8;i++) coor1[i]=coor1bis[i];
                }
                //      cout << "Hayward Used h0= " << h0 << endl;
                //*** Switch principal momenta -> BL: ***

                MakeCoord(coor1,cst,coordout1);
                break;

            } //err>1 if-loop end


        } // while loop end

    //line -> checkPhiTheta(coordout1);
    //for some reason the above line leads to infinite iteration with fixed star

    return 0;
}

int Hayward::CheckCons(const double coor_init[8], const double cst[5], double coor_fin[8]) const {
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
     As the equations of motion (cf diff) are independent of Q,
     it is necessary to check whether this constant is conserved.
     */

    double costh, sinth;
    sincos(mycoor[2], &sinth, &costh);
    double sinthm2=1./(sinth*sinth), costh2=costh*costh;
    double mu=cst[0], EE=cst[1], LL=cst[2], QQ=cst[3], QQm1=cst[4];
    double argsqrt, limarg=1e-6*QQ, limargbis=0.1*QQ;
    double mu2=mu*mu, EE2=EE*EE, LL2=LL*LL;
    double Sigma=mycoor[1]*mycoor[1]+a2_*costh2;
    double Sigma2=Sigma*Sigma;
    double Qtest=
    Sigma2*mycoor[6]*mycoor[6]+costh2*(a2_*(mu2-EE2)+LL2*sinthm2);
    //this should be equal to constant QQ
    int thdotpos=1;
# if GYOTO_DEBUG_ENABLED
    GYOTO_DEBUG
    << "mu="<<mu<<", EE="<<EE<<", LL="<<LL<<", QQ="<<QQ<<", QQm1="<<QQm1
    << ", Qtest="<<Qtest<< ", fabs(Qtest-QQ)/QQm1="<< fabs(Qtest-QQ)*QQm1
    << endl;
# endif
    if (fabs(Qtest-QQ)*QQm1 > 1e-6){//Then change thetadot to allow Qtest=QQ
        // Note: if Q==0, QQm1==1.

        if (mycoor[6]<0.) thdotpos=0; //to preserve the sign of thetadot
        argsqrt=QQ-costh2*(a2_*(mu2-EE2)+LL2*sinthm2);//thetadot should be the sqrt of this quantity

        if (argsqrt<0 && fabs(argsqrt)<=limarg) {//if argsqrt is <0, but "small enough", put it to zero
            argsqrt=0.;
        }else if (argsqrt<0 && fabs(argsqrt)>limarg) {//if it is <0 and "big", two cases: either we are too close to z-axis which causes theta coordinate to behave badly (then return 1: "z-axis problem" in RK4 adaptive); or Error
            if (fabs(fmod(coor_init[2]+M_PI/2, M_PI)-M_PI/2) < M_PI/50.) {
                return 1; // z-axis problem
            }else{
                if (fabs(argsqrt)>limargbis)
                    throwError("In Hayward::CheckCons Impossible to determine thetadot; "
                               "maybe try to increase parameter limarg");
                GYOTO_INFO << "Hayward::CheckCons argsqrt= " << argsqrt << " at theta= " << coor_init[2] << ". Putting it to 0..." << endl;
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

void Hayward::Normalize4v(double coord[8], const double part_mass) const {
    /*
     Changes rdot to allow norm conservation.
     */

    double aa=spin_;
    double rhor=1.+sqrt(1.-aa*aa);

    double g[4][4];
    gmunu(g, coord);

    double gtt=g[0][0], gtph=g[0][3], grr=g[1][1], gthth=g[2][2], gphph=g[3][3];
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
                            << "in Hayward::Normalize4v impossible to normalize 0-mass "
                            << "particule next to horizon. Putting argrac to 0. "
                            << "Effective value of argrac= " << argrac << endl
                            << "with coord= ";
                            for (int myind=0;myind<8;myind++) cerr << coord[myind] << " ";
                            cerr << endl;
                        }
                        argrac=0.;
                    }else{//if the quantity is <0 and "big", and we are not close to horizon --> error
                        throwError( "In Hayward::Normalize4v impossible to normalize 0-mass particle outside horizon!" );
                    }
                }

                coord[5]=sqrt(argrac);//update rdot
                if (!valuepos) coord[5]*=-1.;

            }else if (part_mass>0.){ // *** -1 - NORM CASE ***

                double argrac=-(1.+gtt*tdot*tdot+2.*gtph*phdot*tdot+gthth*thdot*thdot+gphph*phdot*phdot)/grr;
                double arglim=1e-7;

                GYOTO_DEBUG_ARRAY(coord, 8);
                GYOTO_DEBUG_EXPR(argrac);
                GYOTO_DEBUG_EXPR(rr/rhor);

                if (argrac<0 && fabs(argrac)<arglim) argrac=0.;
                if (argrac<0 && fabs(argrac)>arglim) {
                    if (rr/rhor < 2.) {//A AFFINER??
                        if (verbose()>=GYOTO_WARNING_VERBOSITY) {
                            cerr << "WARNING -1 - NORM CLOSE TO HORIZON : "
                            << "in Hayward::Normalize4v impossible to normalize massive "
                            << "particle next to horizon. Putting argrac to 0. "
                            << "Effective value of argrac= " << argrac << endl
                            << "with coord= ";
                            for (int myind=0;myind<8;myind++) cerr << coord[myind] << " ";
                            cerr << endl;
                        }
                        argrac=0.;
                    }else{
                        throwError( "In Hayward::Normalize4v impossible to normalize massive particule outside horizon!" );
                    }
                }

                coord[5]=sqrt(argrac);//update rdot
                if (!valuepos) coord[5]*=-1.;

            }else{
                throwError("In Hayward::Normalize4v: negative mass!");
            }

}


void Hayward::MakeCoord(const double coordin[8], const double cst[5], double coord[8]) const {
    double tt=coordin[0], rr = coordin[1], theta=coordin[2], phi=coordin[3],
    pr=coordin[5], ptheta=coordin[6];

    double g[4][4];
    gmunu(g, coordin);

    double gtt=g[0][0], gtp=g[0][3], gpp =g[3][3],
    det=gtp*gtp-gtt*gpp, detm1=1./det,
    guprr=gmunu_up(coordin,1,1),
    gupthth=gmunu_up(coordin,2,2);

    double EE=cst[1], LL=cst[2];

    double rdot=guprr*pr,
    thetadot = gupthth*ptheta,
    phidot=-(gtt*LL+gtp*EE)*detm1,
    tdot=(gtp*LL+gpp*EE)*detm1;

    coord[0]=tt;coord[1]=rr;coord[2]=theta;coord[3]=phi;coord[4]=tdot;
    coord[5]=rdot;coord[6]=thetadot;coord[7]=phidot;

}

void Hayward::MakeMomentum(const double coord[8], const double cst[5], double coordout[8]) const{
    double EE=cst[1], LL=cst[2];

    double tt=coord[0], rr = coord[1], theta=coord[2], phi=coord[3],
    rdot=coord[5], thetadot=coord[6];

    double g[4][4];
    gmunu(g, coord);

    double pr=g[1][1]*rdot, ptheta=g[2][2]*thetadot;

    coordout[0]=tt;coordout[1]=rr;coordout[2]=theta;
    coordout[3]=phi;coordout[4]=-EE;coordout[5]=pr;
    coordout[6]=ptheta;coordout[7]=LL;
}

void Hayward::nullifyCoord(double coord[8]) const {
    double tdot2;
    nullifyCoord(coord,tdot2);
}

void Hayward::nullifyCoord(double coord[4], double & tdot2) const {

    int i;
    double a, b=0., c=0.;

    double g[4][4];
    gmunu(g, coord);

    a=g[0][0];
    b=g[0][3]*coord[7];

    for (i=1;i<=3;++i){
        c+=g[i][i]*coord[4+i]*coord[4+i];
    }

    double sDelta=sqrt(b*b-a*c), am1=1./a;
    tdot2=(-b+sDelta)*am1;
    coord[4]=(-b-sDelta)*am1;
}

void Hayward::computeCst(const double coord[8], double cst[5]) const{
    double norm=ScalarProd(coord, coord+4, coord+4);
    double mu;//Particule mass: 0 or 1
    if (fabs(norm)<fabs(norm+1.)){
        mu=0.;
    }else{
        mu=1.;
    }

    double g[4][4];
    gmunu(g, coord);

    // EE and LL are generic for any axisymmetric spacetime
    double tdot=coord[4], phidot=coord[7];
    double EE=-g[0][0]*tdot-g[0][3]*phidot, LL=g[3][3]*phidot+g[0][3]*tdot;

    // Carter constant is Kerr-specific:
    double sinth, costh;
    double theta=coord[2], thetadot=coord[6];
    sincos(theta, &sinth, &costh);
    double costheta2=costh*costh, sintheta2=sinth*sinth;
    double gthth = g[2][2];
    double QQ=gthth*thetadot*gthth*thetadot
    +costheta2*(a2_*(mu*mu-EE*EE)+LL*LL/sintheta2);

    cst[0]=mu;cst[1]=EE;cst[2]=LL;cst[3]=QQ;
    cst[4]= cst[3]==0. ? 1. : 1./cst[3]; // cache 1/Q or 1. if Q==0
}

void Hayward::setParticleProperties(Worldline * line, const double* coord) const
{
    double cst[5];
    computeCst(coord,cst);
    line -> setCst(cst,5);
}

void Hayward::observerTetrad(string const obskind,
                             double const pos[4], double fourvel[4],
                             double screen1[4], double screen2[4],
                             double screen3[4]) const{

    double gtt = gmunu(pos,0,0),
    grr      = gmunu(pos,1,1),
    gthth    = gmunu(pos,2,2),
    gpp      = gmunu(pos,3,3),
    gtp      = gmunu(pos,0,3);

    if (obskind=="ZAMO"){
        double fact = gtt*gpp - gtp*gtp;
        if (fact == 0.){
            throwError("In Hayward::observerTetrad: "
                       "bad values");
        }
        double ut2 = -gpp/fact;
        if (grr==0. || gthth==0. || gpp==0. || ut2<=0.){
            throwError("In Hayward::observerTetrad: "
                       "bad values");
        }
        double omega = -gtp/gpp;
        double ut = sqrt(ut2);
        double fourv[4]={ut,0.,0.,omega*ut};
        double e3[4] = {0.,-1./sqrt(grr),0.,0.};
        double e2[4] = {0.,0.,-1./sqrt(gthth),0.};
        double e1[4] = {0.,0.,0.,-1/sqrt(gpp)};

        for (int ii=0;ii<4;ii++){
            fourvel[ii]=fourv[ii];
            screen1[ii]=e1[ii];
            screen2[ii]=e2[ii];
            screen3[ii]=e3[ii];
        }

    }else if (obskind=="KeplerianObserver"){
        double omega = 1./(pow(pos[1],1.5)+spin_);
        double ut2 = -1/(gtt+gpp*omega*omega+2.*gtp*omega);
        if (ut2 <= 0. || grr<=0. || gthth <=0.) {
            throwError("In Hayward::observerTetrad: "
                       "bad values");
        }
        double ut = sqrt(ut2);
        double fourv[4]={ut,0.,0.,omega*ut};
        double e3[4] = {0.,-1./sqrt(grr),0.,0.};
        double e2[4] = {0.,0.,-1./sqrt(gthth),0.};

        double fact1 = (gtp+gpp*omega)/(gtt+gtp*omega),
        fact2 = gtt*fact1*fact1+gpp-2.*gtp*fact1;
        if (fact2 <= 0.) throwError("In Hayward::observerTetrad: "
                                    "bad values");
        double a2 = 1./sqrt(fact2), a1 = -a2*fact1;
        double e1[4] = {-a1,0.,0.,-a2};

        for (int ii=0;ii<4;ii++){
            fourvel[ii]=fourv[ii];
            screen1[ii]=e1[ii];
            screen2[ii]=e2[ii];
            screen3[ii]=e3[ii];
        }
    }else{
        throwError("In Hayward::observerTetrad "
                   "unknown observer kind");
    }
    Generic::observerTetrad(obskind,pos,fourvel,screen1,screen2,screen3);
}
