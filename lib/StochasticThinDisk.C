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

#include "GyotoPhoton.h"
#include "GyotoStochasticThinDisk.h"
#include "GyotoProperty.h"
#include "GyotoUtils.h"
#include "GyotoFactoryMessenger.h"
#include "GyotoKerrBL.h"
#include "GyotoKerrKS.h"
#include "GyotoRezzollaZhidenko.h"
#include "GyotoKonoplyaRezzollaZhidenko.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <cmath>
#include <limits>
#include <string>
#include <random>
#include <boost/math/special_functions/bessel.hpp>
#include <vector>

using namespace std;
using namespace Gyoto;
using namespace Gyoto::Astrobj;
using namespace Gyoto::Metric;

GYOTO_PROPERTY_START(StochasticThinDisk)
GYOTO_PROPERTY_VECTOR_DOUBLE(StochasticThinDisk, Model_param, model_param,
			     "Parameters useful for the disk, max number NPAR_MAX")
GYOTO_PROPERTY_STRING(StochasticThinDisk, EquationKind, equationKind,
		      "Propagation equation of a scalar field which must be analytically solved")
GYOTO_PROPERTY_STRING(StochasticThinDisk, MotionKind, motionKind,
		      "Keplerian or Radial or Mixed (subkeplerian) velocity of the disk")
GYOTO_PROPERTY_END(StochasticThinDisk, ThinDisk::properties)

#define ath 3
#define ar 3

#define Mmax 64
#define Nmax 64

#define NPAR_MAX 10 // Max allowed number of parameters

#define WAVE 0
#define ADVECTIONDIFFUSION 1

#define STATIC 0
#define CIRCULAR 1
#define RADIAL 2

void StochasticThinDisk::equationKind(string const &kind) {
  if (kind == "Wave")
    equationkind_ = WAVE;
  else if (kind == "AdvectionDiffusion")
    equationkind_ = ADVECTIONDIFFUSION;
  else
    throwError("StochasticThinDisk: Unknown velocity kind");
}
string StochasticThinDisk::equationKind() const {
  switch (equationkind_) {
  case WAVE:
    return "Wave";
  case ADVECTIONDIFFUSION:
    return "AdvectionDiffusion";
  default:
    throwError("StochasticThinDisk: Unknown equation kind tag");
  }
  return "will not reach here, this line to avoid compiler warning"; 
}

void StochasticThinDisk::motionKind(string const &kind) {
  if (kind == "Static")
    motionkind_ = STATIC;
  if (kind == "Circular")
    motionkind_ = CIRCULAR;
  else if (kind == "Radial")
    motionkind_ = RADIAL;
  else
    throwError("StochasticThinDisk: Unknown velocity kind");
}
string StochasticThinDisk::motionKind() const {
  switch (motionkind_) {
  case STATIC:
    return "Static";
  case CIRCULAR:
    return "Circular";
  case RADIAL:
    return "Radial";
  default:
    throwError("StochasticThinDisk: Unknown velocity kind tag");
  }
  return "will not reach here, this line to avoid compiler warning"; 
}

void StochasticThinDisk::model_param(std::vector<double> const &v) {
  size_t n = v.size();
  if (n>NPAR_MAX) throwError("Too many parameters in model_param");
  for (size_t i=0; i<n; ++i) model_param_[i]=v[i];
}
std::vector<double> StochasticThinDisk::model_param() const {
  std::vector<double> v(NPAR_MAX, 0.);
  for (size_t i=0; i<NPAR_MAX; ++i) v[i]=model_param_[i];
  return v;
}

StochasticThinDisk::StochasticThinDisk() :
  ThinDisk("StochasticThinDisk"),
  equationkind_(WAVE),
  motionkind_(STATIC),
  model_param_(NULL)
{
  GYOTO_DEBUG << endl;

  model_param_ = new double[NPAR_MAX];
  for (int ii = 0; ii < NPAR_MAX; ii++) model_param_[ii] = 0.;

  const int size = (2*Mmax+1)*Nmax;
  Cmn_.resize(size, 1.);
  Nmn_.resize(size, 1.);
  Phimn_.resize(size, 0.);
  Psimn_.resize(size, 0.);
  Kmn_.resize(size, 0.);
}

StochasticThinDisk::StochasticThinDisk(const StochasticThinDisk& o) :
  ThinDisk(o),
  equationkind_(o.equationkind_),
  motionkind_(o.motionkind_),
  model_param_(NULL),
  Cmn_(o.Cmn_),
  Nmn_(o.Nmn_),
  Phimn_(o.Phimn_),
  Psimn_(o.Psimn_),
  Kmn_(o.Kmn_)
{
  if (o.gg_()) gg_=o.gg_->clone();
  Generic::gg_=gg_;

  model_param_ = new double[NPAR_MAX];
  for (int ii=0;ii<NPAR_MAX;ii++) model_param_[ii]=o.model_param_[ii];
}
StochasticThinDisk* StochasticThinDisk::clone() const
{ return new StochasticThinDisk(*this); }

bool StochasticThinDisk::isThreadSafe() const {
  return ThinDisk::isThreadSafe();
}

StochasticThinDisk::~StochasticThinDisk() {
  GYOTO_DEBUG << endl;
  delete [] model_param_;
}

double StochasticThinDisk::spectrum(double const alpha_r, double const alpha_theta, 
                                int const m, int const n) const{
  
  double S, S_r, S_theta;
  
  S_theta = pow(1.+m*m,-alpha_theta/2.);
  S_r = pow(1.+n*n,-alpha_r/2.);
  S = S_theta*S_r;
  
  return S;
}

void StochasticThinDisk::modalQuantities(){
  /*std::default_random_engine generator_amp(42); 
  std::default_random_engine generator_phi(16);
  std::default_random_engine generator_psi(54);*/
  std::default_random_engine generator(42);
    
  std::normal_distribution<double> gaussian(0.,1.); 
  std::uniform_real_distribution<double> uniform(0.,1.);
    
  double rout = ThinDisk::outerRadius();
  Cmn_.resize((2*Mmax+1)*Nmax);
  Nmn_.resize((2*Mmax+1)*Nmax);
  Phimn_.resize((2*Mmax+1)*Nmax);
  Psimn_.resize((2*Mmax+1)*Nmax);
  Kmn_.resize((2*Mmax+1)*Nmax);
    
  std::vector<double> zeros(Nmax);
    
  for (int m = -Mmax; m <= Mmax; m++) {
    int m_abs = abs(m);
    boost::math::cyl_bessel_j_zero(double(m_abs), 1, Nmax, zeros.begin());
    for (int n = 1; n <= Nmax; n++) {
      int idx = (m+Mmax)*Nmax+(n-1);
      if (equationkind_==WAVE) {
        double Smn = spectrum(ar, ath, m, n);
        double sigma = std::sqrt(Smn);
        Cmn_[idx] = sigma*gaussian(generator);
        Nmn_[idx] = 1.;
        Phimn_[idx] = 2.*M_PI*uniform(generator);
        Psimn_[idx] = 2.*M_PI*uniform(generator);
        Kmn_[idx] = zeros[n-1]/rout;
      }
    }
  }
}

double StochasticThinDisk::solution(double const coord_obj[8]) const{
  
  double t = coord_obj[0], r = coord_obj[1], phi = coord_obj[3];
  double kmn, omegamn, Jm;
  double Amn, phimn, psimn;
  double u=0;
  
  for (int m = -Mmax; m <= Mmax; m++) {
    int m_abs = abs(m);
    for (int n = 1; n <= Nmax; n++) {
      int idx = (m+Mmax)*Nmax+(n-1);
      
      Amn = Cmn_[idx]/Nmn_[idx];
      phimn = Phimn_[idx];
      psimn = Psimn_[idx];
      kmn = Kmn_[idx];
      omegamn = kmn;
      Jm = std::cyl_bessel_j(m_abs, kmn*r);
      
      if (equationkind_==WAVE) {
      
        u += Amn*Jm*cos(m*phi+phimn)*cos(omegamn*t+psimn);
        
        }
    }
  }
  
  return u;
}

double StochasticThinDisk::envelope(double nu, state_t const &coord_ph, double const coord_obj[8]) const{
			    
  string emission_model = "Constant"; // should be in: "Constant", "Gralla_et_al", "Thermal_Synchrotron"  
  double rr = coord_obj[1];
  double env=0.;	
  
  if (emission_model == "Constant"){
    // Uniform emission model equal to 1
    env = 1.;
  }
  else if (emission_model == "Gralla_et_al"){
    // Emission from Gralla et al 2020
    // ****************************** //
    // Here model_param must contain:
    // model_param = [gamma, mu, sigG]
    
    // This model is Kerr specific
    string kin = gg_->kind();      
    if (kin != "KerrBL")
      GYOTO_ERROR("StochasticThinDisk: KerrBL needed!");
    double SPIN = static_cast<SmartPointer<Metric::KerrBL> >(gg_) -> spin(),
      a2 = SPIN*SPIN;
      
    double rhor=1.+sqrt(1.-a2), rminus=1.-sqrt(1.-a2);
    
    // Choose profile here:
    double gamma=model_param_[0],
      mu=model_param_[1],
      sigG=model_param_[2];
    double norm =1.e-5; // adjust the normalisation to get a reasonable flux for your source
    // e.g n=1.e-5 for M87
    
    double tmp = gamma+asinh((rr-mu)/sigG);
    env = norm*exp(-0.5*tmp*tmp)/sqrt((rr-mu)*(rr-mu)+sigG*sigG);
  }
  
  if (emission_model == "Thermal_Synchrotron"){
    // Emission from Vincent+22 thick disk paper synchrotron formula
    // Equation B.7, Appendix B (nu_em=cst=230GHz)
    // ****************************** //
    // Here model_param must contain:
    //model_param = [zeta]
    
    double rin = ThinDisk::innerRadius();
    
    // Choose profile here:
     double zeta = model_param_[0];
     double norm=1.e-3;  // adjust the normalisation to get a reasonable flux for your source
     // the 1e-3 for zeta=3 is just there to get a reasonable flux for M87
   
     env = norm*exp(-zeta*rr/rin);
  }
  	    
  return env;
}

double StochasticThinDisk::emission(double nu, state_t const &coord_ph, double const coord_obj[8]) const{

  double emiss=0.; 
  emiss = envelope(nu, coord_ph, coord_obj)*solution(coord_obj);
  
  return emiss;
}

void StochasticThinDisk::getVelocity(double const pos[4], double vel[4])
{  
    string kin = gg_->kind();
    double gtt = gg_->gmunu(pos,0,0),
           grr = gg_->gmunu(pos,1,1),
           gpp = gg_->gmunu(pos,3,3),
           gtp = gg_->gmunu(pos,0,3),
           guptt = gg_->gmunu_up(pos,0,0),
           guptp = gg_->gmunu_up(pos,0,3),
           guppp = gg_->gmunu_up(pos,3,3),
           guprr = gg_->gmunu_up(pos,1,1);
    double rr = pos[1];
    
    if (motionkind_==STATIC) {
      // u^\mu = (u^t,0,0,0)
      vel[1] = 0.;
      vel[2] = 0.;
      vel[3] = 0.;
      vel[0] = sqrt(-1./gtt);
    } else if (motionkind_==RADIAL) {
      // u_\mu = (u_t,u_r,0,0)
      // with u_t=-1 (null radial velocity at infinity) and u_r obtained by unit normalisation
      vel[0] = -guptt;
      vel[1] = -sqrt((-1.-guptt)/grr);
      vel[2] = 0.;
      vel[3] = -guptp;
    } else if (motionkind_==CIRCULAR) {
      // u_\mu = (u_t,u_r,0,u_phi) equatorial
      // u_t=-E, u_phi=L
      double risco = 0.;
      if (kin!="Minkowski" && kin!="Hayward")
      risco=gg_->getRms(); // prograde Kerr ISCO
      //cout << "In StochasticThinDisk::getVelocity: r, isco= " << pos[1] << ", " << risco << endl;
      if (rr > risco){
        // KERPLERIAN MOTION
        gg_->circularVelocity(pos, vel, 1);
      } else {
        // RADIAL INFALL
        vel[0] = -guptt;
        vel[1] = -sqrt((-1.-guptt)/grr);
        vel[2] = 0.;
        vel[3] = -guptp;
      }
    }
    
  // Check the normalisation of the 4-velocity
  double tol=1e-5;
  double u2 = gg_->ScalarProd(pos,vel,vel);
  if (fabs(u2+1.)>tol or u2!=u2) { 
    cerr << " *** 4-velocity squared norm = " << u2 << endl;
    throwError("In StochasticThinDisk: 4vel is not properly normalized!");
  }
  //cout << "4vel = " << vel[0] << " " << vel[1] << " " << vel[2] << " " << vel[3]<< endl;
}

void StochasticThinDisk::processHitQuantities(Photon* ph,
      					   state_t const &coord_ph_hit,
 					   double const *coord_obj_hit,
 					   double dt,
 					   Properties* data) const {
 #if GYOTO_DEBUG_ENABLED
   GYOTO_DEBUG << endl;
 #endif
   SmartPointer<Spectrometer::Generic> spr = ph -> spectrometer();
   size_t nbnuobs = spr() ? spr -> nSamples() : 0 ;
   //cout << "nbnuobs " << nbnuobs << endl;
   if (nbnuobs!=1) GYOTO_ERROR("nbnuobs should be 1"); //spectro.set("NSamples", 1)
   double const * const nuobs = nbnuobs ? spr -> getMidpoints() : NULL;
   double dlambda = dt/coord_ph_hit[4]; //dlambda = dt/tdot
   double ggredm1 = -gg_->ScalarProd(&coord_ph_hit[0],coord_obj_hit+4,
 				    &coord_ph_hit[4]);// / 1.; 
                                        //this is nu_em/nu_obs
   if (noredshift_) ggredm1=1.;
   double ggred = 1./ggredm1;           //this is nu_obs/nu_em
   double dsem = dlambda*ggredm1; // *1.
   double inc =0.;

   if (data){ // this check is necessary as process can be called
     // with data replaced by NULL (see ComplexAstrobj and
     // Photon nb_cross_eqplane)
     if (data->user4) {
       // *** CAREFUL!! ***
       /*
 	I have to include here the "fudge factor" of Gralla+.
 	Do not forget to remove it to consider some other disk profile.
       */
       double max_cross_eqplane = ph->maxCrossEqplane();
       if (max_cross_eqplane==DBL_MAX) cout << "WARNING: in StochasticThinDisk::process: max_cross_eqplane is DBL_MAX and probably should not be" << endl;
       int nb_cross_eqplane = ph->nb_cross_eqplane();
       double fudge_Gralla=1.;
       if (nb_cross_eqplane>0) fudge_Gralla=1.;
       inc = fudge_Gralla
 	* (emission(ggredm1, coord_ph_hit, coord_obj_hit))
 	* (ph -> getTransmission(size_t(-1)))
 	* ggred*ggred*ggred*ggred; // I/nu^4 invariant
       *data->user4 += inc;
 #if GYOTO_DEBUG_ENABLED
       GYOTO_DEBUG_EXPR(*data->user4);
 #endif
      
     } else if (data->spectrum) {
       // *** CAREFUL!! ***
       /*
 	I have to include here the "fudge factor" of Gralla+.
 	Do not forget to remove it to consider some other disk profile.
       */
       double * nuem         = new double[nbnuobs];

	for (size_t ii=0; ii<nbnuobs; ++ii) {
	  nuem[ii]=nuobs[ii]*ggredm1;
	}
       double max_cross_eqplane = ph->maxCrossEqplane();
       if (max_cross_eqplane==DBL_MAX) cout << "WARNING: in StochasticThinDisk::process: max_cross_eqplane is DBL_MAX and probably should not be" << endl;
       int nb_cross_eqplane = ph->nb_cross_eqplane();
       double fudge_Gralla=1.;
       
       if (nb_cross_eqplane>0) fudge_Gralla=1.;
       for (size_t ii=0; ii<nbnuobs; ++ii) {
         inc = fudge_Gralla
 	  * (emission(nuem[ii], coord_ph_hit, coord_obj_hit))
 	  * (ph -> getTransmission(size_t(-1)))
 	  * ggred*ggred*ggred; // I/nu^3 invariant
         *data->spectrum += inc; // data->spectrum[ii*data->offset] += inc  !Error
       }
       #if GYOTO_DEBUG_ENABLED
         GYOTO_DEBUG_EXPR(*data->spectrum);
       #endif
       
       
     } else if (data->redshift) {
       *data->redshift=ggred;
       #if GYOTO_DEBUG_ENABLED
         GYOTO_DEBUG_EXPR(*data->redshift);
       #endif
    } else if (data->impactcoords ) { // && data->impactcoords[0]==DBL_MAX to get only the n=0 coordinates
        if (coord_ph_hit.size() > 8) GYOTO_ERROR("ImpactCoords is incompatible with parallel transport");
        memcpy(data->impactcoords, coord_obj_hit, 8 * sizeof(double));
        memcpy(data->impactcoords+8, &coord_ph_hit[0], 8 * sizeof(double));
        //cout << setprecision(18);
       //cout << endl << "t,r,theta,phi: " << coord_obj_hit[0] << " " << coord_obj_hit[1] << " " << coord_obj_hit[2] << " " << coord_obj_hit[3] << endl; 
    } else GYOTO_ERROR("unimplemented data");
 }
}
