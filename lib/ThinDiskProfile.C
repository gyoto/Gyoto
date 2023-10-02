/*
    Copyright 2020 Frederic Vincent, Thibaut Paumard

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
#include "GyotoThinDiskProfile.h"
#include "GyotoProperty.h"
#include "GyotoUtils.h"
#include "GyotoFactoryMessenger.h"
#include "GyotoKerrBL.h"
#include "GyotoKerrKS.h"


#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <cmath>
#include <limits>
#include <string>

using namespace std;
using namespace Gyoto;
using namespace Gyoto::Astrobj;

GYOTO_PROPERTY_START(ThinDiskProfile)
GYOTO_PROPERTY_BOOL(ThinDiskProfile, CircularMotion, NoCircularMotion,
		    circularMotion)
GYOTO_PROPERTY_VECTOR_DOUBLE(ThinDiskProfile, Model_param, model_param,
			     "Parameters useful for the disk, max number NPAR_MAX")
GYOTO_PROPERTY_END(ThinDiskProfile, ThinDisk::properties)

//#define SPIN 0.94 // Kerr spin parameter for Kerr-specific formulas...
#define NPAR_MAX 10 // Max allowed number of parameters

bool ThinDiskProfile::circularMotion() const {return circular_motion_;}
void ThinDiskProfile::circularMotion(bool circ) {circular_motion_=circ;}

void ThinDiskProfile::model_param(std::vector<double> const &v) {
  size_t n = v.size();
  if (n>NPAR_MAX) throwError("Too many parameters in model_param");
  for (size_t i=0; i<n; ++i) model_param_[i]=v[i];
}
std::vector<double> ThinDiskProfile::model_param() const {
  std::vector<double> v(NPAR_MAX, 0.);
  for (size_t i=0; i<NPAR_MAX; ++i) v[i]=model_param_[i];
  return v;
}

ThinDiskProfile::ThinDiskProfile() :
  ThinDisk("ThinDiskProfile"),
  circular_motion_(1),
  model_param_(NULL)
{
  GYOTO_DEBUG << endl;
  model_param_ = new double[NPAR_MAX];
  for (int ii=0;ii<NPAR_MAX;ii++) model_param_[ii]=0.;
}

ThinDiskProfile::ThinDiskProfile(const ThinDiskProfile& o) :
  ThinDisk(o),
  circular_motion_(o.circular_motion_),
  model_param_(NULL)
{
  if (o.gg_()) gg_=o.gg_->clone();
  Generic::gg_=gg_;

  model_param_ = new double[NPAR_MAX];
  for (int ii=0;ii<NPAR_MAX;ii++) model_param_[ii]=o.model_param_[ii];
}
ThinDiskProfile* ThinDiskProfile::clone() const
{ return new ThinDiskProfile(*this); }

bool ThinDiskProfile::isThreadSafe() const {
  return ThinDisk::isThreadSafe();
}

ThinDiskProfile::~ThinDiskProfile() {
  GYOTO_DEBUG << endl;
  delete [] model_param_;
}

double ThinDiskProfile::emission(double nu, double,
			    state_t const &,
			    double const coord_obj[8]) const{
  string emission_model = "Thermal_Synchrotron"; // should be in: "Gralla_et_al", "Thermal_Synchrotron"

  double rr = coord_obj[1],
    emiss=0.; // model-dependent emission defined below

  if (emission_model == "Gralla_et_al"){
    // Emission from Gralla et al 2020
    // ****************************** //
    // Here model_param must contain:
    // model_param = [gamma, mu, sigG]
    
    // This model is Kerr specific
    string kin = gg_->kind();      
    if (kin != "KerrBL")
      GYOTO_ERROR("ThinDiskProfile: KerrBL needed!");
    double SPIN = static_cast<SmartPointer<Metric::KerrBL> >(gg_) -> spin(),
      a2 = SPIN*SPIN;
      
    double rhor=1.+sqrt(1.-a2), rminus=1.-sqrt(1.-a2),
      risco=gg_->getRms();
    
    // Choose profile here:
    double gamma=model_param_[0],
      mu=model_param_[1],
      sigG=model_param_[2];
    //double gamma=-3./2., mu=rminus, sigG=1./2.;
    //double gamma=-3., mu=risco-0.33, sigG=0.25;
    
    double tmp = gamma+asinh((rr-mu)/sigG);
    emiss = 1e-5*exp(-0.5*tmp*tmp)/sqrt((rr-mu)*(rr-mu)+sigG*sigG);
    // the 1e-5 is just there to get a reasonable flux for M87
  }
  
  if (emission_model == "Thermal_Synchrotron"){
    // Emission from Vincent+22 thick disk paper synchrotron formula
    // ****************************** //
    // Here model_param must contain:
    // model_param = [zeta, rin]
    double zeta = model_param_[0],
      rin = model_param_[1]; //1.+sqrt(1-SPIN*SPIN);
    emiss = 1e-3*exp(-zeta*rr/rin);
    // the 1e-3 is just there to get a reasonable flux for M87
  }
  
  return emiss;
}

void ThinDiskProfile::getVelocity(double const pos[4], double vel[4])
{

  string kin = gg_->kind();
  if (kin != "KerrBL")
    GYOTO_ERROR("ThinDiskProfile: KerrBL needed!");
  double SPIN = static_cast<SmartPointer<Metric::KerrBL> >(gg_) -> spin();
  
  double risco = 0.;
  if (gg_->kind()!="Minkowski" && gg_->kind()!="Hayward")
    risco=gg_->getRms(); // prograde Kerr ISCO
  // let risco=0 if metric is Minko; then ISCO not defined
  // let also risco=0 for Hayward as we would need to
  // compute it numerically and give it in xml Metric field,
  // not implemented so far

  //cout << "in velo, r isco= " << pos[1] << " " << risco << endl;
  double rr = pos[1];
  //cout << "circ=" << circular_motion_ <<endl;
  if (circular_motion_) {
    // CIRCULAR ROTATION
      if (rr > risco){
	// Keplerian velocity above ISCO
	gg_ -> circularVelocity(pos, vel, 1);
      }else{
	// See formulas in Gralla, Lupsasca & Marrone 2020, Eqs B8-B14
	// initally from Cunnigham 1975
	double lambda_ms = (risco*risco - 2.*SPIN*sqrt(risco) + SPIN*SPIN)/(pow(risco,1.5) - 2.*sqrt(risco) + SPIN),
	gamma_ms = sqrt(1.-2./(3.*risco)),
	delta = rr*rr - 2.*rr + SPIN*SPIN,
	hh = (2.*rr - SPIN*lambda_ms)/delta;
	
	vel[0] = gamma_ms*(1.+2./rr*(1.+hh)); // this is: -Ems*g^{tt} + Lms*g^{tp}
	vel[1] = -sqrt(2./(3.*risco))*pow(risco/rr-1.,1.5); // this is: -sqrt{(-1 - g_{tt}*u^t - g_{pp}*u^p - 2*g_{tp}*u^t*u^p)/grr}
	vel[2] = 0.;
	vel[3] = gamma_ms/(rr*rr)*(lambda_ms+SPIN*hh);
	
	//cout << "u2 = " << gg_->ScalarProd(pos,vel,vel) << endl;
      }
    }else{
    // RADIAL FALL
    double gtt = gg_->gmunu(pos,0,0),
      grr = gg_->gmunu(pos,1,1),
      guptt = gg_->gmunu_up(pos,0,0),
      guptp = gg_->gmunu_up(pos,0,3),
      guprr = gg_->gmunu_up(pos,1,1);
    
    // 4-vel obtained by imposing: u_t=-1, u_phi=0, u^theta=0
    // see FV notes SphericalVelocity.pdf for details
    vel[0] = -guptt;
    vel[1] = -sqrt((-1.-guptt)*guprr);
    vel[2] = 0;
    vel[3] = -guptp;

    double tol=1e-5;
    double u2 = gg_->ScalarProd(pos,vel,vel);
    //cout << "4vel,u2= " << rr << " " << pos[2] << " " << gtt << " " << grr << " " << vel[0] << " " << vel[1] << " " << vel[2] << " " << vel[3] << " " << u2 << endl;
    
    if (fabs(u2+1.)>tol or u2!=u2) {
      cerr << " *** 4-velocity squared norm= " << u2 << endl;
      throwError("In ThinDiskProfile: 4vel "
		 "is not properly normalized!");
    }
  }
  
  //cout << "4vel= " << vel[0] << " " << vel[1] << " " << vel[2] << " " << vel[3]<< endl;
}

void ThinDiskProfile::processHitQuantities(Photon* ph,
					   state_t const &coord_ph_hit,
					   double const *coord_obj_hit,
					   double dt,
					   Properties* data) const {
#if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << endl;
#endif
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
      if (max_cross_eqplane==DBL_MAX) cout << "WARNING: in ThinDiskProfile::process: max_cross_eqplane is DBL_MAX and probably should not be" << endl;
      int nb_cross_eqplane = ph->nb_cross_eqplane();
      double fudge_Gralla=1.;
      //if (nb_cross_eqplane>0) fudge_Gralla=1.5;
      inc = fudge_Gralla
	* (emission(ggredm1, dsem, coord_ph_hit, coord_obj_hit))
	* (ph -> getTransmission(size_t(-1)))
	* ggred*ggred*ggred*ggred; // I/nu^4 invariant
      *data->user4 += inc;
#if GYOTO_DEBUG_ENABLED
      GYOTO_DEBUG_EXPR(*data->user4);
#endif
      
    }else GYOTO_ERROR("unimplemented data");
  }
}
