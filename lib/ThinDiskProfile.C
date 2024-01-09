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
#include "GyotoRezzollaZhidenko.h"

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
using namespace Gyoto::Metric;

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
  
  // Emission radius
  //cerr << "Emission radius"; //cout
  //std::ofstream outfile;
  //outfile.open("./Bug.txt", std::ios_base::app); // append instead of overwrite
  //outfile << "Radius: " << rr << "\n"; 
  //outfile.close();

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
     //model_param = [zeta, rin, norm, indx1=alpha-2beta, indx2=gamma+2beta], with n_e propto r^-alpha, Theta_e propto r^-beta, B propto r^-gamma
     double zeta = model_param_[0],
      rin = model_param_[1], //1.+sqrt(1-SPIN*SPIN);
      norm = model_param_[2],
      indx1 = model_param_[3],
      indx2 = model_param_[4];
      //cout << "zeta= " << zeta << endl;
      //cout << "nu_em= " << nu << endl;    
     
     // Equation B.7, Appendix B (nu_em=cst=230GHz)
     //emiss = norm*exp(-zeta*rr/rin);
     // the 1e-3 for zeta=3 is just there to get a reasonable flux for M87
     
     // Equation B.4 with nu_em dependence and free indices of power law
     // Decomment the following line (with indx1 = 0, indx2 = 3) to check that in this case the simplified emission is retrieved
     //nu = 230e9; 
     // Power laws for n_e, Theta_e, B propto r^-2, r^-1, r^-1
     //emiss = norm*nu/230*exp(-zeta*pow(230,-1./3)*pow(nu,1./3.)*rr/rin);
     // Power laws for n_e, Theta_e, B propto r^-alpha, r^-beta, r^-gamma
     emiss = norm*nu*1e-9/230*pow(rr,-indx1)*exp(-zeta*pow(230,-1./3)*pow(nu*1e-9,1./3.)*pow(rr/rin,indx2/3.));
     //cout << setprecision(16);
     //cout << "emiss " << emiss << endl; 
    
     // Check
     //std::ofstream outfile;
     //outfile.open("./Check/Check_Emission_B4_nu_general_new.txt", std::ios_base::app); // append instead of overwrite
     //outfile << "Emission: " << emiss << "\n"; 
     //outfile.close();
  }
  
  return emiss;
}

void ThinDiskProfile::getVelocity(double const pos[4], double vel[4])
{

  string kin = gg_->kind();
  
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
       if (kin == "KerrBL"){
         double SPIN = static_cast<SmartPointer<Metric::KerrBL> >(gg_) -> spin();
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
	} else if (kin == "RezzollaZhidenko") {
	// See formulas in Cárdenas, Godfrey, Yunes & Lohfink 2019, Eqs 11 and 12 for E and L
	// initally from Cunnigham 1975
	  double N2ms = static_cast<SmartPointer<Metric::RezzollaZhidenko> >(gg_) -> N2(risco);
	  double Nprimems = static_cast<SmartPointer<Metric::RezzollaZhidenko> >(gg_) -> Nprime(risco);
	  double g00 = static_cast<SmartPointer<Metric::RezzollaZhidenko> >(gg_) -> gmunu(pos,0,0);
	  double grr = static_cast<SmartPointer<Metric::RezzollaZhidenko> >(gg_) -> gmunu(pos,1,1);
	  double gpp = static_cast<SmartPointer<Metric::RezzollaZhidenko> >(gg_) -> gmunu(pos,3,3);
	  double NNms = sqrt(N2ms); 
	  double Ems = sqrt(pow(NNms,3)/(NNms-risco*Nprimems));
	  double Lms = sqrt(pow(risco,3)*Nprimems/(NNms-risco*Nprimems));
	  
	  vel[0] = -Ems/g00; // this is: -Ems*(g_{00}^{-1})
	  vel[2] = 0.; // this is: 0.
	  vel[3] = Lms/pow(rr,2); // this is: Lms/r**2
	  vel[1] = -sqrt(-pow(grr,-1)*(1.+g00*pow(vel[0],2)+gpp*pow(vel[3],2))); 
	  // this is: -sqrt{-(g_{rr}^{-1})*(1 + g_{00}*(u^t)^2 + g_{pp}*(u^p)^2)}
	  
	  // Velocity to see the discontinuity for n=0
	  //vel[1] = 0.; 
	  //vel[2] = 0.;
	  //vel[3] = 0.;
	  //vel[0] = sqrt(-pow(g00,-1));

	  // Check the normalisation of the 4-velocity
	  double tol=1e-5;
	  double u2 = gg_->ScalarProd(pos,vel,vel);
	  if (fabs(u2+1.)>tol or u2!=u2) {
	  cerr << " *** 4-velocity squared norm= " << u2 << endl;
	  throwError("In ThinDiskProfile: 4vel RezzollaZhidenko is not properly normalized!");
	  }
	} else {
	  GYOTO_ERROR("ThinDiskProfile: KerrBL or RezzollaZhidenko needed!");
	}
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
       if (max_cross_eqplane==DBL_MAX) cout << "WARNING: in ThinDiskProfile::process: max_cross_eqplane is DBL_MAX and probably should not be" << endl;
       int nb_cross_eqplane = ph->nb_cross_eqplane();
       double fudge_Gralla=1.;
       if (nb_cross_eqplane>0) fudge_Gralla=1.5;
       inc = fudge_Gralla
 	* (emission(ggredm1, dsem, coord_ph_hit, coord_obj_hit))
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
       if (max_cross_eqplane==DBL_MAX) cout << "WARNING: in ThinDiskProfile::process: max_cross_eqplane is DBL_MAX and probably should not be" << endl;
       int nb_cross_eqplane = ph->nb_cross_eqplane();
       double fudge_Gralla=1.;
       //cout << "nb cross fudge " << nb_cross_eqplane << endl;
       if (nb_cross_eqplane>0) fudge_Gralla=1.5;
       for (size_t ii=0; ii<nbnuobs; ++ii) {
         inc = fudge_Gralla
 	  * (emission(nuem[ii], dsem, coord_ph_hit, coord_obj_hit))
 	  * (ph -> getTransmission(size_t(-1)))
 	  * ggred*ggred*ggred; // I/nu^3 invariant
         *data->spectrum += inc; // data->spectrum[ii*data->offset] += inc  !Error
       }
 #if GYOTO_DEBUG_ENABLED
       GYOTO_DEBUG_EXPR(*data->spectrum);
 #endif
     }
       else GYOTO_ERROR("unimplemented data");
   }
 }
