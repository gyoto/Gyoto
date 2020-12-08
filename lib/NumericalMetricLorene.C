/*
    Copyright 2014-2020 Frederic Vincent, Thibaut Paumard

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

// Lorene headers
#include "metric.h"
#include "nbr_spx.h"
#include "utilitaires.h"
#include "graphique.h"
#include "cmp.h"

using namespace Lorene;

//Gyoto headers
#include "GyotoUtils.h"
#include "GyotoNumericalMetricLorene.h"
#include "GyotoError.h"
#include "GyotoFactoryMessenger.h"
#include "GyotoProperty.h"

//Std headers
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <dirent.h>
#include <ctime>

using namespace Gyoto;
using namespace Gyoto::Metric;

GYOTO_PROPERTY_START(NumericalMetricLorene)
GYOTO_PROPERTY_BOOL(NumericalMetricLorene, MapEt, MapAf, mapEt)
GYOTO_PROPERTY_BOOL(NumericalMetricLorene, AxisymCirc, NoAxisymCirc, axisymCirc)
GYOTO_PROPERTY_BOOL(NumericalMetricLorene,
		    SpecifyMarginalOrbits, DontSpecifyMarginalOrbits,
		    specifyMarginalOrbits)
GYOTO_PROPERTY_BOOL(NumericalMetricLorene,
		    HasSurface, HasNoSurface, hasSurface)
GYOTO_PROPERTY_BOOL(NumericalMetricLorene,
		    HasAccelerationVector, HasNoAccelerationVector,
		    hasAccelerationVector)
GYOTO_PROPERTY_BOOL(NumericalMetricLorene, BosonStarCircular, NonBosonStarCircular, bosonstarcircular)
GYOTO_PROPERTY_DOUBLE(NumericalMetricLorene, Horizon, horizon)
GYOTO_PROPERTY_DOUBLE(NumericalMetricLorene, Time, initialTime)
GYOTO_PROPERTY_DOUBLE(NumericalMetricLorene, Rico, rico)
// Keep File last here, so it is processed last in fillElement() 
// (just before the generic Properties, that is
GYOTO_PROPERTY_FILENAME(NumericalMetricLorene, File, directory)
GYOTO_PROPERTY_END(NumericalMetricLorene, Generic::properties)

// Lorene Metrics are not thread-safe
GYOTO_PROPERTY_THREAD_UNSAFE(NumericalMetricLorene)

#define GYOTO_NML_PPHI_TOL 5 // tolerance on p_phi drift, percentage

NumericalMetricLorene::NumericalMetricLorene() :
  Generic(GYOTO_COORDKIND_SPHERICAL, "NumericalMetricLorene"),
  filename_(NULL),
  mapet_(true),
  axisymCirc_(false),
  bosonstarcircular_(false),
  has_surface_(0),
  has_acceleration_vector_(0),
  specify_marginalorbits_(0),
  horizon_(0.),
  initial_time_(0.),
  lapse_tab_(NULL),
  shift_tab_(NULL),
  gamcov_tab_(NULL),
  gamcon_tab_(NULL),
  kij_tab_(NULL),
  times_(NULL),
  nb_times_(0),
  nssurf_tab_(NULL),
  vsurf_tab_(NULL),
  accel_tab_(NULL),
  lorentz_tab_(NULL),
  hor_tab_(NULL),
  risco_(0.),
  rico_(0.),
  rmb_(0.)
{
  GYOTO_DEBUG << endl;
}

NumericalMetricLorene::NumericalMetricLorene(const NumericalMetricLorene&o) :
  Generic(GYOTO_COORDKIND_SPHERICAL,"NumericalMetricLorene"),
  filename_(NULL),
  mapet_(o.mapet_),
  axisymCirc_(o.axisymCirc_),
  bosonstarcircular_(o.bosonstarcircular_),
  has_surface_(o.has_surface_),
  has_acceleration_vector_(o.has_acceleration_vector_),
  specify_marginalorbits_(o.specify_marginalorbits_),
  horizon_(o.horizon_),
  initial_time_(o.initial_time_),
  lapse_tab_(NULL),
  shift_tab_(NULL),
  gamcov_tab_(NULL),
  gamcon_tab_(NULL),
  kij_tab_(NULL),
  times_(NULL),
  nb_times_(0),
  nssurf_tab_(NULL),
  vsurf_tab_(NULL),
  accel_tab_(NULL),
  lorentz_tab_(NULL),
  hor_tab_(NULL),
  risco_(o.risco_),
  rico_(o.rico_),
  rmb_(o.rmb_)
{
  GYOTO_DEBUG << endl;
  if (o.filename_) directory(o.filename_);
}

NumericalMetricLorene* NumericalMetricLorene::clone() const{
  GYOTO_DEBUG << endl;
  return new NumericalMetricLorene(*this);
}

NumericalMetricLorene::~NumericalMetricLorene() 
{
  GYOTO_DEBUG<< endl;
  free();
}

void NumericalMetricLorene::free() {
  GYOTO_DEBUG << "freeing memory\n";
  if (filename_)   { delete [] filename_;   filename_=NULL;  }
  if (lapse_tab_)  { delete [] lapse_tab_;  lapse_tab_=NULL; }
  if (shift_tab_)  { delete [] shift_tab_;  shift_tab_=NULL; }
  if (gamcov_tab_) { delete [] gamcov_tab_; gamcov_tab_=NULL;}
  if (gamcon_tab_) { delete [] gamcon_tab_; gamcon_tab_=NULL;}
  if (kij_tab_)    { delete [] kij_tab_;    kij_tab_=NULL;   }
  if (times_)      { delete [] times_;      times_=NULL;     }
  if (nssurf_tab_) { delete [] nssurf_tab_;      nssurf_tab_=NULL;     }
  if (vsurf_tab_) { delete [] vsurf_tab_;      vsurf_tab_=NULL;     }
  if (accel_tab_) { delete [] accel_tab_;      accel_tab_=NULL;     }
  if (lorentz_tab_) { delete [] lorentz_tab_;      lorentz_tab_=NULL;     }
  if (hor_tab_) { delete [] hor_tab_;      hor_tab_=NULL;     }
}

void NumericalMetricLorene::setMetricSource() {
  GYOTO_DEBUG << endl;
  DIR *dp;
  struct dirent *dirp;
  if((dp  = opendir(filename_)) == NULL) {
    GYOTO_ERROR("In NumericalMetricLorene.C constructor : bad filename_");
  }

  nb_times_=0;
  while ((dirp = readdir(dp)) != NULL) {
    nb_times_++;
  }
  nb_times_-=2; //for directories . and .. 

  GYOTO_DEBUG << "Nb of metric files= " << nb_times_ << endl;
  /*
    NB: ***Caution***, here it is assumed that filename_ contains ONLY the .d 
    Lorene result files, nothing else.
  */
  closedir(dp);

  if (nb_times_<1) 
    GYOTO_ERROR("In NumericalMetricLorene.C: bad nb_times_ value");

  lapse_tab_ = new Scalar*[nb_times_] ; 
  shift_tab_ = new Vector*[nb_times_] ;
  gamcov_tab_ = new Sym_tensor*[nb_times_] ;  
  gamcon_tab_ = new Sym_tensor*[nb_times_] ;  
  kij_tab_ = new Sym_tensor*[nb_times_] ; 
  times_ = new double[nb_times_];

  if (has_surface_){
    nssurf_tab_ = new Valeur*[nb_times_] ;
    vsurf_tab_ = new Vector*[nb_times_] ;
    lorentz_tab_ = new Scalar*[nb_times_] ;
    hor_tab_ = new Valeur*[nb_times_] ;
    if (has_acceleration_vector_)
      accel_tab_ = new Vector*[nb_times_] ;
  }

  if (debug()) {
    cout << "In NumericalMetricLorene" << endl;
    cout << "File name=" << filename_ << endl;
    cout << "Number of time slices=" << nb_times_ << endl;
  }

  if (debug()) cout << "NumericalMetricLorene.C: "
		 "initializing geometrical quantities..." << endl;
  
  for (int i=1; i<=nb_times_; i++) {
    ostringstream stream_name ;
    stream_name << filename_ << "metric" << setw(6) << setfill('0') 
		<< i << ".d" ;

    if (debug()) cout << "Reading file: " << stream_name.str() << endl ;
    FILE* resu = fopen(stream_name.str().data(), "r") ;
    if (resu == 0x0) {
      cerr << "With file name: " << stream_name.str() << endl ;
      GYOTO_ERROR("NumericalMetricLorene.C: Problem opening file!");
    }
    if (debug()) cout << "File read normally." << endl ;
    double cLor = GYOTO_C*1e-3*1e-4;
    /*
      this is c in Lorene units, allows to translate between 
      Lorene times and Gyoto times (Lorene speaks in ms, 10^4 m; 
      Gyoto speaks in natural units): t(Gyoto) = t(Lorene)*cLor
     */
    
    double time ;
    fread_be(&time, sizeof(double), 1, resu) ;
    //cout << "Lorene time= " << setprecision(15) << i << " " << time << endl;
    //cout << "Gyoto time = " << i-1 << " " << initial_time_+time*cLor << endl;
    setTimes(initial_time_+time*cLor,i-1); // ***COLLAPSE TIME A GERER
    
    Mg3d* grid = new Mg3d(resu) ;
    Map* map;
    /* Use Map_af for collapse + Kerr + BS, Map_et for star imaging */
    if (mapet_) { // Map_et case
      //Map_et* map = new Map_et(*grid, resu) ;
      map = new Map_et(*grid, resu) ;
    } else {      // Map_af case
      //      Map_af* map = new Map_af(*grid, resu) ;
      map = new Map_af(*grid, resu) ;
    }
    Scalar* lapse = new Scalar(*map, *grid, resu) ;
    (*lapse).std_spectral_base() ;
    setLapse_tab(lapse,i-1);
    Vector* shift = new Vector(*map, (*map).get_bvect_spher(), resu) ;
    setShift_tab(shift,i-1);
    Sym_tensor* g_ij = new Sym_tensor(*map, (*map).get_bvect_spher(), resu) ;
    Sym_tensor* g_up_ij = new Sym_tensor(*map, (*map).get_bvect_spher(), resu) ;
    setGamcov_tab(g_ij,i-1);
    setGamcon_tab(g_up_ij,i-1);
    Sym_tensor* kij = new Sym_tensor(*map, (*map).get_bvect_spher(), resu) ;
    
    if (has_surface_){
      // This seems to be only necessary for collapsing or not collapsing star
      // --> F.V. October 2015: seems outdated, now produces a bug on dzpuis
      for (int l=1; l<=3; l++)
	for (int c=l; c<=3; c++)
	  (*kij).set(l,c).dec_dzpuis(2) ;
    }

    setKij_tab(kij,i-1);

    if (has_surface_){
      Scalar* lorentz_factor = new Scalar(*map, *grid, resu) ;
      lorentz_tab_[i-1] = lorentz_factor;
      Vector* v_i = new Vector(*map, (*map).get_bvect_spher(), resu) ;
      vsurf_tab_[i-1] = v_i;
      Mg3d* grid_surf = new Mg3d(resu) ;
      Valeur* ns_surf = new Valeur(*grid_surf, resu) ;
      nssurf_tab_[i-1] = ns_surf;
      if (has_acceleration_vector_){
	Vector* a_i = new Vector(*map, (*map).get_bvect_spher(), resu) ;
	accel_tab_[i-1] = a_i ;
      }
      Mg3d* grid_ah = new Mg3d(resu) ;
      Valeur* horizon = new Valeur(*grid_ah, resu) ;
      hor_tab_[i-1] = horizon;
    }

    if (specify_marginalorbits_){
      double r_isco ;
      fread_be(&r_isco, sizeof(double), 1, resu) ;
      risco_ = r_isco ;
      
      if (debug()) cout << "DEBUG: READ Risco = " << r_isco << endl ;
      
      double r_mb ;
      fread_be(&r_mb, sizeof(double), 1, resu) ;
      rmb_ = r_mb;
      
      if (debug()) cout << "DEBUG: READ Rmb = " << r_mb << endl ;
    }

    fclose(resu) ;
  }

  if (debug()) cout << "NumericalMetricLorene.C constructor: "
		 "geometrical quantities initialized." << endl;

}

Sym_tensor** NumericalMetricLorene::getGamcon_tab() const {
  GYOTO_DEBUG << endl;
  return gamcon_tab_;}
Sym_tensor** NumericalMetricLorene::getGamcov_tab() const {
  GYOTO_DEBUG << endl;
  return gamcov_tab_;}
Vector** NumericalMetricLorene::getShift_tab() const {
  GYOTO_DEBUG << endl;
  return shift_tab_;}
Scalar** NumericalMetricLorene::getLapse_tab() const {
  GYOTO_DEBUG << endl;
  return lapse_tab_;}
double* NumericalMetricLorene::getTimes() const {
  GYOTO_DEBUG << endl;
  return times_;}
int NumericalMetricLorene::getNbtimes() const {
  GYOTO_DEBUG << endl;
  return nb_times_;}
Valeur** NumericalMetricLorene::getNssurf_tab() const {
  GYOTO_DEBUG << endl;
  return nssurf_tab_;}
Vector** NumericalMetricLorene::getVsurf_tab() const {  
  GYOTO_DEBUG << endl;
  return vsurf_tab_;}
Vector** NumericalMetricLorene::getAccel_tab() const {  
  GYOTO_DEBUG << endl;
  return accel_tab_;}
Scalar** NumericalMetricLorene::getLorentz_tab() const {
  GYOTO_DEBUG << endl;
  return lorentz_tab_;}
Valeur** NumericalMetricLorene::getHor_tab() const {
  GYOTO_DEBUG << endl;
  return hor_tab_;}
double NumericalMetricLorene::getRms() const {
  GYOTO_DEBUG << endl;
  if (rico()!=0.) return rico();
  else return risco_;
}  
double NumericalMetricLorene::getRmb() const {
  GYOTO_DEBUG << endl;
  return rmb_;}

double NumericalMetricLorene::getSpecificAngularMomentum(double rr) const {
  // Computes the Keplerian specific angular momentum \ell = -u_phi / u_t
  // for circular geodesics,
  // for a general axisym metric in the equatorial plane

  /*
    // Standard formula: leads to problem at least in Kerr where
    // gphph_dr goes to zero outside horizon.

    double pos[4]={0., rr, M_PI/2., 0.};

    double gtt_dr  =gmunu_up_dr(pos, 0, 0),
    gtph_dr =gmunu_up_dr(pos, 0, 3), 
    gphph_dr=gmunu_up_dr(pos, 3, 3);
    
    double lKep = gtph_dr/gphph_dr + 
    sqrt(gtph_dr/gphph_dr * gtph_dr/gphph_dr - gtt_dr/gphph_dr);
    
    if (lKep!=lKep || lKep==lKep+1.){
    cerr << "At r= " << rr << endl;
    GYOTO_ERROR("In NML::getSpecificAngMom: lKep not defined here!"
    " You are probably below the innermost circular orbit.");
    }
  */

  if (nb_times_>1) GYOTO_ERROR("In NML::getSpecificAngularMomentum:"
			      "so far only stationary metric implemented");

  int indice_time=0;
  double th=M_PI/2., ph=0.; // in equatorial plane
  double rm1 = 1./rr, rsm1 = rm1, rm2 = rm1*rm1, sm1 = 1.; // NB: sinth=1
  const Sym_tensor& g_ij = *(gamcov_tab_[indice_time]) ;
  double B2 = g_ij(3,3).val_point(rr,th,ph); // no mistake here, B2 is g_pp_Lorene, but g_pp_Gyoto/(r2sinth2)
  if (B2<=0.) GYOTO_ERROR("In NML::getSpecificAngMom: bad B2");
  double BB = sqrt(B2);
  double Br = g_ij(3,3).dsdr().val_point(rr,th,ph)/(2.*BB);
  const Vector& shift = *(shift_tab_[indice_time]);
  double beta_p = rsm1*shift(3).val_point(rr,th,ph),
    beta_p_r = rsm1*shift(3).dsdr().val_point(rr,th,ph)
    -rm2*sm1*shift(3).val_point(rr,th,ph);
  Scalar* lapse = lapse_tab_[indice_time];
  double NN = lapse -> val_point(rr,th,ph);
  if (NN==0.) GYOTO_ERROR("In NML::getSpecificAngMom: bad N");
  double Nr = lapse->dsdr().val_point(rr,th,ph);
  double DD = B2*rr*rr/(NN*NN)*beta_p_r*beta_p_r
    + 4.*Nr/NN*(Br/BB+rm1);
  if (DD<0.) GYOTO_ERROR("In NML::getSpecificAngMom: bad D");
  double Vzamo = 0.5*(-BB*rr/NN*beta_p_r+sqrt(DD))/(rm1+Br/BB);
  
  // 3+1 l_Kep for any circular QI-coord spacetime:
  double lKep = BB*rr*Vzamo/(NN-beta_p*BB*rr*Vzamo); 

  return lKep;
}

double NumericalMetricLorene::getPotential(double const pos[4], double l_cst) const {
  // returns W= -log(abs(u_t)), so that PD::operator, returning Wsurf-W,
  // is negative inside doughnut

  double gtt=gmunu(pos,0,0);
  double gtph=gmunu(pos,0,3); 
  double gphph=gmunu(pos,3,3);

  double u_t_squared = (gtph*gtph-gtt*gphph)
    /(gtt*l_cst*l_cst+2.*l_cst*gtph+gphph);

  if (u_t_squared<0.) return -DBL_MAX; // so Wsurf-W is always >0, not hit,
                                       // 4-velocity not defined
  
  return  -log(sqrt(u_t_squared)) ;
}


void NumericalMetricLorene::setLapse_tab(Scalar* lapse, int ii) {
  GYOTO_DEBUG << endl;
  lapse_tab_[ii]=lapse;}
void NumericalMetricLorene::setShift_tab(Vector* shift, int ii) {
  GYOTO_DEBUG << endl;
  shift_tab_[ii]=shift;}
void NumericalMetricLorene::setGamcov_tab(Sym_tensor* gamcov, int ii) {
  GYOTO_DEBUG << endl;
  gamcov_tab_[ii]=gamcov;}
void NumericalMetricLorene::setGamcon_tab(Sym_tensor* gamcon, int ii) {  
  GYOTO_DEBUG << endl;
  gamcon_tab_[ii]=gamcon;}
void NumericalMetricLorene::setKij_tab(Sym_tensor* kij, int ii) {
  GYOTO_DEBUG << endl;
  kij_tab_[ii]=kij;}
void NumericalMetricLorene::setTimes(double time, int ii) {
  GYOTO_DEBUG << endl;
  times_[ii]=time;}

int NumericalMetricLorene::diff(state_t const &coord,
				state_t &res,
				double mass) const{
  double rhor=computeHorizon(&coord[0]);
  if (coord[1]<rhor && rhor>0.) {
    GYOTO_DEBUG << "rr, rhor= " << coord[1] << " " << rhor << endl;
    GYOTO_DEBUG << "Sub-horizon r, stop" << endl;
    return 1; 
  }

  return Generic::diff(coord, res, mass);
}

int NumericalMetricLorene::diff31(const state_t &x,
			  state_t &dxdt,
			  double /* mass */) const {
  return diff(0.,&x[0],&dxdt[0]); // first slot should be time!!
}

int NumericalMetricLorene::diff(double tt, 
				const double y[7], double res[7]) const
{
  GYOTO_DEBUG << endl;
  /*
    3+1 diff called by RK4 WITH ENERGY INTEG, that itself calls the correct 
    diff(double*,double*,int indice_time)
    with indice_time indicating which metric to consider
    This diff then returns the simple linear interpolation between 
    the 2 metrics such that t(metric1)>t>=t(metric2)
   */
  double rr=y[1], th=y[2], phi=y[3];
  double pos[4]={tt,rr,th,phi};
  double rhor=computeHorizon(pos);

  //cout << endl;
  //cout << "current t,r,rhoriz= " << setprecision(10) << tt << " " << rr << " " << th << " " << rhor << endl;

  if (rr<rhor && rhor>0.) {
    //horizon: stop integration
    if (debug()){
      cout << "In NumericalMetricLorene::diff() ";
      cout << "rr, rhor= " << rr << " " << rhor << endl;
      cout << "Sub-horizon r, stop"
	   << endl;
    }
    return 1; 
  }

  int it=nb_times_-1;
  while(tt<times_[it] && it>=0){ //ASSUMES backward integration, to generalize
    it--;
  }

  //if (it==0) it=-1; //TEST!!!
  //  if (rr<0.187) it=1305; // TEST!!!! 
  // if (rr<1.) return 1; // TEST!!!! 
  //if (it!=4056) return 1; // TEST
  
  if (debug()){
    cout << "**** metric number= " << it << endl;
  }

  if (it==nb_times_-1) {
    return diff(y,res,nb_times_-1); //use metric nb nb_times_-1
                                    //for all times > max(metric times)
  }

  if (it==-1) {
    return diff(y,res,0); //use metric nb 0 for all times < min(metric times)
  }

  if (it==nb_times_-2 || it==0){ // Linear interp for extremal-1 points
    double res1[7], res2[7];
    double t1=times_[it], t2=times_[it+1];
    if (diff(y,res1,it) || diff(y,res2,it+1)) return 1;
    for (int ii=0;ii<7;ii++) res[ii] = 
			       (res2[ii]-res1[ii])/(t2-t1)*(tt-t1) 
			       + res1[ii];
    return 0;
  }
    
  //Else: use 3rd order interp
  double res1[7], res2[7], res3[7], res4[7];
  if (diff(y,res1,it-1) || diff(y,res2,it) || 
      diff(y,res3,it+1) || diff(y,res4,it+2)) return 1;

  double values[4];
  for (int ii=0;ii<7;ii++) {
    values[0]=res1[ii];
    values[1]=res2[ii];
    values[2]=res3[ii];
    values[3]=res4[ii];
    res[ii] = Interpol3rdOrder(tt,it,values);
  }
  
  return 0;
}

int NumericalMetricLorene::diff(const double y[7], 
				double res[7], int indice_time) const
{
  GYOTO_DEBUG << endl;
  /*
    3+1 diff function WITH ENERGY INTEG, computing the derivatives of
    the 7-vector : y=[E,r,th,ph,Vr,Vth,Vph] (see CQG 3+1 paper for
    definition) Calls metric nb indice_time among the various Lorene
    metric at hand.

    This diff is the most general possible, there is no symmetry assumption.
   */

  /*
    clock_t time1, time2;
    double diftime, clocks = CLOCKS_PER_SEC;
    time1 = clock();
  */

  if (indice_time<0 || indice_time>nb_times_-1) {
    GYOTO_ERROR("NumericalMetricLorene::diff: incoherent value of indice_time");
  }

  //NB: here t=theta, not time!
  double EE=y[0], rr=y[1], th=y[2], phi=y[3], sth=0, cth=0;
  sincos(th, &sth, &cth);
  double rsinth = rr*sth,
    r2sinth2 = rsinth*rsinth, sinth2 = sth*sth;
  if (rr==0.) GYOTO_ERROR("In NumericalMetricLorene.C::diff r is 0!");
  if (rsinth==0.) GYOTO_ERROR("In NumericalMetricLorene.C::diff on z axis!");
  double rm1 = 1./rr, rm2 = rm1*rm1, r2 = rr*rr, rsm1 = 1./rsinth, 
    sm1 = 1./sth, sm2 = sm1*sm1;
  double Vr=y[4], Vth=y[5], Vph=y[6];

  /*
    Important remark!  Lorene uses the orthonormal spherical tetrad
    d/dr, 1/r d/dtheta, 1/(r*sin(theta)) d/dphi, and not the natural
    basis d/dr, d/dtheta, d/dphi There's thus a change of basis to
    come back to the natural basis of spherical coordinates.  This is
    why there are plenty of factors r2, 1/r2, 1/rsinth2, ... everywhere
   */
  //clock_t t1=clock();

  if (!axisymCirc_){
    // GENERAL SPACETIME WITH NO SYMMETRY
    
    //LAPSE
    Scalar* lapse = lapse_tab_[indice_time];
    double NN = lapse -> val_point(rr,th,phi), NNm1 = 1./NN,
      N_dr = lapse->dsdr().val_point(rr,th,phi), 
      N_dt = lapse->dsdt().val_point(rr,th,phi),
      N_dp = sth*lapse->stdsdp().val_point(rr,th,phi); // stdsdp is 1/sin(th)*d/dphi
    //cout << "at r th ph= " << rr << " " << th << " " << phi << endl; 
    //cout << "lapse= " << NN << " " << N_dr << " " << N_dt << " " << N_dp << endl;
    //SHIFT
    const Vector& shift = *(shift_tab_[indice_time]);
    double betar = shift(1).val_point(rr,th,phi),
      betar_dr = shift(1).dsdr().val_point(rr,th,phi),
      betar_dt = shift(1).dsdt().val_point(rr,th,phi),
      betar_dp = sth*shift(1).stdsdp().val_point(rr,th,phi),
      betat_Lorene = shift(2).val_point(rr,th,phi),
      betat = rm1*betat_Lorene,
      betat_dr = rm1*shift(2).dsdr().val_point(rr,th,phi)
      -rm2*betat_Lorene,
      betat_dt = rm1*shift(2).dsdt().val_point(rr,th,phi),
      betat_dp = sth*rm1*shift(2).stdsdp().val_point(rr,th,phi),
      betap_Lorene = shift(3).val_point(rr,th,phi),
      betap = rsm1*betap_Lorene,
      betap_dr = rsm1*shift(3).dsdr().val_point(rr,th,phi)
      -rm2*sm1*betap_Lorene,
      betap_dt = rsm1*shift(3).dsdt().val_point(rr,th,phi)
      -cth*sm2*rm1*betap_Lorene,
      betap_dp = rm1*shift(3).stdsdp().val_point(rr,th,phi);
    //cout << "shift= " << betar << " " << betat << " " << betap << endl;
    
    //3-METRIC DERIVATIVES (for Christo)
    const Sym_tensor& g_ij = *(gamcov_tab_[indice_time]);
    double g_rr_dr = g_ij(1,1).dsdr().val_point(rr,th,phi),
      g_rr_dt = g_ij(1,1).dsdt().val_point(rr,th,phi),
      g_rr_dp = sth*g_ij(1,1).stdsdp().val_point(rr,th,phi),
      g_tt_dr = r2*g_ij(2,2).dsdr().val_point(rr,th,phi)
      +2.*rr*g_ij(2,2).val_point(rr,th,phi),
      g_tt_dt = r2*g_ij(2,2).dsdt().val_point(rr,th,phi),
      g_tt_dp = sth*r2*g_ij(2,2).stdsdp().val_point(rr,th,phi),
      g_pp_Lorene  = g_ij(3,3).val_point(rr,th,phi),
      g_pp_dr = r2sinth2*g_ij(3,3).dsdr().val_point(rr,th,phi)
      +2.*rr*sinth2*g_pp_Lorene,
      g_pp_dt = r2sinth2*g_ij(3,3).dsdt().val_point(rr,th,phi)
      +2.*cth*sth*r2*g_pp_Lorene,
      g_pp_dp = sth*r2sinth2*g_ij(3,3).stdsdp().val_point(rr,th,phi),
      g_rt_dr = rr*g_ij(1,2).dsdr().val_point(rr,th,phi)
      +g_ij(1,2).val_point(rr,th,phi),
      g_rt_dt = rr*g_ij(1,2).dsdt().val_point(rr,th,phi),
      g_rt_dp = rsinth*g_ij(1,2).stdsdp().val_point(rr,th,phi),
      g_rp_Lorene = g_ij(1,3).val_point(rr,th,phi),
      g_rp_dr = rsinth*g_ij(1,3).dsdr().val_point(rr,th,phi)
      + sth*g_rp_Lorene,
      g_rp_dt = rsinth*g_ij(1,3).dsdt().val_point(rr,th,phi)
      + rr*cth*g_rp_Lorene,
      g_rp_dp = sth*rsinth*g_ij(1,3).stdsdp().val_point(rr,th,phi),
      g_tp_Lorene = g_ij(2,3).val_point(rr,th,phi), 
      g_tp_dr = rr*rsinth*g_ij(2,3).dsdr().val_point(rr,th,phi)
      + 2.*rsinth*g_tp_Lorene,
      g_tp_dt = rr*rsinth*g_ij(2,3).dsdt().val_point(rr,th,phi)
      + r2*cth*g_tp_Lorene,
      g_tp_dp = r2sinth2*g_ij(2,3).stdsdp().val_point(rr,th,phi);

    //cout << "gamma Lorene= " << g_ij(1,1).val_point(rr,th,phi) << " " << g_ij(2,2).val_point(rr,th,phi) << " " << g_ij(3,3).val_point(rr,th,phi) << " " << g_ij(1,2).val_point(rr,th,phi) << " " << g_ij(1,3).val_point(rr,th,phi) << " " << g_ij(2,3).val_point(rr,th,phi) << endl;
    //cout << "gamma= " << g_ij(1,1).val_point(rr,th,phi) << " " << r2*g_ij(2,2).val_point(rr,th,phi) << " " << r2sinth2*g_ij(3,3).val_point(rr,th,phi) << " " << rr*g_ij(1,2).val_point(rr,th,phi) << " " << rsinth*g_ij(1,3).val_point(rr,th,phi) << " " << rr*rsinth*g_ij(2,3).val_point(rr,th,phi) << endl;
    
    //INVERSE 3-METRIC
    //NB: these are gamma^ij, not g^ij
    const Sym_tensor& g_up_ij = *(gamcon_tab_[indice_time]);
    double grr=g_up_ij(1,1).val_point(rr,th,phi), 
      gtt=rm2*g_up_ij(2,2).val_point(rr,th,phi),
      gpp=rm2*sm2*g_up_ij(3,3).val_point(rr,th,phi),
      grt=rm1*g_up_ij(1,2).val_point(rr,th,phi),
      grp=rsm1*g_up_ij(1,3).val_point(rr,th,phi),
      gtp=rm1*rsm1*g_up_ij(2,3).val_point(rr,th,phi);
    //cout << "gamma= " << grr << " " << gtt << " " << gpp << " " << grt << " " << gtp << " " << grp << endl;
    
    //EXTRINSIC CURVATURE
    const Sym_tensor& kij = *(kij_tab_[indice_time]);
    double K_rr = kij(1,1).val_point(rr,th,phi),
      K_tt = r2*kij(2,2).val_point(rr,th,phi),
      K_pp = rsinth*rsinth*kij(3,3).val_point(rr,th,phi),
      K_rt = rr*kij(1,2).val_point(rr,th,phi),
      K_rp = rsinth*kij(1,3).val_point(rr,th,phi), 
      K_tp = rr*rsinth*kij(2,3).val_point(rr,th,phi);
    //cout << "Kij= " << K_rr << " " << K_tt << " " << K_pp << " " << K_rp << " " << K_rt << " " << K_tp << endl;
    
    //3-CHRISTOFFELS
    double Grrr = 0.5*grr*g_rr_dr+0.5*grt*(2.*g_rt_dr-g_rr_dt)
      + 0.5*grp*(2.*g_rp_dr-g_rr_dp),
      Grrt = 0.5*grr*g_rr_dt+0.5*grt*g_tt_dr
      +0.5*grp*(g_tp_dr+g_rp_dt-g_rt_dp),
      Grtt = 0.5*grr*(2.*g_rt_dt-g_tt_dr)
      +0.5*grt*g_tt_dt+0.5*grp*(2.*g_tp_dt-g_tt_dp),
      Grpp = 0.5*grr*(2.*g_rp_dp-g_pp_dr)+0.5*grt*(2.*g_tp_dp-g_pp_dt)
      +0.5*grp*g_pp_dp,
      Grrp = 0.5*grr*g_rr_dp+0.5*grt*(g_rt_dp+g_tp_dr-g_rp_dt)+0.5*grp*g_pp_dr,
      Grtp = 0.5*grr*(g_rt_dp + g_rp_dt - g_tp_dr) + 0.5*grt*g_tt_dp
      + 0.5*grp*g_pp_dt,
      Gtrr = 0.5*grt*g_rr_dr+0.5*gtt*(2.*g_rt_dr-g_rr_dt)
      + 0.5*gtp*(2.*g_rp_dr-g_rr_dp),
      Gtrt = 0.5*grt*g_rr_dt+0.5*gtt*g_tt_dr
      +0.5*gtp*(g_tp_dr+g_rp_dt-g_rt_dp),
      Gttt = 0.5*grt*(2.*g_rt_dt-g_tt_dr)
      +0.5*gtt*g_tt_dt+0.5*gtp*(2.*g_tp_dt-g_tt_dp),
      Gtpp = 0.5*grt*(2.*g_rp_dp-g_pp_dr)+0.5*gtt*(2.*g_tp_dp-g_pp_dt)
      +0.5*gtp*g_pp_dp,
      Gtrp = 0.5*grt*g_rr_dp+0.5*gtt*(g_rt_dp+g_tp_dr-g_rp_dt)+0.5*gtp*g_pp_dr,
      Gttp = 0.5*grt*(g_rt_dp + g_rp_dt - g_tp_dr) + 0.5*gtt*g_tt_dp
      + 0.5*gtp*g_pp_dt,
      Gprr = 0.5*grp*g_rr_dr+0.5*gtp*(2.*g_rt_dr-g_rr_dt)
      + 0.5*gpp*(2.*g_rp_dr-g_rr_dp),
      Gprt = 0.5*grp*g_rr_dt+0.5*gtp*g_tt_dr
      +0.5*gpp*(g_tp_dr+g_rp_dt-g_rt_dp),
      Gptt = 0.5*grp*(2.*g_rt_dt-g_tt_dr)
      +0.5*gtp*g_tt_dt+0.5*gpp*(2.*g_tp_dt-g_tt_dp),
      Gppp = 0.5*grp*(2.*g_rp_dp-g_pp_dr)+0.5*gtp*(2.*g_tp_dp-g_pp_dt)
      +0.5*gpp*g_pp_dp,
      Gprp = 0.5*grp*g_rr_dp+0.5*gtp*(g_rt_dp+g_tp_dr-g_rp_dt)+0.5*gpp*g_pp_dr,
      Gptp = 0.5*grp*(g_rt_dp + g_rp_dt - g_tp_dr) + 0.5*gtp*g_tt_dp
      + 0.5*gpp*g_pp_dt;
    
    // 3+1 GEODESIC EQUATION
    double VN_d = Vr*N_dr+Vth*N_dt+Vph*N_dp,
      KV2 = K_rr*Vr*Vr+K_tt*Vth*Vth+K_pp*Vph*Vph
      +2.*(K_rt*Vr*Vth+K_rp*Vr*Vph+K_tp*Vth*Vph),
      factor = NNm1*VN_d - KV2,
      KV_r = K_rr*Vr+K_rt*Vth+K_rp*Vph,
      KV_t = K_rt*Vr+K_tt*Vth+K_tp*Vph,
      KV_p = K_rp*Vr+K_tp*Vth+K_pp*Vph;
    
    res[0] = EE*(NN*KV2 - VN_d);
    res[1] = NN*Vr-betar;
    res[2] = NN*Vth-betat;
    res[3] = NN*Vph-betap;
    res[4] = NN*(Vr*factor
		 + 2.*grr*KV_r + 2.*grt*KV_t + 2.*grp*KV_p 
		 - Grrr*Vr*Vr-2.*Grrt*Vr*Vth-Grtt*Vth*Vth-Grpp*Vph*Vph
		 - 2.*Grrp*Vr*Vph - 2.*Grtp*Vth*Vph) 
      - grr*N_dr - grt*N_dt - grp*N_dp
      - Vr*betar_dr - Vth*betar_dt - Vph*betar_dp;
    res[5] = NN*(Vth*factor
		 +2.*grt*KV_r + 2.*gtt*KV_t + 2.*gtp*KV_p
		 - Gttt*Vth*Vth-Gtpp*Vph*Vph-Gtrr*Vr*Vr-2.*Gtrt*Vr*Vth
		 - 2.*Gtrp*Vr*Vph - 2.*Gttp*Vth*Vph) 
      - grt*N_dr - gtt*N_dt - gtp*N_dp
      - Vr*betat_dr - Vth*betat_dt - Vph*betat_dp;
    res[6] = NN*(Vph*factor
		 +2.*grp*KV_r+2.*gtp*KV_t + 2.*gpp*KV_p
		 - Gprr*Vr*Vr - 2.*Gprt*Vr*Vth - Gptt*Vth*Vth - Gppp*Vph*Vph
		 - 2*Vph*(Vr*Gprp+Vth*Gptp))
      - grp*N_dr - gtp*N_dt - gpp*N_dp
      - Vr*betap_dr - Vth*betap_dt - Vph*betap_dp;
  }else{
    // AXISYMMETRIC CIRCULAR SPACETIME

    //LAPSE
    Scalar* lapse = lapse_tab_[indice_time];
    double NN = lapse -> val_point(rr,th,phi), NNm1 = 1./NN,
      N_dr = lapse->dsdr().val_point(rr,th,phi), 
      N_dt = lapse->dsdt().val_point(rr,th,phi);
    
    //SHIFT
    const Vector& shift = *(shift_tab_[indice_time]);
    double betap = rsm1*shift(3).val_point(rr,th,phi),
      betap_dr = rsm1*shift(3).dsdr().val_point(rr,th,phi)
      -rm2*sm1*shift(3).val_point(rr,th,phi),
      betap_dt = rsm1*shift(3).dsdt().val_point(rr,th,phi)
      -cth*sm2*rm1*shift(3).val_point(rr,th,phi);
    
    //3-METRIC DERIVATIVES (for Christo)
    const Sym_tensor& g_ij = *(gamcov_tab_[indice_time]);
    double g_rr_dr = g_ij(1,1).dsdr().val_point(rr,th,phi),
      g_rr_dt = g_ij(1,1).dsdt().val_point(rr,th,phi),
      g_tt_dr = r2*g_ij(2,2).dsdr().val_point(rr,th,phi)
      +2.*rr*g_ij(2,2).val_point(rr,th,phi),
      g_tt_dt = r2*g_ij(2,2).dsdt().val_point(rr,th,phi),
      g_pp_Lorene  = g_ij(3,3).val_point(rr,th,phi),
      g_pp_dr = r2sinth2*g_ij(3,3).dsdr().val_point(rr,th,phi)
      +2.*rr*sinth2*g_pp_Lorene,
      g_pp_dt = r2sinth2*g_ij(3,3).dsdt().val_point(rr,th,phi)
      +2.*cth*sth*r2*g_pp_Lorene;
    
    //INVERSE 3-METRIC
    //NB: these are gamma^ij, not g^ij
    const Sym_tensor& g_up_ij = *(gamcon_tab_[indice_time]);
    double grr=g_up_ij(1,1).val_point(rr,th,phi), 
      gtt=rm2*g_up_ij(2,2).val_point(rr,th,phi),
      gpp=rm2*sm2*g_up_ij(3,3).val_point(rr,th,phi);
    
    //EXTRINSIC CURVATURE
    const Sym_tensor& kij = *(kij_tab_[indice_time]);
    double K_rp = rsinth*kij(1,3).val_point(rr,th,phi), 
      K_tp = rr*rsinth*kij(2,3).val_point(rr,th,phi);
    
    //3-CHRISTOFFELS
    double Grrr = 0.5*grr*g_rr_dr,
      Grrt = 0.5*grr*g_rr_dt,
      Grtt = 0.5*grr*(-g_tt_dr),
      Grpp = 0.5*grr*(-g_pp_dr),
      Gtrr = 0.5*gtt*(-g_rr_dt),
      Gtrt = 0.5*gtt*g_tt_dr,
      Gttt = 0.5*gtt*g_tt_dt,
      Gtpp = 0.5*gtt*(-g_pp_dt),
      Gprp = 0.5*gpp*g_pp_dr,
      Gptp = 0.5*gpp*g_pp_dt;

    double factor = NNm1*(Vr*N_dr+Vth*N_dt) - 2.*K_rp*Vr*Vph-2.*K_tp*Vth*Vph;
    
    res[0] = EE*(2.*NN*(K_rp*Vr*Vph+K_tp*Vth*Vph) - (Vr*N_dr+Vth*N_dt));
    res[1] = NN*Vr;
    res[2] = NN*Vth;
    res[3] = NN*Vph-betap;
    res[4] = NN*(Vr*factor+2.*grr*(K_rp*Vph) 
		 - Grrr*Vr*Vr-2.*Grrt*Vr*Vth-Grtt*Vth*Vth
		 -Grpp*Vph*Vph) - grr*N_dr;
    res[5] = NN*(Vth*factor+2.*gtt*(K_tp*Vph) 
		 - Gttt*Vth*Vth-Gtpp*Vph*Vph
		 -Gtrr*Vr*Vr-2.*Gtrt*Vr*Vth) - gtt*N_dt;
    res[6] = NN*(Vph*factor+2.*gpp*(K_rp*Vr+K_tp*Vth) 
		 - 2*Vph*(Vr*Gprp+Vth*Gptp)) - Vr*betap_dr-Vth*betap_dt;
  }

  for (int ii=0;ii<7;ii++){
    if (res[ii]!=res[ii]){
      cout << "i, res[i]= " << ii << " " << res[ii] << endl;
      GYOTO_ERROR("In NumericalMetricLorene::diff(): "
		 "derivatives are nan");
    }
    if (res[ii]==res[ii]+1.){
      cout << "i, res[i]= " << ii << " " << res[ii] << endl;
      GYOTO_ERROR("In NumericalMetricLorene::diff(): "
		 "derivatives are inf");
    }
  }
  
  /*
    time2 = clock();
    diftime = time2 - time1;
    cout << "Time elapsed in 3+1 diff (in sec)= " << setprecision(GYOTO_PREC) << setw(GYOTO_WIDTH) << diftime/clocks << endl;
  */

  /*
    Homogeneize equations by multiplying by c in Lorene units
    eg: dX^i/dt = N*V^i-beta^i -> dX^i/cdt = N*V^i-beta^i which is
    homogeneous (no units on both sides -- in Lorene units,
    V^i homogeneous to dX^i/cdt, see def of 4-vel with c!=1).
    Thus, just multiply res by c in Lorene units.

    Lorene units : time in 10^-3 s, length in 10^4 m
    So: c_lorene = c_SI*1e-3*1e-4

    NB: this is not necessary if time has been recast to Gyoto units
    in NumericalMetricLorene.C (for collapsing NS, for other metrics, no pb)
    so that t(Gyoto) = t(Lorene)*c_lorene.
    Then time is in Gyoto units (natural) and r is in Lorene
    units (10^4 m), and this mix of units is OK.
    
    for (int ii=0;ii<=6;ii++)
    res[ii]*=GYOTO_C*1e-3*1e-4;
  */

  return 0;
}

void NumericalMetricLorene::computeNBeta(const double coord[4],
					 double &NN,double beta[3]) const
{
  GYOTO_DEBUG << endl;
  double tt=coord[0], rr=coord[1],th=coord[2],rsinth = rr*sin(th),ph=coord[3];
  if (rr==0.) GYOTO_ERROR("In NumericalMetricLorene.C::computeNBeta r is 0!");
  if (rsinth==0.) GYOTO_ERROR("In NumericalMetricLorene.C::computeNBeta "
			     "on z axis!");
  double rm1 = 1./rr, rsm1 = 1./rsinth;

  int it=nb_times_-1;
  while(tt<times_[it] && it>=0){ //ASSUMES backward integration, to generalize
    it--;
  }

  int ind_nb=it;
  if (it==-1) ind_nb=0;

  Scalar* lapse = (lapse_tab_[ind_nb]);
  NN = lapse->val_point(rr,th,ph);
  const Vector& shift = *(shift_tab_[ind_nb]);
  double beta_r = shift(1).val_point(rr,th,ph), 
    beta_t = rm1*shift(2).val_point(rr,th,ph),
    beta_p = rsm1*shift(3).val_point(rr,th,ph);

  if ((it==nb_times_-2 && it!=-1) || (it==0 && it!=nb_times_-1)){ 
    // linear interpol
    //NB: the && conditions are made for the very special case of
    //nb_times_=1...

    //Lapse:
    Scalar* lapse1 = (lapse_tab_[it]);
    double NN1 = lapse1->val_point(rr,th,ph);
    Scalar* lapse2 = (lapse_tab_[it+1]);
    double NN2 = lapse2->val_point(rr,th,ph);
    double t1=times_[it], t2=times_[it+1];
    NN = (NN2-NN1)/(t2-t1)*(tt-t1)+NN1;
    //Shift:
    const Vector& shift1 = *(shift_tab_[it]);
    double beta_r1 = shift1(1).val_point(rr,th,ph), 
      beta_t1 = rm1*shift1(2).val_point(rr,th,ph),
      beta_p1 = rsm1*shift1(3).val_point(rr,th,ph);
    const Vector& shift2 = *(shift_tab_[it+1]);
    double beta_r2 = shift2(1).val_point(rr,th,ph), 
      beta_t2 = rm1*shift2(2).val_point(rr,th,ph),
      beta_p2 = rsm1*shift2(3).val_point(rr,th,ph);
    beta_r = (beta_r2-beta_r1)/(t2-t1)*(tt-t1)+beta_r1;
    beta_t = (beta_t2-beta_t1)/(t2-t1)*(tt-t1)+beta_t1;
    beta_p = (beta_p2-beta_p1)/(t2-t1)*(tt-t1)+beta_p1;
  }else if (it > 0 && it < nb_times_-2){ // 3rd order interpo
    // Lapse:
    Scalar* lapse1 = (lapse_tab_[it-1]);
    double NN1 = lapse1->val_point(rr,th,ph);
    Scalar* lapse2 = (lapse_tab_[it]);
    double NN2 = lapse2->val_point(rr,th,ph);
    Scalar* lapse3 = (lapse_tab_[it+1]);
    double NN3 = lapse3->val_point(rr,th,ph);
    Scalar* lapse4 = (lapse_tab_[it+2]);
    double NN4 = lapse4->val_point(rr,th,ph);
    double values[4]={NN1,NN2,NN3,NN4};
    NN = Interpol3rdOrder(tt,it,values);
    // Shift:
    const Vector& shift1 = *(shift_tab_[it-1]);
    const Scalar& betar1=shift1(1), betat1=shift1(2), betap1=shift1(3) ;
    double beta_r1 = shift1(1).val_point(rr,th,ph), 
      beta_t1 = rm1*shift1(2).val_point(rr,th,ph),
      beta_p1 = rsm1*shift1(3).val_point(rr,th,ph);
    const Vector& shift2 = *(shift_tab_[it]);
    double beta_r2 = shift2(1).val_point(rr,th,ph), 
      beta_t2 = rm1*shift2(2).val_point(rr,th,ph),
      beta_p2 = rsm1*shift2(3).val_point(rr,th,ph);
    const Vector& shift3 = *(shift_tab_[it+1]);
    double beta_r3 = shift3(1).val_point(rr,th,ph), 
      beta_t3 = rm1*shift3(2).val_point(rr,th,ph),
      beta_p3 = rsm1*shift3(3).val_point(rr,th,ph);
    const Vector& shift4 = *(shift_tab_[it+2]);
    double beta_r4 = shift4(1).val_point(rr,th,ph), 
      beta_t4 = rm1*shift4(2).val_point(rr,th,ph),
      beta_p4 = rsm1*shift4(3).val_point(rr,th,ph);
    double values_r[4]={beta_r1,beta_r2,beta_r3,beta_r4};
    double values_t[4]={beta_t1,beta_t2,beta_t3,beta_t4};
    double values_p[4]={beta_p1,beta_p2,beta_p3,beta_p4};
    beta_r = Interpol3rdOrder(tt,it,values_r);
    beta_t = Interpol3rdOrder(tt,it,values_t);
    beta_p = Interpol3rdOrder(tt,it,values_p);
  }

  beta[0]=beta_r;
  beta[1]=beta_t;
  beta[2]=beta_p;

}

void NumericalMetricLorene::jacobian(double jac[4][4][4],
				     const double x0[4]
				     ) const
{
  // Compute jac[alpha][mu][nu] = \partial_alpha(g_{mu,nu})]

  // Special case alpha=0, derivative wrt t, TO BE CODED

  // d/dr and d/dth
  double gmunudr[4][4], gmunudth[4][4];
  gmunu_di(x0,gmunudr,gmunudth);

  // alpha=1 (d/dr) redirect to gmunu_dr
  // alpha=1 (d/dtheta) redirect to gmunu_dth
  // alpha=3 (d/dphi) is zero
  for (int mu=0; mu<4; ++mu){
    for (int nu=0; nu<4; ++nu){
      jac[1][mu][nu]=gmunudr[mu][nu];
      jac[2][mu][nu]=gmunudth[mu][nu];
      jac[3][mu][nu]=0.;
    }
  }

}

void NumericalMetricLorene::gmunu_di(const double pos[4],
				     double gmunudr[4][4],
				     double gmunudth[4][4]
				     ) const
{
  GYOTO_DEBUG << endl;
  double tt=pos[0];
  int it=nb_times_-1;
  while(tt<times_[it] && it>=0){ //ASSUMES backward integration, to generalize
    it--;
  }

  double pos3[3]={pos[1],pos[2],pos[3]};
  if (it==nb_times_-1) {
    double gmunudrnbt[4][4], gmunudthnbt[4][4];
    gmunu_di(pos3,nb_times_-1,gmunudrnbt,gmunudthnbt);
    for (int mu=0;mu<4;++mu){
      for (int nu=0;nu<4;++nu){
	gmunudr[mu][nu]=gmunudrnbt[mu][nu];
	gmunudth[mu][nu]=gmunudthnbt[mu][nu];
	//use metric nb nb_times_-1
	//for all times > max(metric times)
      }
    }

  }
  if (it==-1) {
    double gmunudr0[4][4], gmunudth0[4][4];
    gmunu_di(pos3,0,gmunudr0,gmunudth0);
    for (int mu=0;mu<4;++mu){
      for (int nu=0;nu<4;++nu){
	gmunudr[mu][nu]=gmunudr0[mu][nu];
	gmunudth[mu][nu]=gmunudth0[mu][nu];
	//use metric nb 0
	//for all times > max(metric times)
      }
    }
  }
  if (it==nb_times_-2 || it==0){ //linear inter for extremal-1 points
    double t1=times_[it], t2=times_[it+1];
    double gmunudr1[4][4], gmunudr2[4][4],
      gmunudth1[4][4], gmunudth2[4][4];
    gmunu_di(pos3,it,gmunudr1,gmunudth1);
    gmunu_di(pos3,it+1,gmunudr2,gmunudth2);
    for (int mu=0;mu<4;++mu){
      for (int nu=0;nu<4;++nu){
	gmunudr[mu][nu]=(gmunudr1[mu][nu] - gmunudr2[mu][nu])/(t1-t2)*(tt-t1)
	  + gmunudr1[mu][nu];
	gmunudth[mu][nu]=(gmunudth1[mu][nu] - gmunudth2[mu][nu])/(t1-t2)*(tt-t1)
	  + gmunudth1[mu][nu];
      }
    }    
  }
  //Else : use 3rd order interp
  double gmunudr1[4][4], gmunudr2[4][4], gmunudr3[4][4], gmunudr4[4][4],
    gmunudth1[4][4], gmunudth2[4][4], gmunudth3[4][4], gmunudth4[4][4];
  gmunu_di(pos3,it-1,gmunudr1,gmunudth1);
  gmunu_di(pos3,it,gmunudr2,gmunudth2);
  gmunu_di(pos3,it+1,gmunudr3,gmunudth3);
  gmunu_di(pos3,it+2,gmunudr4,gmunudth4);

  for (int mu=0;mu<4;++mu){
    for (int nu=0;nu<4;++nu){
      double y1=gmunudr1[mu][nu],
	y2=gmunudr2[mu][nu],
	y3=gmunudr3[mu][nu],
	y4=gmunudr4[mu][nu];
      double values[4]={y1,y2,y3,y4};
      gmunudr[mu][nu] = Interpol3rdOrder(tt,it,values);

      y1=gmunudth1[mu][nu];
      y2=gmunudth2[mu][nu];
      y3=gmunudth3[mu][nu];
      y4=gmunudth4[mu][nu];
      double values2[4]={y1,y2,y3,y4};
      gmunudth[mu][nu] = Interpol3rdOrder(tt,it,values2);
    }
  }    
}

void NumericalMetricLorene::gmunu_di(const double pos[4],
				     int indice_time,
				     double gmunudr[4][4],
				     double gmunudth[4][4]
				     ) const
{
  // This provides gmunu_dr and gmunu_dth for all mu,nu
  if (indice_time<0 || indice_time>nb_times_-1) 
    GYOTO_ERROR("NumericalMetricLorene::gmunu_di: "
		"incoherent value of indice_time");

  double rr=pos[0], r2=rr*rr, th=pos[1], costh=cos(th), sinth=sin(th),
    sinth2=sinth*sinth, rsinth=rr*sin(th), ph=pos[2];
  Scalar* lapse = (lapse_tab_[indice_time]);
  double lapse_val = lapse->val_point(rr,th,ph), // lapse value
    lapsedr = lapse->dsdr().val_point(rr,th,ph), // d(lapse)/dr
    lapsedth = lapse->dsdt().val_point(rr,th,ph); // d(lapse)/dtheta

  const Vector& shift = *(shift_tab_[indice_time]);
  double betap = shift(3).val_point(rr,th,ph), // omega
    betapdr = shift(3).dsdr().val_point(rr,th,ph), // domega/dr
    betapdth = shift(3).dsdt().val_point(rr,th,ph); // domega/dtheta

  const Sym_tensor& g_ij = *(gamcov_tab_[indice_time]);
  double g_rr = g_ij(1,1).val_point(rr,th,ph), // gamma_rr
    g_rrdr = g_ij(1,1).dsdr().val_point(rr,th,ph), // d(gamma_rr)/dr
    g_rrdth = g_ij(1,1).dsdt().val_point(rr,th,ph), // d(gamma_rr)/dtheta
    g_thth = g_ij(2,2).val_point(rr,th,ph), // idem for gamma_thth
    g_ththdr = g_ij(2,2).dsdr().val_point(rr,th,ph),
    g_ththdth = g_ij(2,2).dsdt().val_point(rr,th,ph),
    g_pp = g_ij(3,3).val_point(rr,th,ph), // idem for gamma_pp
    g_ppdr = g_ij(3,3).dsdr().val_point(rr,th,ph),
    g_ppdth = g_ij(3,3).dsdt().val_point(rr,th,ph);

  // Okay here we have all possible metric quantities and derivatives

  // g_{t,mu}
  gmunudr[0][0] = -2.*lapsedr*lapse_val + 2.*(betapdr - betap/rr)*betap*g_pp
    + betap*betap*g_ppdr + 2.*betap*betap/rsinth*g_pp;
  gmunudr[1][0] = gmunudr[0][1] = 0.;
  gmunudr[2][0] = gmunudr[0][2] = 0.;
  gmunudr[0][3] = gmunudr[3][0] = (betapdr - betap/rr)*g_pp*rsinth
    + betap*g_ppdr*rsinth + 2.*betap*g_pp*sinth;

  gmunudth[0][0] = -2.*lapsedth*lapse_val
    + 2.*(-costh/sinth*betap + betapdth)*betap*g_pp
    + betap*betap*g_ppdth + 2.*betap*betap*g_pp*costh/sinth;
  gmunudth[1][0] = gmunudth[0][1] = 0.;
  gmunudth[2][0] = gmunudth[0][2] = 0.;
  gmunudth[0][3] = gmunudth[3][0] = (-costh/sinth*betap + betapdth)*g_pp*rsinth
    + betap*g_ppdth*rsinth + 2.*betap*g_pp*rr*costh;

  // g_{i,mu}
  
  gmunudr[1][1] = g_rrdr;
  gmunudr[1][0] = gmunudr[0][1] = 0.;
  gmunudr[1][2] = gmunudr[2][1] = 0.;
  gmunudr[1][3] = gmunudr[3][1] = 0.;

  gmunudth[1][1] = g_rrdth;
  gmunudth[1][0] = gmunudth[0][1] = 0.;
  gmunudth[1][2] = gmunudth[2][1] = 0.;
  gmunudth[1][3] = gmunudth[3][1] = 0.;

  gmunudr[2][2] = r2*g_ththdr+ 2.*rr*g_thth;
  gmunudr[2][0] = gmunudr[0][2] = 0.;
  gmunudr[2][3] = gmunudr[3][2] = 0.;

  gmunudth[2][2] = r2*g_ththdth;
  gmunudth[2][0] = gmunudth[0][2] = 0.;
  gmunudth[2][3] = gmunudth[3][2] = 0.;

  gmunudr[3][3] = (g_ppdr*rr + 2.*g_pp)*rr*sinth2;

  gmunudth[3][3] = r2*sinth*(g_ppdth*sinth + 2.*g_pp*costh);
      
}


void NumericalMetricLorene::gmunu_up(double gup[4][4], const double x[4]
				     ) const
{
  GYOTO_DEBUG << endl;
  double tt=x[0];
  int it=nb_times_-1;
  while(tt<times_[it] && it>=0){ //ASSUMES backward integration, to generalize
    it--;
  }

  double pos3[3]={x[1],x[2],x[3]};
  if (it==nb_times_-1) {
    double gupnbt[4][4];
    gmunu_up(gupnbt,pos3,nb_times_-1);
    for (int mu=0;mu<4;++mu){
      for (int nu=0;nu<4;++nu){
	gup[mu][nu]=gupnbt[mu][nu];
	//use metric nb nb_times_-1
	//for all times > max(metric times)
      }
    }

  }
  if (it==-1) {
    double gup0[4][4];
    gmunu_up(gup0,pos3,0);
    for (int mu=0;mu<4;++mu){
      for (int nu=0;nu<4;++nu){
	gup[mu][nu]=gup0[mu][nu];
	//use metric nb 0
	//for all times > max(metric times)
      }
    }

  }
  if (it==nb_times_-2 || it==0){ //linear inter for extremal-1 points
    double t1=times_[it], t2=times_[it+1];
    double gup1[4][4], gup2[4][4];
    gmunu_up(gup1,pos3,it);
    gmunu_up(gup2,pos3,it+1);
    for (int mu=0;mu<4;++mu){
      for (int nu=0;nu<4;++nu){
	gup[mu][nu]=(gup1[mu][nu] - gup2[mu][nu])/(t1-t2)*(tt-t1)
	  + gup1[mu][nu];
      }
    }    
  }
  //Else : use 3rd order interp
  double gup1[4][4], gup2[4][4], gup3[4][4], gup4[4][4];
  gmunu_up(gup1,pos3,it-1);
  gmunu_up(gup2,pos3,it);
  gmunu_up(gup3,pos3,it+1);
  gmunu_up(gup4,pos3,it+2);

  for (int mu=0;mu<4;++mu){
    for (int nu=0;nu<4;++nu){
      double y1=gup1[mu][nu],
	y2=gup2[mu][nu],
	y3=gup3[mu][nu],
	y4=gup4[mu][nu];
      double values[4]={y1,y2,y3,y4};
      gup[mu][nu] = Interpol3rdOrder(tt,it,values);
    }
  }    
}

void NumericalMetricLorene::gmunu_up(double gup[4][4], const double x[4],
				     int indice_time
				     ) const
{
  // Contravariant metric coefs
  if (indice_time<0 || indice_time>nb_times_-1) 
    GYOTO_ERROR("NumericalMetricLorene::gmunu_up: "
		"incoherent value of indice_time");

  double rr=x[0], r2=rr*rr, th=x[1], costh=cos(th), sinth=sin(th),
    sinth2=sinth*sinth, rsinth=rr*sin(th), ph=x[2];
  Scalar* lapse = (lapse_tab_[indice_time]);
  double lapse_val = lapse->val_point(rr,th,ph); // lapse value
  double lapsem2 = 1./(lapse_val*lapse_val);

  const Vector& shift = *(shift_tab_[indice_time]);
  double betap = shift(3).val_point(rr,th,ph); // omega

  const Sym_tensor& g_up_ij = *(gamcov_tab_[indice_time]);
  double g_up_rr = g_up_ij(1,1).val_point(rr,th,ph), // gamma^rr
    g_up_thth = g_up_ij(2,2).val_point(rr,th,ph), // idem for gamma^thth
    g_up_pp = g_up_ij(3,3).val_point(rr,th,ph); // idem for gamma^pp

  // Okay here we have all possible metric quantities

  // g^{t,mu}

  gup[0][0] = -lapsem2;
  gup[0][1] = gup[1][0] = 0.;
  gup[0][2] = gup[2][0] = 0.;
  gup[0][3] = gup[3][0] = betap/rsinth*lapsem2;
  // remember that betap is in Lorene orthonormal basis
  // and that betap=-omega

  // g^{i,mu}

  gup[1][1] = g_up_rr;
  gup[1][2] = gup[2][1] = 0.;
  gup[1][3] = gup[3][1] = 0.;

  gup[2][2] = g_up_thth*1./r2;
  gup[2][3] = gup[3][2] = 0.;

  gup[3][3] = 1./(r2*sinth2) * (-betap*betap*lapsem2 + g_up_pp);
  
}
  

double NumericalMetricLorene::gmunu(const double pos[4], 
				     int mu, 
				     int nu) const
{
  GYOTO_DEBUG << endl;
  double tt=pos[0];
  int it=nb_times_-1;
  while(tt<times_[it] && it>=0){ //ASSUMES backward integration, to generalize
    it--;
  }
  //  it=-1; //DEBUGIT
  //  if (pos[1]<0.187) it=1305; // TEST!!!! 
  double pos3[3]={pos[1],pos[2],pos[3]};
  if (it==nb_times_-1) return gmunu(pos3,nb_times_-1,mu,nu); 
                       //use metric nb nb_times_-1
                       //for all times > max(metric times)
  if (it==-1) return gmunu(pos3,0,mu,nu); //use metric nb 0 for all
					  //times < min(metric times)
  if (it==nb_times_-2 || it==0){ //linear inter for extremal-1 points
    double t1=times_[it], t2=times_[it+1];
    return (gmunu(pos3,it,mu,nu)
	    -gmunu(pos3,it+1,mu,nu))/(t1-t2)*(tt-t1) 
      + gmunu(pos3,it,mu,nu);
  }
  //Else : use 3rd order interp
  double y1=gmunu(pos3,it-1,mu,nu),
    y2=gmunu(pos3,it,mu,nu),
    y3=gmunu(pos3,it+1,mu,nu),
    y4=gmunu(pos3,it+2,mu,nu);
  double values[4]={y1,y2,y3,y4};
  return Interpol3rdOrder(tt,it,values);
}

double NumericalMetricLorene::gmunu(const double pos[3], 
				    int indice_time, int mu, int nu) const
{
  GYOTO_DEBUG << endl;
  /*
    4D metric.
    NB: 3-metric supposed to be conformally flat
   */
  if (indice_time<0 || indice_time>nb_times_-1) 
    GYOTO_ERROR("NumericalMetricLorene::gmunu: incoherent value of indice_time");
  
  if ( mu<0 || mu>3 || nu<0 || nu>3)
       GYOTO_ERROR("In NumericalMetricLorene::gmunu bad indice value");

  double rr=pos[0], r2=rr*rr, th=pos[1], rsinth=rr*sin(th);
  if (rr==0.) GYOTO_ERROR("In NumericalMetricLorene.C::gmunu r is 0!");
  if (rsinth==0.) GYOTO_ERROR("In NumericalMetricLorene.C::gmunu on z axis!");
  double rm1=1./rr, r2sinth2=r2*sin(th)*sin(th), ph=pos[2];

  Scalar* lapse = (lapse_tab_[indice_time]);

  double lapse_val = lapse->val_point(rr,th,ph);

  const Vector& shift = *(shift_tab_[indice_time]);

  const Sym_tensor& g_ij = *(gamcov_tab_[indice_time]);

  double res=0.;

  if ((mu==0) && (nu==0)) 
    {
      res = -lapse_val*lapse_val // NB: no correction factors due to
	                         // change of basis because they
	                         // cancel each other
	+g_ij(1,1).val_point(rr,th,ph)
	*shift(1).val_point(rr,th,ph)
	*shift(1).val_point(rr,th,ph)
	+g_ij(2,2).val_point(rr,th,ph)
	*shift(2).val_point(rr,th,ph)
	*shift(2).val_point(rr,th,ph)
	+g_ij(3,3).val_point(rr,th,ph)
	*shift(3).val_point(rr,th,ph)
	*shift(3).val_point(rr,th,ph);
    }else if ((mu==1) && (nu==1))
    {
      res = g_ij(1,1).val_point(rr,th,ph);
    }else if ((mu==2) && (nu==2))
    {
      res = r2*g_ij(2,2).val_point(rr,th,ph);
    }else if ((mu==3) && (nu==3))
    {
      res = r2sinth2*g_ij(3,3).val_point(rr,th,ph);

    }else if (((mu==0) && (nu==1)) || ((mu==1) && (nu==0)) )
    {
      res = g_ij(1,1).val_point(rr,th,ph)*shift(1).val_point(rr,th,ph);
    }else if (((mu==0) && (nu==2)) || ((mu==2) && (nu==0)) )
    {
      res = r2*g_ij(2,2).val_point(rr,th,ph)*rm1*shift(2).val_point(rr,th,ph);
    }else if (((mu==0) && (nu==3)) || ((mu==3) && (nu==0)) )
    {
      res = rsinth*g_ij(3,3).val_point(rr,th,ph)*shift(3).val_point(rr,th,ph);
    }
  if (res!=res) GYOTO_ERROR("NumericalMetricLorene::gmunu is nan!");
  if (res==res+1.) {
    //cout << "r,th,ph= " << rr << " " << th << " " << ph << endl;
    //cout << "mu,nu,res= " << mu << " " << nu << " " << res << endl;
    //GYOTO_ERROR("NumericalMetricLorene::gmunu is inf!");
  }
  return res; 
}

double NumericalMetricLorene::gmunu_up_dr(const double pos[4], 
				      int mu, 
				      int nu) const
{
  GYOTO_DEBUG << endl;
  double tt=pos[0];
  int it=nb_times_-1;
  while(tt<times_[it] && it>=0){ //ASSUMES backward integration, to generalize
    it--;
  }

  double pos3[3]={pos[1],pos[2],pos[3]};
  if (it==nb_times_-1) return gmunu_up_dr(pos3,nb_times_-1,mu,nu); 
                       //use metric nb nb_times_-1
                       //for all times > max(metric times)
  if (it==-1) return gmunu_up_dr(pos3,0,mu,nu); //use metric nb 0 for all
					  //times < min(metric times)
  if (it==nb_times_-2 || it==0){ //linear inter for extremal-1 points
    double t1=times_[it], t2=times_[it+1];
    return (gmunu_up_dr(pos3,it,mu,nu)
	    -gmunu_up_dr(pos3,it+1,mu,nu))/(t1-t2)*(tt-t1) 
      + gmunu_up_dr(pos3,it,mu,nu);
  }
  //Else : use 3rd order interp
  double y1=gmunu_up_dr(pos3,it-1,mu,nu),
    y2=gmunu_up_dr(pos3,it,mu,nu),
    y3=gmunu_up_dr(pos3,it+1,mu,nu),
    y4=gmunu_up_dr(pos3,it+2,mu,nu);
  double values[4]={y1,y2,y3,y4};
  return Interpol3rdOrder(tt,it,values);
}

double NumericalMetricLorene::gmunu_up_dr(const double pos[3], 
				      int indice_time, int mu, int nu) const
{
  GYOTO_DEBUG << endl;
  /*
    gmunu contravariant, derived wrt r.
    NB: 3-metric supposed to be conformally flat
   */
  if (indice_time<0 || indice_time>nb_times_-1) 
    GYOTO_ERROR("NumericalMetricLorene::gmunu_up_dr: "
	       "incoherent value of indice_time");
  
  if ( (mu!=0 && mu!=3) || (nu!=0 && nu!=3))
       GYOTO_ERROR("In NumericalMetricLorene::gmunu_up_dr bad indice value");
  // NB: so far only t and phi components are coded

  double rr=pos[0], th=pos[1], rsinth=rr*sin(th), ph=pos[2];
  if (rr==0.) GYOTO_ERROR("In NumericalMetricLorene.C::gmunu_up_dr r is 0!");
  if (rsinth==0.) GYOTO_ERROR("In NumericalMetricLorene.C::gmunu_up_dr "
			     "on z axis!");
  double rm1=1./rr, rsinthm1=1./(rsinth), rsinthm2=rsinthm1*rsinthm1;

  Scalar* lapse = (lapse_tab_[indice_time]);
  double lapse_val = lapse->val_point(rr,th,ph),
    lapse_valm1=1./lapse_val,
    lapse_valm2=lapse_valm1*lapse_valm1,
    lapse_valm3=lapse_valm2*lapse_valm1,
    lapsedr = lapse->dsdr().val_point(rr,th,ph);
 
  const Vector& shift = *(shift_tab_[indice_time]);
  double betap = shift(3).val_point(rr,th,ph),
    betapdr = shift(3).dsdr().val_point(rr,th,ph);

  const Sym_tensor& g_ij = *(gamcov_tab_[indice_time]);
  double g_pp = g_ij(3,3).val_point(rr,th,ph),
    g_ppdr = g_ij(3,3).dsdr().val_point(rr,th,ph);

  double res=0.;

  if ((mu==0) && (nu==0)) 
    {
      res = 2.*lapsedr*lapse_valm3;
    }else if ((mu==1) && (nu==1))
    {
      res = 0.;
    }else if ((mu==2) && (nu==2))
    {
      res = 0.;
    }else if ((mu==3) && (nu==3))
    {
      res = rsinthm2*(-2.*rm1*(1./g_pp-betap*betap*lapse_valm2)
		      -(g_ppdr/(g_pp*g_pp)
			+2.*betap*(betapdr*lapse_valm2
				   -betap*lapsedr*lapse_valm3)));

    }else if (((mu==0) && (nu==1)) || ((mu==1) && (nu==0)) )
    {
      res = 0.;
    }else if (((mu==0) && (nu==2)) || ((mu==2) && (nu==0)) )
    {
      res = 0.;
    }else if (((mu==0) && (nu==3)) || ((mu==3) && (nu==0)) )
    {
      res = lapse_valm2*rsinthm1*(-rm1*betap+betapdr-2.*betap*lapsedr*lapse_valm1);
    }
  if (res!=res) GYOTO_ERROR("NumericalMetricLorene::gmunu_up_dr is nan!");
  if (res==res+1.) GYOTO_ERROR("NumericalMetricLorene::gmunu_up_dr is inf!");
  return res; 
}

/////////

double NumericalMetricLorene::christoffel(const double coord[4],
					  const int alpha,
					  const int mu, const int nu) const
{
  // 4D christoffels: time interpolation
  GYOTO_DEBUG << endl;
  
  if (nb_times_>1) GYOTO_ERROR("In NML::christoffel:"
			      "so far only stationary metric implemented");

  double tt = coord[0];

  int it=nb_times_-1;
  while(tt<times_[it] && it>=0) it--;

  if (it==nb_times_-1) {
    return christoffel(coord,alpha,mu,nu,nb_times_-1); 
    //use metric nb nb_times_-1 for all times > max(metric times)
  }
  if (it==-1) {
    return christoffel(coord,alpha,mu,nu,0); 
    //use metric nb 0 for all times < min(metric times)
  }
  if (it==nb_times_-2 || it==0){ // LINEAR interp for extremal-1 points
    double t1=times_[it], t2=times_[it+1];
    double chris1 = christoffel(coord,alpha,mu,nu,it), 
      chris2 = christoffel(coord,alpha,mu,nu,it+1); 
    double chris = (chris2-chris1)/(t2-t1)*(tt-t1)+chris1;
    return chris;
  }
    
  //Else: use THIRD ORDER interp
  double chris1 = christoffel(coord,alpha,mu,nu,it-1), 
    chris2 = christoffel(coord,alpha,mu,nu,it), 
    chris3 = christoffel(coord,alpha,mu,nu,it+1), 
    chris4 = christoffel(coord,alpha,mu,nu,it+2);
  double values[4]={chris1,chris2,chris3,chris4};
  double chris = Interpol3rdOrder(tt,it,values);
  return chris;
}

double NumericalMetricLorene::christoffel(const double coord[4],
					  const int alpha,
					  const int mu, const int nu,
					  const int indice_time) const
{
  // 4D christoffels: actual computation on a given time slice

  /* 
                     !!! *** CAUTION *** !!!

     Here (and here only in this class) it is assumed that the metric
     is STATIONARY, AXISYM, and that the spacetime is CIRCULAR
     (typically, rotating relativistic stars framework).

     For other spacetimes, the 3+1 integration should be used (which
     is completely general), not the 4D one.
  */
  
  if (coord[1]==0. || sin(coord[2])==0.) GYOTO_ERROR("NML::christoffel:"
						    " bad location");

  double rr=coord[1], th=coord[2], ph=coord[3], sinth=0., costh=0.;
  sincos(th, &sinth, &costh);
  double r2=rr*rr, rsinth=rr*sinth, sm1=1./sinth, rm1=1./rr, rsm1 = rm1*sm1,
    sinth2=sinth*sinth, r2sinth2=r2*sinth2, rm2 = rm1*rm1;
  if ((alpha==0 && mu==0 && nu==1) || (alpha==0 && mu==1 && nu==0)){
    Scalar* lapse = lapse_tab_[indice_time];
    double NN = lapse -> val_point(rr,th,ph),
      Nr = lapse->dsdr().val_point(rr,th,ph);
    if (NN==0.) GYOTO_ERROR("In NML::christoffel: bad lapse value");
    const Vector& shift = *(shift_tab_[indice_time]);
    double beta_p = rsm1*shift(3).val_point(rr,th,ph);
    const Sym_tensor& kij = *(kij_tab_[indice_time]);
    double Krp = rsinth*kij(1,3).val_point(rr,th,ph);
    return 1./NN*(Nr-Krp*beta_p);
  }else if ((alpha==0 && mu==0 && nu==2) || (alpha==0 && mu==2 && nu==0)){
    Scalar* lapse = lapse_tab_[indice_time];
    double NN = lapse -> val_point(rr,th,ph),
      Nt = lapse->dsdt().val_point(rr,th,ph);
    if (NN==0.) GYOTO_ERROR("In NML::christoffel: bad lapse value");
    const Vector& shift = *(shift_tab_[indice_time]);
    double beta_p = rsm1*shift(3).val_point(rr,th,ph);
    const Sym_tensor& kij = *(kij_tab_[indice_time]);
    double Ktp = rr*rsinth*kij(2,3).val_point(rr,th,ph);
    return 1./NN*(Nt-Ktp*beta_p);
  } else if ((alpha==0 && mu==1 && nu==3) || (alpha==0 && mu==3 && nu==1)) {
    Scalar* lapse = lapse_tab_[indice_time];
    double NN = lapse -> val_point(rr,th,ph);
    if (NN==0.) GYOTO_ERROR("In NML::christoffel: bad laspe value");
    const Sym_tensor& kij = *(kij_tab_[indice_time]);
    double Krp = rsinth*kij(1,3).val_point(rr,th,ph);
    return -Krp/NN;
  } else if ((alpha==0 && mu==2 && nu==3) || (alpha==0 && mu==3 && nu==2)) {
    Scalar* lapse = lapse_tab_[indice_time];
    double NN = lapse -> val_point(rr,th,ph);
    if (NN==0.) GYOTO_ERROR("In NML::christoffel: bad lapse value");
    const Sym_tensor& kij = *(kij_tab_[indice_time]);
    double Ktp = rr*rsinth*kij(2,3).val_point(rr,th,ph);
    return -Ktp/NN;
  } else if (alpha==1 && mu==0 && nu==0) {
    Scalar* lapse = lapse_tab_[indice_time];
    double NN = lapse -> val_point(rr,th,ph),
      Nr = lapse->dsdr().val_point(rr,th,ph); 
    const Vector& shift = *(shift_tab_[indice_time]);
    double beta_p = rsm1*shift(3).val_point(rr,th,ph);
    const Sym_tensor& kij = *(kij_tab_[indice_time]);
    double Krp = rsinth*kij(1,3).val_point(rr,th,ph);
    const Sym_tensor& g_up_ij = *(gamcon_tab_[indice_time]);
    double grr=g_up_ij(1,1).val_point(rr,th,ph);
    return NN*grr*(Nr - 2.*Krp*beta_p);
  } else if (alpha==2 && mu==0 && nu==0) {
    Scalar* lapse = lapse_tab_[indice_time];
    double NN = lapse -> val_point(rr,th,ph),
      Nt = lapse->dsdt().val_point(rr,th,ph); 
    const Vector& shift = *(shift_tab_[indice_time]);
    double beta_p = rsm1*shift(3).val_point(rr,th,ph);
    const Sym_tensor& kij = *(kij_tab_[indice_time]);
    double Ktp = rr*rsinth*kij(1,3).val_point(rr,th,ph);
    const Sym_tensor& g_up_ij = *(gamcon_tab_[indice_time]);
    double gtt=1./r2*g_up_ij(1,1).val_point(rr,th,ph);
    return NN*gtt*(Nt - 2.*Ktp*beta_p);
  } else if ((alpha==1 && mu==0 && nu==3) || (alpha==1 && mu==3 && nu==0)) {
    Scalar* lapse = lapse_tab_[indice_time];
    double NN = lapse -> val_point(rr,th,ph);
    const Sym_tensor& kij = *(kij_tab_[indice_time]);
    double Krp = rsinth*kij(1,3).val_point(rr,th,ph);
    const Sym_tensor& g_up_ij = *(gamcon_tab_[indice_time]);
    double grr=g_up_ij(1,1).val_point(rr,th,ph);
    return -NN*grr*Krp;
  } else if ((alpha==2 && mu==0 && nu==3) || (alpha==2 && mu==3 && nu==0)) {
    Scalar* lapse = lapse_tab_[indice_time];
    double NN = lapse -> val_point(rr,th,ph);
    const Sym_tensor& kij = *(kij_tab_[indice_time]);
    double Ktp = rr*rsinth*kij(2,3).val_point(rr,th,ph);
    const Sym_tensor& g_up_ij = *(gamcon_tab_[indice_time]);
    double gtt=1./r2*g_up_ij(2,2).val_point(rr,th,ph);
    return -NN*gtt*Ktp;
  } else if ((alpha==3 && mu==1 && nu==0) || (alpha==3 && mu==0 && nu==1)) {
    Scalar* lapse = lapse_tab_[indice_time];
    double NN = lapse -> val_point(rr,th,ph),
      Nr = lapse->dsdr().val_point(rr,th,ph); 
    const Vector& shift = *(shift_tab_[indice_time]);
    double beta_p = rsm1*shift(3).val_point(rr,th,ph);
    double beta_pr = rsm1*shift(3).dsdr().val_point(rr,th,ph)
      -1./(r2*sinth)*beta_p;
    const Sym_tensor& kij = *(kij_tab_[indice_time]);
    double Krp = rsinth*kij(1,3).val_point(rr,th,ph);
    const Sym_tensor& g_up_ij = *(gamcon_tab_[indice_time]);
    double gpp=rsm1*rsm1*g_up_ij(3,3).val_point(rr,th,ph);
    const Sym_tensor& g_ij = *(gamcov_tab_[indice_time]) ;
    double gamma_prp = 0.5*gpp
      *(r2sinth2*g_ij(3,3).dsdr().val_point(rr,th,ph)
	+2.*rr*sinth2*g_ij(3,3).val_point(rr,th,ph));
    if (NN==0.) GYOTO_ERROR("In NML::christoffel: bad lapse value");
    return beta_pr + gamma_prp*beta_p
      -NN*gpp*Krp+beta_p/NN*(Krp*beta_p-Nr);
  } else if ((alpha==3 && mu==2 && nu==0) || (alpha==3 && mu==0 && nu==2)) {
    Scalar* lapse = lapse_tab_[indice_time];
    double NN = lapse -> val_point(rr,th,ph),
      Nt = lapse->dsdt().val_point(rr,th,ph); 
    const Vector& shift = *(shift_tab_[indice_time]);
    double beta_p = rsm1*shift(3).val_point(rr,th,ph);
    double beta_pt = rsm1*shift(3).dsdt().val_point(rr,th,ph)
      -costh*rsm1*sm1*beta_p;
    const Sym_tensor& kij = *(kij_tab_[indice_time]);
    double Ktp = rr*rsinth*kij(2,3).val_point(rr,th,ph);
    const Sym_tensor& g_up_ij = *(gamcon_tab_[indice_time]);
    double gpp=rsm1*rsm1*g_up_ij(3,3).val_point(rr,th,ph);
    const Sym_tensor& g_ij = *(gamcov_tab_[indice_time]) ;
    double gamma_ptp = 0.5*gpp
      *(r2sinth2*g_ij(3,3).dsdt().val_point(rr,th,ph)
	+2.*costh*sinth*r2*g_ij(3,3).val_point(rr,th,ph));
    if (NN==0.) GYOTO_ERROR("In NML::christoffel: bad lapse value");
    return beta_pt + gamma_ptp*beta_p
      -NN*gpp*Ktp+beta_p/NN*(Ktp*beta_p-Nt);
  } else if (alpha==1 && mu==1 && nu==1) {
    const Sym_tensor& g_up_ij = *(gamcon_tab_[indice_time]);
    const Sym_tensor& g_ij = *(gamcov_tab_[indice_time]) ;
    return 0.5*g_up_ij(1,1).val_point(rr,th,ph)
      *g_ij(1,1).dsdr().val_point(rr,th,ph);
  } else if (alpha==1 && mu==3 && nu==3) {
    const Sym_tensor& g_up_ij = *(gamcon_tab_[indice_time]);
    const Sym_tensor& g_ij = *(gamcov_tab_[indice_time]) ;
    return -0.5*g_up_ij(1,1).val_point(rr,th,ph)
      *(r2sinth2*g_ij(3,3).dsdr().val_point(rr,th,ph)
	+2.*rr*sinth2*g_ij(3,3).val_point(rr,th,ph));
  } else if ((alpha==1 && mu==1 && nu==2) || (alpha==1 && mu==2 && nu==1)) {
    const Sym_tensor& g_up_ij = *(gamcon_tab_[indice_time]);
    const Sym_tensor& g_ij = *(gamcov_tab_[indice_time]) ;
    return 0.5*g_up_ij(1,1).val_point(rr,th,ph)
      *g_ij(1,1).dsdt().val_point(rr,th,ph); 
  } else if (alpha==1 && mu==2 && nu==2) {
    const Sym_tensor& g_up_ij = *(gamcon_tab_[indice_time]);
    const Sym_tensor& g_ij = *(gamcov_tab_[indice_time]) ;
    return -0.5*g_up_ij(1,1).val_point(rr,th,ph)
      *(r2*g_ij(2,2).dsdr().val_point(rr,th,ph) 
	+ 2.*rr*g_ij(2,2).val_point(rr,th,ph));
  } else if (alpha==2 && mu==3 && nu==3) {
    const Sym_tensor& g_up_ij = *(gamcon_tab_[indice_time]);
    const Sym_tensor& g_ij = *(gamcov_tab_[indice_time]) ;
    return -0.5*rm2*g_up_ij(2,2).val_point(rr,th,ph)
      *(r2sinth2*g_ij(3,3).dsdt().val_point(rr,th,ph)
	+2.*costh*sinth*r2*g_ij(3,3).val_point(rr,th,ph));
  } else if (alpha==2 && mu==1 && nu==1) {
    const Sym_tensor& g_up_ij = *(gamcon_tab_[indice_time]);
    const Sym_tensor& g_ij = *(gamcov_tab_[indice_time]) ;
    return -0.5*rm2*g_up_ij(2,2).val_point(rr,th,ph)
      *g_ij(1,1).dsdt().val_point(rr,th,ph);
  } else if (alpha==2 && mu==2 && nu==2) {
    const Sym_tensor& g_up_ij = *(gamcon_tab_[indice_time]);
    const Sym_tensor& g_ij = *(gamcov_tab_[indice_time]) ;
    return 0.5*rm2*g_up_ij(2,2).val_point(rr,th,ph)
      *(r2*g_ij(2,2).dsdt().val_point(rr,th,ph));
  } else if ((alpha==2 && mu==1 && nu==2) || (alpha==2 && mu==2 && nu==1)) {
    const Sym_tensor& g_up_ij = *(gamcon_tab_[indice_time]);
    const Sym_tensor& g_ij = *(gamcov_tab_[indice_time]) ;
    return 0.5*rm2*g_up_ij(2,2).val_point(rr,th,ph)
      *(r2*g_ij(2,2).dsdr().val_point(rr,th,ph)
	+2.*rr*g_ij(2,2).val_point(rr,th,ph));
  } else if ((alpha==3 && mu==1 && nu==3) || (alpha==3 && mu==3 && nu==1)) {
    const Sym_tensor& g_up_ij = *(gamcon_tab_[indice_time]);
    const Sym_tensor& g_ij = *(gamcov_tab_[indice_time]) ;
    Scalar* lapse = lapse_tab_[indice_time];
    double NN = lapse -> val_point(rr,th,ph);
    if (NN==0.) GYOTO_ERROR("In NML::christoffel: bad lapse value");
    const Vector& shift = *(shift_tab_[indice_time]);
    double beta_p = rsm1*shift(3).val_point(rr,th,ph);
    const Sym_tensor& kij = *(kij_tab_[indice_time]);
    double Krp = rsinth*kij(1,3).val_point(rr,th,ph);
    return 0.5*rsm1*rsm1*g_up_ij(3,3).val_point(rr,th,ph)
      *(r2sinth2*g_ij(3,3).dsdr().val_point(rr,th,ph)
	+2.*rr*sinth2*g_ij(3,3).val_point(rr,th,ph)) + beta_p/NN*Krp;
  } else if ((alpha==3 && mu==2 && nu==3) || (alpha==3 && mu==3 && nu==2)) {
    const Sym_tensor& g_up_ij = *(gamcon_tab_[indice_time]);
    const Sym_tensor& g_ij = *(gamcov_tab_[indice_time]) ;
    Scalar* lapse = lapse_tab_[indice_time];
    double NN = lapse -> val_point(rr,th,ph);
    if (NN==0.) GYOTO_ERROR("In NML::christoffel: bad lapse value");
    const Vector& shift = *(shift_tab_[indice_time]);
    double beta_p = rsm1*shift(3).val_point(rr,th,ph);
    const Sym_tensor& kij = *(kij_tab_[indice_time]);
    double Ktp = rr*rsinth*kij(1,3).val_point(rr,th,ph);
    return 0.5*rsm1*rsm1*g_up_ij(3,3).val_point(rr,th,ph)
      *(r2sinth2*g_ij(3,3).dsdt().val_point(rr,th,ph)
	+2.*costh*sinth*r2*g_ij(3,3).val_point(rr,th,ph))
      + beta_p/NN*Ktp;
  }
  // Other christo are zero
  return 0.;
}

int NumericalMetricLorene::christoffel(double dst[4][4][4], 
				       const double coord[4]) const {
  // all at once computation of christoffel 4D: time interpolation
  GYOTO_DEBUG << endl;

  double tt = coord[0];

  if (nb_times_>1) GYOTO_ERROR("In NML::christoffel all at once:"
			      "so far only stationary metric implemented");

  int it=nb_times_-1;
  while(tt<times_[it] && it>=0) it--;

  if (it==nb_times_-1) {
    return christoffel(dst,coord,nb_times_-1); 
    //use metric nb nb_times_-1 for all times > max(metric times)
  }
  if (it==-1) {
    return christoffel(dst,coord,0); 
    //use metric nb 0 for all times < min(metric times)
  }
  if (it==nb_times_-2 || it==0){ // LINEAR interp for extremal-1 points
    double t1=times_[it], t2=times_[it+1];
    double dst1[4][4][4], dst2[4][4][4];
    if (christoffel(dst1,coord,it) || christoffel(dst2,coord,it+1)) return 1;
    int alpha, mu, nu;
    for (alpha=0; alpha<4; ++alpha) {
      for (mu=0; mu<4; ++mu) {
	double dst1c = dst1[alpha][mu][mu], dst2c = dst2[alpha][mu][mu];
	dst[alpha][mu][mu]=(dst2c-dst1c)/(t2-t1)*(tt-t1)+dst1c;
	for (nu=mu+1; nu<4; ++nu){
	  double dst1c = dst1[alpha][mu][nu], dst2c = dst2[alpha][mu][nu];
	  dst[alpha][mu][nu]=dst[alpha][nu][mu]=
	    (dst2c-dst1c)/(t2-t1)*(tt-t1)+dst1c;
	}
      }
    }
    return 0;
  }
    
  //Else: use THIRD ORDER interp
  double dst1[4][4][4], dst2[4][4][4], dst3[4][4][4], dst4[4][4][4];
  if (christoffel(dst1,coord,it-1) || christoffel(dst2,coord,it)
      || christoffel(dst3,coord,it+1) || christoffel(dst4,coord,it+2)) return 1;
  int alpha, mu, nu;
  for (alpha=0; alpha<4; ++alpha) {
    for (mu=0; mu<4; ++mu) {
      double values[4]={dst1[alpha][mu][mu],dst2[alpha][mu][mu],
			dst3[alpha][mu][mu],dst4[alpha][mu][mu]};
      dst[alpha][mu][mu]=Interpol3rdOrder(tt,it,values);
      for (nu=mu+1; nu<4; ++nu){
	double values[4]={dst1[alpha][mu][nu],dst2[alpha][mu][nu],
			  dst3[alpha][mu][nu],dst4[alpha][mu][nu]};
	dst[alpha][mu][nu]=dst[alpha][nu][mu]=Interpol3rdOrder(tt,it,values);
      }
    }
  }
  return 0;
}

int NumericalMetricLorene::christoffel(double dst[4][4][4], 
				       const double coord[4],
				       const int indice_time) const {
  // all at once computation of christoffel 4D: actual computation
  GYOTO_DEBUG << endl;
  double sinth=0., costh=0, rr=coord[1], th=coord[2], ph=coord[3];
  sincos(th, &sinth, &costh);
  if (rr==0. || sinth==0.) GYOTO_ERROR("NML::christoffel:"
						    " bad location");
  double
    r2=rr*rr, rsinth=rr*sinth, rm1=1./rr, sm1=1./sinth,
    rsm1 = rm1*sm1, sinth2=sinth*sinth, r2sinth2=r2*sinth2,
    rm2 = rm1*rm1;

    Scalar* lapse = lapse_tab_[indice_time];
    double NN = lapse -> val_point(rr,th,ph),
      Nr = lapse->dsdr().val_point(rr,th,ph),
      Nt = lapse->dsdt().val_point(rr,th,ph);
    if (NN==0.) GYOTO_ERROR("In NML::christoffel: bad laspe value");
    const Vector& shift = *(shift_tab_[indice_time]);
    double beta_p = rsm1*shift(3).val_point(rr,th,ph);
    double beta_pr = rsm1*shift(3).dsdr().val_point(rr,th,ph)
      -rm1*rsm1*shift(3).val_point(rr,th,ph);
    double beta_pt = rsm1*shift(3).dsdt().val_point(rr,th,ph)
      -costh*rsm1*sm1*shift(3).val_point(rr,th,ph);
    const Sym_tensor& kij = *(kij_tab_[indice_time]);
    double Krp = rsinth*kij(1,3).val_point(rr,th,ph);
    double Ktp = rr*rsinth*kij(2,3).val_point(rr,th,ph);
    const Sym_tensor& g_up_ij = *(gamcon_tab_[indice_time]);
    // contravariant 3-metric
    double grr=g_up_ij(1,1).val_point(rr,th,ph),
      gtt=rm2*g_up_ij(2,2).val_point(rr,th,ph),
      gpp=rsm1*rsm1*g_up_ij(3,3).val_point(rr,th,ph);
    const Sym_tensor& g_ij = *(gamcov_tab_[indice_time]) ;
    // derivation of covariant 3-metric
    double g_rr_r = g_ij(1,1).dsdr().val_point(rr,th,ph),
      g_rr_t = g_ij(1,1).dsdt().val_point(rr,th,ph),
      g_tt_r = r2*g_ij(2,2).dsdr().val_point(rr,th,ph) 
      + 2.*rr*g_ij(2,2).val_point(rr,th,ph),
      g_tt_t = r2*g_ij(2,2).dsdt().val_point(rr,th,ph),
      g_pp_r = r2sinth2*g_ij(3,3).dsdr().val_point(rr,th,ph)
      +2.*rr*sinth2*g_ij(3,3).val_point(rr,th,ph),
      g_pp_t = r2sinth2*g_ij(3,3).dsdt().val_point(rr,th,ph)
      +2.*costh*sinth*r2*g_ij(3,3).val_point(rr,th,ph);
    
    dst[0][0][1]=dst[0][1][0]=1./NN*(Nr-Krp*beta_p); //checked
    dst[0][0][2]=dst[0][2][0]=1./NN*(Nt-Ktp*beta_p); //checked
    dst[0][1][3]=dst[0][3][1]=-Krp/NN; //checked
    dst[0][2][3]=dst[0][3][2]=-Ktp/NN; //checked
    dst[1][0][0]=NN*grr*(Nr - 2.*Krp*beta_p - beta_p*beta_p/(2.*NN)*g_pp_r); //checked
    dst[2][0][0]=NN*gtt*(Nt - 2.*Ktp*beta_p - beta_p*beta_p/(2.*NN)*g_pp_t); //checked
    dst[1][0][3]=dst[1][3][0]=-grr*(NN*Krp+0.5*beta_p*g_pp_r); //checked
    dst[2][0][3]=dst[2][3][0]=-gtt*(NN*Ktp+0.5*beta_p*g_pp_t); //checked
    dst[3][1][0]=dst[3][0][1]=beta_pr + 0.5*gpp*g_pp_r*beta_p
    -NN*gpp*Krp+beta_p/NN*(Krp*beta_p-Nr); //checked
    dst[3][2][0]=dst[3][0][2]=beta_pt + 0.5*gpp*g_pp_t*beta_p
    -NN*gpp*Ktp+beta_p/NN*(Ktp*beta_p-Nt); //checked
    dst[1][1][1]=0.5*grr*g_rr_r; //checked
    dst[1][3][3]=-0.5*grr*g_pp_r; //checked
    dst[1][1][2]=dst[1][2][1]=0.5*grr*g_rr_t; //checked
    dst[1][2][2]=-0.5*grr*g_tt_r; //checked
    dst[2][3][3]=-0.5*gtt*g_pp_t; //checked
    dst[2][1][1]=-0.5*gtt*g_rr_t; //checked
    dst[2][2][2]=0.5*gtt*g_tt_t; //checked
    dst[2][1][2]=dst[2][2][1]=0.5*gtt*g_tt_r; //checked
    dst[3][1][3]=dst[3][3][1]=0.5*gpp*g_pp_r + beta_p/NN*Krp; //checked
    dst[3][2][3]=dst[3][3][2]=0.5*gpp*g_pp_t + beta_p/NN*Ktp; //checked
    dst[0][0][0]=0.;
    dst[0][0][3]=0.;
    dst[0][3][0]=0.;
    dst[0][1][1]=0.;
    dst[0][2][2]=0.;
    dst[0][3][3]=0.;
    dst[0][1][2]=0.;
    dst[0][2][1]=0.;
    dst[3][0][0]=0.;
    dst[1][0][1]=0.;
    dst[1][1][0]=0.;
    dst[1][0][2]=0.;
    dst[1][2][0]=0.;
    dst[2][0][1]=0.;
    dst[2][1][0]=0.;
    dst[2][0][2]=0.;
    dst[2][2][0]=0.;
    dst[3][0][3]=0.;
    dst[3][3][0]=0.;
    dst[1][1][3]=0.;
    dst[1][3][1]=0.;
    dst[1][2][3]=0.;
    dst[1][3][2]=0.;
    dst[2][1][3]=0.;
    dst[2][3][1]=0.;
    dst[2][2][3]=0.;
    dst[2][3][2]=0.;
    dst[3][1][1]=0.;
    dst[3][2][2]=0.;
    dst[3][3][3]=0.;
    dst[3][1][2]=0.;
    dst[3][2][1]=0.;
    // Oh yeah, this makes 64 christoffels
    
    return 0;
}

double NumericalMetricLorene::computeHorizon(const double* pos) const{
  GYOTO_DEBUG << endl;
  if (!hor_tab_ && !horizon_)
    return 0.;

  if (horizon_ && !hor_tab_)
    return horizon_;

  if (hor_tab_ && !horizon_){
    int it=nb_times_-1;
    double tt=pos[0];
    double* times=getTimes();

    while(tt<times[it] && it>=0){ //ASSUMES backward integration, to generalize
      it--;
    }

    double rhor;
    if (it==nb_times_-1){
      return computeHorizon(pos,nb_times_-1);
    }

    if (it==-1)
      return computeHorizon(pos,0);

    if (it==nb_times_-2 || it==0){ 
      double t1=times[it], t2=times[it+1];
      double rhor1=computeHorizon(pos,it),
	rhor2=computeHorizon(pos,it+1);
      rhor=(rhor2-rhor1)/(t2-t1)*(tt-t1)+rhor1;
      return rhor;
    }

    double rhor1=computeHorizon(pos,it-1),
      rhor2=computeHorizon(pos,it),
      rhor3=computeHorizon(pos,it+1),
      rhor4=computeHorizon(pos,it+2);
    double values[4]={rhor1,rhor2,rhor3,rhor4};
    rhor=Interpol3rdOrder(tt,it,values);
    return rhor;
  }

  GYOTO_ERROR("In NumericalMetricLorene::computeHorizon: "
	     "impossible case");
  return 0.;
  
}

double NumericalMetricLorene::computeHorizon(const double* pos, 
					     int indice_time) const{  
  GYOTO_DEBUG << endl;
  if (indice_time<0 || indice_time>nb_times_-1){
    GYOTO_ERROR("NumericalMetricLorene::computeHorizon"
	       ": incoherent value of indice_time");
  }
  
  double th=pos[2], phi=pos[3];
  Valeur* horizon = (hor_tab_[indice_time]);
  horizon->std_base_scal();
  return horizon->val_point(0,0.,th,phi);
}

double NumericalMetricLorene::Interpol3rdOrder(double tt, 
					       int indice_time, 
					       double values[4]) const {
  GYOTO_DEBUG << endl;
  //Interpolation at order 3 at point tt, the considered function 
  //taking the values "values" at time indices: indice_time-1,+0,+1,+2.
  double t1=times_[indice_time-1], t2=times_[indice_time],
    t3=times_[indice_time+1], t4=times_[indice_time+2];
  double y1=values[0], y2=values[1],
    y3=values[2], y4=values[3];

  //Neuville's algorithm
  //3 first order poly
  double P12=((tt-t2)*y1+(t1-tt)*y2)/(t1-t2),
    P23=((tt-t3)*y2+(t2-tt)*y3)/(t2-t3),
    P34=((tt-t4)*y3+(t3-tt)*y4)/(t3-t4);
  //2 second order poly
  double P123=((tt-t3)*P12+(t1-tt)*P23)/(t1-t3),
    P234=((tt-t4)*P23+(t2-tt)*P34)/(t2-t4);
  //1 third order poly : the solution
  double P1234=((tt-t4)*P123+(t1-tt)*P234)/(t1-t4);
  
  return P1234;
}

void NumericalMetricLorene::directory(std::string const &dir) {
  char const * const cdir=dir.c_str();
  filename_ = new char[strlen(cdir)+1];
  strcpy(filename_, cdir);
  setMetricSource();
}

std::string NumericalMetricLorene::directory() const {
  return filename_?string(filename_):string("");
}

bool NumericalMetricLorene::hasSurface() const {return  has_surface_;}
void NumericalMetricLorene::hasSurface(bool s) {
  has_surface_ = s;
  if (filename_!=NULL){
    GYOTO_ERROR("In NumericalMetricLorene::hasSurface "
	       "please provide Surface information before File in XML");
  }
}

bool NumericalMetricLorene::hasAccelerationVector() const {return  has_acceleration_vector_;}
void NumericalMetricLorene::hasAccelerationVector(bool aa) {
  has_acceleration_vector_ = aa;
  if (filename_!=NULL){
    GYOTO_ERROR("In NumericalMetricLorene::hasAccelerationVector "
	       "please provide Acceleration vector info before File in XML");
  }
}

bool NumericalMetricLorene::bosonstarcircular() const {
  return bosonstarcircular_;}
void NumericalMetricLorene::bosonstarcircular(bool t) {bosonstarcircular_=t;}

bool NumericalMetricLorene::specifyMarginalOrbits() const {
  return specify_marginalorbits_;
}
void NumericalMetricLorene::specifyMarginalOrbits(bool s) {
  specify_marginalorbits_=s;
  if (filename_!=NULL){
    GYOTO_ERROR("In NumericalMetricLorene::specifyMarginalOrbits "
	       "please provide Marginal orbits information "
	       "before File in XML");
  }
}

double NumericalMetricLorene::rico() const {return rico_;}
void NumericalMetricLorene::rico(double r0) {rico_=r0;}

bool NumericalMetricLorene::mapEt() const {return  mapet_;}
void NumericalMetricLorene::mapEt(bool s) {
  mapet_ = s;
  if (filename_!=NULL){
    GYOTO_ERROR("In NumericalMetricLorene::mapEt "
	       "please provide MapET/MapAF information before File in XML");
  }
}

bool NumericalMetricLorene::axisymCirc() const {return  axisymCirc_;}
void NumericalMetricLorene::axisymCirc(bool s) {
  axisymCirc_ = s;
}

double NumericalMetricLorene::initialTime() const {return initial_time_;}
void NumericalMetricLorene::initialTime(double t0) {initial_time_=t0;}

double NumericalMetricLorene::horizon() const {return horizon_;}
void NumericalMetricLorene::horizon(double r0) {horizon_=r0;}

void NumericalMetricLorene::circularVelocity(double const * coord, 
					     double* vel,
					     double dir) const {
  GYOTO_DEBUG << endl;
  //  return Generic::circularVelocity(coord,vel,dir); // TEST!!

  double tt = coord[0];

  int it=nb_times_-1;
  while(tt<times_[it] && it>=0) it--;

  if (it==nb_times_-1) {
    return circularVelocity(coord,vel,dir,nb_times_-1); 
  }
  if (it==-1) {
    return circularVelocity(coord,vel,dir,0);
  }

  // Linear interp for extremal-1 points
  if (it==nb_times_-2 || it==0){ 
    double vel1[4], vel2[4];
    double t1=times_[it], t2=times_[it+1];
    circularVelocity(coord,vel1,dir,it);
    circularVelocity(coord,vel2,dir,it+1);
    for (int ii=0;ii<4;ii++) vel[ii] = 
			       (vel2[ii]-vel1[ii])/(t2-t1)*(tt-t1) 
			       + vel1[ii];
    return;
  }

  //Else: use 3rd order interp
  double vel1[4], vel2[4], vel3[4], vel4[4];
  circularVelocity(coord,vel1,dir,it-1);
  circularVelocity(coord,vel2,dir,it);
  circularVelocity(coord,vel3,dir,it+1);
  circularVelocity(coord,vel4,dir,it+2);
  double values[4];
  for (int ii=0;ii<4;ii++) {
    values[0]=vel1[ii];
    values[1]=vel2[ii];
    values[2]=vel3[ii];
    values[3]=vel4[ii];
    vel[ii] = Interpol3rdOrder(tt,it,values);
  }
  return;

}

void NumericalMetricLorene::circularVelocity(double const * coor, 
					     double* vel,
					     double dir, 
					     int indice_time) const {
  //cout << "IN CIRCULAR" << endl;
  if (bosonstarcircular_){
    // This expression is related to the ZAMO 3-velocity derived
    // in Grandclement+14 boson star paper
    //for (int ii=1;ii<=5000;ii++){
    //double coorbis[4]={0.,double(ii*0.01),M_PI/2.,0.};
    double rr=coor[1], th=coor[2], sinth=sin(th), ph=coor[3];
    //double rr=coorbis[1], th=coorbis[2], sinth=sin(th), ph=coorbis[3];
    if (rr<=0. || sinth==0.) GYOTO_ERROR("In NML::circularv: bad coor");
    double rsm1 = 1./(rr*sinth), rm2 = 1/(rr*rr), sm1 = 1./sinth;
    const Sym_tensor& g_ij = *(gamcov_tab_[indice_time]) ;
    double B2 = g_ij(3,3).val_point(rr,th,ph);
    if (B2<=0.) GYOTO_ERROR("In NML::circularv: bad B2");
    double BB = sqrt(B2);
    double Br = g_ij(3,3).dsdr().val_point(rr,th,ph)/(2.*BB);
    const Vector& shift = *(shift_tab_[indice_time]);
    double beta_p = rsm1*shift(3).val_point(rr,th,ph),
      beta_p_r = rsm1*shift(3).dsdr().val_point(rr,th,ph)
      -rm2*sm1*shift(3).val_point(rr,th,ph);
    Scalar* lapse = lapse_tab_[indice_time];
    double NN = lapse -> val_point(rr,th,ph);
    if (NN==0.) GYOTO_ERROR("In NML::circularv: bad N");
    double Nr = lapse->dsdr().val_point(rr,th,ph);
    double DD = B2*rr*rr/(NN*NN)*beta_p_r*beta_p_r
      + 4.*Nr/NN*(Br/BB+1./rr);
    if (DD<0.) GYOTO_ERROR("In NML::circularv: bad D");
    //    double g_tt = gmunu(coor,0,0), g_tp = gmunu(coor,0,3);
    double g_pp = gmunu(coor,3,3);
    //double g_tt = gmunu(coorbis,0,0), g_tp = gmunu(coorbis,0,3),
    //g_pp = gmunu(coorbis,3,3);
    if (g_pp<=0.) GYOTO_ERROR("In NML::circularv: bad g_pp");
    double Vzamo = 0.5*(-BB*rr/NN*beta_p_r+sqrt(DD))/(1./rr+Br/BB);
    double Omega = NN*Vzamo/sqrt(g_pp) - beta_p;
    double ut = 1./(NN*sqrt(1.-Vzamo*Vzamo));
    vel[0] = ut; vel[1] = 0.; vel[2] = 0.; vel[3] = Omega*ut;
    
    //    double ell=2.5;
    //    double pot = 0.5*log((g_tp*g_tp-g_tt*g_pp)
    //			 /(g_tt*ell*ell+2.*ell*g_tp+g_pp));
    //cout << rr << " " << g_tp*g_tp-g_tt*g_pp << " " << g_tt*ell*ell+2.*ell*g_tp+g_pp << " " << pot << endl;
    //}
    double normtol = 1e-6;
    double norm = ScalarProd(coor,vel,vel);
    if (fabs(norm+1.)>normtol) {
      cerr << "At rr=" << coor[1] << endl;
      GYOTO_ERROR("In NML::circularv: bad norm");
    }
    return;
  }
  
  GYOTO_ERROR("In NML::circularVelocity: circular velocity not implemented"
	     " for this particular metric");
}
