
// Lorene headers
#include "metric.h"
#include "nbr_spx.h"
#include "utilitaires.h"
#include "graphique.h"
#include "cmp.h"

//Gyoto headers
#include "GyotoUtils.h"
#include "GyotoNumericalMetricLorene.h"
#include "GyotoError.h"
#include "GyotoFactoryMessenger.h"

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

NumericalMetricLorene::NumericalMetricLorene() :
  WIP("Metric::NumericalMetricLorene"),
  Generic(GYOTO_COORDKIND_SPHERICAL, "NumericalMetricLorene"),
  mapkind_(NULL),
  initial_time_(0.),
  has_surface_(0),
  refine_(0),
  r_refine_(0.),
  h0_refine_(0.),
  filename_(NULL),
  lapse_tab_(NULL),
  shift_tab_(NULL),
  gamcov_tab_(NULL),
  gamcon_tab_(NULL),
  kij_tab_(NULL),
  times_(NULL),
  nssurf_tab_(NULL),
  vsurf_tab_(NULL),
  lorentz_tab_(NULL),
  hor_tab_(NULL),
  horizon_(0.)
{
  GYOTO_DEBUG << endl;
}

NumericalMetricLorene::NumericalMetricLorene(const NumericalMetricLorene&o) :
  Generic(GYOTO_COORDKIND_SPHERICAL,"NumericalMetricLorene"),
  mapkind_(NULL),
  initial_time_(o.initial_time_),
  has_surface_(o.has_surface_),
  refine_(o.refine_),
  r_refine_(o.r_refine_),
  h0_refine_(o.h0_refine_),
  filename_(NULL),
  lapse_tab_(NULL),
  shift_tab_(NULL),
  gamcov_tab_(NULL),
  gamcon_tab_(NULL),
  kij_tab_(NULL),
  times_(NULL),
  nssurf_tab_(NULL),
  vsurf_tab_(NULL),
  lorentz_tab_(NULL),
  hor_tab_(NULL),
  horizon_(o.horizon_)
{
  GYOTO_DEBUG << endl;
  NumericalMetricLorene::setMetricSource();
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
  if (lorentz_tab_) { delete [] lorentz_tab_;      lorentz_tab_=NULL;     }
  if (hor_tab_) { delete [] hor_tab_;      hor_tab_=NULL;     }
}

void NumericalMetricLorene::setMetricSource() {
  GYOTO_DEBUG << endl;
  DIR *dp;
  struct dirent *dirp;
  if((dp  = opendir(filename_)) == NULL) {
    throwError("In NumericalMetricLorene.C constructor : bad filename_");
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
    throwError("In NumericalMetricLorene.C: bad nb_times_ value");

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
      throwError("NumericalMetricLorene.C: Problem opening file!");
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
    if (mapkind_=="MapEt"){
      //Map_et* map = new Map_et(*grid, resu) ;
      map = new Map_et(*grid, resu) ;
    }else if (mapkind_=="MapAf"){
      //      Map_af* map = new Map_af(*grid, resu) ;
      map = new Map_af(*grid, resu) ;
    }else{
      cout << "MapKind= " << mapkind_ << endl;
      throwError("In NumericalMetricLorene: bad mapping kind.");
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
      Mg3d* grid_ah = new Mg3d(resu) ;
      Valeur* horizon = new Valeur(*grid_ah, resu) ;
      hor_tab_[i-1] = horizon;
    }

    fclose(resu) ;
  }

  if (debug()) cout << "NumericalMetricLorene.C constructor: "
		 "geometrical quantities initialized." << endl;

}

char const *  NumericalMetricLorene::getFileName() const 
{
  GYOTO_DEBUG << endl; 
  return filename_; }

Sym_tensor** NumericalMetricLorene::getGamcon_tab() const {
  GYOTO_DEBUG << endl;
  return gamcon_tab_;}
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
Scalar** NumericalMetricLorene::getLorentz_tab() const {
  GYOTO_DEBUG << endl;
  return lorentz_tab_;}
Valeur** NumericalMetricLorene::getHor_tab() const {
  GYOTO_DEBUG << endl;
  return hor_tab_;}

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
    metric at hand

    ASSUMED: AXISYMMETRY
   */

  /*
    clock_t time1, time2;
    double diftime, clocks = CLOCKS_PER_SEC;
    time1 = clock();
  */

  if (indice_time<0 || indice_time>nb_times_-1) {
    throwError("NumericalMetricLorene::diff: incoherent value of indice_time");
  }

  //NB: here t=theta, not time!
  double EE=y[0], rr=y[1], th=y[2], phi=y[3], sth = sin(th), rsinth = rr*sth,
    r2sinth2 = rsinth*rsinth, sinth2 = sth*sth;
  if (rr==0.) throwError("In NumericalMetricLorene.C::diff r is 0!");
  if (rsinth==0.) throwError("In NumericalMetricLorene.C::diff on z axis!");
  double rm1 = 1./rr, rm2 = rm1*rm1, r2 = rr*rr, rsm1 = 1./rsinth, 
    sm1 = 1./sth, sm2 = sm1*sm1, Vr=y[4], Vth=y[5], Vph=y[6];

  /*
    Important remark!  Lorene uses the orthonormal spherical tetrad
    d/dr, 1/r d/dtheta, 1/(r*sin(theta)) d/dphi, and not the natural
    basis d/dr, d/dtheta, d/dphi There's thus a change of basis to
    come back to the natural basis of spherical coordinates.  This is
    why there are plenty of factors r2, 1/r2, 1/rsinth2, ... everywhere
   */
  //clock_t t1=clock();

  double myth=th; // TEST
  //LAPSE
  Scalar* lapse = lapse_tab_[indice_time];
  //double NN = lapse -> val_point(rr,th,phi), NNm1 = 1./NN,
  double NN = lapse -> val_point(rr,myth,phi), NNm1 = 1./NN,
    Nr = lapse->dsdr().val_point(rr,th,phi), 
    Nt = lapse->dsdt().val_point(rr,th,phi);
  //SHIFT
  const Vector& shift = *(shift_tab_[indice_time]);
  //double beta_r = shift(1).val_point(rr,th,phi),
  double beta_r = shift(1).val_point(rr,myth,phi),
    beta_t = rm1*shift(2).val_point(rr,th,phi),
    beta_r_r = shift(1).dsdr().val_point(rr,th,phi),
    beta_r_t = shift(1).dsdt().val_point(rr,th,phi),
    beta_t_r = rm1*shift(2).dsdr().val_point(rr,th,phi)
    -rm2*shift(2).val_point(rr,th,phi),
    beta_t_t = rm1*shift(2).dsdt().val_point(rr,th,phi),
    beta_p = rsm1*shift(3).val_point(rr,th,phi),
    beta_p_r = rsm1*shift(3).dsdr().val_point(rr,th,phi)
    -rm2*sm1*shift(3).val_point(rr,th,phi),
    beta_p_t = rsm1*shift(3).dsdt().val_point(rr,th,phi)
    -cos(th)*sm2*rm1*shift(3).val_point(rr,th,phi);

  //cout << "betar= " << beta_r << endl;

  //3-METRIC (assumed diagonal)
  const Sym_tensor& g_ij = *(gamcov_tab_[indice_time]);
  double g_rt=g_ij(1,2).val_point(rr,th,phi),
    g_rp=g_ij(1,3).val_point(rr,th,phi),
    g_tp=g_ij(2,3).val_point(rr,th,phi);
  if (g_rt!=0. || g_rp!=0. || g_tp!=0.){
    throwError("In NumericalMetricLorene.C: 3-metric should be diagonal");
  }
  //cout << "r= " << rr << endl;
  //cout << "metric= " << NN << " " << beta_r << " " << beta_t << " " << beta_p << endl;
  //cout << "betar= " << beta_r << " " << NN << " " << Vr << " " <<  NN*Vr-beta_r << endl;
  //cout << rr << " " << NN*Vr-beta_r << endl;
  //cout << g_ij(1,1).val_point(rr,th,phi) << " " << g_ij(2,2).val_point(rr,th,phi) << " " << g_ij(3,3).val_point(rr,th,phi) << endl;

  // Metrics for christoffels
  double g_rrr = g_ij(1,1).dsdr().val_point(rr,th,phi),
    g_rrt = g_ij(1,1).dsdt().val_point(rr,th,phi),
    g_tt  = g_ij(2,2).val_point(rr,th,phi),
    g_ttr = r2*g_ij(2,2).dsdr().val_point(rr,th,phi)+2.*rr*g_tt,
    g_ttt = r2*g_ij(2,2).dsdt().val_point(rr,th,phi),
    g_pp  = g_ij(3,3).val_point(rr,th,phi),
    g_ppr = r2sinth2*g_ij(3,3).dsdr().val_point(rr,th,phi)+2.*rr*sinth2*g_pp,
    g_ppt = r2sinth2*g_ij(3,3).dsdt().val_point(rr,th,phi)
    +2.*cos(th)*sin(th)*r2*g_pp;
  
  //INVERSE 3-METRIC
  const Sym_tensor& g_up_ij = *(gamcon_tab_[indice_time]);
  //  double grr=g_up_ij(1,1).val_point(rr,th,phi), 
  double grr=g_up_ij(1,1).val_point(rr,myth,phi), 
    gtt=rm2*g_up_ij(2,2).val_point(rr,th,phi),
    gpp=rm2*sm2*g_up_ij(3,3).val_point(rr,th,phi); //NB: these are gamma^ij, not g^ij

  //rr=0.19;th=M_PI/2.;phi=0.;
  //cout << setprecision(10) << "coef met= " << " " << lapse -> val_point(rr,th,phi) << " " << shift(1).val_point(rr,th,phi) << " " << g_up_ij(1,1).val_point(rr,th,phi) << " " << shift(3).val_point(rr,th,phi) << endl;

  //EXTRINSIC CURVATURE
  const Sym_tensor& kij = *(kij_tab_[indice_time]);
  double Krr = kij(1,1).val_point(rr,th,phi),
    Ktt = r2*kij(2,2).val_point(rr,th,phi),
    Kpp = rsinth*rsinth*kij(3,3).val_point(rr,th,phi),
    Krt = rr*kij(1,2).val_point(rr,th,phi),
    Krp = rsinth*kij(1,3).val_point(rr,th,phi), 
    Ktp = rr*rsinth*kij(2,3).val_point(rr,th,phi);
  //  cout << Krr << " " << Ktt << " " << Kpp << " " << Krt << " " << Krp << " " << Ktp << endl;
  //3-CHRISTOFFELS
  double Grrr = 0.5*grr*g_rrr,
    Grrt = 0.5*grr*g_rrt,
    Grtt = -0.5*grr*g_ttr,
    Grpp = -0.5*grr*g_ppr,
    Gttt = 0.5*gtt*g_ttt,
    Gtpp = -0.5*gtt*g_ppt,
    Gtrr = -0.5*gtt*g_rrt,
    Gtrt = 0.5*gtt*g_ttr,
    Gprp = 0.5*gpp*g_ppr,
    Gptp = 0.5*gpp*g_ppt;

  double factor = NNm1*(Vr*Nr+Vth*Nt)-Krr*Vr*Vr-Ktt*Vth*Vth-Kpp*Vph*Vph 
    - 2.*Krt*Vr*Vth-2.*Krp*Vr*Vph-2.*Ktp*Vth*Vph;

  res[0] = EE*NN*(Krr*Vr*Vr+Ktt*Vth*Vth+Kpp*Vph*Vph
		  +2.*(Krt*Vr*Vth+Krp*Vr*Vph+Ktp*Vth*Vph)) 
    - EE*(Vr*Nr+Vth*Nt);
  res[1] = NN*Vr-beta_r;
  res[2] = NN*Vth-beta_t;
  res[3] = NN*Vph-beta_p;
  res[4] = NN*(Vr*factor+2.*grr*(Krr*Vr+Krt*Vth+Krp*Vph) 
	       - Grrr*Vr*Vr-2.*Grrt*Vr*Vth-Grtt*Vth*Vth-Grpp*Vph*Vph) 
    - grr*Nr-Vr*beta_r_r-Vth*beta_r_t;
  res[5] = NN*(Vth*factor+2.*gtt*(Krt*Vr+Ktt*Vth+Ktp*Vph) 
	       - Gttt*Vth*Vth-Gtpp*Vph*Vph-Gtrr*Vr*Vr-2.*Gtrt*Vr*Vth) 
    - gtt*Nt-Vr*beta_t_r-Vth*beta_t_t;
  res[6] = NN*(Vph*factor+2.*gpp*(Krp*Vr+Ktp*Vth+Kpp*Vph) 
	       - 2*Vph*(Vr*Gprp+Vth*Gptp)) 
    - Vr*beta_p_r-Vth*beta_p_t;

  //  cout << beta_p << endl;
  for (int ii=0;ii<7;ii++){
    if (res[ii]!=res[ii]){
      cout << "i, res[i]= " << ii << " " << res[ii] << endl;
      throwError("In NumericalMetricLorene::diff(): "
		 "derivatives are nan");
    }
    if (res[ii]==res[ii]+1.){
      cout << "i, res[i]= " << ii << " " << res[ii] << endl;
      throwError("In NumericalMetricLorene::diff(): "
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

int NumericalMetricLorene::myrk4(double tt, const double coorin[7], 
				 double h, double res[7]) const
{
  GYOTO_DEBUG << endl;

  /*
    3+1 RK4, WITH ENERGY INTEG. Called by the 3+1 RK4_ada.
   */

  double k1[7], k2[7], k3[7], k4[7], coor_plus_halfk1[7], sixth_k1[7], 
    coor_plus_halfk2[7], third_k2[7], coor_plus_k3[7], 
    third_k3[7], sixth_k4[7]; 

  if(diff(tt,coorin,k1)) return 1;
  /*  cout << "res diff= " ;
  for (int ii=0;ii<7;ii++) cout << k1[ii] << " ";
  cout << endl;*/
	  
  for (int i=0;i<7;i++) {
    k1[i]=h*k1[i];
    coor_plus_halfk1[i]=coorin[i]+0.5*k1[i];
    reverseR(tt+0.5*h,coor_plus_halfk1);
    sixth_k1[i]=1./6.*k1[i];
  }

  if(diff(tt+h/2.,coor_plus_halfk1,k2)) return 1;

  for (int i=0;i<7;i++) {
    k2[i]=h*k2[i];
    coor_plus_halfk2[i]=coorin[i]+0.5*k2[i];
    reverseR(tt+0.5*h,coor_plus_halfk2);
    third_k2[i]=1./3.*k2[i];
  }

  if(diff(tt+h/2.,coor_plus_halfk2,k3)) return 1;

  for (int i=0;i<7;i++) {
    k3[i]=h*k3[i];
    coor_plus_k3[i]=coorin[i]+k3[i];
    reverseR(tt+h,coor_plus_k3);
    third_k3[i]=1./3.*k3[i];
  }

  if(diff(tt+h,coor_plus_k3,k4)) return 1;

  for (int i=0;i<7;i++) {
    k4[i]=h*k4[i];
    sixth_k4[i]=1./6.*k4[i];
  }

  for (int i=0;i<7;i++) {
    res[i]=coorin[i]+sixth_k1[i]+third_k2[i]+third_k3[i]+sixth_k4[i];
  }
  reverseR(tt+h,res);
  
  return 0;
}

//Non adaptive Runge Kutta (called by WorldlineIntegState if fixed step)
int NumericalMetricLorene::myrk4(Worldline * line, const double coord[8], 
				 double h, double res[8]) const{
  GYOTO_DEBUG << endl;
  double tt=coord[0], rr=coord[1], r2=rr*rr, rm1=1./rr, 
    th=coord[2],rsinth = rr*sin(th), rsm1=1./rsinth, ph=coord[3],
    tdot=coord[4], rdot=coord[5],thdot=coord[6],phdot=coord[7];

  //Check p_phi conservation:
  double const * const cst = line -> getCst();
  double cst_p_ph = cst[1];
  double pphi_err=fabs(cst_p_ph-(gmunu(coord,0,3)*tdot + gmunu(coord,3,3)
				   *phdot))/fabs(cst_p_ph)*100.; 
  double pphi_err_tol=1.;//in percent
  /*
    Next test done on pphi/tdot as pphi propto
    tdot and tdot can diverge close to horizon.
    Just check that pphi/tdot is close to zero.
   */
  if (pphi_err/fabs(tdot)>pphi_err_tol){
    GYOTO_SEVERE << "tdot: " << fabs(tdot) << endl;
    if (verbose() >= GYOTO_SEVERE_VERBOSITY){
      cerr << "***WARNING: in NumericalMetricLorene::myrk4: p_phi is drifting"
	" - with error p_phi,x1,x2,x3= " << pphi_err << " %, at " <<
	rr << " " << th << " " << ph << endl;
    }
  }

  //Check p_t conservation if stationary
  if (nb_times_==1){
    double cst_p_t = cst[2];
    double pt_err=fabs(cst_p_t-(gmunu(coord,0,0)*tdot + gmunu(coord,0,3)
				*phdot))/fabs(cst_p_t)*100.; 
    double pt_err_tol=1.;//in percent
    if (pt_err>pt_err_tol){
      if (verbose() >= GYOTO_SEVERE_VERBOSITY){
	cout << "***WARNING: in NumericalMetricLorene::myrk4: p_t is drifting"
	  " - with error p_t,x1,x2,x3= " << pt_err << " %, at " <<
	  rr << " " << th << " " << ph << endl;
      }
    }
  }
  
  if (tdot==0.) throwError("In NumericalMetricLorene.C::myrk4_ada tdot is 0!");
  double rprime=rdot/tdot, thprime=thdot/tdot,phprime=phdot/tdot;
  if (rr==0.) throwError("In NumericalMetricLorene.C::myrk4_ada r is 0!");
  if (rsinth==0.) 
    throwError("In NumericalMetricLorene.C::myrk4_ada on z axis!");

  // Lapse and shift at tt:
  double NN, beta[3];
  computeNBeta(coord,NN,beta);
  double beta_r=beta[0], beta_t=beta[1], beta_p=beta[2];

  double Vr = 1./NN*(rprime+beta_r), Vth = 1./NN*(thprime+beta_t),
    Vph = 1./NN*(phprime+beta_p);
  // Photon's energy as measured by Eulerian observer:
  double EE = tdot*NN;

  double coor[7]={EE,rr,th,ph,Vr,Vth,Vph};
  double coornew[7];

  double tdot_used=tdot;//, tdot_bef=tdot;

  if (myrk4(tt,coor,h,coornew))
    return 1;

  double tnow = tt+h, rnow=coornew[1], thnow=coornew[2], phnow=coornew[3];
  double posend[4]={tnow,rnow,thnow,phnow};
  //Lapse and shift at tnow:
  computeNBeta(posend,NN,beta);
  beta_r=beta[0];
  beta_t=beta[1];
  beta_p=beta[2];

  rprime=NN*coornew[4]-beta_r;
  thprime=NN*coornew[5]-beta_t;
  phprime=NN*coornew[6]-beta_p;

  double EEend = coornew[0];

  tdot_used = EEend/NN;

  if (tdot_used<0.) 
    GYOTO_SEVERE << "In NumericalMetricLorene.C: WARNING TDOT IS <0" << endl;
  
  rdot=rprime*tdot_used;
  thdot=thprime*tdot_used;
  phdot=phprime*tdot_used;
  res[0]=tnow;
  res[1]=coornew[1];
  res[2]=coornew[2];
  res[3]=coornew[3];
  res[4]=tdot_used;
  res[5]=rdot;
  res[6]=thdot;
  res[7]=phdot;
  
  return 0;
}

//3+1 SPATIAL ADAPTIVE RK4 WITH ENERGY INTEGRATION
int NumericalMetricLorene::myrk4_adaptive(double tt, const double coord[7], 
					  double, double normref, 
					  double coordnew[7], 
					  const double cst[2], 
					  double& tdot_used, 
					  double h0, double& h1, 
					  double& hused,
					  double h1max) const{
  GYOTO_DEBUG << endl;
  /*
    3+1 RK4_ada, for internal use only 
    (not called directly by WorldlineIntegState).
    coord=[E,r,th,ph,Vr,Vth,Vph]
    This function is called by the 4D RK4_ada below.
  */
  
  double delta0[7], dcoord[7];
  
  //Standard RK4 parameters:
  double delta0min=1e-15;
  double S=0.9;
  double errmin=1e-6;
  double sigh1=1.;

  h1max=deltaMax(coord, h1max);

  /* ***Parameter to fine tuned (mainly for spacetimes with horizon) *** */
  /*
    eps:     defines the quality of RK4 integration
             it must be low enough to ensure proper integration
	     but high enough to prevent from too long integration

    stepmin: minimum allowed RK4 step, below which integration is stopped
             important parameter for spacetimes with horizon
	     must be low enough to prevent stopping before hitting
	     and high enough to prevent from too long integration

    rhormax: is above the largest value of apparent horizon radius
             for all envolved metrics
	     the stepmin condition is only read for r<rhormax

    NB: don't forget that the Coconut resolutions (time res, spectral res)
        are also crucial parameters for tracing horizon properly
   */
  double eps=0.005; //"std" value is 0.0001
  double stepmin=1e-8;//1e-6;
  double rhormax=0.2;
  /* ****************************************************************** */
 
  if (diff(tt,coord,dcoord)) return 1;

  for (int i = 0;i<7;i++) {
    delta0[i]=delta0min+eps*(fabs(h0*dcoord[i]));
  }

  double hbis=0.5*h0;
  double coordhalf[7];
  double coord2[7];
  double delta1[7];
  
  double err;
  int count=0;
  int zaxis=0;
  double thetatol=1e-5; //launch z-axis pb if theta is below

  /* WHILE LOOP FOR STEP DETERMINATION */
  while (1){

    count++;
    if (count > 100){
      throwError("NumericalMetricLorene: too many iterations in RK4");
    }
    err=0.;

    int step1=myrk4(tt,coord,h0,coordnew);
    int step21=myrk4(tt,coord,hbis,coordhalf);
    int step22=myrk4(tt+hbis,coordhalf,hbis,coord2);
    
    while (step1 || step21 || step22) {
      /*
	If here, then NumColStar::myrk4 returned 1, which means that one step
	of integration computed above leads to sub-app. horizon
	location.
	Divide integration step by 10.
	Stop condition if step less than stepmin value.
       */
      h0/=10.;hbis/=10.; 
      //Update delta0 with new h0
      for (int i = 0;i<7;i++) {
	delta0[i]=delta0min+eps*(fabs(h0*dcoord[i]));
      }
      if (debug()){
	cout << "Step divided to " << h0 << endl;
      }
      if (fabs(h0)<stepmin) {
	//if (debug()){
	  cout
	    << "Stop condition: at t,r= " 
	    << tt << " " << coord[1] 
	    << ", due to too small integration step" 
	    << " after dividing step: too close to horizon." << endl;
	  //}
	return 1;
      }
      step1=myrk4(tt,coord,h0,coordnew);
      step21=myrk4(tt,coord,hbis,coordhalf);
      step22=myrk4(tt+hbis,coordhalf,hbis,coord2);
    } // End step computation     
    
 
    /* *** Special zaxis treatment  *** */
    if (fabs(fmod(fabs(coordnew[2])+M_PI/2, M_PI)-M_PI/2) < thetatol){ 
      //checks whether theta is close to 0[pi]
      //launch z-axis special treatment if yes
      zaxis=1;
      h0*=1.1; hbis*=1.1;
      if (myrk4(tt,coord,h0,coordnew) 
	  || myrk4(tt,coord,hbis,coordhalf)
	  || myrk4(tt+hbis,coordhalf,hbis,coord2)) return 1;
# if GYOTO_DEBUG_ENABLED
      GYOTO_DEBUG << endl << "!!!!NOTE!!!!: Passing close to z-axis at theta= "
		  << coord[2] << " and r= " << coord[1]
		  << ", jumping ahead with h0= " << h0 << endl;
# endif
    } //End zaxis

    /* *** Error determination  *** */
    int ierr=10;
    if (!zaxis){
      for (int i = 0;i<7;i++){
	delta1[i]=coord2[i]-coordnew[i];
	if (err<fabs(delta1[i]/delta0[i])) {
	  err=fabs(delta1[i]/delta0[i]);
	  ierr=i;
	}
      }
    }else{ /* z-axis pb case : forget phi coordinate, which is a
	      function of 1/sin(theta) Indeed phidotdot is a function
	      of Gprp and Gptp, Christoffels that are themselves
	      function of 1/sin(theta) terms.  NB: this special
	      treatment is used very rarely (only for i=N/2) 
	      so it's not a pb for the overall RK4
	    */
      for (int i = 0;i<3;i++){
	delta1[i]=coord2[i]-coordnew[i];
	if (err<fabs(delta1[i]/delta0[i])) {
	  err=fabs(delta1[i]/delta0[i]);
	  ierr=i;
	}
      }
      for (int i = 4;i<7;i++){
	delta1[i]=coord2[i]-coordnew[i];
	if (err<fabs(delta1[i]/delta0[i])) {
	  err=fabs(delta1[i]/delta0[i]);
	  ierr=i;
	}
      }
    } //End of err determination

    /* *** Final treatment depending on err *** */
    if (err>1) { 
      h0=S*h0*pow(err,-0.25);
      hbis=0.5*h0;
    }else{
      double rr=coord[1];
      if (fabs(h0)<stepmin && rr<rhormax) {
	/*
	  Stop condition if used step smaller than stepmin,
	  and r less than ~ max value of app. horizon.
	  This means that the geodesic is "accumulating"
	  close to an horizon.
	 */
	GYOTO_DEBUG
	  << "Stop condition: at t,r= " 
	  << tt << " " << rr
	  << ", due to too small integration step. "
	  " Too close to horizon." << endl;
	return 1;
      }
      //pour éviter les explosions:
      h1=(err > errmin ? S*h0*pow(err,-0.2) : 4.*h0);
      if (h1<0.) sigh1=-1.;//why sigh1 and fabs(h1)? because otherwise
			   //if h1<0 (possible here if backwards
			   //integration), h1 is < h1min, so h1 is
			   //always set to h1min...
      if (fabs(h1)<delta_min_) h1=sigh1*delta_min_;
      if (fabs(h1)>h1max) h1=sigh1*h1max;
      hused=h0;
      break;
    }
 
  }
  /* END OF WHILE LOOP FOR STEP DETERMINATION */
  
  return 0;
} //End spatial ada RK with E integ

//Main adaptive RK4
int NumericalMetricLorene::myrk4_adaptive(Worldline* line, 
					  const double coord[8], 
					  double lastnorm, double normref, 
					  double coordnew[8], double h0, 
					  double& h1,
					  double h1max) const
{
  GYOTO_DEBUG << endl;
  double tt=coord[0], rr=coord[1],th=coord[2],rsinth = rr*sin(th),ph=coord[3],
    tdot=coord[4], rdot=coord[5],thdot=coord[6],phdot=coord[7];
  
  //Check p_phi conservation:
  double const * const cst = line -> getCst();
  double cst_p_ph = cst[1];
  double pphi_err=fabs(cst_p_ph-(gmunu(coord,0,3)*tdot + gmunu(coord,3,3)
				 *phdot))/fabs(cst_p_ph)*100.; 
  double pphi_err_tol=1.;//in percent
  /*
    Next test done on pphi/tdot as pphi propto
    tdot and tdot can diverge close to horizon.
    Just check that pphi/tdot is close to zero.
  */
  if (pphi_err/fabs(tdot)>pphi_err_tol){
    if (verbose() >= GYOTO_SEVERE_VERBOSITY){
      cerr << "***WARNING: in NumericalMetricLorene::myrk4_adaptive:"
	" p_phi is drifting"
	" - with error p_phi,x1,x2,x3= " << pphi_err << " %, at " <<
	rr << " " << th << " " << ph << endl;
    }
  }
  
  //Check p_t conservation if stationary
  if (nb_times_==1){
    double cst_p_t = cst[2];
    double pt_err=fabs(cst_p_t-(gmunu(coord,0,0)*tdot + gmunu(coord,0,3)
				*phdot))/fabs(cst_p_t)*100.; 
    double pt_err_tol=1.;//in percent
    if (pt_err>pt_err_tol){
      if (verbose() >= GYOTO_SEVERE_VERBOSITY){
	cerr << "***WARNING: in NumericalMetricLorene::myrk4: p_t is drifting"
	  " - with error p_t,x1,x2,x3= " << pt_err << " %, at " <<
	  rr << " " << th << " " << ph << endl;
      }
    }
  }
  
  if (tdot==0.) throwError("In NumericalMetricLorene.C::myrk4_ada tdot is 0!");
  double rprime=rdot/tdot, thprime=thdot/tdot,phprime=phdot/tdot;
  if (rr==0.) throwError("In NumericalMetricLorene.C::myrk4_ada r is 0!");
  if (rsinth==0.) 
    throwError("In NumericalMetricLorene.C::myrk4_ada on z axis!");
  double rm1 = 1./rr, rsm1 = 1./rsinth;
  
  // Lapse and shift at tt:
  double NN, beta[3];
  computeNBeta(coord,NN,beta);
  double beta_r=beta[0], beta_t=beta[1], beta_p=beta[2];
  
  double Vr = 1./NN*(rprime+beta_r), Vth = 1./NN*(thprime+beta_t),
    Vph = 1./NN*(phprime+beta_p);
  // Photon's energy as measured by Eulerian observer:
  double EE = tdot*NN;
  double coor[7]={EE,rr,th,ph,Vr,Vth,Vph};
  double coornew[7];
  double hused=1000.;
  
  if (tdot<0. && h0>0.) h0*=-1.;//to integrate backwards if tdot<0
  
  double tdot_used=tdot;//, tdot_bef=tdot;
  
  //to debug:
  int it=nb_times_-1;
  while(tt<times_[it] && it>=0){ //ASSUMES backward integration, to generalize
    it--;
  }

  if (refine_){
    /*
      Refined integration:
      if asked in the XML, the integration step is
      imposed below the value h0_refine_.
     */
    double h0tmp=h0;
    if (rr<r_refine_ && fabs(h0)>fabs(h0_refine_)) h0=h0_refine_;
    if (h0*h0tmp<0.) h0*=-1;
  }

  if (myrk4_adaptive(tt,coor,lastnorm,normref,coornew,
		     cst,tdot_used,h0,h1,hused, h1max)) {
    return 1;
  }
  
  double tnow = tt+hused, rnow=coornew[1], thnow=coornew[2], phnow=coornew[3];
  
  rm1 = 1./rnow;rsm1 = 1./(rnow*sin(thnow));
  double posend[4]={tnow,rnow,thnow,phnow};
  
  //Lapse and shift at tnow:
  computeNBeta(posend,NN,beta);
  beta_r=beta[0];
  beta_t=beta[1];
  beta_p=beta[2];
  
  rprime=NN*coornew[4]-beta_r;
  thprime=NN*coornew[5]-beta_t;
  phprime=NN*coornew[6]-beta_p;
  
  double EEend = coornew[0];

  tdot_used = EEend/NN;
  if (tdot_used<0.) 
    GYOTO_SEVERE << "In NumericalMetricLorene.C: WARNING TDOT IS <0" << endl;
  
  rdot=rprime*tdot_used;
  thdot=thprime*tdot_used;
  phdot=phprime*tdot_used;

  coordnew[0]=tnow;
  coordnew[1]=coornew[1];
  coordnew[2]=coornew[2];
  coordnew[3]=coornew[3];
  coordnew[4]=tdot_used;
  coordnew[5]=rdot;
  coordnew[6]=thdot;
  coordnew[7]=phdot;

  return 0;
}

void NumericalMetricLorene::reverseR(double tt, double coord[7]) const{
  GYOTO_DEBUG << endl;
  if (coord[1]<0.) {
    double rhor=computeHorizon(coord);
    if (rhor==0.){
      // then no horizon, r=0 is allowed, e.g. boson star case
      // r vector -> - itself

      //reverse spatial vector
      coord[1]*=-1.;
      coord[2]=M_PI-coord[2];
      coord[3]=M_PI+coord[3];

      //reverse derivative
      double pos[4]={tt,coord[1],coord[2],coord[3]};
      double NN, beta[3];
      computeNBeta(pos,NN,beta);

      double beta_r=beta[0], beta_t=beta[1];
      coord[4]=-coord[4]+2.*beta_r/NN;
      coord[5]=-coord[5]+2.*beta_t/NN;	
    }
  } 
}

void NumericalMetricLorene::computeNBeta(const double coord[4],
					 double &NN,double beta[3]) const
{
  GYOTO_DEBUG << endl;
  double tt=coord[0], rr=coord[1],th=coord[2],rsinth = rr*sin(th),ph=coord[3];
  if (rr==0.) throwError("In NumericalMetricLorene.C::computeNBeta r is 0!");
  if (rsinth==0.) throwError("In NumericalMetricLorene.C::computeNBeta "
			     "on z axis!");
  double rm1 = 1./rr, rsm1 = 1./rsinth;

  int it=nb_times_-1;
  while(tt<times_[it] && it>=0){ //ASSUMES backward integration, to generalize
    it--;
  }

  // if (rr<0.187) it=1305; // TEST!!! 

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
    throwError("NumericalMetricLorene::gmunu: incoherent value of indice_time");
  
  if ( mu<0 || mu>3 || nu<0 || nu>3)
       throwError("In NumericalMetricLorene::gmunu bad indice value");

  double rr=pos[0], r2=rr*rr, th=pos[1], rsinth=rr*sin(th);
  if (rr==0.) throwError("In NumericalMetricLorene.C::gmunu r is 0!");
  if (rsinth==0.) throwError("In NumericalMetricLorene.C::gmunu on z axis!");
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
  if (res!=res) throwError("NumericalMetricLorene::gmunu is nan!");
  if (res==res+1.) throwError("NumericalMetricLorene::gmunu is inf!");
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
    throwError("NumericalMetricLorene::gmunu_up_dr: "
	       "incoherent value of indice_time");
  
  if ( (mu!=0 && mu!=3) || (nu!=0 && nu!=3))
       throwError("In NumericalMetricLorene::gmunu_up_dr bad indice value");
  // NB: so far only t and phi components are coded

  double rr=pos[0], th=pos[1], rsinth=rr*sin(th), ph=pos[2];
  if (rr==0.) throwError("In NumericalMetricLorene.C::gmunu_up_dr r is 0!");
  if (rsinth==0.) throwError("In NumericalMetricLorene.C::gmunu_up_dr "
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
      res = 2.*lapsedr*lapse_valm1;
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
      res = lapse_valm2*rsinthm1*(-rm1+betapdr-2.*betap*lapsedr*lapse_valm1);
    }
  if (res!=res) throwError("NumericalMetricLorene::gmunu_up_dr is nan!");
  if (res==res+1.) throwError("NumericalMetricLorene::gmunu_up_dr is inf!");
  return res; 
}

/////////

double NumericalMetricLorene::christoffel(const double coord[8], 
					  const int alpha, 
					  const int mu, const int nu) const
{
  GYOTO_DEBUG << endl;
  //4D Christoffels
  throwError("In NumericalMetricLorene.C: christoffel not implemented");
  return 0.; //not to have a warning when compiling
}

double NumericalMetricLorene::christoffel3(const double coord[6],
					    const int indice_time, 
					    const int ii, 
					    const int jj, 
					    const int kk) const
{
  GYOTO_DEBUG << endl;
  //Computation of 3D Christoffels from the metric : \Gamma^{ii}_{jj kk}
  //NB: 3-metric supposed to be conformally flat
  if (indice_time<0 || indice_time>nb_times_-1) 
    throwError("NumericalMetricLorene::christoffel3: "
	       "incoherent value of indice_time");

  if ( ii<1 || ii>3 || jj<1 || jj>3 || kk<1 || kk>3 )
       throwError("In NumericalMetricLorene::christoffel3 bad indice value");

  double rr=coord[0], r2=rr*rr, th=coord[1], rsinth=rr*sin(th);
  if (rr==0.) throwError("In NumericalMetricLorene.C::christoffel3 r is 0!");
  if (rsinth==0.) throwError("In NumericalMetricLorene.C::christoffel3 "
			     "on z axis!");
  double rm2=1./r2, rsm1 = 1./rsinth, r2sinth2=r2*sin(th)*sin(th), 
    sinth2=sin(th)*sin(th), ph=coord[2];
  
  Scalar* lapse = (lapse_tab_[indice_time]);

  Scalar lapseder = lapse->dsdr() ;
  lapseder.std_spectral_base() ;
  
  const Vector& shift = *(shift_tab_[indice_time]);
  
  const Sym_tensor& g_ij = *(gamcov_tab_[indice_time]) ;

  const Sym_tensor& g_up_ij = *(gamcon_tab_[indice_time]) ;

  double res=0.;

  if ((ii==1) && (jj==1) && (kk==1)) {
    res = 0.5*g_up_ij(1,1).val_point(rr,th,ph)
      *g_ij(1,1).dsdr().val_point(rr,th,ph);
  }else if(((ii==1) && (jj==1) && (kk==2)) || ((ii==1) && (jj==2) && (kk==1))){
    res = 0.5*g_up_ij(1,1).val_point(rr,th,ph)
      *g_ij(1,1).dsdt().val_point(rr,th,ph);
  }else if((ii==1) && (jj==2) && (kk==2)){
    res = -0.5*g_up_ij(1,1).val_point(rr,th,ph)
      *(r2*g_ij(2,2).dsdr().val_point(rr,th,ph) 
	+ 2.*rr*g_ij(2,2).val_point(rr,th,ph));
  }else if((ii==1) && (jj==3) && (kk==3)){
    res = -0.5*g_up_ij(1,1).val_point(rr,th,ph)
      *(r2sinth2*g_ij(3,3).dsdr().val_point(rr,th,ph)
	+2.*rr*sinth2*g_ij(3,3).val_point(rr,th,ph));
  }else if((ii==2) && (jj==2) && (kk==2)){
    res = 0.5*rm2*g_up_ij(2,2).val_point(rr,th,ph)
      *(r2*g_ij(2,2).dsdt().val_point(rr,th,ph));
  }else if((ii==2) && (jj==3) && (kk==3)){
    res = -0.5*rm2*g_up_ij(2,2).val_point(rr,th,ph)
      *(r2sinth2*g_ij(3,3).dsdt().val_point(rr,th,ph)
	+2.*cos(th)*sin(th)*r2*g_ij(3,3).val_point(rr,th,ph));
  }else if((ii==2) && (jj==1) && (kk==1)){
    res = -0.5*rm2*g_up_ij(2,2).val_point(rr,th,ph)
      *g_ij(1,1).dsdt().val_point(rr,th,ph);
  }else if(((ii==2) && (jj==1) && (kk==2)) || ((ii==2) && (jj==2) && (kk==1))){
    res = 0.5*rm2*g_up_ij(2,2).val_point(rr,th,ph)
      *(r2*g_ij(2,2).dsdr().val_point(rr,th,ph)
	+2.*rr*g_ij(2,2).val_point(rr,th,ph));
  }else if(((ii==3) && (jj==1) && (kk==3)) || ((ii==3) && (jj==3) && (kk==1))){
    res = 0.5*rsm1*rsm1*g_up_ij(3,3).val_point(rr,th,ph)
      *(r2sinth2*g_ij(3,3).dsdr().val_point(rr,th,ph)
	+2.*rr*sinth2*g_ij(3,3).val_point(rr,th,ph));
  }else if(((ii==3) && (jj==2) && (kk==3)) || ((ii==3) && (jj==3) && (kk==2))){
    res = 0.5*rsm1*rsm1*g_up_ij(3,3).val_point(rr,th,ph)
      *(r2sinth2*g_ij(3,3).dsdt().val_point(rr,th,ph)
	+2.*cos(th)*sin(th)*r2*g_ij(3,3).val_point(rr,th,ph));
  }
  //Other Christoffels are 0

  if (res!=res) {
    throwError("NumericalMetricLorene::christoffel3 is nan!");
  }
  if (res==res+1.) throwError("NumericalMetricLorene::christoffel3 is inf!");
  return res; 
}

double NumericalMetricLorene::computeHorizon(const double* pos) const{
  GYOTO_DEBUG << endl;
  if (!hor_tab_ && !horizon_)
    return 0.;

  if (horizon_ && !hor_tab_)
    return horizon_;

  if (hor_tab_ && !horizon_){
    int it=nb_times_-1;
    double tt=pos[0], th=pos[2], phi=pos[3];
    double* times=getTimes();

    while(tt<times[it] && it>=0){ //ASSUMES backward integration, to generalize
      it--;
    }

    //    if (pos[1]<0.187) it=1305; // TEST!!!! 

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

    double t1=times[it-1], t2=times[it], t3=times[it+1], t4=times[it+2];
    double rhor1=computeHorizon(pos,it-1),
      rhor2=computeHorizon(pos,it),
      rhor3=computeHorizon(pos,it+1),
      rhor4=computeHorizon(pos,it+2);
    double values[4]={rhor1,rhor2,rhor3,rhor4};
    rhor=Interpol3rdOrder(tt,it,values);
    return rhor;
  }

  throwError("In NumericalMetricLorene::computeHorizon: "
	     "impossible case");
  
}

double NumericalMetricLorene::computeHorizon(const double* pos, 
					     int indice_time) const{  
  GYOTO_DEBUG << endl;
  if (indice_time<0 || indice_time>nb_times_-1){
    throwError("NumericalMetricLorene::computeHorizon"
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

void NumericalMetricLorene::setParticleProperties(Worldline * line, 
						  const double* coord) const
{ 
  GYOTO_DEBUG << endl;
  double cst[3];
  //norm (ALWAYS conserved here)
  cst[0]=ScalarProd(coord,coord+4,coord+4);
  double tdot=coord[4], phdot=coord[7];
  //p_phi (ALWAYS conserved here)
  cst[1]=gmunu(coord,0,3)*tdot + gmunu(coord,3,3)*phdot;
  if (nb_times_==1){
    //p_t (ONLY conserved for stationary case)
    cst[2]=gmunu(coord,0,0)*tdot + gmunu(coord,0,3)*phdot;
  }
  else{
    // p_t is never 0, put it to 0 when not conserved
    cst[2]=0.;
  }
  line -> setCst(cst,3);
}

void  NumericalMetricLorene::setParameter(string name,
					  string content,
					  string unit) {
  GYOTO_DEBUG << endl;
  if (name=="File") {
    if (mapkind_!="MapEt" && mapkind_!="MapAf"){
      throwError("In NumericalMetricLorene::setParameter "
		 "please provide Map kind before File in XML");
    }
    filename_ = new char[strlen(content.c_str())+1];
    strcpy(filename_,content.c_str());
    setMetricSource();
  }
  else if (name=="MapEt") {
    mapkind_="MapEt";
  }
  else if (name=="MapAf") {
    mapkind_="MapAf";
  }
  else if (name=="Time") initial_time_=atof(content.c_str());
  else if (name=="HasSurface") {
    if (filename_!=NULL){
      throwError("In NumericalMetricLorene::setParameter "
		 "please provide Surface information before File in XML");
    }
    has_surface_=1;
  }
  else if (name=="Horizon") horizon_=atof(content.c_str());
  else if (name=="RefineIntegStep"){
    refine_=1;
    char * tc = const_cast<char*>(content.c_str());
    r_refine_  = strtod(tc, &tc);
    h0_refine_ = strtod(tc, &tc);
  }
  else  Generic::setParameter(name, content, unit);
}

#ifdef GYOTO_USE_XERCES
void NumericalMetricLorene::fillElement(Gyoto::FactoryMessenger *fmp) {
  GYOTO_DEBUG << endl;
  fmp -> setParameter("File", filename_);
  //fmp -> setParameter("Time", collapse_time_);
  Generic::fillElement(fmp);
}

void NumericalMetricLorene::setParameters(FactoryMessenger* fmp) {
  GYOTO_DEBUG << endl;
  string name="", content="", unit="";
  if (fmp) {
    while (fmp->getNextParameter(&name, &content, &unit)) {
      if (name == "File") content = fmp -> fullPath( content );
      GYOTO_DEBUG << "Setting \"" << name << "\" to \"" << content
		  << "\" (in unit \"" << unit << "\")\n";
      setParameter(name, content, unit);
    }
  }
}

#endif
