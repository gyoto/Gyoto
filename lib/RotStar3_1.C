/*
    Copyright 2011-2014, 2016, 2018-2020 Frederic Vincent & Thibaut Paumard

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
#include "star_rot.h"
#include "eos.h"
#include "utilitaires.h"
#include "graphique.h"
#include "nbr_spx.h"
#include "unites.h"

using namespace Lorene;

#include "GyotoUtils.h"
#include "GyotoRotStar3_1.h"
#include "GyotoError.h"
#include "GyotoFactoryMessenger.h"
#include "GyotoProperty.h"

#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cstring>
#include <ctime>

using namespace Gyoto;
using namespace Gyoto::Metric;

GYOTO_PROPERTY_START(RotStar3_1)
GYOTO_PROPERTY_BOOL(RotStar3_1, GenericIntegrator, SpecificIntegrator,
		    genericIntegrator)
GYOTO_PROPERTY_FILENAME(RotStar3_1, File, file)
GYOTO_PROPERTY_END(RotStar3_1, Generic::properties)

// Lorene Metrics are not thread-safe
GYOTO_PROPERTY_THREAD_UNSAFE(RotStar3_1)

RotStar3_1::RotStar3_1() : 
Generic(GYOTO_COORDKIND_SPHERICAL, "RotStar3_1"),
  filename_(NULL),
  star_(NULL),
  integ_kind_(1)
{}

RotStar3_1::RotStar3_1(const RotStar3_1& o) : 
  Generic(o),
  filename_(NULL),
  star_(NULL),
  integ_kind_(o.integ_kind_)
{
  kind("RotStar3_1");
  fileName(o.fileName());
}

RotStar3_1* RotStar3_1::clone() const {
  return new RotStar3_1(*this);
}

RotStar3_1::~RotStar3_1() 
{
  if (star_) {
    const Map& mp=star_ -> get_mp();
    const Mg3d* mg=mp.get_mg();
    const Map* mpp=&mp;
    delete star_;
    delete mpp;
    delete mg;
  }
  
  if (filename_) delete [] filename_;

  if (debug()) cout << "RotStar3_1 Destruction" << endl;
}

void RotStar3_1::file(std::string const &fname) {
  cerr << "Setting file name to '" << fname << "'" << endl;
  fileName(fname.c_str());
}

std::string RotStar3_1::file() const {
  if (!filename_) return "";
  return filename_;
}

void RotStar3_1::fileName(char const * lorene_res) {
  if (filename_) { delete[] filename_; filename_=NULL; }
  if (star_) {
    const Map& mp=star_ -> get_mp();
    const Mg3d* mg=mp.get_mg();
    const Map* mpp=&mp;
    delete star_; star_=NULL;
    delete mpp;
    delete mg;
  }
  if (!lorene_res) return;

  filename_ = new char[strlen(lorene_res)+1];
  strcpy(filename_,lorene_res);
  FILE* resfile=fopen(lorene_res,"r");
  if (!resfile) GYOTO_ERROR(string("No such file or directory: ")+lorene_res);
  Mg3d* mg = new Mg3d(resfile);
  Map_et* mps = new Map_et(*mg,resfile);
  Eos* p_eos = Eos::eos_from_file(resfile);
  star_ = new Lorene::Star_rot(*mps,*p_eos,resfile);
  star_ -> equation_of_state();
  star_ -> update_metric();
  star_ -> hydro_euler();

  tellListeners();
}

char const * RotStar3_1::fileName() const { return filename_; }

void RotStar3_1::integKind(int ik) { integ_kind_ = ik; }
int RotStar3_1::integKind() const { return integ_kind_; }

bool RotStar3_1::genericIntegrator() const {return !integ_kind_;}
void RotStar3_1::genericIntegrator(bool t) {integ_kind_=!t;}

int RotStar3_1::diff(state_t const &coord, state_t &res, double /* mass */) const
{
  //4-DIMENSIONAL INTEGRATION
  //NB: this diff is only called by Generic::RK4

  //if (debug()) cout << "In 4D RotStar diff [8]..." << endl;
  //clock_t time1, time2;
  //double diftime, clocks = CLOCKS_PER_SEC;
  //time1 = clock();

  double rr=coord[1],r2=rr*rr,th=coord[2],sinth2=sin(th)*sin(th),ph=coord[3];
  //LAPSE
  const Scalar & NNscal=star_ -> get_nn();
  double NN=NNscal.val_point(rr,th,ph), N2=NN*NN;
  double N_r=NNscal.dsdr().val_point(rr,th,ph);
  double N_th=NNscal.dsdt().val_point(rr,th,ph);

  //SHIFT (OMEGA)
  const Scalar & omega_scal=star_ -> get_nphi();//NB: avoid using star_ -> get_beta(), it is very time-consuming
  double omega=omega_scal.val_point(rr,th,ph), omega2=omega*omega;
  double omega_r=omega_scal.dsdr().val_point(rr,th,ph);
  double omega_th=omega_scal.dsdt().val_point(rr,th,ph);

  //METRIC POTENTIALS
  const Scalar & A2scal=star_ -> get_a_car();
  const Scalar & B2scal=star_ -> get_b_car();
  double A2=A2scal.val_point(rr,th,ph), B2=B2scal.val_point(rr,th,ph);
  double A2_r=A2scal.dsdr().val_point(rr,th,ph), B2_r=B2scal.dsdr().val_point(rr,th,ph);
  double A2_th=A2scal.dsdt().val_point(rr,th,ph), B2_th=B2scal.dsdt().val_point(rr,th,ph);

  /*  time2 = clock();
  diftime = time2 - time1;
  if (debug()) cout << "Time elapsed in temp 4D diff (in sec)= " << setprecision(GYOTO_PREC) << setw(GYOTO_WIDTH) << diftime/clocks << endl;*/

  //METRIC COEF
  double gtt=-1./N2, grr=1./A2, gthth=1./(A2*r2), gpp=1./(B2*r2*sinth2)-omega2/N2, gtp=-omega/N2;
  double g_ttr=-2.*NN*N_r+B2_r*omega2*r2*sinth2+2.*omega*omega_r*B2*r2*sinth2+2.*rr*B2*omega2*sinth2, g_ttth=-2.*NN*N_th+B2_th*omega2*r2*sinth2+2.*omega*omega_th*B2*r2*sinth2+2.*cos(th)*sin(th)*r2*B2*omega2;
  double g_rrr=A2_r, g_rrth=A2_th;
  double g_ththr=r2*A2_r+2.*rr*A2, g_ththth=r2*A2_th;
  double g_ppr=sinth2*(r2*B2_r+2.*rr*B2), g_ppth=r2*(sinth2*B2_th+2.*cos(th)*sin(th)*B2);
  double g_tpr=-omega_r*B2*r2*sinth2-omega*B2_r*r2*sinth2-2.*rr*omega*B2*sinth2, g_tpth=-omega_th*B2*r2*sinth2-omega*B2_th*r2*sinth2-2.*cos(th)*sin(th)*omega*B2*r2;

  //4D CHRISTOFFELS in stationnary axisym spacetime
  double ch001=1./2.*gtt*g_ttr+1./2.*gtp*g_tpr, ch002=1./2.*gtt*g_ttth+1./2.*gtp*g_tpth, ch031=1./2.*gtt*g_tpr+1./2.*gtp*g_ppr, ch032=1./2.*gtt*g_tpth+1./2.*gtp*g_ppth, ch331=1./2.*gpp*g_ppr+1./2.*gtp*g_tpr, ch332=1./2.*gpp*g_ppth+1./2.*gtp*g_tpth, ch301=1./2.*gpp*g_tpr+1./2.*gtp*g_ttr, ch302=1./2.*gpp*g_tpth+1./2.*gtp*g_ttth, ch122=-1./2.*grr*g_ththr, ch133=-1./2.*grr*g_ppr, ch100=-1./2.*grr*g_ttr, ch111=1./2.*grr*g_rrr, ch112=1./2.*grr*g_rrth, ch103=-1./2.*grr*g_tpr, ch200=-1./2.*gthth*g_ttth, ch211=-1./2.*gthth*g_rrth, ch233=-1./2.*gthth*g_ppth, ch222=1./2.*gthth*g_ththth, ch221=1./2.*gthth*g_ththr, ch203=-1./2.*gthth*g_tpth;

  //RESULT: Derivatives
  res[0]=coord[4];
  res[1]=coord[5];
  res[2]=coord[6];
  res[3]=coord[7];
        //From the 4D equation of geodesics:
  res[4]=-2.*ch001*coord[4]*coord[5]-2.*ch002*coord[4]*coord[6]-2.*ch031*coord[7]*coord[5] - 2.*ch032*coord[7]*coord[6];//added 2011-02-27 ch032 term forgotten! to be checked
  res[5]=-ch111*coord[5]*coord[5]-ch122*coord[6]*coord[6]-ch133*coord[7]*coord[7]-ch100*coord[4]*coord[4]-2.*ch112*coord[5]*coord[6]-2.*ch103*coord[4]*coord[7];
  res[6]=-ch211*coord[5]*coord[5]-ch222*coord[6]*coord[6]-ch233*coord[7]*coord[7]-ch200*coord[4]*coord[4]-2.*ch221*coord[5]*coord[6]-2.*ch203*coord[4]*coord[7];
  res[7]=-2.*ch301*coord[4]*coord[5]-2.*ch302*coord[4]*coord[6]-2.*ch331*coord[7]*coord[5]-2.*ch332*coord[7]*coord[6];

  /*  time2 = clock();
  diftime = time2 - time1;
  if (debug()) cout << "TOTAL Time elapsed in 4D diff (in sec)= " << setprecision(GYOTO_PREC) << setw(GYOTO_WIDTH) << diftime/clocks << endl;*/

  return 0;
  
}

int RotStar3_1::diff(const double y[6], double res[6], int) const
{
  //3+1 INTEGRATION
  //NB: this diff is only called by RotStar::RK4
  //NBB: here t=theta, not time!

  //if (debug()) cout << "In 3+1 D RotStar::diff" << endl;
  // clock_t time1, time2;
  //double diftime, clocks = CLOCKS_PER_SEC;
  //time1 = clock();
  //NB: all the computations used here are detailed in my 3+1 notes.

  //NB: Lorene coordinates are spherical-like
  //Variables: y=[r,theta,phi,Vr,Vtheta,Vphi] see definitions 3+1 geod paper
  double rr=y[0], r2=rr*rr, th=y[1], sinth2=sin(th)*sin(th), phi=y[2];
  /*
    Important remark! [actually it is not used because the metric coef are known in function of A, B, N... ; but interesting to keep this in mind]
    Lorene gcon and gcov function return metric coef expressed in the orthonormal spherical tetrad d/dr, 1/r d/dtheta, 1/(r*sin(theta)) d/dphi
    There's thus a change of basis to come back to the natural basis of spherical coordinates.
   */

  //LAPSE
  const Scalar & NNscal=star_ -> get_nn();
  double NN=NNscal.val_point(rr,th,phi);//, NN2=NN*NN;
  if (NN == 0.) GYOTO_ERROR("In RotStar3_1.C: NN==0!!");
  double Nr=NNscal.dsdr().val_point(rr,th,phi);
  double Nt=NNscal.dsdt().val_point(rr,th,phi);

  //SHIFT (OMEGA)
  const Scalar & omega_scal=star_ -> get_nphi();
  double omega=omega_scal.val_point(rr,th,phi);
  double omega_r=omega_scal.dsdr().val_point(rr,th,phi);
  double omega_t=omega_scal.dsdt().val_point(rr,th,phi);

  //METRIC POTENTIALS
  const Scalar & A2scal=star_ -> get_a_car();
  const Scalar & B2scal=star_ -> get_b_car();
  double A2=A2scal.val_point(rr,th,phi), B2=B2scal.val_point(rr,th,phi);
  double A2_r=A2scal.dsdr().val_point(rr,th,phi), B2_r=B2scal.dsdr().val_point(rr,th,phi);
  double A2_th=A2scal.dsdt().val_point(rr,th,phi),B2_th=B2scal.dsdt().val_point(rr,th,phi);

  /*  time2 = clock();
  diftime = time2 - time1;
  if (debug()) cout << "Time elapsed in 3+1 diff (in sec)= " << setprecision(GYOTO_PREC) << setw(GYOTO_WIDTH) << diftime/clocks << endl;*/

  //METRIC COEF
  double grr=1./A2, gtt=1./(A2*r2), gpp=1./(B2*r2*sinth2), g_rrr=A2_r, g_rrt=A2_th, g_ttr=r2*A2_r+2.*rr*A2, g_ttt=r2*A2_th, g_ppr=r2*sinth2*B2_r+2.*rr*B2*sinth2, g_ppt=r2*sinth2*B2_th+2.*sin(th)*cos(th)*r2*B2;
  
  /* //Change of basis [Not used]
  gtt_=1./r2*gtt_;gpp_=1./(r2*sinth2)*gpp_;//natural basis contrav metric
  g_ttr_=2.*rr*g_tt_+r2*g_ttr;g_ttt_=r2*g_ttt_;g_ppr_=sinth2*(2.*rr*g_pp_+r2*g_ppr_);g_ppt_=r2*(2.*cos(th)*sin(th)*g_pp_+sinth2*g_ppt_);//natural basis derived covar metric (NB: the cov metric ceof that appear here must be in the orthonormal basis)
  */

  //EXTRINSIC CURVATURE
  double Krp=-1./(2.*NN)*B2*r2*sinth2*omega_r;
  double Ktp=-1./(2.*NN)*B2*r2*sinth2*omega_t;

  //3-CHRISTOFFELS
  double Grrr=1./2.*grr*g_rrr,
    Grrt=1./2.*grr*g_rrt,
    Gttr=1./2.*gtt*g_ttr,
    Gttt=1./2.*gtt*g_ttt,
    Gppr=1./2.*gpp*g_ppr,
    Gptp=1./2.*gpp*g_ppt,
    Grtt=-1./2.*grr*g_ttr,
    Grpp=-1./2.*grr*g_ppr,
    Gtrr=-1./2.*gtt*g_rrt,
    Gtpp=-1./2.*gtt*g_ppt;//these are the only non-0 Christo

  /* OBSOLETE VERSION */
//   //Christoffel terms (xdot^j+beta^j)(xdot^k+beta^k)Gamma^i_jk in geodesic eq.:
//   double Christo_r=rdot*rdot*Grrr+2.*rdot*thdot*Grrt+thdot*thdot*Grtt+(phidot-omega)*(phidot-omega)*Grpp, 
//     Christo_th=2.*rdot*thdot*Gttr+thdot*thdot*Gttt+rdot*rdot*Gtrr+(phidot-omega)*(phidot-omega)*Gtpp,
//     Christo_ph=2.*rdot*(phidot-omega)*Gppr+2.*thdot*(phidot-omega)*Gptp;
//   double kappa=1./NN2*(rdot*Nr+thdot*Nt)-2./NN2*rdot*(phidot-omega)*Krp-2./NN2*thdot*(phidot-omega)*Ktp;

//   //RESULT: derivatives (see geodesic equation in 3+1)
//   res[0]=rdot;
//   res[1]=thdot;
//   res[2]=phidot;
//   res[3]=kappa*NN*rdot-NN*grr*Nr+1./NN*(Nr*rdot+Nt*thdot)*rdot+2.*NN*grr*(phidot-omega)*Krp-Christo_r;//-Grtt*thdot*thdot-Grpp*phidot*phidot;
//   res[4]=kappa*NN*thdot-NN*gtt*Nt+1./NN*(Nr*rdot+Nt*thdot)*thdot+2.*NN*gtt*(phidot-omega)*Ktp-Christo_th;//-2.*Gttr*rdot*thdot-Gtpp*phidot*phidot;
//   res[5]=kappa*NN*(phidot-omega)+1./NN*(Nr*rdot+Nt*thdot)*(phidot-omega)+2.*rdot*omega_r+2.*thdot*omega_t+2.*NN*gpp*(rdot*Krp+thdot*Ktp)-Christo_ph;//-2.*Gppr*rdot*phidot-2.*Gptp*thdot*phidot;

//   /* NEW VERSION JULY 2011 */

  double Vr = y[3], Vth = y[4], Vph = y[5];
  double Christo_r=Vr*Vr*Grrr+2.*Grrt*Vr*Vth+Grtt*Vth*Vth+Grpp*Vph*Vph,
    Christo_th=Gtrr*Vr*Vr+2.*Gttr*Vr*Vth+Gttt*Vth*Vth+Gtpp*Vph*Vph,
    Christo_ph=2.*Gppr*Vr*Vph+2.*Gptp*Vth*Vph;

  //See 3+1 equation of geodesics
  double prefact=1./NN*Vr*Nr+1./NN*Vth*Nt-2.*Krp*Vr*Vph-2.*Ktp*Vth*Vph;
  double Vrdot = NN*(Vr*prefact+2.*grr*Krp*Vph-Christo_r)-grr*Nr,
    Vthdot = NN*(Vth*prefact+2.*gtt*Ktp*Vph-Christo_th)-gtt*Nt,
    Vphdot = NN*(Vph*prefact+2.*gpp*(Krp*Vr+Ktp*Vth)-Christo_ph)
    +Vr*omega_r+Vth*omega_t; //NB: dot = wrt coordinate time t
    
    
  res[0]=NN*Vr; //dr/dt
  res[1]=NN*Vth; //dtheta/dt
  res[2]=NN*Vph+omega; //dphi/dt
  res[3]=Vrdot; //dVr/dt
  res[4]=Vthdot; //dVtheta/dt
  res[5]=Vphdot; //dVphi/dt

  
  /*  time2 = clock();
  diftime = time2 - time1;
  if (debug()) cout << "TOTAL Time elapsed in 3+1 diff (in sec)= " << setprecision(GYOTO_PREC) << setw(GYOTO_WIDTH) << diftime/clocks << endl;*/
  
  return 0;
}


int RotStar3_1::myrk4(const double coorin[6], double h, double res[6]) const
{
  //if (debug()) cout << "In RotStar::rk4" << endl;

  //Here the integration must be 3+1:
  if (!integ_kind_) GYOTO_ERROR("In RotStar3_1::myrk4: Impossible case");
  
  double k1[6], k2[6], k3[6], k4[6], coor_plus_halfk1[6], sixth_k1[6], 
    coor_plus_halfk2[6], third_k2[6], coor_plus_k3[6], third_k3[6], sixth_k4[6]; 

  if(diff(coorin,k1,1)) return 1;
	  
  for (int i=0;i<6;i++) {
    k1[i]=h*k1[i];
    coor_plus_halfk1[i]=coorin[i]+0.5*k1[i];
    sixth_k1[i]=1./6.*k1[i];
  }

  if(diff(coor_plus_halfk1,k2,1)) return 1;
  
  for (int i=0;i<6;i++) {
    k2[i]=h*k2[i];
    coor_plus_halfk2[i]=coorin[i]+0.5*k2[i];
    third_k2[i]=1./3.*k2[i];
  }

  if(diff(coor_plus_halfk2,k3,1)) return 1;

  for (int i=0;i<6;i++) {
    k3[i]=h*k3[i];
    coor_plus_k3[i]=coorin[i]+k3[i];
    third_k3[i]=1./3.*k3[i];
  }

  if(diff(coor_plus_k3,k4,1)) return 1;

  for (int i=0;i<6;i++) {
    k4[i]=h*k4[i];
    sixth_k4[i]=1./6.*k4[i];
  }

  for (int i=0;i<6;i++) {
    res[i]=coorin[i]+sixth_k1[i]+third_k2[i]+third_k3[i]+sixth_k4[i];
  }
  
  return 0;
}

//int RotStar3_1::myrk4_adaptive(const double coord[8], double lastnorm, double normref, double coordnew[8], double h0, double& h1, int &) const
int RotStar3_1::myrk4_adaptive(const double coord[6], double, double normref, double coordnew[6], double cst[2], double& tdot_used, double h0, double& h1, double h1max, double& hused) const
{

  // if (debug()) cout << "In Rotstar::adaptive [6]" << endl;

  double delta0[6];
  double delta0min=1e-15;
  double dcoord[6];
  double eps=0.0001;
  double S=0.9;
  double errmin=1e-6;
  // double factnorm=2.;
  double sigh1=1.;

  h1max=deltaMax(coord, h1max);

  /*if (debug()) cout << "RotStar.C: coord in rk=";
  for (int ii=0;ii<8;ii++) if (debug()) cout << coord[ii] << " " ;
  if (debug()) cout << endl;*/

  diff(coord,dcoord,1);

  for (int i = 0;i<6;i++) {
    delta0[i]=delta0min+eps*(fabs(h0*dcoord[i]));
    //if (debug()) cout << "Rot: delta0[i]=" << delta0[i] << endl;
  }

  double hbis=0.5*h0;
  double coordhalf[6];
  double coord2[6];
  double delta1[6];
  
  double err;
  int count=0;
  //  double newnorm;

  while (1){
    count++;
    //if (debug()) cout << "count in rk Rot= " << count << endl;
    if (count > 100){
      GYOTO_ERROR("RotStar: bad rk");
    }
    err=0.;
    myrk4(coord,h0,coordnew);
    myrk4(coord,hbis,coordhalf);
    myrk4(coordhalf,hbis,coord2);

    double coordnewbis[6], coordhalfbis[6], coord2bis[6], tdot_unused;

    int normalize=1;
    /*
      NB: for the time being (2011-03-11), normalize=1 imposes to update phprime and tdot to insure conservation of cst of motion;
      In order to insure norm conservation, there is an other flag in Normalize4v: consnorm, to be put to the desired value
      Caution: be sure it is really needed before putting consnorm to 1, it slows down A LOT the computation
     */

    if (normalize){
      Normalize4v(coordnew,coordnewbis,cst,tdot_used);
      Normalize4v(coordhalf,coordhalfbis,cst,tdot_unused);
      Normalize4v(coord2,coord2bis,cst,tdot_unused);

      for (int ii=0;ii<6;ii++){
	coordnew[ii]=coordnewbis[ii];
	coordhalf[ii]=coordhalfbis[ii];
	coord2[ii]=coord2bis[ii];
      }
      
    }else{

      if (fabs(normref)<fabs(normref+1.)){ // NULL GEODESIC
	if (tdot_used<=0.) {//NB: default value of tdot_used has correct sign
	  tdot_used=-1.;//-1 or +1 is chosen by hand ; the norm is proportionnal to tdot_used^2 and is supposed to be =0, so the actual value of tdot_used doesn't matter
	}else{
	  tdot_used=1.;
	}
      }else{ // TIMELIKE GEODESIC
	double rprime=coordnew[3], thprime=coordnew[4], phprime=coordnew[5];
	double pos[4]={0.,coordnew[0],coordnew[1],coordnew[2]};
	double g_tt=gmunu(pos,0,0), g_tp=gmunu(pos,0,3), g_rr=gmunu(pos,1,1), g_thth=gmunu(pos,2,2), g_pp=gmunu(pos,3,3);
	//	if (debug()) cout << "time integ" << endl;
	double ds2=g_tt+2.*g_tp*phprime+g_rr*rprime*rprime+g_thth*thprime*thprime+g_pp*phprime*phprime;
	if (ds2>0) GYOTO_ERROR("In RotStar3_1.C: impossible to compute timelike norm!");
	if (tdot_used<=0.){//NB: default value of tdot_used has correct sign
	  tdot_used=-sqrt(-1./ds2);
	}else{
	  tdot_used=sqrt(-1./ds2);
	}
	
      }
    }


    
    for (int i = 0;i<6;i++){
      delta1[i]=coord2[i]-coordnew[i];
      if (err<fabs(delta1[i]/delta0[i])) err=fabs(delta1[i]/delta0[i]);
    }

    if (err>1) {
      h0=S*h0*pow(err,-0.25);
      hbis=0.5*h0;
    }else{
      h1=(err > errmin ? S*h0*pow(err,-0.2) : 4.*h0);//pour éviter les explosions
      if (h1<0.) sigh1=-1.;//why sigh1 and fabs(h1)? because otherwise if h1<0 (possible here if backwards integration), h1 is < delta_min_, so h1 is always set to delta_min_...
      if (fabs(h1)<delta_min_) h1=sigh1*delta_min_;
      if (fabs(h1)>h1max) h1=sigh1*h1max;
      hused=h0;

      break;
    }
  }

  return 0;
}

int RotStar3_1::myrk4_adaptive(Worldline* line, state_t const &coord, 
			       double lastnorm, double normref, 
			       state_t &coordnew, double h0, 
			       double& h1, double h1max) const
{
  //  if (debug()) cout << "In Rotstar::adaptive [8]" << endl;
  if (coord[1] < 2.5) {//inside rotating star -> a ameliorer
    if (debug()) cout << "In RotStar3_1.C: Particle has reached the rotating star. Stopping integration." << endl;
    return 1;
  }
  
  if (!integ_kind_) {
    /*
      To use 4D integration
      Here, the integration is performed by the most general Generic::myrk4 + Generic::diff functions that only call the christoffels (basic geodesic equation).
      The function christoffel being defined here in RotStar3_1, it is the 4D-christo computed thanks to 3+1 quantities that are used.
    */

    if (Generic::myrk4_adaptive(line,coord,lastnorm,normref,coordnew,h0,h1,h1max)) {
      return 1;
    }else{
      return 0;
    }
  }

  //Below: 3+1 integration
  //Variables used: [r,theta,phi,Vr,Vtheta,Vphi] for definition see 3+1 geod paper

  double rr=coord[1],th=coord[2],ph=coord[3],tdot=coord[4],rdot=coord[5],thdot=coord[6],phdot=coord[7],rprime=rdot/tdot,thprime=thdot/tdot,phprime=phdot/tdot;

  const Scalar & NNscal=star_ -> get_nn();
  double NN=NNscal.val_point(rr,th,ph);//, NN2=NN*NN;
  if (NN == 0.) GYOTO_ERROR("In RotStar3_1.C: NN==0!!");
  const Scalar & omega_scal=star_ -> get_nphi();
  double omega=omega_scal.val_point(rr,th,ph);
  
  double Vr = 1./NN*rprime, Vth=1./NN*thprime, Vph=1./NN*(phprime-omega);

  double g_tt=gmunu(coord.data(),0,0), g_tp=gmunu(coord.data(),0,3), g_pp=gmunu(coord.data(),3,3);
  double cst_p_t=g_tt*tdot+g_tp*phdot, cst_p_ph=g_pp*phdot+g_tp*tdot;//Cst of motion because the vectors d/dt and d/dphi are Killing
  //if (debug()) cout << "Rot: cst= " << cst_p_t << " " << cst_p_ph << endl;
  double cst[2]={cst_p_t,cst_p_ph};
  //cout << "3+1 Cst of motion= " << cst_p_t << " " << cst_p_ph << endl;
  double coor[6]={rr,th,ph,Vr,Vth,Vph};
 
  double coornew[6];
  double hused=1000.;

  if (tdot<0. && h0>0.) h0*=-1.;//to integrate backwards if tdot<0

  double tdot_used=tdot;
  /*if (tdot<=0.){
    tdot_used=-1000.;
  }else{
    tdot_used=1000.;
    }//tdot_used thus has the correct sign*/

  if (myrk4_adaptive(coor,lastnorm,normref,coornew,cst,tdot_used,h0,h1,delta_max_,hused)) return 1;
  //  if (debug()) cout << "tdot_used in rk-8= " << tdot_used << endl;
  
  //phdot=coornew[5]*tdot_used;rdot=coornew[3]*tdot_used;thdot=coornew[4]*tdot_used;
  
  NN=NNscal.val_point(coornew[0],coornew[1],coornew[2]);
  omega=omega_scal.val_point(coornew[0],coornew[1],coornew[2]);
  phdot=(NN*coornew[5]+omega)*tdot_used;rdot=NN*coornew[3]*tdot_used;thdot=NN*coornew[4]*tdot_used;

  coordnew[0]=coord[0]+hused;coordnew[1]=coornew[0];coordnew[2]=coornew[1];coordnew[3]=coornew[2];coordnew[4]=tdot_used;coordnew[5]=rdot;coordnew[6]=thdot;coordnew[7]=phdot;
  //  if (debug()) cout << "Rot: norm at end rk-8= " << ScalarProd(coordnew,coordnew+4,coordnew+4) << endl;

  return 0;
}

void RotStar3_1::Normalize4v(const double coordin[6], double coordout[6], const double cst[2], double& tdot_used) const{
  
  //Here coordin=[r,theta,phi,Vr,Vtheta,Vphi]

  double posin[4]={0.,coordin[0],coordin[1],coordin[2]};//posin={t,rr,th,ph} with t=anything (g_munu are independent of t), {rr,th,ph}=cf coordin
  double g_tt=gmunu(posin,0,0), g_rr=gmunu(posin,1,1), g_thth=gmunu(posin,2,2), g_tp=gmunu(posin,0,3), g_pp=gmunu(posin,3,3), cst_p_t=cst[0], cst_p_ph=cst[1];
  double phdot,phprime;
  //double phprime_init=coordin[5],dphpr=0.01;
  const Scalar & NNscal=star_ -> get_nn();
  double NN=NNscal.val_point(coordin[0],coordin[1],coordin[2]);//, NN2=NN*NN;
  if (NN == 0.) GYOTO_ERROR("In RotStar3_1.C: NN==0!!");
  const Scalar & omega_scal=star_ -> get_nphi();
  double omega=omega_scal.val_point(coordin[0],coordin[1],coordin[2]);
  double phprime_init=NN*coordin[5]+omega,dphpr=0.01;

  // Changing phdot and tdot (thus phprime) to insure conservation of cst of motion
  if (g_tt != 0. && (g_tt*g_pp != g_tp*g_tp)){
    phdot=(cst_p_ph - g_tp/g_tt*cst_p_t)/(g_pp-g_tp*g_tp/g_tt); tdot_used=(cst_p_t - g_tp*phdot)/g_tt; phprime=phdot/tdot_used;//cf constants of motion
    //    if (debug()) cout << "tdot use in Normalize= " << tdot_used << endl;
  }else{
    GYOTO_ERROR("RotStar3_1.C: special case metric coef=0 to handle in Normalize4v...");
  }
  if (fabs(phprime-phprime_init)>dphpr*fabs(phprime_init)){
    if (verbose() >= GYOTO_SEVERE_VERBOSITY)
      cerr << "WARNING (severe):" << endl
	   << "Too big change in phprime: " << phprime_init << " " << phprime << endl;
  }
  
  /*double fourvel[4]={tdot_used,NN*coordin[3]*tdot_used,NN*coordin[4]*tdot_used,phprime*tdot_used}; //4-velocity with correct value of tdot and phdot
  cout << "DEBUT AFFICHE" << endl;
  cout << "Vr= " << coordin[3] << " " << NN << " " << tdot_used << " " << NN*coordin[3]*tdot_used << endl;
  cout << "fourvel= " ;
  for (int ii=0;ii<4;ii++) cout << fourvel[ii] << " " ;
  cout << endl;
  cout << "-------> tdot= " << tdot_used << endl;
  double normaffich=ScalarProd(posin,fourvel,fourvel);
  cout << "cst motion in Normalize= " << normaffich << " " << cst_p_t << " " << cst_p_ph << endl;
  if (fabs(normaffich)>1000. || normaffich!=normaffich) GYOTO_ERROR("testrotstar --");
  */
  //  double rprime=coordin[3], thprime=coordin[4];
  double rprime=NN*coordin[3], thprime=NN*coordin[4];
  double rprime_new=rprime, thprime_new=thprime;
  int consnorm=0;
  if (consnorm){
    // Changing rprime and thprime to insure norm conservation
    
    double norminit=tdot_used*tdot_used*(g_tt+2.*g_tp*phprime+g_rr*rprime*rprime+g_thth*thprime*thprime+g_pp*phprime*phprime);//ScalarProd(tempvec,tempvec+4,tempvec+4);//
    double normref;
    if (fabs(norminit)<fabs(norminit+1.)){
      normref=0.;
    }else{
      normref=-1.;
    }
    
    double aa=tdot_used*tdot_used*(g_tt+2.*g_tp*phprime+g_pp*phprime*phprime), bb=tdot_used*tdot_used*g_rr, cc=tdot_used*tdot_used*g_thth;
    //Aim: find X and Y so that a+b*X^2+c*Y^2=0 (or -1), then rprime_new=X, thprime_new=Y (NB: phprime was already changed above, so its value is fixed, that's why it's put in the cst 'a')
    
    //if (debug()) cout << "aa, bb, cc= " << aa << " " << bb << " " << cc << endl;
    //if (debug()) cout << "rdot, thdot ini= " << rdot << " " << thdot << endl;
    
    
    
    double fact=.001;//allows changes of order 0.1% of original value (don't put it too big!! or the corrected rprime and thprime will depend a lot on the chosen nstep -- no convergence towards a given value when nstep grows...)
    int nstep=100;//nstep search steps on both sides of initial rdot (or thdot) value ; search interval between 2 trials = fact*[init value]/nstep = dr/nstep (or dth/nstep)
    double dth=fabs(fact*thprime), dr=fabs(fact*rprime), newnorm, lastnorm=norminit, rp, thp, intth=dth/double(nstep), intr=dr/double(nstep);
    
    // cout << "Norm init= "<<norminit<<endl;
    
    for (int ii=-nstep/2;ii<=nstep/2;ii++){ // ii = -dr/2,-dr/2+dr/nstep,...,dr/2
      for (int jj=-nstep/2;jj<=nstep/2;jj++){
	rp=rprime+ii*intr;
	thp=thprime+jj*intth;
	newnorm=aa+bb*rp*rp+cc*thp*thp;
	//  cout << "New norm= " << newnorm << endl;
	if (fabs(newnorm-normref)<fabs(lastnorm-normref)){
	  rprime_new=rp;thprime_new=thp;lastnorm=newnorm;
	}
      }
    }
  } // end if consnorm
  
  //cout << "New norm= " << norm_new << endl;
  //coordout[0]=coordin[0];coordout[1]=coordin[1];coordout[2]=coordin[2];coordout[3]=rprime_new;coordout[4]=thprime_new;coordout[5]=phprime;
  coordout[0]=coordin[0];coordout[1]=coordin[1];coordout[2]=coordin[2];coordout[3]=1./NN*rprime_new;coordout[4]=1./NN*thprime_new;coordout[5]=1./NN*(phprime-omega);

}

double RotStar3_1::gmunu(const double * pos, int mu, int nu) const
{
  /*
    4-metric coefficients
    cf Eric's Rotating Stars Notes Eq. 2.32
  */
  
  //if (debug()) cout << "In gmunu Rot" << endl;
  double rr=pos[1],r2=rr*rr,th=pos[2],sinth2=sin(th)*sin(th),ph=pos[3];
  const Scalar & NNtemp=star_ -> get_nn();
  double NN=NNtemp.val_point(rr,th,ph), N2=NN*NN;
  const Scalar & omega_scal=star_ -> get_nphi();
  double omega=omega_scal.val_point(rr,th,ph);
  double B2=(star_ -> get_b_car()).val_point(rr,th,ph);
  double A2=(star_ -> get_a_car()).val_point(rr,th,ph);
  double g_tt=(B2*r2*sinth2*omega*omega-N2), g_tp=-omega*B2*r2*sinth2, g_rr=A2, 
    g_thth=A2*r2, g_pp=B2*r2*sinth2;

  //DEBUG!!!! flat metric ***************
  //g_tt=-1.;g_tp=0.;g_rr=1.;g_thth=r2;g_pp=r2*sinth2;
  //*************************************

  if ((mu==0) && (nu==0)) return g_tt;
  if ((mu==1) && (nu==1)) return g_rr;
  if ((mu==2) && (nu==2)) return g_thth;
  if ((mu==3) && (nu==3)) return g_pp;
  if (((mu==0) && (nu==3)) || ((mu==3) && (nu==0))) return g_tp;

  return 0.;
}

//double RotStar3_1::christoffel(const double coord[8], const int alpha, const int mu, 
//		       const int nu) const
double RotStar3_1::christoffel(const double coord[8], const int alpha, 
			       const int mu, const int nu) const
{
  /*
    The computation of the christo is easy since we know the expression of gmunu as a function of 3+1 quantities, and since Lorene allows to perform derivatives on quantities. So gmunu,sigma is computable. 
   */

  //if (debug()) cout << "In RotStar::christo" << endl;
  //GYOTO_ERROR("RotStar3_1::christoffel: Not implemented yet");

  //24/05/10: stationary axisymmetric version

  double rr=coord[1],r2=rr*rr,th=coord[2],sinth2=sin(th)*sin(th),ph=coord[3];

  const Scalar & NNscal=star_ -> get_nn();
  double NN=NNscal.val_point(rr,th,ph), N2=NN*NN;
  const Scalar & NNrscal = NNscal.dsdr();
  const Scalar & NNthscal = NNscal.dsdt();
  double N_r=NNrscal.val_point(rr,th,ph);
  double N_th=NNthscal.val_point(rr,th,ph);
  const Scalar & omega_scal=star_ -> get_nphi();
  double omega=omega_scal.val_point(rr,th,ph), omega2=omega*omega;
  double omega_r=omega_scal.dsdr().val_point(rr,th,ph);
  double omega_th=omega_scal.dsdt().val_point(rr,th,ph);
  const Scalar & A2scal=star_ -> get_a_car();
  const Scalar & B2scal=star_ -> get_b_car();
  const Scalar & A2rscal = A2scal.dsdr();
  const Scalar & A2thscal = A2scal.dsdt(); 
  const Scalar & B2rscal = B2scal.dsdr();
  const Scalar & B2thscal = B2scal.dsdt(); 
  double B2=B2scal.val_point(rr,th,ph);
  double A2=A2scal.val_point(rr,th,ph);
  double A2_r=A2rscal.val_point(rr,th,ph);
  double A2_th=A2thscal.val_point(rr,th,ph);
  double B2_r=B2rscal.val_point(rr,th,ph);
  double B2_th=B2thscal.val_point(rr,th,ph);

  double gtt=-1./N2, grr=1./A2, gthth=1./(A2*r2), gpp=1./(B2*r2*sinth2)-omega2/N2, gtp=-omega/N2;
  double g_ttr=-2.*NN*N_r+B2_r*omega2*r2*sinth2+2.*omega*omega_r*B2*r2*sinth2+2.*rr*B2*omega2*sinth2, g_ttth=-2.*NN*N_th+B2_th*omega2*r2*sinth2+2.*omega*omega_th*B2*r2*sinth2+2.*cos(th)*sin(th)*r2*B2*omega2;
  double g_rrr=A2_r, g_rrth=A2_th;
  double g_ththr=r2*A2_r+2.*rr*A2, g_ththth=r2*A2_th;
  double g_ppr=sinth2*(r2*B2_r+2.*rr*B2), g_ppth=r2*(sinth2*B2_th+2.*cos(th)*sin(th)*B2);
  double g_tpr=-omega_r*B2*r2*sinth2-omega*B2_r*r2*sinth2-2.*rr*omega*B2*sinth2, g_tpth=-omega_th*B2*r2*sinth2-omega*B2_th*r2*sinth2-2.*cos(th)*sin(th)*omega*B2*r2;

  //  if (debug()) cout << "at r t p= " << rr << " " << th << " " << ph << endl;
  //if (debug()) cout << "gthth, g_ppr= " << gthth << " " << g_ppr << endl;

  //32 non-0 christo
    if ((alpha==0 && mu==0 && nu==1) || (alpha==0 && mu==1 && nu==0)) 
      return 1./2.*gtt*g_ttr+1./2.*gtp*g_tpr;
    if ((alpha==0 && mu==0 && nu==2) || (alpha==0 && mu==2 && nu==0)) 
      return 1./2.*gtt*g_ttth+1./2.*gtp*g_tpth;
    if ((alpha==0 && mu==3 && nu==1) || (alpha==0 && mu==1 && nu==3)) 
      return 1./2.*gtt*g_tpr+1./2.*gtp*g_ppr;
    if ((alpha==0 && mu==3 && nu==2) || (alpha==0 && mu==2 && nu==3)) 
      return 1./2.*gtt*g_tpth+1./2.*gtp*g_ppth;
    if ((alpha==3 && mu==3 && nu==1) || (alpha==3 && mu==1 && nu==3)) 
      return 1./2.*gpp*g_ppr+1./2.*gtp*g_tpr;
    if ((alpha==3 && mu==3 && nu==2) || (alpha==3 && mu==2 && nu==3)) 
      return 1./2.*gpp*g_ppth+1./2.*gtp*g_tpth;
    if ((alpha==3 && mu==0 && nu==1) || (alpha==3 && mu==1 && nu==0)) 
      return 1./2.*gpp*g_tpr+1./2.*gtp*g_ttr;
    if ((alpha==3 && mu==0 && nu==2) || (alpha==3 && mu==2 && nu==0)) 
      return 1./2.*gpp*g_tpth+1./2.*gtp*g_ttth;
    if (alpha==1 && mu==2 && nu==2) 
      return -1./2.*grr*g_ththr;
    if (alpha==1 && mu==3 && nu==3) 
      return -1./2.*grr*g_ppr;
    if (alpha==1 && mu==0 && nu==0) 
      return -1./2.*grr*g_ttr;
    if (alpha==1 && mu==1 && nu==1) 
      return 1./2.*grr*g_rrr;
    if ((alpha==1 && mu==1 && nu==2) || (alpha==1 && mu==2 && nu==1)) 
      return 1./2.*grr*g_rrth;
    if ((alpha==1 && mu==0 && nu==3) || (alpha==1 && mu==3 && nu==0)) 
      return -1./2.*grr*g_tpr;
    if (alpha==2 && mu==0 && nu==0) 
      return -1./2.*gthth*g_ttth;
    if (alpha==2 && mu==1 && nu==1) 
      return -1./2.*gthth*g_rrth;
    if (alpha==2 && mu==3 && nu==3) 
      return -1./2.*gthth*g_ppth;
    if (alpha==2 && mu==2 && nu==2) 
      return 1./2.*gthth*g_ththth;
    if ((alpha==2 && mu==2 && nu==1) || (alpha==2 && mu==1 && nu==2)) 
      return 1./2.*gthth*g_ththr;
    if ((alpha==2 && mu==0 && nu==3) || (alpha==2 && mu==3 && nu==0)) 
      return -1./2.*gthth*g_tpth;

    return 0.;//all other christo are 0

}

double RotStar3_1::ScalarProd(const double pos[4],
			  const double u1[4], const double u2[4]) const {
  //cout << "in RotStar ScalarProd" << endl;
  if (debug()) 
    cout << "u1,u2 in Scal= " ;
  for (int ii=0;ii<4;ii++) {
    if (debug()) 
      cout << u1[ii] << " " << u2[ii] << " ";
  }
  if (debug()) 
    cout << endl;
  double g_tt=gmunu(pos,0,0), g_tp=gmunu(pos,0,3), g_rr=gmunu(pos,1,1), g_thth=gmunu(pos,2,2), g_pp=gmunu(pos,3,3);
  //if (debug()) 
 
  //cout << "metrics in rotstar scalar prod= " << g_tt << " " << g_rr << " " << g_thth << " " << g_pp << " " << g_tp << endl;
  //cout << "norm sub-terms in rotstar scalar prod= " << g_tt*u1[0]*u2[0] << " " <<g_rr*u1[1]*u2[1]<< " " <<g_thth*u1[2]*u2[2]<< " " <<g_pp*u1[3]*u2[3]<< " " <<g_tp*u1[0]*u2[3]<< " " <<g_tp*u1[3]*u2[0]<<endl;
  return g_tt*u1[0]*u2[0]+g_rr*u1[1]*u2[1]+g_thth*u1[2]*u2[2]+g_pp*u1[3]*u2[3]+g_tp*u1[0]*u2[3]+g_tp*u1[3]*u2[0];
}

int RotStar3_1::setParameter(string name, string content, string unit){
  if      (name=="IntegKind") {
    GYOTO_WARNING << "IntegKind is deprecated, please use "
      "<GenericIntegrator/> or <SpecificIntegrator/> instead\n";
    integKind(atoi(content.c_str()));
  } else return Generic::setParameter(name, content, unit);
  return 0;
}
