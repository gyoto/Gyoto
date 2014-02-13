/*
    Copyright (c) 2012 Frederic Vincent, Odele Straub, Thibaut Paumard

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
#include "GyotoPolishDoughnut.h"
#include "GyotoPhoton.h"
#include "GyotoFactoryMessenger.h"
#include "GyotoDefs.h"

#include <iostream>
#include <cmath>
#include <string>
#include <cstdlib>
#include <sstream>

using namespace std;
using namespace Gyoto;
using namespace Gyoto::Astrobj;

#define CST_POLY_INDEX 1.5//polytropic index n (gamma=1+1/n=5/3)
#define CST_POLY_INDEX_M1 0.666666666666666666666666666666666666666667
#define CST_HYDRO_FRAC 0.75//hydrogen fraction
#define CST_Z_1 1.//atomic number
#define CST_Z_2 2.
#define CST_MU_ION 1.2307692307692308375521861 //(4./(1. + 3. * CST_HYDRO_FRAC))
#define CST_MU_ELEC 1.1428571428571427937015414 //(2./(1. + CST_HYDRO_FRAC))

#define GYOTO_C2_CGS 8.98755178736817668096e+20 //c^2 in cgs
#define GYOTO_C2_CGS_M1 1.1126500560536184087938986e-21 // 1./GYOTO_C2_CGS
#define w_tol 1e-9

PolishDoughnut::PolishDoughnut() :
  Standard("PolishDoughnut"),
  gg_(NULL),
  l0_(0.),
  lambda_(0.5),
  W_surface_(0.),
  W_centre_(0.),
  r_cusp_(0.),
  r_centre_(0.),
  central_density_(1.),
  centraltemp_over_virial_(1.),
  beta_(0.),
  use_specific_impact_(0),
  //aa_ and aa2_ are set by setLambda()
  spectral_oversampling_(10),
//intersection doesn't need initilizing
  komissarov_(0)
{  
#ifdef GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << endl;
#endif
  critical_value_=0.; safety_value_=.1; //rmax_=25.;
}

PolishDoughnut::PolishDoughnut(const PolishDoughnut& orig) :
  Standard(orig),
  gg_(NULL),
  l0_(orig.l0_),
  lambda_(orig.lambda_),
  W_surface_(orig.W_surface_),
  W_centre_(orig.W_centre_),
  r_cusp_(orig.r_cusp_),
  r_centre_(orig.r_centre_),
  DeltaWm1_(orig.DeltaWm1_),
  central_density_(orig.central_density_),
  centraltemp_over_virial_(orig.centraltemp_over_virial_),
  beta_(orig.beta_),
  use_specific_impact_(orig.use_specific_impact_),
  aa_(orig.aa_),
  aa2_(orig.aa2_),
		 spectral_oversampling_(orig.spectral_oversampling_),
		 intersection(orig.intersection),
		 komissarov_(orig.komissarov_)
{
  if (orig.gg_()) {
    gg_=orig.gg_->clone();
    Standard::gg_ = gg_;
  }
}
PolishDoughnut* PolishDoughnut::clone() const
{return new PolishDoughnut(*this);}

double PolishDoughnut::getL0() const { return l0_; }
//void   PolishDoughnut::setL0(double l0) { l0_ = l0; }
double PolishDoughnut::getWsurface() const { return W_surface_; }
double PolishDoughnut::getWcentre() const { return W_centre_; }
double PolishDoughnut::getRcusp() const { return r_cusp_; }
double PolishDoughnut::getRcentre() const { return r_centre_; }

double PolishDoughnut::getLambda() const { return lambda_; }
void   PolishDoughnut::setLambda(double lambda) {
  if (!gg_) throwError("Metric but be set before lambda in PolishDoughnut");
  //Computing marginally stable and marginally bound radii:
  lambda_=lambda;  
  aa_=gg_->getSpin(), aa2_=aa_*aa_;

  double  rms = gg_->getRms();//(3. + z2 - pow((3. - z1)*(3. + z1 +
  //                          2.*z2),1./2.));
  double  rmb = gg_->getRmb();//pow(1. + sqrt(1. - aa),2);
 
  // marginally stable & marginally bound keplerian angular momentum
  // lK(rms), lK(rmb): NB: this is Bardeen 74 L/E (2.12-2.13), not L
  // (Polish doughnut l is rescaled)
  double  l_ms
    = (rms*rms - 2.*aa_*sqrt(rms) + aa2_)/(pow(rms,3./2.) - 2.*sqrt(rms) + aa_);
  double  l_mb
    = (rmb*rmb - 2.*aa_*sqrt(rmb) + aa2_)/(pow(rmb,3./2.) - 2.*sqrt(rmb) + aa_);
 
  l0_ = lambda_*(l_mb-l_ms)+l_ms ;//torus angular momentum

  //Computing the potential at the photon position:
  double r1_min = rmb ; 
  double r1_max = rms ;
  double r2_min = rms ; 
  double r2_max = 1000. ;
 
  // update intersection functor:
  intersection.aa_=aa_;
  intersection.aa2_=aa2_;
  intersection.l0_=l0_;

  r_cusp_   = intersection.ridders(r1_min, r1_max) ;
  r_centre_ = intersection.ridders(r2_min, r2_max) ;
  W_surface_ = potential(r_cusp_, M_PI/2.) ;
  W_centre_  = potential(r_centre_, M_PI/2.) ;

  GYOTO_IF_DEBUG
    GYOTO_DEBUG_EXPR(lambda_);
  GYOTO_DEBUG_EXPR(l0_);
  GYOTO_DEBUG_EXPR(aa_);
  GYOTO_DEBUG_EXPR(r_cusp_);
  GYOTO_DEBUG_EXPR(r_centre_);
  GYOTO_DEBUG_EXPR(W_surface_);
  GYOTO_DEBUG_EXPR(W_centre_);
  GYOTO_ENDIF_DEBUG

    DeltaWm1_= 1./(W_centre_ - W_surface_);

}

double PolishDoughnut::getCentralDensity() const {return central_density_;}
double PolishDoughnut::getCentralDensity(string unit) const {
  double dens = getCentralDensity();
  if (unit != "") {
#   ifdef HAVE_UDUNITS
    dens = Units::Converter("kg/L", unit)(dens);
#   else
    GYOTO_WARNING << "Units ignored, please recompile Gyoto with --with-udunits"
		  << endl ;
# endif
  }
  return dens;
}
void   PolishDoughnut::setCentralDensity(double dens) {
  central_density_=dens;
}
void   PolishDoughnut::setCentralDensity(double dens, string unit) {
  if (unit != "") {
# ifdef HAVE_UDUNITS
    dens = Units::Converter(unit, "kg/L")(dens);
# else
    GYOTO_WARNING << "Units ignored, please recompile Gyoto with --with-udunits"
		  << endl ;
# endif
  }
  setCentralDensity(dens);
}

double PolishDoughnut::getCentralTempOverVirial() const
{return centraltemp_over_virial_;}
void   PolishDoughnut::setCentralTempOverVirial(double val)
{centraltemp_over_virial_=val;}

double PolishDoughnut::getBeta() const { return beta_; }
void   PolishDoughnut::setBeta(double beta)   { beta_ = beta; }

size_t PolishDoughnut::getSpectralOversampling() const
{ return spectral_oversampling_; }
void   PolishDoughnut::setSpectralOversampling(size_t val)
{ spectral_oversampling_ = val; }

bool PolishDoughnut::komissarov() const
{return komissarov_;}
void PolishDoughnut::komissarov(bool komis)
{komissarov_=komis;}

PolishDoughnut::~PolishDoughnut() {
  GYOTO_DEBUG << "PolishDoughnut Destruction" << endl;
  if (gg_) gg_ -> unhook(this);
}

Gyoto::SmartPointer<Gyoto::Metric::Generic> PolishDoughnut::getMetric() const {
  return SmartPointer<Metric::Generic>(gg_);
}
void PolishDoughnut::setMetric(Gyoto::SmartPointer<Gyoto::Metric::Generic> met)
{
  if (met->getKind() != "KerrBL")
    throwError("PolishDoughnut::setMetric(): only KerrBL, please");
  if (gg_) gg_ -> unhook(this);
  gg_ = SmartPointer<Metric::KerrBL>(met);
  Generic::gg_ = gg_;
  if (gg_) gg_ -> hook(this);
  GYOTO_DEBUG << "Metric set, calling setLambda\n";
  setLambda(lambda_); // setLambda_ initializes other members
}

void PolishDoughnut::tell(Hook::Teller * met) {
  if (met == gg_) setLambda(lambda_);
  else throwError("BUG: PolishDoughnut::tell(Hook::Teller * met) called with"
		  "wrong metric");
}

int PolishDoughnut::Impact(Photon *ph, size_t index,
			   Astrobj::Properties *data) {
  if (beta_==1.) throwError("Please set beta to != 1.");
  GYOTO_DEBUG_EXPR(use_specific_impact_);
  if (use_specific_impact_)
    return Impact_(ph, index, data);
  return Standard::Impact(ph, index, data);
}

int PolishDoughnut::Impact_(Photon *ph, size_t index,
			    Astrobj::Properties *data) {
  /*
    Should not be used anymore, use Standard::Impact
  */
  GYOTO_DEBUG << endl;

  //  int width=15, prec=12;
  double coord_ph_hit[8], coord_obj_hit[8];
  double dt = 0.1; // seems to be OK if not so little (0.1)
  int hit=0;
  double coord[8], coord2[8];
  ph->getCoord(index, coord2);  // coord2 is the earliest point, t2<t1
  ph->getCoord(index+1, coord); // which must be processed last

  double t1=coord[0], r1=coord[1], theta1=coord[2],
    tprev=t1, rprev=r1, thetaprev=theta1,
    t2=coord2[0], r2=coord2[1], theta2=coord2[2],
    tcur, rcur, thetacur, phicur, tdotcur, rdotcur, thetadotcur, phidotcur;
  // remember : t1 > t2, we integrate backwards in time
  //NB: to keep in mind: inside doughnut 0<=w<=1 with w=0 at the surface, 
  // outside, excluding the funnel w<0, in the funnel w>0 (can be >1) 
  // Funnel = open field lines with projected distance to the axis < r_cusp_
  double wcur, w1 = (potential(r1, theta1) - W_surface_)*DeltaWm1_, wprev = w1 ;
  double w2 = (potential(r2, theta2) - W_surface_)*DeltaWm1_;

  double w_lim=3.;

  if (debug())
    cerr << "DEBUG: PolishDoughnut::Impact():"
	 << "w1=" << w1 << ", w2=" << w2
	 << ", r1*sin(theta1)/r_cusp_="<<r1*sin(theta1)/r_cusp_
	 << ", r2*sin(theta2)/r_cusp_="<<r2*sin(theta2)/r_cusp_
	 << endl;

  //checks whether w > 1 outside funnel (should not happen!)
  if ((w1 > 1. && r1*fabs(sin(theta1)) > r_cusp_)
      || (w2 > 1. && r2*fabs(sin(theta2)) > r_cusp_))
    throwError("I was wrong, w may be > 1 outside the funnel!");
  
  //checks whether both w are far from 0 (surface value) or inside funnel
  //if yes: no impact
  if (   (w1 < 0.-w_lim || r1*fabs(sin(theta1)) < r_cusp_)
	 && (w2 < 0.-w_lim || r2*fabs(sin(theta2)) < r_cusp_) ) {
    return 0;
  }

  coord_obj_hit[5] = 0.;
  coord_obj_hit[6] = 0.;
  double tfirst = t1, tlast=t2, wfirst=w1, wlast=w2;


  // At least one of the ends is close to the doughnut. Refine.  We
  // want t1 inside the object. If entering the object, we want to
  // stay close to the surface.
  if ( r1*fabs(sin(theta1)) < r_cusp_ ) {
    // t1 is in the funnel. So t2 is not (see if loop above).
    // Bring t1 closer to t2, outside funnel.
    tfirst=t1; tlast=t2;
    while ((tfirst-tlast)>GYOTO_T_TOL) {
      tcur = 0.5*(tfirst+tlast);
      ph -> getCoord( &tcur, 1, &rcur, &thetacur, &phicur,
		      &tdotcur, &rdotcur, &thetadotcur, &phidotcur);
      if (rcur*fabs(sin(thetacur))<r_cusp_) tfirst=tcur;
      else tlast=tcur;
    }
    t1=tcur;
    w1 =(potential(rcur, thetacur) - W_surface_)*DeltaWm1_;
    if (debug()) 
      cerr << "DEBUG: PD::Impact: Funnel handling for t1,"
	   << " new w1="<<w1
	   <<", new t1="<<t1
	   << ", new rsinth/rcusp= " << rcur*sin(thetacur)/r_cusp_ <<endl;
  }

  // t1 is outside (not in the funnel) --> bring it inside doughnut
  if ( w1 < 0. ) {
    if (debug())
      cerr << "DEBUG: PD::Impact(): w1<0\n"; 
    // look for entry point at W=0
    wcur=w2; tcur=t2; rcur=r2; thetacur=theta2;
    if ( w2 < 0.) { // t1 and t2 outside, but close to surface --> look for 
                    // a maximum in between, is it inside the object ?
      if (debug()) 
	cout << "DEBUG: PD::Impact():"
	     << "t1 and t2 close to surface"
	     << endl;
      while (wcur < 0.) {
	if ((tfirst-tlast)<GYOTO_T_TOL) return 0;
	tcur = 0.5*(tfirst+tlast);
	ph -> getCoord( &tcur, 1, &rcur, &thetacur, &phicur,
			&tdotcur, &rdotcur, &thetadotcur, &phidotcur);
	wcur = (potential(rcur, thetacur) - W_surface_)*DeltaWm1_;
	if (wlast < wfirst) {
	  tlast=tcur; wlast=wcur;
	} else {
	  tfirst=tcur; wfirst=wcur;
	}
      }
    }
    wlast=wcur; tlast=tcur;

    if (rcur*fabs(sin(thetacur))<r_cusp_) {
      if (debug()) 
	cerr << "DEBUG: PD::Impact():"
	     << "Earliest end of geodesic in funnel: " 
	     << "tcur= " << tcur 
	     << " ,rcur= " << rcur << endl;
      return 0; 
      //In this case t1 is outside, not in the funnel
      //and tcur is in the funnel: just go on with the
      //integration, no impact
    }

    tfirst=t1; wfirst=w1;
    //At this stage, t1 is outside, w1<0, t2 is inside, w2>0
    //Find crossing point
    while ((tfirst-tlast)>GYOTO_T_TOL) {
      tcur = 0.5*(tfirst+tlast);
      ph -> getCoord( &tcur, 1, &rcur, &thetacur, &phicur,
		      &tdotcur, &rdotcur, &thetadotcur, &phidotcur);
      wcur = (potential(rcur, thetacur) - W_surface_)*DeltaWm1_;
      if (debug()) 
	cerr << "DEBUG: PD::Impact: wfirst="<<wfirst
	     <<", wlast="<<wlast
	     <<" wcur="<<wcur<<endl;
      if (wcur < wfirst) throwError("Doughnut fishy heuristic wcur < wfirst");
      /*cout << "In PolishDoughnut::Impact fishy heuristic pb"
	<< " continuing..." << endl;*/
      //Keep in mind: the other heuristic lead to bug 
      //--> check this one if new bug
      if (wcur<0) {tfirst=tcur; wfirst=wcur; }
      else { tlast=tcur; wlast=wcur; }
    }
    // tlast, wlast is inside of the object, close to the surface:
    t1=tlast; w1=wlast;
    if (debug())
      cerr << "DEBUG: PD::Impact(): found entry point: t1="<<t1
	   << ", w1="<<w1<<endl;
  }
  // OK, now t1, w1 is inside the object.  If we are entering, t1 is
  // right "after" entry (latest date inside the object).

  // Same for t2: we look for the earliest date inside the object.
  if ( r2*fabs(sin(theta2)) < r_cusp_ ) {
    // t2 is in the funnel (w2>0). t1 is in the object by now.
    // Bring t2 out of the funnel, just outside doughnut (with w2<0)
    tfirst=t1; tlast=t2;
    while ((tfirst-tlast)>GYOTO_T_TOL) {
      tcur = 0.5*(tfirst+tlast);
      ph -> getCoord( &tcur, 1, &rcur, &thetacur, &phicur,
		      &tdotcur, &rdotcur, &thetadotcur, &phidotcur);
      if (rcur*fabs(sin(thetacur))>r_cusp_) tfirst=tcur;
      else tlast=tcur;
    }
    t2=tcur;
    w2 =(potential(rcur, thetacur) - W_surface_)*DeltaWm1_;
    if (debug()) 
      cerr << "DEBUG: PD::Impact: Funnel handling for t2,"
	   << " new w2="<<w2
	   <<", new t2="<<t2
	   << ", new rsinth/rcusp= " << rcur*sin(thetacur)/r_cusp_ <<endl;
  }

  double funtol=3e-2; // for heuristic special case when t2 near funnel
                      // chosen by hand...
  if (w2 < 0.) {
    // now we want t2 inside the doughnut.
    if (debug())
      cerr << "DEBUG: PD::Impact(): w2<0\n"; 
    tfirst = t1; wfirst=w1;
    tlast = t2; wlast=w2;
    while ((tfirst-tlast)>GYOTO_T_TOL) {
      tcur = 0.5*(tfirst+tlast);
      ph -> getCoord( &tcur, 1, &rcur, &thetacur, &phicur,
		      &tdotcur, &rdotcur, &thetadotcur, &phidotcur);
      wcur = (potential(rcur, thetacur) - W_surface_)*DeltaWm1_;
      if (debug()) 
	cerr << "DEBUG: PD::Impact: wfirst="<<wfirst
	     <<", wlast="<<wlast
	     <<" wcur="<<wcur<<endl;
      if (wcur < wlast && r2*fabs(sin(theta2))/r_cusp_>1+funtol) {
	throwError("Doughnut fishy heuristic wcur < wlast");
	//Condition on r2sinth2 : if t2 was inside funnel, wcur can be < wlast
      }

      if (wcur > 0.) { tfirst=tcur; wfirst=wcur; }
      else { tlast=tcur; wlast=wcur; }
    } 
    t2=tcur; w2=wcur;
  }
  // Phew, now both t1 and t2 are inside the object. Guaranteed.
  // Now we split the line in small bits by interpolation.
  if (debug())
    cerr << "DEBUG: PD::Impact(): ***HIT*** interpolating worldline between "
	 << "t1="<<t1<<" (w1="<<w1<<") & t2="<<t2<<"(w2="<<w2<<")\n";
  hit=1;

  double ut2, gtt, gtph, gphph, Omega;
  //cout << "BEFORE FOR" << endl;

  //NON RADIATIVE TRANSFER CASE
  if (!flag_radtransf_){
    tcur=t1;
    ph -> getCoord( &tcur, 1, &rcur, &thetacur, &phicur,
		    &tdotcur, &rdotcur, &thetadotcur, &phidotcur);

    coord_obj_hit[0] = coord_ph_hit[0] = tcur;
    coord_obj_hit[1] = coord_ph_hit[1] = rcur;
    coord_obj_hit[2] = coord_ph_hit[2] = thetacur;
    coord_obj_hit[3] = coord_ph_hit[3] = phicur;
    coord_ph_hit[4] = tdotcur;	   
    coord_ph_hit[5] = rdotcur;	   
    coord_ph_hit[6] = thetadotcur;
    coord_ph_hit[7] = phidotcur;

    gtt=gg_->gmunu(coord_ph_hit,0,0);
    gtph=gg_->gmunu(coord_ph_hit,0,3); 
    gphph=gg_->gmunu(coord_ph_hit,3,3);
    Omega=-(l0_*gtt+gtph)/(l0_*gtph+gphph);
    ut2=-1./(gtt+2.*gtph*Omega+gphph*Omega*Omega);
    if (ut2 < 0.)
      throwError("In PolishDoughnut.C PROBLEM!!!: ut^2 is negative.");
    coord_obj_hit[4] = sqrt(ut2);
    coord_obj_hit[7] = Omega*sqrt(ut2);

    processHitQuantities(ph, coord_ph_hit, coord_obj_hit, 1., data);
    return hit;
  }

  //RADIATIVE TRANSFER CASE
  for (tprev = t1, tcur=t1-dt;
       tprev>t2;
       tprev=tcur, wprev=wcur, rprev=rcur, thetaprev=thetacur, tcur-=dt) {
    if (tcur < t2) tcur=t2;
    ph -> getCoord( &tcur, 1, &rcur, &thetacur, &phicur,
		    &tdotcur, &rdotcur, &thetadotcur, &phidotcur);
    wcur = (potential(rcur, thetacur) - W_surface_)*DeltaWm1_;
    if (debug())
      cerr << "**********DEBUG: PD::Impact(): tcur="<<tcur<<", wcur="<<wcur
	   << ", rcur="<<rcur<<", thetacur="<<thetacur<<", phicur="<<phicur
	   <<endl;
    if (wcur < 0. || wprev < 0.
	|| rcur*fabs(sin(thetacur)) < r_cusp_ 
	|| rprev*fabs(sin(thetaprev)) < r_cusp_) {
      if (debug())
	cerr << "WARNING: PolishDoughnut::Impact: skipping\n";
      continue;
    }

    coord_obj_hit[0] = coord_ph_hit[0] = tcur;
    coord_obj_hit[1] = coord_ph_hit[1] = rcur;
    coord_obj_hit[2] = coord_ph_hit[2] = thetacur;
    coord_obj_hit[3] = coord_ph_hit[3] = phicur;
    coord_ph_hit[4] = tdotcur;	   
    coord_ph_hit[5] = rdotcur;	   
    coord_ph_hit[6] = thetadotcur;
    coord_ph_hit[7] = phidotcur;

    gtt=gg_->gmunu(coord_ph_hit,0,0);
    gtph=gg_->gmunu(coord_ph_hit,0,3); 
    gphph=gg_->gmunu(coord_ph_hit,3,3);
    Omega=-(l0_*gtt+gtph)/(l0_*gtph+gphph);
    ut2=-1./(gtt+2.*gtph*Omega+gphph*Omega*Omega);
    if (ut2 < 0.)
      throwError("In PolishDoughnut.C PROBLEM!!!: ut^2 is negative.");
    coord_obj_hit[4] = sqrt(ut2);
    coord_obj_hit[7] = Omega*sqrt(ut2);

    /*cout << "coords ph obj at hit= " << endl;
      for (int ii=0;ii<8;ii++)
      cout << coord_ph_hit[ii] << " ";
      cout << endl;
      for (int ii=0;ii<8;ii++)
      cout << coord_obj_hit[ii] << " ";
      cout << endl;*/

    processHitQuantities(ph, coord_ph_hit, coord_obj_hit, tprev-tcur, data);
      
  }
  if (debug()) 
    cerr<< "DEBUG: PD::Impact() returning after HIT" << endl;
  return hit;
}

double PolishDoughnut::operator()(double const coord[4]) {
  // w1 = ((potential(r1, theta1, aa) - W_surface_)
  //      /(W_centre_ - W_surface_));
  //
  // w1 < 0. outside polishdoughnut, anything inside funnel, 0<w<1
  // inside doughnut.
  //
  // so: operator()() < 0. <=> inside PolishDoughnut.

  double tmp =  W_surface_-potential(coord[1], coord[2]);
  double rproj = coord[1] * sin(coord[2]);
  if (rproj<r_cusp_) tmp = fabs(tmp)+(r_cusp_-rproj);
  return tmp;
}

void PolishDoughnut::getVelocity(double const pos[4], double vel[4]) 
{
  double gtt=gg_->gmunu(pos,0,0);
  double gtph=gg_->gmunu(pos,0,3); 
  double gphph=gg_->gmunu(pos,3,3);
  double Omega=-(l0_*gtt+gtph)/(l0_*gtph+gphph);
  double ut2=-1./(gtt+2.*gtph*Omega+gphph*Omega*Omega);
  if (ut2 < 0.) {
    stringstream ss;
    ss << "PolishDoughnut::getVelocity(pos=[";
    for (int i=0; i<3; ++i) ss << pos[i] << ", ";
    ss << pos[3] << "]): ut^2 is negative.";
    throwError(ss.str());
  }
  vel[0] = sqrt(ut2);
  vel[1] = vel[2] = 0.;
  vel[3] = Omega*sqrt(ut2);
}

void PolishDoughnut::integrateEmission
(double * I, double * boundaries,
 size_t * chaninds, size_t nbnu,
 double dsem, double *cph, double *co) const
{
  // The original channels may or may not be contiguous.  We split
  // each original channels into spectral_oversampling_ subchannels.
  // All we know is that each chunk of spectral_oversampling_
  // subchannels are contiguous. Don't try to recover contiguousness
  // in the original channels, it's too hard for now.
  double som1=1./double(spectral_oversampling_);
  size_t onbnu=nbnu*spectral_oversampling_; // number of subchannels
  size_t onbb = onbnu+nbnu; // number of subchannel boundaries : most
			    // are used twice as subchannels are
			    // contiguous in each channel.
  double * Inu = new double[onbnu+1];
  double * bo = new double[onbb];
  size_t * ii = new size_t[2*onbnu]; // two indices for each subchannel
  double dnu;
  size_t k=0;
  for (size_t i=0; i<nbnu; ++i) {
    dnu=(boundaries[chaninds[2*i+1]]-boundaries[chaninds[2*i]])*som1;
    for (size_t j=0; j<spectral_oversampling_; ++j) {
      k=i*spectral_oversampling_+j;
      ii[2*k]=k+i;
      ii[2*k+1]=k+i+1;
      bo[ii[2*k]]=boundaries[chaninds[2*i]]+double(j)*dnu;
    }
    bo[ii[2*(i*spectral_oversampling_+spectral_oversampling_-1)+1]]
      =boundaries[chaninds[2*i+1]];
  }
  emission(Inu, bo, onbb, dsem, cph, co);
  for (size_t i=0; i<nbnu; ++i) {
    I[i]=0.;
    for (size_t j=0; j<spectral_oversampling_; ++j) {
      k=i*spectral_oversampling_+j;
      I[i]+=(Inu[ii[2*k+1]]+Inu[ii[2*k]])*0.5*fabs(bo[ii[2*k+1]]-bo[ii[2*k]]);
    }
  }
  delete [] Inu;
  delete [] bo;
  delete [] ii;
}

double PolishDoughnut::emission(double nu_em, double dsem,
				double *cph, double *co) const
{
  GYOTO_DEBUG << endl;
  double Inu;
  emission(&Inu, &nu_em, 1, dsem, cph, co);
  return Inu;
}

void PolishDoughnut::emission_komissarov(double Inu[], // output
					 double nu_ems[], size_t nbnu, // input
					 double dsem,
					 double coord_ph[8],
					 double coord_obj[8]) const {
  /*
    This is the emission function for the 
    Komissarov model: 
    Komissarov, MNRAS, 368, 993 (2006)
    Only synchrotron radiation is implemented,
    the emissivity coming from:
    Wardzinski & Zdziarski, MNRAS, 314, 183 (2000)
  */

  /* COMPUTING PHYS QUANTITIES */
  double rr = coord_ph[1], theta = coord_ph[2];//NB: rr is units of GM/c^2
  double rcgs = rr * gg_ -> unitLength() * 100.;//rr in cgs
  double r_centre_cgs = r_centre_ * gg_ -> unitLength() * 100.;

  double ww = (potential(rr, theta) - W_surface_)*DeltaWm1_;
  if (ww<=0.){//Will generate nan in computations w must be strictly positive
    if (fabs(ww)<w_tol) {
      if (ww!=0.) ww=fabs(ww);
      else ww=w_tol;//can be the case if w at entrance in doughnut is exactly 0
    }else{
      throwError("In PolishDoughnut::emission() w<0!");
    }
  }

  double Msgr = gg_->getMass()*1e3; // Gyoto speaks in SI --> here cgs

  double enthalpy_c=central_density_; // Warning: central_density_ is here
  // p+rho*c2 (enthalpy), not rho; model is different from std doughnut
  double beta = beta_; // magnetic pressure parameter at torus center

  double g_tt=gg_->gmunu(coord_ph,0,0),
    g_pp=gg_->gmunu(coord_ph,3,3),
    g_tp=gg_->gmunu(coord_ph,0,3),
    LL=g_tp*g_tp-g_tt*g_pp;
  double posc[4]={0.,r_centre_,M_PI/2.,0.};
  double g_ttc=gg_->gmunu(posc,0,0),
    g_ppc=gg_->gmunu(posc,3,3),
    g_tpc=gg_->gmunu(posc,0,3),
    LLc=g_tpc*g_tpc-g_ttc*g_ppc;

  double pc = enthalpy_c*(W_centre_-W_surface_)
    /((CST_POLY_INDEX+1.)*(1+1./beta_)),
    pmc = pc/beta_,
    kappa = pow(enthalpy_c,-CST_POLY_INDEX_M1)*(W_centre_-W_surface_)
    /((CST_POLY_INDEX+1)*(1+1./beta_)),
    kappam = pow(LLc,-CST_POLY_INDEX_M1)/beta_*kappa;

  double enthalpy = enthalpy_c*
    pow(
	ww*
	(kappa+kappam*pow(LLc,CST_POLY_INDEX_M1))
	/(kappa+kappam*pow(LL,CST_POLY_INDEX_M1))
	,CST_POLY_INDEX
	);

  double number_density = (enthalpy-kappa*pow(enthalpy,1.+CST_POLY_INDEX_M1))/(GYOTO_C2_CGS*CST_MU_ELEC*GYOTO_ATOMIC_MASS_UNIT_CGS);
  double number_density_central = (enthalpy_c-kappa*pow(enthalpy_c,1.+CST_POLY_INDEX_M1))/(GYOTO_C2_CGS*CST_MU_ELEC*GYOTO_ATOMIC_MASS_UNIT_CGS);

  //  double gas_pressure = kappa*pow(enthalpy,1.+CST_POLY_INDEX_M1);  

  double magnetic_pressure = kappam*pow(LL,CST_POLY_INDEX_M1)*pow(enthalpy,1.+CST_POLY_INDEX_M1);

  double bphi = sqrt(2*magnetic_pressure/(g_pp+2*l0_*g_tp+l0_*l0_*g_tt));

  double b4vec[4]={bphi*l0_,0,0,bphi}; // B 4-vector in BL frame

  double vel[4]; // 4-velocity of emitter
  const_cast<PolishDoughnut*>(this)->getVelocity(coord_obj, vel);

  double photon_emframe[4]; // photon tgt vector projected in emitter frame
  double b_emframe[4]; // mag. field vector projected in emitter frame
  for (int ii=0;ii<4;ii++){
    photon_emframe[ii]=coord_ph[ii+4]
      +vel[ii]*gg_->ScalarProd(coord_ph,coord_ph+4,vel);
    b_emframe[ii]=b4vec[ii]
      +vel[ii]*gg_->ScalarProd(coord_ph,b4vec,vel);
  }

  double bnorm = gg_->ScalarProd(coord_ph,b_emframe,b_emframe);
  if (bnorm<=0.) throwError("In PolishDoughnut::emission_komissarov"
			   " b_emframe should be spacelike");
  bnorm=sqrt(bnorm);

  double lnorm = gg_->ScalarProd(coord_ph,photon_emframe,photon_emframe);
  if (lnorm<=0.) throwError("In PolishDoughnut::emission_komissarov"
			   " photon_emframe should be spacelike");
  lnorm=sqrt(lnorm);

  double lscalb = gg_->ScalarProd(coord_ph,photon_emframe,b_emframe);
  double theta_mag = acos(lscalb/(lnorm*bnorm)),
    sth = sin(theta_mag), cth = cos(theta_mag);
  
  // The virial temperature at doughnut's centre,
  // derived at energy equipartition from : 3/2*k*T = G*M*mp/r_centre
  double Tvir = 2./3. * GYOTO_G_CGS * Msgr * GYOTO_PROTON_MASS_CGS 
    / (GYOTO_BOLTZMANN_CGS * r_centre_cgs) ;
  // doughnut's central temperature
  /*
    NB: central temperature defined by a fraction of virial temperature
  */
  double T0   = centraltemp_over_virial_*Tvir;

  double kappabis = T0*pow(number_density_central,-CST_POLY_INDEX_M1);
  double T_electron = kappabis*pow(number_density,CST_POLY_INDEX_M1);
  //cout << T_electron << endl;
  if (T_electron < 5e8){
    for (size_t i=0; i<nbnu; ++i) Inu[i]=0.;
    return;
    //for temperatures below, synch is 0, and compton behaves badly
    //Checked: ~1% cases rejected
  }

  double Theta_elec = GYOTO_BOLTZMANN_CGS*T_electron
    /(GYOTO_ELECTRON_MASS_CGS*GYOTO_C2_CGS);

  double nuc = GYOTO_ELEMENTARY_CHARGE_CGS*bnorm
    /(2.*M_PI*GYOTO_ELECTRON_MASS_CGS*GYOTO_C_CGS);
    
  /* SYNCHRO COMPUTATION */

  for (size_t ii=0; ii<nbnu; ++ii){
    // Formuli below are from Wardzinski&Zdziarski 2000:
    double nuem = nu_ems[ii];
    double gamma0=0., chi0=0.;
    if (Theta_elec<=0.08){
      gamma0 = sqrt(1+2.*nuem*Theta_elec/nuc
		    *pow(1.+9.*nuem*Theta_elec*sth*sth/(2.*nuc),-0.3333333333));
      chi0 = sqrt((2.*Theta_elec*(gamma0*gamma0-1.))
		  /(gamma0*(3.*gamma0*gamma0-1.)));
    }else{
      gamma0 = sqrt(1.+pow(4.*nuem*Theta_elec/(3.*nuc*sth),0.6666666666));
      chi0 = sqrt(2.*Theta_elec/(3.*gamma0));
    }
    double tt = sqrt(gamma0*gamma0-1.)*sth,
      nn = nuem*(1.+tt*tt)/(nuc*gamma0);
    double Z0 = pow((tt*exp(1./sqrt(1.+tt*tt)))/(1.+sqrt(1.+tt*tt)),2.*nn);
    double K2 = bessk(2,1./Theta_elec);
    double ne0 = number_density/Theta_elec*gamma0*sqrt(gamma0*gamma0-1.)/K2
      *exp(-gamma0/Theta_elec);
    // this is j_nu synchro:
    double emis_synch = 
      M_PI*GYOTO_ELEMENTARY_CHARGE_CGS*GYOTO_ELEMENTARY_CHARGE_CGS
      /(2.*GYOTO_C_CGS)*sqrt(nuc*nuem)*chi0*ne0
      *(1.+cth*cth/(sth*sth*gamma0*gamma0))
      *pow(1.-(1.-1./(gamma0*gamma0))*cth*cth,0.25)
      *Z0;

    if (emis_synch!=emis_synch) {
      throwError("In PolishDoughnut::emission_komissarov: "
		 "emissivity is nan");
    }
    if (emis_synch==emis_synch+1.) 
      throwError("In PolishDoughnut::emission_komissarov "
		 "emissivity is infinite");

    Inu[ii]=
      emis_synch*dsem*GYOTO_G_CGS*Msgr*GYOTO_C2_CGS_M1 * GYOTO_INU_CGS_TO_SI;
  }
}

void PolishDoughnut::emission(double Inu[], // output
			      double nu_ems[], size_t nbnu, // input
			      double dsem,
			      double coord_ph[8],
			      double coord_obj[8]) const {

  // Beware: all computations are done in cgs, output must be in SI
  GYOTO_DEBUG << "entering emission()\n";

  if (!flag_radtransf_) {//NON RADIATIVE TRANSFER CASE
    for (size_t i=0; i<nbnu; ++i) Inu[i]=1.; //assumes cst I_nu = 1
    return;
  }

  /*NB: to obtain a trivial rad transfer image with constant emissivity:
    put here 'return dsem;' (not 'return 1;')*/
  //return dsem;

  if (komissarov_){
    return emission_komissarov(Inu,nu_ems,nbnu,dsem,coord_ph,coord_obj);
  }

  double emisstot=0., emiss_synch=0., emiss_brems=0.,
    emiss_Cbrems=0., emiss_Csynch=0.;
  int usesynch=1, usebrems=1, usecompton=1;
  //usesynch, usebrems: obvious meaning
  if (usecompton && (!usesynch || !usebrems))
    throwError("In PolishDoughnut::emission() Compton can't be computed"
	       "without Synchrotron"
	       " and Bremsstrahlung");

  /* ***DOUGHNUT PHYSICAL QUANTITIES COMPUTATION*** */

  double vel[4];
  const_cast<PolishDoughnut*>(this)->getVelocity(coord_obj, vel);
  //double Omega=vel[3]/vel[0];

  double Msgr = gg_->getMass()*1e3; // Gyoto speaks in SI --> here we
  // switch to cgs units
  double rr = coord_ph[1], theta = coord_ph[2];//NB: rr is units of GM/c^2

  double rcgs = rr * gg_ -> unitLength() * 100.;//rr in cgs
  //r_centre_ in cgs:
  double r_centre_cgs = r_centre_ * gg_ -> unitLength() * 100.;
  double ww = (potential(rr, theta) - W_surface_)*DeltaWm1_;

  if (ww<=0.){//Will generate nan in computations w must be strictly positive
    if (fabs(ww)<w_tol) {
      if (ww!=0.) ww=fabs(ww);
      else ww=w_tol;//can be the case if w at entrance in doughnut is exactly 0
    }else{
      throwError("In PolishDoughnut::emission() w<0!");
    }
  }

  // The virial temperature at doughnut's centre,
  // derived at energy equipartition from : 3/2*k*T = G*M*mp/r_centre
  double Tvir = 2./3. * GYOTO_G_CGS * Msgr * GYOTO_PROTON_MASS_CGS 
    / (GYOTO_BOLTZMANN_CGS * r_centre_cgs) ;
  // doughnut's central temperature
  /*
    NB: central temperature defined by a fraction of virial temperature
   */
  double T0   = centraltemp_over_virial_*Tvir;

  GYOTO_DEBUG << "T0="<< T0 << endl;

  // gas mass density [g cm^-3]
  double beta = beta_; // magnetic pressure parameter
  GYOTO_DEBUG << "beta=" << beta << endl;

  GYOTO_DEBUG << "kappa=";
  double kappa = (exp((W_centre_-W_surface_)/(CST_POLY_INDEX+1.))-1.)
    *pow(central_density_*GYOTO_C2_CGS,-CST_POLY_INDEX_M1);
  if (debug()) cerr << kappa << endl;

  /*
    NB: 'central_density_' is given in g/cm3 in XML (mass density)
        'density' is as well a mass density (see the 1/c2 factor)
   */
  GYOTO_DEBUG << "density=";
  double density = GYOTO_C2_CGS_M1
    *pow(1./kappa
	 *(pow(1.
	       +kappa*pow(central_density_*GYOTO_C2_CGS,CST_POLY_INDEX_M1)
	       ,ww)
	   -1.)
	 ,CST_POLY_INDEX);
  if (debug()) cerr << density << endl;

  GYOTO_DEBUG_EXPR(central_density_);
  GYOTO_DEBUG_EXPR(GYOTO_C2_CGS);
  GYOTO_DEBUG_EXPR(CST_POLY_INDEX);
  GYOTO_DEBUG_EXPR(ww);

  // gas number density [cm^-3] --> 'density/mass' is a number density
  double n_e = density/(CST_MU_ELEC * GYOTO_ATOMIC_MASS_UNIT_CGS) ;
  double n_i = density/(CST_MU_ION * GYOTO_ATOMIC_MASS_UNIT_CGS) ;  
  double n_j = (CST_Z_1*CST_Z_1 * CST_HYDRO_FRAC) * n_i
    + (CST_Z_2*CST_Z_2 * (1.-CST_HYDRO_FRAC)) * n_i ;
  //adjusted for H and He ion species ; n_j = n_i = n_e if
  //CST_HYDRO_FRAC=1 (only protons in disk)
      
  // pressure
  double PP = kappa*pow(density*GYOTO_C2_CGS,1.+CST_POLY_INDEX_M1);

  // electron temperature
  /* 
     NB: polish doughnut cannot be a perfect gas.
     Following definition of temperature is as close as possible
     to perfect gas law, but it is not perfect gas law.
     T = cst * p/rho, cst defined from chosen central temperature;
     perfect gas: T = cst' * p/rho, cst' is defined from cst of nature.
   */
  double kappabis = T0*pow(central_density_,-CST_POLY_INDEX_M1);
  double T_electron = kappabis*pow(density,CST_POLY_INDEX_M1);
  GYOTO_DEBUG << "T_electron=" << T_electron << endl;

  // dimensionless electron temp
  double temp_e     = GYOTO_BOLTZMANN_CGS * T_electron 
    / (GYOTO_ELECTRON_MASS_CGS*GYOTO_C2_CGS) ;  
  //Frequency of maximum BB emission (see Rybicki-Lightman)
  double numax = 2.82*GYOTO_BOLTZMANN_CGS*T_electron/GYOTO_PLANCK_CGS;

  double BB = sqrt(24.*M_PI*beta*PP); // Pmagn = B^2/24pi
  double nu_0    = GYOTO_ELEMENTARY_CHARGE_CGS* BB 
    / (2. * M_PI * GYOTO_ELECTRON_MASS_CGS * GYOTO_C_CGS) ;
  if (T_electron < 5e8){
    for (size_t i=0; i<nbnu; ++i) Inu[i]=0.;
    return;
    //for temperatures below, synch is 0, and compton behaves badly
    //Checked: ~1% cases rejected
  }

  /* ***PRELIMINARY COMPUTATIONS OF EMISSION-LINKED QUANTITIES*** */

  //Prefactor for synch emission, see NY95 Eq 3.9
  //preff is 1/4pi times NY95 equivalent factor to have a result in str^-1
  double preff=n_e*GYOTO_ELEMENTARY_CHARGE_CGS
    *GYOTO_ELEMENTARY_CHARGE_CGS
    /(sqrt(3)*GYOTO_C_CGS*bessk(2,1./temp_e));

  double amplification=1., Csynch=1., Cbrems=1.;
  
  //See Mahadevan 96 Table 1
  double temp1,temp2;
  double alpha1=1., alpha2=1., alpha3=1.;
  if (T_electron>3.2e10){ //see Mahadevan 96 Table 1
    alpha1 = 1.;
    alpha2 = 1.;
    alpha3 = 1.;
  }else if (T_electron<=(temp2=3.2e10) && T_electron>(temp1=1.6e10)){
    alpha1 = 0.9768+(0.9768-0.9788)/(temp1-temp2)*(T_electron-temp1);
    alpha2 = 1.095+(1.095-1.021)/(temp1-temp2)*(T_electron-temp1);
    alpha3 = 0.8332+(0.8332-1.031)/(temp1-temp2)*(T_electron-temp1);
  }else if (T_electron<=(temp2=1.6e10) && T_electron>(temp1=8e9)){
    alpha1 = 0.9774+(0.9774-0.9768)/(temp1-temp2)*(T_electron-temp1);
    alpha2 = 1.16+(1.16-1.095)/(temp1-temp2)*(T_electron-temp1);
    alpha3 = 0.2641+(0.2641-0.8332)/(temp1-temp2)*(T_electron-temp1);
  }else if (T_electron<=(temp2=8e9) && T_electron>(temp1=4e9)){
    alpha1 = 1.045+(1.045-0.9774)/(temp1-temp2)*(T_electron-temp1);
    alpha2 = -0.1897+(-0.1897-1.16)/(temp1-temp2)*(T_electron-temp1);
    alpha3 = 0.0595+(0.0595-0.2641)/(temp1-temp2)*(T_electron-temp1);
  }else if (T_electron<=(temp2=4e9) && T_electron>(temp1=2e9)){
    alpha1 = 1.18+(1.18-1.045)/(temp1-temp2)*(T_electron-temp1);
    alpha2 = -4.008+(-4.008+0.1897)/(temp1-temp2)*(T_electron-temp1);
    alpha3 = 1.559+(1.559-0.0595)/(temp1-temp2)*(T_electron-temp1);
  }else if (T_electron<=(temp2=2e9) && T_electron>(temp1=1e9)){
    alpha1 = 1.121+(1.121-1.18)/(temp1-temp2)*(T_electron-temp1);
    alpha2 = -10.65+(-10.65+4.008)/(temp1-temp2)*(T_electron-temp1);
    alpha3 = 9.169+(9.169-1.559)/(temp1-temp2)*(T_electron-temp1);
  }else if (T_electron<=(temp2=1e9) && T_electron>(temp1=5e8)){
    alpha1 = 0.0431+(0.0431-1.121)/(temp1-temp2)*(T_electron-temp1);
    alpha2 = 10.44+(10.44+10.65)/(temp1-temp2)*(T_electron-temp1);
    alpha3 = 16.61+(16.61-9.169)/(temp1-temp2)*(T_electron-temp1);
  }else{ 
    GYOTO_DEBUG << "Too low temperature for synchrotron emission" << endl;
    //see Mahadevan 96, no fit given for T_e<5e8K
  }

  double param[7]={rcgs,n_e,BB,T_electron,alpha1,alpha2,alpha3};
  double nu_test = 1e17;
  //NB: see mma Transcendental.nb, nu_crit is well below 1e17
  transcendental_t transcendental;
  transcendental.par = param;
  double xM_crit = transcendental.secant(50, 5000) ;
  if (transcendental.status) { // maxiter reached in secant method
    double xmin = 2./3.*nu_test/(nu_0*temp_e*temp_e), xmax;
    while (transcendental(xmin)<0.){
      xmin*=0.1;
    }
    xmax=10.*xmin;
    while (transcendental(xmax)>0.){
      xmax*=10.;
    }
    xM_crit = transcendental.ridders(xmin, xmax) ;
  }
  //  cout << "xmin="<<xmin<<", xmax="<<xmax<<", xM_crit="<<xM_crit<<endl;
  //  cout << "secant: xM_crit="<<xM_crit<<endl;
  double nu_crit = 3./2. * nu_0 * temp_e * temp_e * xM_crit ;
  
  if (usecompton){
    // nu_em-independent part
    /*  DOUGHNUT HEIGHT */

    //compute H(r), doughnut's height at given value of r
    int inside=1;//will be 0 when getting out of doughnut
    double newr=rr,newth=theta,rsth=rr*sin(theta),neww,wbef=ww,rbef=rr;
    double dr=1e-2;//crude, but there will be a linear interpolation later
    //if (debug()) cout << "***ww= " << ww << endl;

    // TO CHANGE!!!! put in the non-nu part above
    while (inside){
      newr+=dr;
      newth=asin(rsth/newr);//NB: no pb with asin, even if it's
			    //defined between -pi/2 and pi/2 and theta
			    //is between 0 and pi; we just want one of
			    //the 2 points at the surface of the
			    //doughnut with same value of x and y as
			    //the rr,theta point and one of them will
			    //be found by this asin (which one we
			    //don't care)
      neww=(potential(newr, newth) - W_surface_)*DeltaWm1_;
      if (neww>1. || neww<0.) inside=0;
      else {wbef=neww;rbef=newr;}
    }
    double lina=(neww-wbef)/(newr-rbef), linb=neww-lina*newr;
    //linear interpolation to find surface
    if (neww>1.) newr=(1-linb)/lina;
    else newr=-linb/lina;
    newth=asin(rsth/newr);
    double Hofr = newr*cos(newth);

    /* COMPTON ENHANCEMENT */

    double x_crit        = GYOTO_PLANCK_CGS * nu_crit 
      / (GYOTO_ELECTRON_MASS_CGS * GYOTO_C2_CGS) ;
    double tau_es        = 2.*(n_e * GYOTO_THOMSON_CGS)
      *Hofr*GYOTO_G_CGS * Msgr*GYOTO_C2_CGS_M1;
    // See Narayan&Yi 95 eq. 2.15
    double probability   = 1. - exp(-tau_es) ; 
    amplification = 1. + 4. * temp_e + 16. * temp_e * temp_e ;
	
    double eta1          = probability * (amplification - 1.)
      /(1. - probability * amplification) ;
    double eta3          = -1. - log(probability) / log(amplification) ;
    double eta2          = pow(1./3.,eta3) * eta1 ;
    
    //Formula for NY95 Eq. 3.23 corrected
    Cbrems  = 3. * eta1 * temp_e
      * ((1./3. - x_crit/(3.*temp_e)) - 1./(eta3 + 1.)
	 * (pow(1./3.,eta3 + 1.) - pow(x_crit/(3.*temp_e),eta3 + 1.))) ;
    //NY95 Eq. 3.24
    Csynch  = (eta1 - eta2 * pow(x_crit/temp_e,eta3)) ;

  }

  double nu_em;
  for (size_t i=0; i<nbnu; ++i) {
    nu_em=nu_ems[i];
    /* ***SYNCHROTRON COMPUTATION*** */
  
    if (usesynch)
      emiss_synch=emissionSynch(nu_em,nu_crit,numax,nu_0,
				T_electron,1.,1.,
				alpha1,alpha2,alpha3,
				preff,
				0); 
  
    /* ***BREMSSTRAHLUNG COMPUTATION*** */ 

    if (usebrems)
      emiss_brems = emissionBrems(nu_em, nu_crit, numax, T_electron,
				  n_e, n_j,
				  1., 1., 0);
  
    /* ***COMPTONIZATION COMPUTATION*** */ 

    if (usecompton){
      // nu_em-dependent part

      /* COMPTON BREMSSTRAHLUNG */

      emiss_Cbrems=emissionBrems(nu_em,nu_crit,numax,T_electron,
				 n_e, n_j,
				 amplification,Cbrems,1);//first order
      emiss_Cbrems+=emissionBrems(nu_em,nu_crit,numax,T_electron,
				  n_e, n_j,
				  amplification,Cbrems,2);//second order

      /* COMPTON SYNCHROTRON */

      emiss_Csynch=emissionSynch(nu_em,nu_crit,numax,nu_0,
				 T_electron,amplification, Csynch,
				 alpha1,alpha2,alpha3,preff,1);//first order
      emiss_Csynch+=emissionSynch(nu_em,nu_crit,numax,nu_0,
				  T_electron,amplification, Csynch,
				  alpha1,alpha2,alpha3,preff,2);//second order

      if (emiss_Cbrems<0. || emiss_Csynch<0. || emiss_synch<0. || emiss_brems<0.)
	throwError("In PolishDoughnut::emission emissivity<0!!");
  
    }//end compton
    //  cout << "emiss= " << emiss_synch << " " << emiss_brems<< " " << emiss_Csynch << " " << emiss_Cbrems << endl;
    emisstot=emiss_synch+emiss_brems+emiss_Csynch+emiss_Cbrems;
    GYOTO_DEBUG << "emiss_synch: " << emiss_synch
		<< ", emiss_brems: " << emiss_brems
		<< ", emiss_Csynch: " << emiss_synch
		<< ", emiss_Cbrems: " << emiss_brems
		<< ", emisstot: " << emisstot 
		<< endl;

    if (emisstot!=emisstot) throwError("In PolishDoughnut.C: emissivity is nan");
    if (emisstot==emisstot+1.) throwError("In PolishDoughnut.C: "
					  "emissivity is infinite");

    Inu[i]=emisstot*dsem*GYOTO_G_CGS*Msgr*GYOTO_C2_CGS_M1 * GYOTO_INU_CGS_TO_SI;
    // returned is j_nu*dsem, homogeneous to I_nu
    // Remember PolishDoughnut thinks in c.g.s., output *must* be SI

    GYOTO_DEBUG << "i="<< i
		<< ", nu_em[i]=" << nu_em
		<< ", Inu[i]=" << Inu[i]
		<< endl;
  }
}

double PolishDoughnut::emissionBrems(double nu_em, double nu_crit, 
				     double numax, double T_electron,
				     double n_e, double n_j,
				     double amplification,
				     double Cbrems,
				     int comptonorder) const{

  double amplinu=nu_em;
  if (comptonorder>0){
    amplinu/=pow(amplification,comptonorder);
    Cbrems=pow(Cbrems,comptonorder);
  }else if (Cbrems!=1.)
    throwError("In PolishDoughnut::emissionBrems: Cbrems should be 1"
	       "if no Compton amplification");
  
  /*
    Cases: 
    nu_em<nu_crit -> 0
    nu_crit<nu_em<numax -> shifted brems (NB not shifted if 
    comptonorder=0, in this case amplinu=nu)
    nu_em>numax -> wien tail smoothly connected to previous brems
  */

  double temp_e = GYOTO_BOLTZMANN_CGS * T_electron 
    / (GYOTO_ELECTRON_MASS_CGS*GYOTO_C2_CGS) ;  
  double Fee, Fei;
  
  if (temp_e < 1.) {
    Fee = 20./(9.*sqrt(M_PI))*(44.- 3.*M_PI*M_PI)*pow(temp_e,3./2.)
      *(1.+1.1*temp_e+temp_e*temp_e - 1.25*pow(temp_e,5./2.));
    Fei = 4.*sqrt(2.*temp_e / (M_PI*M_PI*M_PI))*(1.+1.781*pow(temp_e,1.34));
  }else{
    Fee = 24.*temp_e*(log(2.*exp(-GYOTO_EULER_MASCHERONI)*temp_e)+1.28);
    Fei = 9.*temp_e / (2.*M_PI)*(log(1.123*temp_e + 0.48) + 1.5);
  }
  
  double fee = n_e*n_e*GYOTO_C_CGS*GYOTO_ELECTRON_CLASSICAL_RADIUS_CGS
    *GYOTO_ELECTRON_CLASSICAL_RADIUS_CGS
    *GYOTO_ELECTRON_MASS_CGS*GYOTO_C2_CGS*GYOTO_ALPHA_F*Fee;
  double fei = n_e*n_j*GYOTO_THOMSON_CGS*GYOTO_C_CGS
    *GYOTO_ALPHA_F*GYOTO_ELECTRON_MASS_CGS
    *GYOTO_C2_CGS*Fei;
  double f_brems = fee + fei ;
  
  if (nu_em > nu_crit){
    if (nu_em < numax) {
      /* Brems at shifted freq between nu_crit and numax */
      double Gaunt; // Gaunt factor at shifted frequency
      if ( GYOTO_BOLTZMANN_CGS*T_electron/
	   (GYOTO_PLANCK_CGS*nu_em) < 1.){
	Gaunt = sqrt(
		     3./M_PI*GYOTO_BOLTZMANN_CGS*T_electron
		     /(GYOTO_PLANCK_CGS*amplinu)
		     ) ;
      }else{
	Gaunt = sqrt(3.) / M_PI 
	  * log(
		4./GYOTO_EULER_MASCHERONI * GYOTO_BOLTZMANN_CGS*T_electron
		/(GYOTO_PLANCK_CGS*amplinu)
		) ;
      }
      //Brems emission at shifted freq:
      double emissbrems = 1./(4.*M_PI)*f_brems*Gaunt
	*exp(-GYOTO_PLANCK_CGS*amplinu/(GYOTO_BOLTZMANN_CGS*T_electron)) 
	* GYOTO_PLANCK_CGS/(GYOTO_BOLTZMANN_CGS*T_electron);
      //NB: the 1/4pi factor is just to have the per steradian
      //dependency, the emission being assumed isotropic in
      //emitter's frame
      return emissbrems*Cbrems;
    } else {
      /* Wien tail above numax */
      double MaxGaunt; // Gaunt factor with max frequency
      if ( GYOTO_BOLTZMANN_CGS*T_electron/
	   (GYOTO_PLANCK_CGS*numax) < 1.){
	MaxGaunt = sqrt(
			3./M_PI*GYOTO_BOLTZMANN_CGS*T_electron
			/(GYOTO_PLANCK_CGS*numax)
			) ;
      }else{
	MaxGaunt = sqrt(3.) / M_PI 
	  * log(
		4./GYOTO_EULER_MASCHERONI * GYOTO_BOLTZMANN_CGS*T_electron
		/(GYOTO_PLANCK_CGS*numax)
		) ;
      }
      //Wien tail smoothly connected to brems emission below numax:
      double WienEm=BBapprox(nu_em,T_electron);
      double WienEmMax = BBapprox(numax,T_electron);
      double emissbremsmax = 1./(4.*M_PI)*f_brems*MaxGaunt
	*exp(-GYOTO_PLANCK_CGS*numax/(GYOTO_BOLTZMANN_CGS*T_electron)) 
	* GYOTO_PLANCK_CGS/(GYOTO_BOLTZMANN_CGS*T_electron);
      double emiss_bremsmax = Cbrems*emissbremsmax;
      double fact = WienEmMax/emiss_bremsmax;
      
      return WienEm/fact;
    }
  }else{
    return 0.;
  }
}

double PolishDoughnut::emissionSynch(double nu_em, double nu_crit,
				     double numax, double nu_0,
				     double T_electron,
				     double amplification,
				     double Csynch, 
				     double alpha1, double alpha2,
				     double alpha3, double preff,
				     int comptonorder) const{

  double amplinu=nu_em;
  if (comptonorder>0){
    amplinu/=pow(amplification,comptonorder);
    Csynch=pow(Csynch,comptonorder);
  }else if (Csynch!=1.)
    throwError("In PolishDoughnut::emissionSynch: Csynch should be 1"
	       "if no Compton amplification");

  /*
    Cases:
    nu_em<nu_crit -> Rayleigh-Jeans emission smoothly connected to 
    following synch emission
    nu_crit<nu_em<numax -> shifted synchrotron (NB not shifted if 
    comptonorder=0, in this case amplinu=nu)
    nu_em>numax -> 0
  */
  
  double temp_e     = GYOTO_BOLTZMANN_CGS * T_electron 
    / (GYOTO_ELECTRON_MASS_CGS*GYOTO_C2_CGS) ;  
  if (nu_em<nu_crit){
    //Rayleigh-Jeans emission below nu_crit, smoothly connected 
    //to normal comptonized synch above nu_crit 
    double RJEm = BBapprox(amplinu,T_electron);
    double RJEmCrit = BBapprox(nu_crit,T_electron);
    double pref=preff*nu_crit;
    double xM=2.*nu_crit/(3.*nu_0*temp_e*temp_e);
    double func_xM = funcxM(alpha1,alpha2,alpha3,xM);
    double emisssynch_crit = pref*func_xM;
    double fact = RJEmCrit/(Csynch*emisssynch_crit);
    return RJEm/fact;
  }else{
    if (nu_em<numax) {
      //Synchrotron emission at shifted freq. between nu_crit and numax
      double pref = preff*amplinu;
      double xM = 2.*amplinu/(3.*nu_0*temp_e*temp_e);
      double func_xM = funcxM(alpha1,alpha2,alpha3,xM);
      double emisssynch_ampli=pref*func_xM;
      return Csynch*emisssynch_ampli;
    }
    else return 0.; // 0 above numax
  }
  
}

double PolishDoughnut::transmission(double , double , 
				    double *) const {

  if (!flag_radtransf_) return 0.; //Complete absorption for optically thick

  return 1.;//NO ABSORPTION FOR OPTICALLY THIN

  //Following function returns abs=jnu/Bnu
  //To check: is it always 0?

  /* double Msgr = gg_->getMass()*1e3; // Gyoto speaks in SI --> here we */
  /* 			             // switch to cgs units */
  /* double rr = coord_ph[1], theta = coord_ph[2];//NB: rr is units of GM/c^2 */
  /* //  double rcgs = rr * gg_ -> unitLength() * 100.;//rr in cgs */
  /* //r_centre_ in cgs: */
  /* double r_centre_cgs = r_centre_ * gg_ -> unitLength() * 100.; */
  /* double ww = (potential(rr, theta) - W_surface_)*DeltaWm1_; */
  /* if (ww<=0.){//Will generate nan in computations w must be strictly positive */
  /*   if (fabs(ww)<w_tol) { */
  /*     if (ww!=0.) ww=fabs(ww); */
  /*     else ww=w_tol;//can be the case if w at entrance in doughnut is exactly 0 */
  /*   }else{ */
  /*     throwError("In PolishDoughnut::transmission() w<0!"); */
  /*   } */
  /* } */
  /* double mycst_e0=central_density_; */
  /* double Tvir = 2./3. * GYOTO_G_CGS * Msgr * GYOTO_PROTON_MASS_CGS  */
  /*   / (GYOTO_BOLTZMANN_CGS * r_centre_cgs) ; */
  /* // doughnut's central temperature */
  /* double T0   = centraltemp_over_virial_*Tvir; */
  /* double beta = beta_; */
  /* double mycst_xi0=temperature_ratio_; */
  /* double mrond = CST_MU_ION/(CST_MU_ELEC+CST_MU_ION),  */
  /*   mrondxi = CST_MU_ION*mycst_xi0/(CST_MU_ELEC+CST_MU_ION*mycst_xi0); */
  /* double kappa = GYOTO_BOLTZMANN_CGS*T0 */
  /*   /((1.-beta)*pow(mycst_e0,CST_POLY_INDEX_M1) */
  /*     *GYOTO_ATOMIC_MASS_UNIT_CGS*CST_MU_ELEC*mrondxi */
  /*     *pow(GYOTO_C2_CGS,1.+CST_POLY_INDEX_M1)); */
  /* double density = GYOTO_C2_CGS_M1 */
  /*   *pow(1./kappa*(pow(1.+kappa */
  /* 		       *pow(mycst_e0*GYOTO_C2_CGS,CST_POLY_INDEX_M1),ww) */
  /* 		   -1.),CST_POLY_INDEX); */
  /* double PP = kappa*pow(density*GYOTO_C2_CGS,1.+CST_POLY_INDEX_M1); */
  /* // temperature */
  /* double T_electron = (mrond*(1.-ww)+mrondxi*ww)*CST_MU_ELEC*(1.-beta) */
  /*   *PP*GYOTO_ATOMIC_MASS_UNIT_CGS/(density*GYOTO_BOLTZMANN_CGS); */

  /* double Bnu = 2.*GYOTO_PLANCK_CGS*nuem*nuem*nuem*GYOTO_C2_CGS_M1 */
  /* 		    / (exp(GYOTO_PLANCK_CGS*nuem */
  /* 			   /(GYOTO_BOLTZMANN_CGS*T_electron))-1.); */
  /* double jnu = emission(nuem,dsem,coord_ph,coord_ph); */
  /* if (jnu==0.) return 0.; */
  /* double opacity_dsem = jnu/Bnu; */
  /* if (opacity_dsem!=opacity_dsem) */
  /*   throwError("In PD::transmission opacity is nan"); */
  /* if (opacity_dsem==opacity_dsem+1) */
  /*   throwError("In PD::transmission opacity is inf"); */
  /* return exp(-opacity_dsem);  */
  /* //NB : no multiplication by dsem here: */
  /* //     emission has a dsem factor, and opacity is just prop. to emission */
  /* //     thus the multiplication by dsem is already done. */
}

double PolishDoughnut::BBapprox(double nuem, double Te) const{

  double HnuOverKt=GYOTO_PLANCK_CGS*nuem/(GYOTO_BOLTZMANN_CGS*Te);
  double tol=1e-2;
  if (HnuOverKt<tol) //Rayleigh-Jeans
    return 2.*nuem*nuem*GYOTO_C2_CGS_M1*GYOTO_BOLTZMANN_CGS*Te;
  else if (HnuOverKt>1./tol) //Wien
    return 2.*GYOTO_PLANCK_CGS*nuem*nuem*nuem*GYOTO_C2_CGS_M1
      *exp(-HnuOverKt);
  else //Planck
    return 2.*GYOTO_PLANCK_CGS*nuem*nuem*nuem*GYOTO_C2_CGS_M1
      *1./(exp(HnuOverKt)-1.);
}

double PolishDoughnut::funcxM(double alpha1, double alpha2, 
			      double alpha3, double xM) {
  //Mahadevan 96 fit function
  return 4.0505*alpha1/pow(xM,1./6.)
    *(1.+0.4*alpha2/pow(xM,1./4.)+0.5316*alpha3/sqrt(xM))
    *exp(-1.8899*pow(xM,1./3.));
}

// Intersection of the constant angular momentum l0 with the Keplerian one
//double PolishDoughnut::intersection(double rr) const
double PolishDoughnut::intersection_t::operator()(double rr) const
  
{
  double y = ((rr*rr - 2.*aa_*sqrt(rr) + aa2_)/(pow(rr,3./2.) 
						- 2.*sqrt(rr) + aa_)) - l0_ ;
  
  return y ;   // y = 0 gives 2 intersections, 
  //the cusp and the central radius of the torus
}

// Solving the transcendental equation for "xM" numerically at each r
double PolishDoughnut::transcendental_t::operator()(double xM) const
{
  double       rr = par[0] ;
  double      n_e = par[1] ;
  double       BB = par[2] ;
  double       Te = par[3] ;
  double   alpha1 = par[4] ;
  double   alpha2 = par[5] ;
  double   alpha3 = par[6] ;


  double temp_e = GYOTO_BOLTZMANN_CGS*Te
    /(GYOTO_ELECTRON_MASS_CGS*GYOTO_C_CGS*GYOTO_C_CGS);
  double  nu_0 = GYOTO_ELEMENTARY_CHARGE_CGS* BB 
    / (2. * M_PI * GYOTO_ELECTRON_MASS_CGS * GYOTO_C_CGS);
  double func_xM = PolishDoughnut::funcxM(alpha1,alpha2,alpha3,xM);
  double BesselK2 = bessk(2,1./temp_e);
  double nuem = 3./2.*xM*nu_0*temp_e*temp_e;

  //   double y = GYOTO_ELEMENTARY_CHARGE_CGS*GYOTO_ELEMENTARY_CHARGE_CGS
  //     /(sqrt(3.)*GYOTO_C_CGS) * (4.*M_PI*n_e)/BesselK2
  //     * nuem * func_xM * 4./3.*M_PI*rr*rr*rr
  //     - M_PI * 2.*GYOTO_PLANCK_CGS/(GYOTO_C_CGS*GYOTO_C_CGS)*nuem*nuem*nuem
  //     / (exp(GYOTO_PLANCK_CGS*nuem/(GYOTO_BOLTZMANN_CGS*Te))-1.)
  //     * 4.*M_PI*rr*rr;
  double Bnu=0.;
  //if (exp(GYOTO_PLANCK_CGS*nuem/(GYOTO_BOLTZMANN_CGS*Te))==1.){
  //Rayleigh Jeans
  //AJOUTER TEST hnu<<kT
  Bnu = 2.*nuem*nuem/(GYOTO_C_CGS*GYOTO_C_CGS)*GYOTO_BOLTZMANN_CGS*Te;
  /*}else{
  //Planck
  Bnu = 2.*GYOTO_PLANCK_CGS*nuem*nuem*nuem/(GYOTO_C_CGS*GYOTO_C_CGS)
  / (exp(GYOTO_PLANCK_CGS*nuem/(GYOTO_BOLTZMANN_CGS*Te))
  -1.);
  }*/
  double y = GYOTO_ELEMENTARY_CHARGE_CGS*GYOTO_ELEMENTARY_CHARGE_CGS
    /(sqrt(3.)*GYOTO_C_CGS) * (4.*M_PI*n_e)/BesselK2
    * nuem * func_xM * 4./3.*M_PI*rr*rr*rr
    - M_PI * Bnu
    * 4.*M_PI*rr*rr;
  //This is j_nu_synch * 4/3pi*r3 - pi*B_Planck(Te)*4pi*r2

  return y ;   // y = 0 gives 1 intersection for each radius
}


// Potential W = -ln(ut)
double PolishDoughnut::potential(double rr, double theta) const
{
  double sinth = sin(theta), sinth2 = sinth*sinth;
  // the KerrBL metric parameters
  double  sigma = rr*rr + aa2_*cos(theta)*cos(theta);
  double  gtt = -(1. - 2.*rr/sigma);
  double  gtp = -2.*rr*aa_*sinth2/sigma;
  double  gpp = (rr*rr + aa2_ + 2.*rr*aa2_*sinth2/sigma)*sinth2;
  double  Omega = -(gtp + l0_ * gtt)/(gpp + l0_ * gtp) ;
  
  double  W = 0.5 * log(abs(gtt + 2. * Omega * gtp + Omega*Omega * gpp)) 
    - log(abs(gtt + Omega * gtp)) ;
  return  W ;
}

double PolishDoughnut::bessi0(double xx) {
  double ax,ans,y;
  if((ax=fabs(xx))< 3.75){ 
    y=xx/3.75;
    y*=y;
    ans=1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492
					 +y*(0.2659732
					     +y*(0.360768e-1
						 +y*0.45813e-2)))));
  }else{
    y=3.75/ax;
    ans=(exp(ax)/sqrt(ax))
      *(0.39894228
	+y*(0.1328592e-1
	    +y*(0.225319e-2
		+y*(-0.157565e-2
		    +y*(0.916281e-2
			+y*(-0.2057706e-1
			    +y*(0.2635537e-1
				+y*(-0.1647633e-1
				    +y*0.392377e-2))))))));
  }

  return ans;

}

double PolishDoughnut::bessk0(double xx) {
  double ans,y;
  if(xx<=2.0){
    y=xx*xx/4.0;
    ans=(-log(xx/2.0)*bessi0(xx))
      +(-0.57721566
	+y*(0.42278420
	    +y*(0.23069756+y*(0.3488590e-1
			      +y*(0.262698e-2
				  +y*(0.10750e-3+y*0.74e-5))))));
  }else{
    y=2.0/xx;
    ans=(exp(-xx)/sqrt(xx))*(1.25331414
			     +y*(-0.7832358e-1
				 +y*(0.2189568e-1
				     +y*(-0.1062446e-1
					 +y*(0.587872e-2
					     +y*(-0.251540e-2
						 +y*0.53208e-3))))));
  }

  return ans;

}

double PolishDoughnut::bessi1(double xx) {
  double ax,ans,y;
  if((ax=fabs(xx))< 3.75){ 
    y=xx/3.75;
    y*=y;
    ans=ax*(0.5+y*(0.87890594
		   +y*(0.51498869
		       +y*(0.15084934
			   +y*(0.2658733e-1
			       +y*(0.301532e-2+y*0.32411e-3))))));    
  }else{
    y=3.75/ax;
    ans=0.2282967e-1+y*(-0.2895312e-1+y*(0.1787654e-1
					 -y*0.420059e-2));
    ans=0.39894228+y*(-0.3988024e-1+y*(-0.362018e-2
				       +y*(0.163801e-2
					   +y*(-0.1031555e-1+y*ans))));
    ans*=(exp(ax)/sqrt(ax));
  }

  return xx<0.0 ? -ans : ans;


}

double PolishDoughnut::bessk1(double xx) {
  double yy,ans;
  if(xx<=2.0){
    yy=xx*xx/4.0;
    ans=(log(xx/2.0)*bessi1(xx))
      +(1.0/xx)*(1.0+yy*(0.15443144
			 +yy*(-0.67278579
			      +yy*(-0.18156897
				   +yy*(-0.1919402e-1
					+yy*(-0.110404e-2
					     +yy*(-0.4686e-4)))))));
  }else{
    yy=2.0/xx;
    ans=(exp(-xx)/sqrt(xx))*(1.25331414
			     +yy*(0.23498619
				  +yy*(-0.3655620e-1
				       +yy*(0.1504268e-1
					    +yy*(-0.780353e-2
						 +yy*(0.325614e-2
						      +yy*(-0.68245e-3)))))));
  }
  return ans;
}

double PolishDoughnut::bessk(int nn,double xx) {
  double bk,bkm,bkp,tox;
  if(nn< 2) throwError("PolishDoughnut::besselk n>2!");
  tox=2.0/xx;
  bkm=bessk0(xx);
  bk=bessk1(xx);
  for(int j=1;j<nn;j++){
    bkp=bkm+j*tox*bk;
    bkm=bk;
    bk=bkp;
  }
  return bk;
}

void PolishDoughnut::useSpecificImpact(int yes) {
  use_specific_impact_ = (yes != 0);
  cout << "use_specific_impact_==" << use_specific_impact_ << endl;
}

int PolishDoughnut::setParameter(string name, string content, string unit) {
  if      (name=="Lambda") setLambda(atof(content.c_str()));
  else if (name=="CentralDensity")
    setCentralDensity(atof(content.c_str()), unit);
  else if (name=="CentralTempOverVirial")
    centraltemp_over_virial_=atof(content.c_str());
  else if (name=="Beta") beta_=atof(content.c_str());
  else if (name=="UseSpecificImpact") useSpecificImpact();
  else if (name=="SpectralOversampling")
    spectral_oversampling_=atoi(content.c_str());
  else if (name=="Komissarov") komissarov_=1;
  else return Standard::setParameter(name, content, unit);
  return 0;
}

#ifdef GYOTO_USE_XERCES
void PolishDoughnut::fillElement(FactoryMessenger *fmp) const {
  fmp->setMetric(gg_);
  fmp->setParameter("Lambda", lambda_);
  fmp->setParameter("CentralDensity", central_density_);
  fmp->setParameter("CentralTempOverVirial", centraltemp_over_virial_);
  fmp->setParameter("Beta", beta_);
  if (use_specific_impact_) fmp->setParameter("UseSpecificImpact");
  fmp->setParameter("SpectralOversampling", spectral_oversampling_);
  fmp->setParameter("Komissarov", komissarov_);
  Standard::fillElement(fmp);
}
#endif
