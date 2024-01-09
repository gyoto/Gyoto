/*
    Copyright 2011-2015, 2018-2020 Frederic Vincent, Thibaut Paumard

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
#include <iostream>
#include <cstdlib>
#include <GyotoWorldline.h>
#include <cmath>
#include <string>
#include <cstring>
#include <ctime>

using namespace std ; 
using namespace Gyoto;

#ifdef GYOTO_HAVE_BOOST_INTEGRATORS
# include <boost/version.hpp>
# if BOOST_VERSION >= 106400 
#  include <boost/serialization/array_wrapper.hpp>
# endif // BOOST_VERSION >= 106400 
#include <boost/numeric/odeint/stepper/generation.hpp>
using namespace boost::numeric::odeint;

#if defined HAVE_FENV_H
# include <fenv.h>
# pragma STDC FENV_ACCESS ON
# define DISABLE_SIGFPE							\
  fenv_t envp;								\
  if (feholdexcept(&envp)) GYOTO_ERROR("failed holding FPE")
# define REENABLE_SIGFPE						\
  if (feclearexcept(FE_ALL_EXCEPT)) GYOTO_ERROR("failed clearing FPE");	\
  if (fesetenv(&envp)) GYOTO_ERROR("failed setting back FPE")
#else
# define DISABLE_SIGFPE
# define REENABLE_SIGFPE
#endif

#define GYOTO_TRY_BOOST_CONTROLLED_STEPPER(a)				\
  if (kind_==Kind::a) {							\
    typedef boost::numeric::odeint::a<state_t> error_stepper_type;	\
    DISABLE_SIGFPE;							\
    auto controlled=							\
      make_controlled< error_stepper_type >				\
           ( line->absTol() , line->relTol() );				\
    REENABLE_SIGFPE;							\
    try_step_ =								\
      [controlled, system]						\
      (state_t &inout, double &t, double &h)				\
      mutable								\
      -> controlled_step_result						\
    {									\
      return controlled.try_step(system, inout, t, h);			\
    };									\
    do_step_ =								\
      [controlled, system]						\
      (state_t &inout, double h)					\
      mutable								\
    {									\
      controlled.stepper().do_step(system, inout, 0., h);		\
    };									\
  }

#endif // GYOTO_HAVE_BOOST_INTEGRATORS


/// Generic
Worldline::IntegState::Generic::~Generic() {};
Worldline::IntegState::Generic::Generic(Worldline *parent) :
  SmartPointee(), line_(parent), gg_(NULL), integ_31_(false) {};
void
Worldline::IntegState::Generic::init(){
  if (!line_) return;
  adaptive_=line_->adaptive();
  parallel_transport_=line_->parallelTransport();
  gg_=line_->metric();
}
void
Worldline::IntegState::Generic::init(Worldline * line,
				     const state_t &coord,
				     const double delta)
{
  line_=line;
  init();
  delta_=delta;
  if (line_->getImin() <= line_->getImax() && gg_) norm_=normref_= gg_->ScalarProd(&coord[0],&coord[4],&coord[4]);
}

void Worldline::IntegState::Generic::checkNorm(double coord[8])
{
  norm_=gg_ -> ScalarProd(coord,coord+4,coord+4);

  double normtol=.001;
  /* 
     NB: following test done for norm/tdot
     as tdot can diverge close to horizon (it's the case for
     NS integration eg where geodesic can come close to horizon)
     Then just check that norm/tdot does not diverge.
   */
  if (fabs(norm_-normref_)/(coord[4]*coord[4])>normtol) {
    GYOTO_SEVERE << 
      "in Worldline::IntegState.C: "
      "norm is drifting"
      " - with norm,normref= " << norm_ << " " 
  		 << normref_ << " -- x1,x2,x3= " << coord[1] 
  		 << " " << coord[2] << " " << coord[3] << " " << endl;
  }
}

void Worldline::IntegState::Generic::integ31(bool integ) {
  integ_31_=integ;
}

bool Worldline::IntegState::Generic::integ31() const{
  return integ_31_;
}

/// Legacy

Worldline::IntegState::Legacy::Legacy(Worldline *parent) : Generic(parent)
{}


Worldline::IntegState::Legacy *
Worldline::IntegState::Legacy::clone(Worldline *newparent) const
{ return new Legacy(newparent); }

void
Worldline::IntegState::Legacy::init(Worldline * line,
				    const state_t &coord, const double delta) {
  static bool need_warning = true;
  if (need_warning) {
    GYOTO_WARNING << "The Legacy integrator is deprecated and will be removed soon. "
		  << "Please update your code to use the Boost integrators." << endl;
    need_warning=false;
  }
  Generic::init(line, coord, delta);
  coord_ = coord;

}

int Worldline::IntegState::Legacy::nextStep(state_t &coord, double &tau, double h1max) {
  if (!gg_) init();
  static bool need_warning=true;
  if (need_warning) {
    GYOTO_WARNING << "The Legacy integrator does not compute proper time." << endl;
    need_warning=false;
  }
  tau=0.;
  if (parallel_transport_) GYOTO_ERROR("TODO: implement parallel transport");
  GYOTO_DEBUG << h1max << endl;
  int j;
  double h1;
  
  if (adaptive_){
    if (gg_ -> myrk4_adaptive(line_,coord_,norm_,normref_,coord,delta_,h1, h1max)) return 1;
    delta_ = h1;
  }else{
    if (gg_ -> myrk4(line_,coord_,delta_,coord)) return 1;
  }
  for (j=0;j<8;j++) coord_[j] = coord[j];

  checkNorm(&coord[0]);
# if GYOTO_DEBUG_ENABLED
  GYOTO_IF_DEBUG
  GYOTO_DEBUG_ARRAY(coord,8);
  GYOTO_DEBUG_EXPR(delta_);
  GYOTO_ENDIF_DEBUG
# endif

  if (delta_==delta_+1) return 1; // delta == Infinity : stop condition  
  return 0;
}

void Worldline::IntegState::Legacy::doStep(state_t const &coordin,
					   double step, 
					   state_t &coordout) {
  if (!gg_) init();
  gg_ -> myrk4(line_, coordin, step, coordout); 
}

std::string Worldline::IntegState::Legacy::kind() { return "Legacy"; } 

Worldline::IntegState::Legacy::~Legacy() {}

/// Boost
#ifdef GYOTO_HAVE_BOOST_INTEGRATORS
Worldline::IntegState::Boost::~Boost() {};
Worldline::IntegState::Boost::Boost(Worldline*line, std::string type) :
  Generic(line)
{
  if (type=="runge_kutta_cash_karp54") kind_=runge_kutta_cash_karp54;
  else if (type=="runge_kutta_fehlberg78") kind_=runge_kutta_fehlberg78;
  else if (type=="runge_kutta_dopri5") kind_=runge_kutta_dopri5;
  else if (type=="runge_kutta_cash_karp54_classic") kind_=runge_kutta_cash_karp54_classic;
  else GYOTO_ERROR("unknown integrator kind");
}

Worldline::IntegState::Boost::Boost(Worldline*line, Kind type) :
  Generic(line), kind_(type)
{}

void Worldline::IntegState::Boost::init()
{
  Generic::init();
  Worldline* line=line_;
  Metric::Generic* met=line->metric();
  system_t system;
  double mass=line->getMass();

  if (!met)
    system=[](const state_t &/*x*/,
	      state_t & /*dxdt*/,
	      const double /* t*/ ){
      GYOTO_ERROR("Metric not set");
    };
  else{
    if (integ_31_==false){
      system=[this, line, met, mass](const state_t &x,
				     state_t &dxdt,
				     const double t)
	{
	  line->stopcond=met->diff(x, dxdt, mass);
	};
    }else{
      system=[this, line, met, mass](const state_t &x,
				     state_t &dxdt,
				     const double t)
	{
	  line->stopcond=met->diff31(x, dxdt, mass); // time must be passed
	};
    }
  }

  if (line->getImin() > line->getImax() || !met) return;

  GYOTO_TRY_BOOST_CONTROLLED_STEPPER(runge_kutta_cash_karp54)
  else GYOTO_TRY_BOOST_CONTROLLED_STEPPER(runge_kutta_fehlberg78)
  else GYOTO_TRY_BOOST_CONTROLLED_STEPPER(runge_kutta_dopri5)
  else GYOTO_TRY_BOOST_CONTROLLED_STEPPER(runge_kutta_cash_karp54_classic)
	 //else GYOTO_TRY_BOOST_CONTROLLED_STEPPER(rosenbrock4)
  else GYOTO_ERROR("unknown stepper type");
};

Worldline::IntegState::Boost *
Worldline::IntegState::Boost::clone(Worldline*newparent) const
{ return new Boost(newparent, kind_); }



void
Worldline::IntegState::Boost::init(Worldline * line,
				   const state_t &coord, const double delta) {
  Generic::init(line, coord, delta);

}

int Worldline::IntegState::Boost::nextStep(state_t &coord, double& tau, double h1max) {
  if (!gg_) init();
  GYOTO_DEBUG << h1max << endl;
  double dt=0, dtau=0;

  // Transform to proper vector depending on integration kind (4D/3+1)
  state_t xx;
  double told = coord[0];
  if (integ_31_==false) xx = coord;
  else{
    double rr=coord[1], th=coord[2], ph=coord[3],
      tdot=coord[4], rdot=coord[5], thdot=coord[6], phdot=coord[7];
    if (tdot==0.) GYOTO_ERROR("In WlI::nextStep tdot is 0!");
    double rprime=rdot/tdot, thprime=thdot/tdot, phprime=phdot/tdot;
    double NN, beta[3];
    gg_->computeNBeta(&coord[0],NN,beta);
    double betar=beta[0], betat=beta[1], betap=beta[2];
    
    double Vr = 1./NN*(rprime+betar), Vth = 1./NN*(thprime+betat),
      Vph = 1./NN*(phprime+betap);
    // Photon's energy as measured by Eulerian observer:
    double EE = tdot*NN;
    xx = {EE,rr,th,ph,Vr,Vth,Vph};
  }

  if (adaptive_) {
    double h1=delta_;
    double sgn=h1>0?1.:-1.;
    h1max=line_->deltaMax(&coord[0], h1max);
    double delta_min=line_->deltaMin();

    if (abs(h1)>h1max) h1=sgn*h1max;
    if (abs(h1)<delta_min) h1=sgn*delta_min;
    controlled_step_result cres;
    GYOTO_DEBUG << h1 << endl;

    do {
      // try_step_ is a lambda function encapsulating
      // the actual adaptive-step integrator from boost
      cres=try_step_(xx, dt, h1);
    } while (abs(h1)>=delta_min &&
	     cres==controlled_step_result::fail &&
	     abs(h1)<h1max);

    // At this point, if a successful step was found, xx is updated,
    // dt is the increment of the integration variable (typically
    // affine parameter, proper time, or coordinate time) over the
    // successful step, and h1 is the proposed step size for
    // the next iteration.
  
    // Check and report two possible error conditions (possible bugs)
    if (sgn*h1<0) GYOTO_ERROR("h1 changed sign!");
    if (abs(dt)>h1max) GYOTO_ERROR("used step larger than provided");

    // cres is still fail, redo with delta_min using the fixed-step integrator
    if (cres==controlled_step_result::fail) {
      GYOTO_SEVERE << "delta_min is too large: " << delta_min << endl;
      dt=sgn*delta_min;
      do_step_(xx, dt);
    }
    // update adaptive step
    delta_=h1;

  } else {
    // non adaptive case
    // do_Step_ is a lambda function encapsulating a fixed-step integrator
    // from Boost
    dt=delta_;
    do_step_(xx, dt);
  }

  // Transform back to original vector coord
  if (integ_31_==false) {
    coord = xx;
    dtau = dt;
  }else{
    double tnew = told+dt, EE = xx[0], rr=xx[1], th=xx[2], ph=xx[3],
      Vr=xx[4], Vth=xx[5], Vph=xx[6];
    double NN, beta[3];
    double newpos[4]={tnew,rr,th,ph};
    gg_->computeNBeta(newpos,NN,beta);
    double beta_r=beta[0], beta_t=beta[1], beta_p=beta[2];
    
    double rprime=NN*Vr-beta_r,
      thprime=NN*Vth-beta_t,
      phprime=NN*Vph-beta_p,
      tdotnew = EE/NN,
      rdot = rprime*tdotnew,
      thdot = thprime*tdotnew,
      phdot = phprime*tdotnew;
    
    coord = {tnew,rr,th,ph,tdotnew,rdot,thdot,phdot};
    dtau = NN/EE*dt; // affine parameter increment
  }

  tau += dtau;
  checkNorm(&coord[0]);
  if (gg_ -> coordKind() == GYOTO_COORDKIND_SPHERICAL){
     line_->checkPhiTheta(&coord[0]);
  }

  return line_->stopcond;
}

void Worldline::IntegState::Boost::doStep(state_t const &coordin,
					  double step, 
					  state_t &coordout) {
  if (!gg_) init();
  coordout = coordin;

  // Transform to proper vector depending on integration kind (4D/3+1)
  state_t xx;
  double told = coordout[0];

  if (integ_31_==false) xx = coordout;
  else{
    double rr=coordout[1], th=coordout[2], ph=coordout[3],
      tdot=coordout[4], rdot=coordout[5], thdot=coordout[6], phdot=coordout[7];
    if (tdot==0.) GYOTO_ERROR("In WlI::nextStep tdot is 0!");
    double rprime=rdot/tdot, thprime=thdot/tdot, phprime=phdot/tdot;
    double NN, beta[3];
    gg_->computeNBeta(&coordout[0],NN,beta);
    double betar=beta[0], betat=beta[1], betap=beta[2];
    
    double Vr = 1./NN*(rprime+betar), Vth = 1./NN*(thprime+betat),
      Vph = 1./NN*(phprime+betap);
    // Photon's energy as measured by Eulerian observer:
    double EE = tdot*NN;
    xx = {EE,rr,th,ph,Vr,Vth,Vph};
  }

  // We call the Boost stepper
  do_step_(xx, step);

  // Transform back to original vector coord
  if (integ_31_==false) coordout = xx;
  else{
    double tnew = told+step, EE = xx[0], rr=xx[1], th=xx[2], ph=xx[3],
      Vr=xx[4], Vth=xx[5], Vph=xx[6];
    double NN, beta[3];
    double newpos[4]={tnew,rr,th,ph};
    gg_->computeNBeta(newpos,NN,beta);
    double beta_r=beta[0], beta_t=beta[1], beta_p=beta[2];
    
    double rprime=NN*Vr-beta_r,
      thprime=NN*Vth-beta_t,
      phprime=NN*Vph-beta_p,
      tdotnew = EE/NN,
      rdot = rprime*tdotnew,
      thdot = thprime*tdotnew,
      phdot = phprime*tdotnew;
    
    coordout = {tnew,rr,th,ph,tdotnew,rdot,thdot,phdot};
  }
}

std::string Worldline::IntegState::Boost::kind() {
  if (kind_== Kind::runge_kutta_cash_karp54) return "runge_kutta_cash_karp54";
  if (kind_== Kind::runge_kutta_fehlberg78) return "runge_kutta_fehlberg78";
  if (kind_== Kind::runge_kutta_dopri5) return "runge_kutta_dopri5";
  if (kind_== Kind::runge_kutta_cash_karp54_classic) return "runge_kutta_cash_karp5";
  GYOTO_ERROR("unknown enum value");
  return "error";
} 
#endif // GYOTO_HAVE_BOOST_INTEGRATORS
