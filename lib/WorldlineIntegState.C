/*
    Copyright 2011-2015, 2018 Frederic Vincent, Thibaut Paumard

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
  SmartPointee(), line_(parent) {};
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

/// Legacy

Worldline::IntegState::Legacy::Legacy(Worldline *parent) : Generic(parent)
{
  Legacy::init();
}


Worldline::IntegState::Legacy *
Worldline::IntegState::Legacy::clone(Worldline *newparent) const
{ return new Legacy(newparent); }

void
Worldline::IntegState::Legacy::init(Worldline * line,
				    const state_t &coord, const double delta) {
  Generic::init(line, coord, delta);
  coord_ = coord;

}

int Worldline::IntegState::Legacy::nextStep(state_t &coord, double h1max) {
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
  Boost::init();
}

Worldline::IntegState::Boost::Boost(Worldline*line, Kind type) :
  Generic(line), kind_(type)
{
  Boost::init();
}

void Worldline::IntegState::Boost::init()
{
  Generic::init();
  Worldline* line=line_;
  Metric::Generic* met=line->metric();
  system_t system;

  if (!met)
    system=[](const state_t &/*x*/,
	      state_t & /*dxdt*/,
	      const double /* t*/ ){
      GYOTO_ERROR("Metric not set");
    };
  else
    system=[this, line, met](const state_t &x,
			     state_t &dxdt,
			     const double t)
      {
	line->stopcond=met->diff(x, dxdt);
      };

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

int Worldline::IntegState::Boost::nextStep(state_t &coord, double h1max) {
  GYOTO_DEBUG << h1max << endl;
  
  if (adaptive_) {
    double h1=delta_;
    double sgn=h1>0?1.:-1.;
    h1max=line_->deltaMax(&coord[0], h1max);
    double delta_min=line_->deltaMin();
    double dt=0.;

    if (abs(h1)>h1max) h1=sgn*h1max;
    if (abs(h1)<delta_min) h1=sgn*delta_min;
    controlled_step_result cres;
    GYOTO_DEBUG << h1 << endl;

    do {
      // try_step_ is a lambda function encapsulating
      // the actual adaptive-step integrator from boost
      cres=try_step_(coord, dt, h1);
    } while (abs(h1)>=delta_min &&
	     cres==controlled_step_result::fail &&
	     abs(h1)<h1max);
  
    // Check and report two possible error conditions (possible bugs)
    if (sgn*h1<0) GYOTO_ERROR("h1 changed sign!");
    if (abs(dt)>h1max) GYOTO_ERROR("used step larger than provided");

    // cres is still fail, redo with delta_min using the fixed-step integrator
    if (cres==controlled_step_result::fail) {
      GYOTO_SEVERE << "delta_min is too large: " << delta_min << endl;
      do_step_(coord, sgn*delta_min);
    }
    // update adaptive step
    delta_=h1;
  } else {
    // non adaptive case
    // do_Step_ is a lambda function encapsulating a fixed-step integrator
    // from Boost
    do_step_(coord, delta_);
  }

  checkNorm(&coord[0]);

  return line_->stopcond;
}

void Worldline::IntegState::Boost::doStep(state_t const &coordin,
					  double step, 
					  state_t &coordout) {
  coordout = coordin;

  // We call the Boost stepper
  do_step_(coordout, step);
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
