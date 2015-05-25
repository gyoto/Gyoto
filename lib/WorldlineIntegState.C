/*
    Copyright 2011-2015 Frederic Vincent, Thibaut Paumard

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
#include <boost/numeric/odeint/stepper/generation.hpp>
using namespace boost::numeric::odeint;

#define GYOTO_TRY_BOOST_CONTROLLED_STEPPER(a)				\
  if (kind_==Kind::a) {							\
    typedef boost::numeric::odeint::a<state_type> error_stepper_type;	\
    auto controlled=							\
      make_controlled< error_stepper_type >				\
           ( line->absTol() , line->relTol() );				\
    try_step_ =								\
      [controlled, system]						\
      (std::array<double, 8> &inout, double &t, double &h)		\
      mutable								\
      -> controlled_step_result						\
    {									\
      return controlled.try_step(system, inout, t, h);			\
    };									\
    do_step_ =								\
      [controlled, system]						\
      (std::array<double, 8> &inout, double h)		\
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
  gg_=line_->metric();
}
void
Worldline::IntegState::Generic::init(Worldline * line,
				     const double coord[8],
				     const double delta)
{
  line_=line;
  delta_=delta;
  gg_=line->metric();
  norm_=normref_= gg_->ScalarProd(coord,coord+4,coord+4);
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
				    const double coord[8], const double delta) {
  Generic::init(line, coord, delta);

  short i;
  for (i=0;i<8;++i) coord_[i]=coord[i];

}

int Worldline::IntegState::Legacy::nextStep(double coord[8], double h1max) {
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

  checkNorm(coord);
# if GYOTO_DEBUG_ENABLED
  GYOTO_IF_DEBUG
  GYOTO_DEBUG_ARRAY(coord,8);
  GYOTO_DEBUG_EXPR(delta_);
  GYOTO_ENDIF_DEBUG
# endif

  if (delta_==delta_+1) return 1; // delta == Infinity : stop condition  
  return 0;
}

void Worldline::IntegState::Legacy::doStep(double const coordin[8],
					   double step, 
					   double coordout[8]) {
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
  else throwError("unknown integrator kind");
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
  std::function<void(const std::array<double, 8> &/*x*/,
		  std::array<double, 8> & /*dxdt*/,
		     const double /* t*/ )> system;

  if (!met) 
    system=[](const std::array<double, 8> &/*x*/,
		  std::array<double, 8> & /*dxdt*/,
		  const double /* t*/ ){
      throwError("Metric not set");
    };
  else
    system=[this, line, met](const std::array<double, 8> &x,
				  std::array<double, 8> &dxdt,
				  const double t)
      {
	double xx[8]={x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7]};
	double res[8];
	//line->stopcond=met->Generic::diff(xx, res);
	line->stopcond=met->diff(xx, res);
	for (size_t i=0; i<=7; ++i) dxdt[i]=res[i];
      };

  typedef std::array< double, 8 > state_type;
  GYOTO_TRY_BOOST_CONTROLLED_STEPPER(runge_kutta_cash_karp54)
  else GYOTO_TRY_BOOST_CONTROLLED_STEPPER(runge_kutta_fehlberg78)
  else GYOTO_TRY_BOOST_CONTROLLED_STEPPER(runge_kutta_dopri5)
  else GYOTO_TRY_BOOST_CONTROLLED_STEPPER(runge_kutta_cash_karp54_classic)
	 //else GYOTO_TRY_BOOST_CONTROLLED_STEPPER(rosenbrock4)
  else throwError("unknown stepper type");
};

Worldline::IntegState::Boost *
Worldline::IntegState::Boost::clone(Worldline*newparent) const
{ return new Boost(newparent, kind_); }



void
Worldline::IntegState::Boost::init(Worldline * line,
				   const double coord[8], const double delta) {
  Generic::init(line, coord, delta);

}

int Worldline::IntegState::Boost::nextStep(double coord[8], double h1max) {
  GYOTO_DEBUG << h1max << endl;
  
  // We first make a C++ std::array out of the bare C array:
  std::array<double, 8> inout = {
    coord[0], coord[1], coord[2], coord[3],
    coord[4], coord[5], coord[6], coord[7]};

  if (adaptive_) {
    double h1=delta_;
    double sgn=h1>0?1.:-1.;
    h1max=line_->deltaMax(coord, h1max);
    double delta_min=line_->deltaMin();
    double dt=0.;

    if (abs(h1)>h1max) h1=sgn*h1max;
    if (abs(h1)<delta_min) h1=sgn*delta_min;
    controlled_step_result cres;
    GYOTO_DEBUG << h1 << endl;

    do {
      // try_step_ is a lambda function encapsulating
      // the actual adaptive-step integrator from boost
      cres=try_step_(inout, dt, h1);
    } while (abs(h1)>=delta_min &&
	     cres==controlled_step_result::fail &&
	     abs(h1)<h1max);
  
    // Check and report two possible error conditions (possible bugs)
    if (sgn*h1<0) throwError("h1 changed sign!");
    if (abs(dt)>h1max) throwError("used step larger than provided");

    // cres is still fail, redo with delta_min using the fixed-step integrator
    if (cres==controlled_step_result::fail) {
      GYOTO_SEVERE << "delta_min is too large: " << delta_min << endl;
      for (size_t i=0; i<=7; ++i) inout[i]=coord[i];
      do_step_(inout, sgn*delta_min);
    }
    // update adaptive step
    delta_=h1;
  } else {
    // non adaptive case
    // do_Step_ is a lambda function encapsulating a fixed-step integrator
    // from Boost
    do_step_(inout, delta_);
  }

  for (size_t i=0; i<=7; ++i) coord[i]=inout[i];

  checkNorm(coord);

  return line_->stopcond;
}

void Worldline::IntegState::Boost::doStep(double const coordin[8],
					 double step, 
					 double coordout[8]) {
  // We first make a C++ std::array out of the bare C array:
  std::array<double, 8> inout = {
    coordin[0], coordin[1], coordin[2], coordin[3],
    coordin[4], coordin[5], coordin[6], coordin[7]};

  // We call the Boost stepper
  do_step_(inout, step);

  // Copy the result
  for (size_t i=0; i<=7; ++i) coordout[i]=inout[i];
}

std::string Worldline::IntegState::Boost::kind() {
  if (kind_== Kind::runge_kutta_cash_karp54) return "runge_kutta_cash_karp54";
  if (kind_== Kind::runge_kutta_fehlberg78) return "runge_kutta_fehlberg78";
  if (kind_== Kind::runge_kutta_dopri5) return "runge_kutta_dopri5";
  if (kind_== Kind::runge_kutta_cash_karp54_classic) return "runge_kutta_cash_karp5";
  throwError("unknown enum value");
  return "error";
} 
#endif // GYOTO_HAVE_BOOST_INTEGRATORS
