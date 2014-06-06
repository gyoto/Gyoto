/*
    Copyright 2011 Frederic Vincent, Thibaut Paumard

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

#include <boost/numeric/odeint/stepper/generation.hpp>

using namespace std ; 
using namespace Gyoto;
using namespace boost::numeric::odeint;

#define GYOTO_TRY_BOOST_CONTROLLED_STEPPER(a)				\
  if (kind_==#a) {							\
    typedef a<state_type> error_stepper_type;				\
    auto controlled=make_controlled< error_stepper_type >( 1.0e-10 , 1.0e-6 ); \
									\
    double delta_min=met->deltaMin();					\
									\
    stepper_=								\
      [this, controlled, delta_min, system, met, line]			\
      (double coord[8], double h1max)					\
      mutable								\
      {									\
	double h1=this->delta_;                                         \
	double sgn=h1>0?1.:-1.;						\
        h1max=met->deltaMax(coord, h1max);				\
	if (abs(h1)>h1max) h1=sgn*h1max;                                \
	if (abs(h1)<delta_min) h1=sgn*delta_min;			\
	controlled_step_result cres;					\
	std::array<double, 8> inout = {					\
	  coord[0], coord[1], coord[2], coord[3],			\
	  coord[4], coord[5], coord[6], coord[7]};			\
									\
	double dt=0.;							\
									\
	do {								\
	  cres=controlled.try_step(system, inout, dt, h1);		\
	} while (abs(h1)>=delta_min &&					\
		 cres==controlled_step_result::fail &&			\
		 abs(h1)<h1max);					\
									\
	if (sgn*h1<0) throwError("h1 changed sign!");			\
	if (abs(dt)>h1max) throwError("used step larger than provided");\
	if (cres==controlled_step_result::fail) {			\
	  GYOTO_SEVERE << "delta_min is too large: " << delta_min << endl;\
	  for (size_t i=0; i<=7; ++i) inout[i]=coord[i];		\
	  controlled.stepper().do_step(system, inout, dt, sgn*delta_min); \
	}								\
									\
	for (size_t i=0; i<=7; ++i) coord[i]=inout[i];			\
	this->delta_=h1;						\
      };								\
   }

/// Generic
Worldline::IntegState::Generic::~Generic() {};
Worldline::IntegState::Generic::Generic() : SmartPointee() {};
void
Worldline::IntegState::Generic::init(Worldline * line,
				     const double coord[8],
				     const double delta)
{
  line_=line;
  delta_=delta;
}

/// Legacy

Worldline::IntegState::Legacy::Legacy() : Generic() {}


Worldline::IntegState::Legacy *
Worldline::IntegState::Legacy::clone() const
{ return new Legacy(*this); }


void
Worldline::IntegState::Legacy::init(Worldline * line,
				    const double coord[8], const double delta) {
  Generic::init(line, coord, delta);
  gg_=line->metric();
  adaptive_=line->adaptive();

  short i;
  for (i=0;i<8;++i) coord_[i]=coord[i];

  norm_=normref_= gg_->ScalarProd(coord,coord+4,coord+4);
}

int Worldline::IntegState::Legacy::nextStep(double coord[8], double h1max) {
  int j;
  double h1;

  if (adaptive_){
    if (gg_ -> myrk4_adaptive(line_,coord_,norm_,normref_,coord,delta_,h1, h1max)) return 1;
    delta_ = h1;
  }else{
    if (gg_ -> myrk4(line_,coord_,delta_,coord)) return 1; 
  }
  for (j=0;j<8;j++) coord_[j] = coord[j];

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

# if GYOTO_DEBUG_ENABLED
  GYOTO_IF_DEBUG
  GYOTO_DEBUG_ARRAY(coord,8);
  GYOTO_DEBUG_EXPR(delta_);
  GYOTO_ENDIF_DEBUG
# endif

  if (delta_==delta_+1) return 1; // delta == Infinity : stop condition  
  return 0;
}


std::string Worldline::IntegState::Legacy::kind() { return "Legacy"; } 

Worldline::IntegState::Legacy::~Legacy() {}

/// Boost
Worldline::IntegState::Boost::~Boost() {};
Worldline::IntegState::Boost::Boost(std::string type) : Generic(), kind_(type) {
};

Worldline::IntegState::Boost *
Worldline::IntegState::Boost::clone() const
{ return new Boost(*this); }



void
Worldline::IntegState::Boost::init(Worldline * line,
				   const double coord[8], const double delta) {
  line_=line;
  Metric::Generic* met=line_->metric();
  delta_=delta;

  // This, below, is called a lambda function.
  auto system=[this, line, met](const std::array<double, 8> &x,
		 std::array<double, 8> &dxdt,
		 const double t)
    {
      double xx[8]={x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7]};
      double res[8];
      line->stopcond=met->Generic::diff(xx, res);
      for (size_t i=0; i<=7; ++i) dxdt[i]=res[i];
    };

  typedef std::array< double, 8 > state_type;
  GYOTO_TRY_BOOST_CONTROLLED_STEPPER(runge_kutta_cash_karp54)
  else GYOTO_TRY_BOOST_CONTROLLED_STEPPER(runge_kutta_fehlberg78)
  else GYOTO_TRY_BOOST_CONTROLLED_STEPPER(runge_kutta_dopri5)
  else GYOTO_TRY_BOOST_CONTROLLED_STEPPER(runge_kutta_cash_karp54_classic)
	 //else GYOTO_TRY_BOOST_CONTROLLED_STEPPER(rosenbrock4)
  else throwError("unknown stepper type");


}

int Worldline::IntegState::Boost::nextStep(double coord[8], double h1max) {
  stepper_(coord, h1max);
  return line_->stopcond;
}

std::string Worldline::IntegState::Boost::kind() { return kind_; } 
