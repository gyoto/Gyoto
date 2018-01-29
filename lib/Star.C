/*
    Copyright 2011 Thibaut Paumard, Frederic Vincent

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
#include "GyotoWorldline.h"
#include "GyotoStar.h"
#include "GyotoProperty.h"
#include "GyotoPhoton.h"
#include "GyotoPowerLawSpectrum.h"
#include "GyotoBlackBodySpectrum.h"
#include "GyotoFactoryMessenger.h"

#include <iostream>
#include <cmath>
#include <string>
#include <cstdlib>
#include <float.h>
#include <sstream>
#include <string.h>

using namespace std;
using namespace Gyoto;
using namespace Gyoto::Astrobj;

/// Properties
GYOTO_PROPERTY_START(Gyoto::Astrobj::Star,
 "UniformSphere following a time-like Gyoto::Worldline.")
// Star only need to implement the Worldline interface on top of the 
// UniformSphere interface, which is trivially tone with this macro:
GYOTO_WORLDLINE_PROPERTY_END(Star, UniformSphere::properties)

// XML I/O
// We also need to parse and write Position+Velocity in addition to
// InitCoord, which is done by overriding setParameter(), setParameters()
// and fillProperty()
int Star::setParameter(std::string name,
			    std::string content,
			    std::string unit) {
  double coord[8];
  if (name=="InitialCoordinate") {
    name="InitCoord";
    return UniformSphere::setParameter(name, content, unit);
  } else if (name=="Position") {
    if (FactoryMessenger::parseArray(content, coord, 4) != 4)
      throwError("Worldline \"Position\" requires exactly 4 tokens");
    if (init_vel_) {
      setInitCoord(coord, init_vel_);
      delete[] init_vel_; init_vel_=NULL;
    } else setPosition(coord);
    wait_pos_ = 0;
  } else if (name=="Velocity") {
    if (FactoryMessenger::parseArray(content, coord, 3) != 3)
      throwError("Worldline \"Velocity\" requires exactly 3 tokens");
    if (wait_pos_) {
      if (init_vel_) delete [] init_vel_;
      init_vel_ = new double[3];
      memcpy(init_vel_, coord, 3*sizeof(double));
    } else setVelocity(coord);
  }
  else return UniformSphere::setParameter(name, content, unit);
  return 0;
}

#ifdef GYOTO_USE_XERCES
void Star::fillProperty(Gyoto::FactoryMessenger *fmp, Property const &p) const {
  if (p.name == "InitCoord") {
    if (imin_ <= imax_) {
      state_t coord;
      getInitialCoord(coord);
      // For massive particule, express initial condition with 3-velocity
      double vel[3] = {coord[5]/coord[4], coord[6]/coord[4], coord[7]/coord[4]};
      fmp -> setParameter ("Position", &coord[0], 4);
      fmp -> setParameter ("Velocity", vel, 3);
    }
    return;
  }
  UniformSphere::fillProperty(fmp, p);
}

void Star::setParameters(FactoryMessenger* fmp) {
  wait_pos_ = 1;
  metric(fmp->metric());
  UniformSphere::setParameters(fmp);
  wait_pos_ = 0;
  if (init_vel_) {
    delete[] init_vel_; init_vel_=NULL;
    throwError("Worldline::setParameters(): "
	       "Velocity was found but not Position");
  }
}
#endif
///

Star::Star() :
  UniformSphere("Star"),
  Worldline()
{
# ifdef GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << "done." << endl;
# endif
}

Star::Star(SmartPointer<Metric::Generic> met, double rad,
	   double const pos[4],
	   double const v[3]) :
  UniformSphere("Star"),
  Worldline()
{
  if (debug()) {
    cerr << "DEBUG: Star Construction " << endl
	 << "       POS=[" << pos[0];
    for (int i=1; i<4; ++i) cerr << ", " << pos[i];
    cerr << "]\n       VEL=[" << v[0] ;
    for (int i=1; i<3; ++i) cerr << ", " << v[i];
    cerr << "]\n       RADIUS=" << rad << endl;

  }

  metric(met);
  setInitCoord(pos, v);
  radius(rad);
}

Star::Star(const Star& orig) :
  UniformSphere(orig), Worldline(orig)
{
  GYOTO_DEBUG << endl;
  // we have two distinct clones of the metric, not good...
  Worldline::metric(UniformSphere::metric());
}

Star* Star::clone() const { return new Star(*this); }

Star::~Star() {
  if (debug()) cerr << "DEBUG: Star::~Star()\n";
}

string Star::className() const { return  string("Star"); }
string Star::className_l() const { return  string("star"); }

SmartPointer<Metric::Generic> Star::metric() const { return gg_; }
void Star::metric(SmartPointer<Metric::Generic> gg) {
  UniformSphere::metric(gg);
  Worldline::metric(gg);
}

void Star::setInitialCondition(double const coord[8]) {
  if (!metric_) throwError("Please set metric before calling Star::setInitialCondition(double*)");
  Worldline::setInitialCondition(metric_, coord, 0);
}

double Star::getMass() const {return 1. ;}

void Star::getVelocity(double const pos[4], double vel[4]) {
  getCoord(pos, 1, NULL, NULL, NULL, vel, vel+1, vel+2, vel+3);
}

void Star::getCartesian(double const * const t,
			size_t const n,
			double* const x, double*const y, double*const z,
			double*const xp, double*const yp, double*const zp) {
  Worldline::getCartesian(t, n, x, y, z, xp, yp, zp);
}


double Star::rMax() {
  if (rmax_==DBL_MAX && i0_>=imin_ && i0_<=imax_) {
    size_t i;
    rmax_=x1_[i0_];
    int ck=gg_->coordKind();
    for (i=imin_;i<=imax_;++i) {
      if (x1_[i]>rmax_) rmax_=x1_[i];
      if (ck==GYOTO_COORDKIND_CARTESIAN) {
	if (x2_[i]>rmax_) rmax_=x2_[i];
	if (x3_[i]>rmax_) rmax_=x3_[i];
      }
    }
    rmax_ *= 3.;
  }
  return rmax_;
}
  
  
