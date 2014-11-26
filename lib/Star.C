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
#include "GyotoStar.h"
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

Star::Star() :
  UniformSphere("Star"),
  Worldline()
{
# ifdef GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << "done." << endl;
# endif
}

Star::Star(SmartPointer<Metric::Generic> met, double rad,
	   double pos[4],
	   double v[3]) :
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

void Star::setInitialCondition(double coord[8]) {
  if (!metric_) throwError("Please set metric before calling Star::setInitialCondition(double*)");
  Worldline::setInitialCondition(metric_, coord, 1);
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
  if (rmax_==DBL_MAX) {
    size_t i;
    for (i=imin_;i<=imax_;++i) if (x1_[i]>rmax_) rmax_=x1_[i];
    rmax_ *= 3.;
  }
  return rmax_;
}

int Star::setParameter(string name, string content, string unit) {
  if        (!UniformSphere::setParameter(name, content, unit)) ; // if found
  else if   (!Worldline    ::setParameter(name, content, unit)) ; // do nothing
  else return 1;
  return 0;
}

#ifdef GYOTO_USE_XERCES
void Star::fillElement(FactoryMessenger *fmp) const {
  Worldline::fillElement(fmp);
  Astrobj::UniformSphere::fillElement(fmp);
}

void Star::setParameters(FactoryMessenger* fmp) {
  wait_pos_ = 1;
  UniformSphere::setParameters(fmp);
  wait_pos_ = 0;
  if (init_vel_) {
    delete[] init_vel_; init_vel_=NULL;
    throwError("Worldline::setParameters(): "
	       "Velocity was found but not Position");
  }
}
#endif
