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
  Worldline(),
  wait_pos_(0), init_vel_(NULL)
{
  if (debug())
    cerr << "DEBUG: in Star::Star()" << endl;
}

Star::Star(SmartPointer<Metric::Generic> met, double rad,
	   double pos[4],
	   double v[3]) :
  UniformSphere("Star"),
  Worldline(),
  wait_pos_(0), init_vel_(NULL)
{
  if (debug()) {
    cerr << "DEBUG: Star Construction " << endl
	 << "       POS=[" << pos[0];
    for (int i=1; i<4; ++i) cerr << ", " << pos[i];
    cerr << "]\n       VEL=[" << v[0] ;
    for (int i=1; i<3; ++i) cerr << ", " << v[i];
    cerr << "]\n       RADIUS=" << rad << endl;

  }

  setMetric(met);
  setInitCoord(pos, v);
  setRadius(rad);
}

Star::Star(const Star& orig) :
  UniformSphere(orig), Worldline(orig),
  wait_pos_(orig.wait_pos_), init_vel_(NULL)
{
  if (debug()) cerr << "Star copy" << endl;
  if (orig.init_vel_) {
    init_vel_ = new double [3];
    memcpy(init_vel_, orig.init_vel_, 3*sizeof(double));
  }
  gg_ = metric_; // we have two distinct clones of the metric, not good...
}

Star* Star::clone() const { return new Star(*this); }

Star::~Star() {
  if (debug()) cerr << "DEBUG: Star::~Star()\n";
  if (init_vel_) delete[] init_vel_;
}

string Star::className() const { return  string("Star"); }
string Star::className_l() const { return  string("star"); }

SmartPointer<Metric::Generic> Star::getMetric() const { return gg_; }
void Star::setMetric(SmartPointer<Metric::Generic> gg) {
  UniformSphere::setMetric(gg);
  Worldline::setMetric(gg);
}

void Star::setInitCoord(double pos[4], double v[3], int dir) {
  if (!metric_) throwError("Please set metric before calling Star::setInitCoord(double pos[4], double vel[3])");
  double tdot0=metric_->SysPrimeToTdot(pos, v);
  if (debug()) cerr << "DEBUG: Star::setInitCoord(): TDOT0=" << tdot0 << endl;
  double coord[8]={pos[0], pos[1], pos[2], pos[3],
		   tdot0, v[0]*tdot0, v[1]*tdot0, v[2]*tdot0};
  setInitCoord(coord, dir);
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


double Star::getRmax() {
  // rmax may not be up to date. Assume it's OK if it's non-0 warning:
  // if line is extended, need to either update rmax or reset it to 0
  // if (debug()) cerr << "DEBUG: Star::getRmax(): rmax_set_==" 
  //                   << rmax_set_ << endl;
  if (!rmax_set_ && !rmax_) {
    size_t i;
    for (i=imin_;i<=imax_;++i) if (x1_[i]>rmax_) rmax_=x1_[i];
    rmax_*=3.;
  }
  return rmax_;
}

void Star::unsetRmax() {
  rmax_set_=0;
  rmax_=DBL_MAX;
}

void Star::setPosition(double pos[4]) {
  double vel[] = {0., 0., 0.};
  setInitCoord(pos, vel);
}

void Star::setVelocity(double vel[3]) {
  double coord[8];
  getInitialCoord(coord);
  setInitCoord(coord, vel);
}

int Star::setParameter(std::string name, std::string content) {
  double coord[8];
  char* tc = const_cast<char*>(content.c_str());
  if (name=="InitialCoordinate") {
    for (int i=0;i<8;++i) coord[i] = strtod(tc, &tc);
    setInitCoord(coord);
  } else if (name=="Position") {
    for (int i=0;i<4;++i) coord[i] = strtod(tc, &tc);
    if (init_vel_) {
      setInitCoord(coord, init_vel_);
      delete[] init_vel_; init_vel_=NULL;
    } else setPosition(coord);
    wait_pos_ = 0;
  } else if (name=="Velocity") {
    for (int i=0;i<3;++i) coord[i] = strtod(tc, &tc);
    if (wait_pos_) {
      if (init_vel_) delete [] init_vel_;
      init_vel_ = new double[3];
      memcpy(init_vel_, coord, 3*sizeof(double));
    } else setVelocity(coord);
  } else return UniformSphere::setParameter(name, content);
  return 0;
}

#ifdef GYOTO_USE_XERCES
void Star::fillElement(FactoryMessenger *fmp) const {

  if (imin_ <= imax_) {
    double coord[8];
    getInitialCoord(coord);
    fmp -> setParameter ("Position", coord, 4);
    double vel[3] = {coord[5]/coord[4], coord[6]/coord[4], coord[7]/coord[4]};
    fmp -> setParameter ("Velocity", vel, 3);
  }

  Astrobj::UniformSphere::fillElement(fmp);
}

void Star::setParameters(FactoryMessenger* fmp) {
  wait_pos_ = 1;
  UniformSphere::setParameters(fmp);
  if (init_vel_) {
    delete[] init_vel_; init_vel_=NULL;
    throwError("Star::setParameters(): Velocity was found but not Position");
  }
  wait_pos_ = 0;
}
#endif
