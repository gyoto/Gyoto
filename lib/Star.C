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
  if (debug())
    cerr << "DEBUG: in Star::Star()" << endl;
}

Star::Star(SmartPointer<Metric::Generic> met, double rad,
	   double pos[4],
	   double v[3]) :
  UniformSphere("Star", met, rad),
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

  metric_=met;
  gg_=met;

  double tdot0=metric_->SysPrimeToTdot(pos, v);

  if (debug()) cerr << "       TDOT0=" << tdot0 << endl;

  double coord[8]={pos[0], pos[1], pos[2], pos[3],
		   tdot0, v[0]*tdot0, v[1]*tdot0, v[2]*tdot0};

  Worldline::setInitialCondition(metric_, coord, 1);
    //last number : direction of integration + or -1
}

Star::Star(const Star& orig) :
  UniformSphere(orig), Worldline(orig)
{
  if (debug()) cerr << "Star copy" << endl;
  gg_ = metric_; // we have two distinct clones of the metric, not good...
}

Star* Star::clone() const { return new Star(*this); }

Star::~Star() {
  if (debug()) cerr << "DEBUG: Star::~Star()\n";
}

string Star::className() const { return  string("Star"); }
string Star::className_l() const { return  string("star"); }

SmartPointer<Metric::Generic> Star::getMetric() const { return gg_; }
void Star::setMetric(SmartPointer<Metric::Generic> gg) {gg_=gg; metric_=gg;}

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

SmartPointer<Astrobj::Generic> Gyoto::Astrobj::Star::Subcontractor(FactoryMessenger* fmp) {

  string name="", content="";
  int pos_found=0, vel_found=0;
  double pos[4], v[3];
  SmartPointer<Metric::Generic> gg = NULL;
  SmartPointer<Spectrum::Generic> sp = NULL, op = NULL;
  FactoryMessenger * child = NULL;

  gg = fmp->getMetric();

  while (fmp->getNextParameter(&name, &content)) {
    char* tc = const_cast<char*>(content.c_str());
    if      (name=="Position") {
      pos_found=1;
      for (int i=0;i<4;++i) pos[i] = strtod(tc, &tc);
    }
    else if (name=="Velocity") {
      vel_found=1;
      for (int i=0;i<3;++i) v[i] = strtod(tc, &tc);
    }
  }
  if (!pos_found) throwError("Position MUST be set in Star definition");
  if (!vel_found) throwError("Velocity MUST be set in Star definition");

  SmartPointer<Star> st = new Star(gg, 0., pos, v);

  fmp->reset();
  st -> setGenericParameters(fmp);

  return st;
}

void Gyoto::Astrobj::Star::Init() {
  Gyoto::Astrobj::Register("Star", &Gyoto::Astrobj::Star::Subcontractor);
}
#endif
