/*
    Copyright 2013 Thibaut Paumard, Frederic Vincent

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
#include "GyotoStarTrace.h"
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

StarTrace::StarTrace() : Star() { GYOTO_DEBUG << endl; }

StarTrace::StarTrace(SmartPointer<Metric::Generic> met, double rad,
		double pos[4],
		double v[3]) :
  Star(met, rad, pos, v)
{}

StarTrace::StarTrace(const StarTrace& o) : Star(o), tmin_(o.tmin_), tmax_(o.tmax_) {}

StarTrace* StarTrace::clone() const { return new StarTrace(*this); }

StarTrace::~StarTrace() { GYOTO_DEBUG << endl; }

string StarTrace::className() const { return  string("StarTrace"); }
string StarTrace::className_l() const { return  string("startrace"); }

void StarTrace::setInitialCondition(double coord[8]) {
  Star::setInitialCondition(coord);
}

int StarTrace::setParameter(string name, string content, string unit) {
  if (name=="TMin") tmin_=atof(content.c_str());
  else if (name=="TMax") tmax_=atof(content.c_str());
  else return Star::setParameter(name, content, unit);
  return 0;
}

#ifdef GYOTO_USE_XERCES
void StarTrace::fillElement(FactoryMessenger *fmp) const {
  Star::fillElement(fmp);
  fmp->setParameter("TMin", tmin_);
  fmp->setParameter("TMax", tmax_);
}
#endif

double StarTrace::TMin() { return tmin_; }
void StarTrace::TMin(double t) { tmin_=t; }
double StarTrace::TMax() { return tmax_; }
void StarTrace::TMax(double t) { tmax_=t; }

double StarTrace::operator()(double const coord[]) {
  double d2 = DBL_MAX, tmp;
  double ncoord[4];
  memcpy(ncoord, coord, 4*sizeof(double));
  xFill(tmin_);
  xFill(tmax_);
  for (size_t i=imin_; i<=imax_; ++i) {
    if (x0_[i]<tmin_ || x0_[i]>tmax_) continue;
    ncoord[0]=x0_[i];
    if ( (tmp=Star::operator()(ncoord)) < d2) d2=tmp;
  }
  return d2;
}
