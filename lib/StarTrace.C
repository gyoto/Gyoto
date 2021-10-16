/*
    Copyright 2013, 2018 Thibaut Paumard, Frederic Vincent

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
#include "GyotoProperty.h"
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

GYOTO_PROPERTY_START(StarTrace,
		     "All the points that would be inside a Star at any date between TMin and TMax.")
GYOTO_PROPERTY_DOUBLE(StarTrace, TMin, TMin,
		      "Date defining start of the trace (geometrical_time).")
GYOTO_PROPERTY_DOUBLE(StarTrace, TMax, TMax,
		      "Date defining end of the trace (geometrical_time).")
GYOTO_PROPERTY_END(StarTrace, Star::properties)

StarTrace::StarTrace() : Star(), tmin_(0.), tmax_(0.)
{
  Generic::kind_="StarTrace";
  xAllocateXYZ();
# ifdef GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << "done." << endl;
# endif
}

StarTrace::StarTrace(SmartPointer<Metric::Generic> met, double rad,
		double const pos[4],
		double const v[3]) :
  Star(met, rad, pos, v)
{
  Generic::kind_="StarTrace";
  xAllocateXYZ();
  computeXYZ(i0_);
}

StarTrace::StarTrace(const StarTrace& o) : Star(o), tmin_(o.tmin_), tmax_(o.tmax_)
{
  Generic::kind_="StarTrace";
  xAllocateXYZ();
  size_t sz = get_nelements()*sizeof(double);
  memcpy(x_+imin_, o.x_+imin_, sz);
  memcpy(y_+imin_, o.y_+imin_, sz);
  memcpy(z_+imin_, o.z_+imin_, sz);
}

StarTrace::StarTrace(const Star& o, double tmin, double tmax) :
  Star(o), tmin_(tmin), tmax_(tmax)
{
  Generic::kind_="StarTrace";
  xAllocateXYZ();
  computeXYZ();
}

StarTrace* StarTrace::clone() const { return new StarTrace(*this); }

StarTrace::~StarTrace()
{
  GYOTO_DEBUG << endl;
  delete[] x_;
  delete[] y_;
  delete[] z_;
}

void StarTrace::xAllocateXYZ() {
  x_ = new double[x_size_];
  y_ = new double[x_size_];
  z_ = new double[x_size_];
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG_EXPR(x_size_);
# endif
}

void StarTrace::xAllocate(size_t sz)
{
  Star::xAllocate(sz);
  xAllocateXYZ();
}

size_t StarTrace::xExpand(int dir) {
  
  xExpand(x_, dir);
  xExpand(y_, dir);
  xExpand(z_, dir);
  return Star::xExpand(dir);
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG_EXPR(x_size_);
# endif
}

void StarTrace::computeXYZ(size_t i)
{
  if (!gg_) GYOTO_ERROR("Please set metric before calling computeXYZ");
  switch (gg_->coordKind()) {
  case GYOTO_COORDKIND_SPHERICAL: 
    x_[i]=x1_[i]*sin(x2_[i])*cos(x3_[i]);
    y_[i]=x1_[i]*sin(x2_[i])*sin(x3_[i]);
    z_[i]=x1_[i]*cos(x2_[i]);
    break;
  case GYOTO_COORDKIND_CARTESIAN:
    x_[i]=x1_[i];
    y_[i]=x2_[i];
    z_[i]=x3_[i];
    break;
  default: GYOTO_ERROR("in StarTrace::computeXYZ: Incompatible coordinate kind");
  }
}

void StarTrace::computeXYZ()
{
  size_t n;
  int coordkind = gg_ -> coordKind();
  switch(coordkind) {
 case GYOTO_COORDKIND_SPHERICAL: 
    for (n=imin_;n<=imax_;++n) {
      x_[n]=x1_[n]*sin(x2_[n])*cos(x3_[n]);
      y_[n]=x1_[n]*sin(x2_[n])*sin(x3_[n]);
      z_[n]=x1_[n]*cos(x2_[n]);
    }
    break;
  case GYOTO_COORDKIND_CARTESIAN:
    for (n=imin_;n<=imax_;++n) {
      x_[n]=x1_[n];
      y_[n]=x2_[n];
      z_[n]=x3_[n];
    }
    break;
  default: GYOTO_ERROR("in StarTrace::computeXYZ(): Incompatible coordinate kind");
  }
}

void StarTrace::xStore(size_t ind, state_t const &coord, double tau)
{
  Star::xStore(ind, coord, tau);
  computeXYZ(ind);
}

void StarTrace::setInitCoord(const double coord[8], int dir)
{
  Star::setInitCoord(coord, dir);
  computeXYZ();
}

void StarTrace::metric(SmartPointer<Metric::Generic> gg)
{
  Star::metric(gg);
  computeXYZ();
}

string StarTrace::className() const { return  string("StarTrace"); }
string StarTrace::className_l() const { return  string("startrace"); }

void StarTrace::setInitialCondition(double const coord[8]) {
  Star::setInitialCondition(coord);
}

double StarTrace::TMin() const { return tmin_; }
void StarTrace::TMin(double t)
{
  if (t>tmax_) {
    tmin_=tmax_;
    tmax_=t;
  } else tmin_=t;
  GYOTO_DEBUG_EXPR(tmin_);
  GYOTO_DEBUG_EXPR(tmax_);
}

double StarTrace::TMax() const { return tmax_; }
void StarTrace::TMax(double t)
{
  if (t<tmin_) {
    tmax_=tmin_;
    tmin_=t;
  } else tmax_=t;
  GYOTO_DEBUG_EXPR(tmin_);
  GYOTO_DEBUG_EXPR(tmax_);
}

double StarTrace::operator()(double const coord[]) {
  double d2 = DBL_MAX, tmp;
  double ncoord[4];
  memcpy(ncoord, coord, 4*sizeof(double));
  xFill(tmin_, false);
  xFill(tmax_, false);

  double x=0., y=0., z=0.;
  switch (gg_->coordKind()) {
  case GYOTO_COORDKIND_SPHERICAL: 
    x=coord[1]*sin(coord[2])*cos(coord[3]);
    y=coord[1]*sin(coord[2])*sin(coord[3]);
    z=coord[1]*cos(coord[2]);
    break;
  case GYOTO_COORDKIND_CARTESIAN:
    x=coord[1];
    y=coord[2];
    z=coord[3];
    break;
  default: GYOTO_ERROR("in StarTrace::operator()(): Incompatible coordinate kind");
  }
  
  double tmp1;

  for (size_t i=imin_; i<=imax_; ++i) {
    if (x0_[i]<tmin_ || x0_[i]>tmax_) continue;
    tmp1 = x-x_[i];
    tmp  = tmp1 * tmp1;
    tmp1 = y-y_[i];
    tmp += tmp1 * tmp1;
    tmp1 = z-z_[i];
    tmp += tmp1 * tmp1;
    if (tmp < d2) d2=tmp;
  }
  return d2;
}
