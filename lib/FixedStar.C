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
#include "GyotoPhoton.h"
#include "GyotoFixedStar.h"
#include "GyotoProperty.h"
#include "GyotoFactoryMessenger.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <cmath>
#include <limits>
#include <cstring>
#include <float.h>

using namespace std;
using namespace Gyoto;
using namespace Gyoto::Astrobj;

GYOTO_PROPERTY_START(FixedStar,
		     "Coordinate-spherical blob with fixed centre coordinates.")
GYOTO_PROPERTY_VECTOR_DOUBLE(FixedStar, Position, position,
			     "Space coordinates (3 components).")
GYOTO_PROPERTY_BOOL(FixedStar, Rotating, NonRotating, rotating,
		    "Is fluid at rest or in circular rotation in coordinate system.")
GYOTO_PROPERTY_END(FixedStar, UniformSphere::properties)

FixedStar::FixedStar() : UniformSphere("FixedStar"), rotating_(false)
{
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << endl;
# endif
  for (int i=0;i<3;++i) pos_[i]=0.;
}

FixedStar::FixedStar(SmartPointer<Gyoto::Metric::Generic> gg, double StPsn[3],
		     double rad) :
  UniformSphere("FixedStar", gg, rad), rotating_(false)
{
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << "(metric, pos, rad)" << endl;
# endif
  for (int i=0;i<3;++i) pos_[i] = StPsn[i]; 
  radius(rad);
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << "done" << endl;
# endif
}

FixedStar::FixedStar(const FixedStar& orig) :
  UniformSphere(orig), rotating_(orig.rotating_)
{
  for (int i=0; i<3; ++i) pos_[i] = orig.pos_[i];
}
FixedStar* FixedStar::clone() const { return new FixedStar(*this); }

FixedStar::~FixedStar() {
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << endl;
# endif
}

void FixedStar::getCartesian(double const * const , size_t const n_dates,
			     double * const x, double * const y,
			     double * const z, double * const xprime,
			     double * const yprime, 
			     double * const zprime) {
  double xs, ys, zs;
  switch (gg_->coordKind()) {
  case GYOTO_COORDKIND_CARTESIAN:
    xs= pos_[0];
    ys= pos_[1];
    zs= pos_[2];
    break;
  case GYOTO_COORDKIND_SPHERICAL:
    {
      double rs=pos_[0];
      double ths=pos_[1];
      double phs=pos_[2];
      double st, ct, sp, cp;
      sincos(ths, &st, &ct);
      sincos(phs, &sp, &cp);
      xs= rs*st*cp;
      ys= rs*st*sp;
      zs= rs*ct;
    }
    break;
  default:
    GYOTO_ERROR("unsupported coordkind");
    xs=ys=zs=0.;
  }
  for (size_t i=0; i<n_dates; ++i) {
    if (x) x[i] = xs;
    if (y) y[i] = ys;
    if (z) z[i] = zs;
    if (xprime) xprime[i] = 0.;
    if (yprime) yprime[i] = 0.;
    if (zprime) zprime[i] = 0.;
  }
}

void FixedStar::getVelocity(double const pos[4], double vel[4]) {
  if (rotating_) gg_->circularVelocity(pos, vel);
  else {
    for (size_t i=0; i<4; ++i) vel[i]=0.;
    vel[0]=gg_->SysPrimeToTdot(pos, vel+1);
  }
}

void FixedStar::rotating(bool rot) {rotating_=rot;}
bool FixedStar::rotating() const { return rotating_; }

double const * FixedStar::getPos() const { return pos_; }

void FixedStar::getPos(double dst[3]) const
{ for (int i=0; i<3;++i) dst[i]=pos_[i]; }

void FixedStar::metric(SmartPointer<Metric::Generic> gg) {
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << endl;
# endif
 Generic::metric(gg);
 radius(radius_);
}

void FixedStar::radius(double r) {
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG_EXPR(r) ;
# endif
  UniformSphere::radius(r);
  if (!gg_()) {
#   if GYOTO_DEBUG_ENABLED
    GYOTO_DEBUG << "metric is not set yet" << endl;
#   endif
    return;
  }
}

double FixedStar::rMax() {
  if (rmax_==DBL_MAX) {
    switch (gg_ -> coordKind()) {
    case GYOTO_COORDKIND_SPHERICAL:
      rmax_=3.*(pos_[0]+radius_);
      break;
    case GYOTO_COORDKIND_CARTESIAN:
      rmax_=3.*(sqrt(pos_[0]*pos_[0]+pos_[1]*pos_[1]+pos_[2]*pos_[2])+radius_);
      break;
    default:
      GYOTO_ERROR("unimplemented coordinate system in FixedStar");
    } 
  }
  return rmax_;
}

void FixedStar::setPos(const double p[3])
{ for (int i=0; i<3; ++i) pos_[i]=p[i]; radius(radius_);}

void FixedStar::position(std::vector<double> const &v) {
  GYOTO_DEBUG_EXPR(v.size());
  if (v.size() !=3)
    GYOTO_ERROR("FixedStar position needs exactly 3 tokens"); 
  for (int i=0; i<3; ++i) pos_[i]=v[i];
  radius(radius_);
}

std::vector<double> FixedStar::position() const {
  std::vector<double> res(3, 0.);
  for (int i=0; i<3; ++i) res[i]=pos_[i];
  return res;
}
