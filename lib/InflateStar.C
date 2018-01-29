/*
    Copyright 2011, 2018 Thibaut Paumard, Frederic Vincent

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
#include "GyotoInflateStar.h"
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
#include "GyotoProperty.h"
GYOTO_PROPERTY_START(InflateStar, "Star with inflation")
GYOTO_PROPERTY_DOUBLE_UNIT(InflateStar, TimeInflateInit, timeInflateInit,
			   "Start time of inflation (geometrical units)")
GYOTO_PROPERTY_DOUBLE_UNIT(InflateStar, TimeInflateStop, timeInflateStop,
			   "End time of inflation (geometrical units)")
GYOTO_PROPERTY_DOUBLE_UNIT(InflateStar, RadiusStop, radiusStop,
			   "End radius (geometrical units)")
GYOTO_PROPERTY_END(InflateStar, Star::properties)

GYOTO_PROPERTY_ACCESSORS(InflateStar, double,
			 timeinflateinit_, timeInflateInit)
double InflateStar::timeInflateInit(const string &unit) const {
  return Units::FromGeometricalTime(timeInflateInit(), unit, gg_);
}
void InflateStar::timeInflateInit(double t, const string &unit) {
  timeInflateInit(Units::ToGeometricalTime(t, unit, gg_));
}

GYOTO_PROPERTY_ACCESSORS(InflateStar, double,
			 timeinflatestop_, timeInflateStop)
double InflateStar::timeInflateStop(const string &unit) const {
  return Units::FromGeometricalTime(timeInflateStop(), unit, gg_);
}
void InflateStar::timeInflateStop(double t, const string &unit) {
  timeInflateStop(Units::ToGeometricalTime(t, unit, gg_));
}

GYOTO_PROPERTY_ACCESSORS_GEOMETRICAL(InflateStar, radiusstop_, radiusStop, gg_)

double InflateStar::radiusAt(double time) const {
  double radinit = radius();
  double radcur=radinit;
  if (time>=timeinflatestop_) radcur=radiusstop_;
  else if (time>timeinflateinit_){
    // linear increase of radius in between extreme times
    radcur = radinit+(time-timeinflateinit_)/(timeinflatestop_-timeinflateinit_)
      *(radiusstop_-radinit);
  }
  return radcur;
}

double InflateStar::radiusAt(double time, const string &t_unit) const {
  return radiusAt(Units::ToGeometricalTime(time, t_unit, gg_));
}

double InflateStar::radiusAt(double time, const string &t_unit, const string &r_unit) const {
  return Units::FromGeometrical(radiusAt(time, t_unit), r_unit, gg_);
}

InflateStar::InflateStar() :
  Star(),
  timeinflateinit_(0.), timeinflatestop_(0.), radiusstop_(DBL_MAX)
{
  kind_="InflateStar";
# ifdef GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << "done." << endl;
# endif
}

InflateStar::InflateStar(const InflateStar& orig) :
  Star(orig),
  timeinflateinit_(orig.timeinflateinit_), timeinflatestop_(orig.timeinflatestop_),
  radiusstop_(orig.radiusstop_)
{
}

InflateStar* InflateStar::clone() const { return new InflateStar(*this); }

InflateStar::~InflateStar() {
  if (debug()) cerr << "DEBUG: InflateStar::~InflateStar()\n";
}

string InflateStar::className() const { return  string("InflateStar"); }
string InflateStar::className_l() const { return  string("inflate_star"); }

int InflateStar::Impact(Gyoto::Photon* ph, size_t index,
			 Astrobj::Properties *data) {
  state_t p1;
  ph->getCoord(index, p1);
  double time = p1[0];
  double radinit = radius();
  double radcur=radiusAt(time);

  critical_value_=radcur*radcur;

  return UniformSphere::Impact(ph,index,data);
}

double InflateStar::emission(double nu_em, double dsem, 
			     state_t const &coord_ph, double const coord_obj[8]) const {
  double time = coord_ph[0];
  double radinit = radius();
  double radcur=radiusAt(time);
  double volume=radcur*radcur*radcur;
  double volumei=radinit*radinit*radinit;

  return volumei/volume * UniformSphere::emission(nu_em, dsem, coord_ph, coord_obj);

}
