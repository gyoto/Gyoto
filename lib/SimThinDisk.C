/*
    Copyright 2024 Irene Urso


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
#include "GyotoSimThinDisk.h"

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
GYOTO_PROPERTY_START(SimThinDisk, "Equatorial thin disk painted with exterior intensity")
GYOTO_PROPERTY_DOUBLE_UNIT(SimThinDisk, InnerRadius, innerRadius, "Inner radius (geometrical units, 0).")
GYOTO_PROPERTY_DOUBLE_UNIT(SimThinDisk, OuterRadius, outerRadius, "Outer radius (geometrical units, DBL_MAX).")
GYOTO_PROPERTY_DOUBLE_UNIT(SimThinDisk, Thickness, thickness, "Geometrical thickness (geometrical units, 1e-3, for optical depth).")
GYOTO_PROPERTY_END(SimThinDisk, SimBridge::properties)

SimThinDisk::SimThinDisk() : 
  Astrobj::SimBridge(), 
  rin_(0.), 
  rout_(DBL_MAX), 
  thickness_(1e-3)
{
  kind_="SimThinDisk";
# ifdef GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << "done." << endl;
# endif
}

SimThinDisk::SimThinDisk(const SimThinDisk& orig) : 
  Astrobj::SimBridge(orig), 
  rin_(orig.rin_), 
  rout_(orig.rout_), 
  thickness_(orig.thickness_)
{
}

SimThinDisk* SimThinDisk::clone() const { return new SimThinDisk(*this); }

SimThinDisk::~SimThinDisk() {
  if (debug()) cerr << "DEBUG: SimThinDisk::~SimThinDisk()\n";
}

string SimThinDisk::className() const { return  string("SimThinDisk"); }
string SimThinDisk::className_l() const { return  string("simThindisk"); }


double SimThinDisk::innerRadius() const   { return rin_; }
double SimThinDisk::innerRadius(string const &unit) const   {
  return Units::FromGeometrical(innerRadius(), unit, gg_);
}
void   SimThinDisk::innerRadius(double r) { rin_ = r;    }
void   SimThinDisk::innerRadius(double r, string const &unit) {
  innerRadius(Units::ToGeometrical(r, unit, gg_));
}


double SimThinDisk::outerRadius() const   { return rout_;}
double SimThinDisk::outerRadius(string const &unit) const   {
  return Units::FromGeometrical(outerRadius(), unit, gg_);
}
void   SimThinDisk::outerRadius(double r) { rout_ = r;   }
void   SimThinDisk::outerRadius(double r, string const &unit) {
  outerRadius(Units::ToGeometrical(r, unit, gg_));
}


double SimThinDisk::thickness() const     { return thickness_;}
double SimThinDisk::thickness(string const &unit) const   {
  return Units::FromGeometrical(thickness(), unit, gg_);
}
void   SimThinDisk::thickness(double h)   { thickness_ = h;   }
void   SimThinDisk::thickness(double h, string const &unit)   {
  thickness(Units::ToGeometrical(h, unit, gg_));
}


double SimThinDisk::operator()(double const coord[4]) {
  double theta;
  switch (gg_ -> coordKind()) {
  case GYOTO_COORDKIND_SPHERICAL:
    theta = coord[2];
    theta -= M_PI*0.5;
    while (theta < -M_PI) theta += 2.*M_PI;
    while (theta >= M_PI) theta -= 2.*M_PI;
    return theta;
  case GYOTO_COORDKIND_CARTESIAN:
    return coord[3];
  default:
    GYOTO_ERROR("SimThinDisk::Impact(): unknown COORDKIND");
    return 0.;
  }
}


double SimThinDisk::projectedRadius(double const coord[4]) const {
  switch (gg_ -> coordKind()) {
  case GYOTO_COORDKIND_SPHERICAL:
    return coord[1];
  case GYOTO_COORDKIND_CARTESIAN:
    return sqrt(coord[1]*coord[1]+coord[2]*coord[2]);
  default:
    GYOTO_ERROR("ThinDisk::projectedRadius(): unknown COORDKIND");
    return 0.;
  }
}


int SimThinDisk::Impact(Photon *ph, size_t index,
			       Astrobj::Properties *data) {
  state_t coord_ph_hit;
  double coord_obj_hit[8];
  double rcross;
  state_t coord1, coord2;
  double dt=0.;
  ph->getCoord(index, coord1);
  ph->getCoord(index+1, coord2);

  if (gg_ -> coordKind() == GYOTO_COORDKIND_SPHERICAL){
    //Allows theta and phi to be in the correct range
    ph->checkPhiTheta(&coord1[0]);
    ph->checkPhiTheta(&coord2[0]);
  }
  

  if (gg_ -> coordKind() == GYOTO_COORDKIND_SPHERICAL &&
      fabs(coord2[2]-coord1[2]) > M_PI)
    GYOTO_ERROR ("ThinDisk::Impact: fishy heuristic");

  double h1=operator()(&coord1[0]), h2=operator()(&coord2[0]);
  double r1=projectedRadius(&coord1[0]), r2=projectedRadius(&coord2[0]);

  if ( 0.5*r1 > rout_ && 0.5*r2 > rout_) return 0;
  if ( h1 == h2 && h2 != 0 ) return 0;
  if ( (h1 > 0.) == (h2 > 0.) && h1 != 0. && h2 != 0. ) return 0;
  
  double tlow, thigh;
  if (h1 < h2) {
    tlow = coord1[0]; thigh = coord2[0];
  } else {
    tlow = coord2[0]; thigh = coord1[0];
  }
  ph -> findValue(this, 0., tlow, thigh);

  ph -> getCoord(thigh, coord_ph_hit);

  if ((rcross=projectedRadius(&coord_ph_hit[0])) < rin_ ||
      rcross > rout_) return 0;

  for (int i=0;i<4;i++) coord_obj_hit[i]=coord_ph_hit[i];
  getVelocity(coord_obj_hit, coord_obj_hit+4);

  if (flag_radtransf_) {
    double vel[3];
    gg_->cartesianVelocity(&coord_ph_hit[0], vel);
    dt = (vel[2]==0.)
      ? (coord2[0] - coord1[0]) 
      : sqrt(1.+(vel[0]*vel[0]+vel[1]*vel[1])/(vel[2]*vel[2]))*thickness_;
  }

  int store_impact_coord=0; // put to 1 if storing coord is needed
  if (data && store_impact_coord==1) {
    //Store impact time in user1
    if (data->user1) *data->user1=coord_ph_hit[0];
    if (data->user2) *data->user2=coord_ph_hit[1];
    if (data->user3) *data->user3=coord_ph_hit[3];
  }

  processHitQuantities(ph, coord_ph_hit, coord_obj_hit, dt, data);

  return 1;
}
