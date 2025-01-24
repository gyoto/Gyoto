/*
    Copyright 2024 Nicolas Aimar


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
#include "GyotoSimThickDisk.h"

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
GYOTO_PROPERTY_START(SimThickDisk, "Synchrotron-emitting orbiting blob of plasma")
GYOTO_PROPERTY_DOUBLE(SimThickDisk, HoverR, HoverR)
GYOTO_PROPERTY_END(SimThickDisk, SimBridge::properties)

SimThickDisk::SimThickDisk() : 
  Astrobj::SimBridge(),
  HoverR_(0.)
{
  kind_="SimThickDisk";
# ifdef GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << "done." << endl;
# endif
}

SimThickDisk::SimThickDisk(const SimThickDisk& orig) : 
  Astrobj::SimBridge(orig),
  HoverR_(orig.HoverR_)
{
}

SimThickDisk* SimThickDisk::clone() const { return new SimThickDisk(*this); }

SimThickDisk::~SimThickDisk() {
  if (debug()) cerr << "DEBUG: SimThickDisk::~SimThickDisk()\n";
}

string SimThickDisk::className() const { return  string("SimThickDisk"); }
string SimThickDisk::className_l() const { return  string("SimThickDisk"); }

void Gyoto::Astrobj::SimThickDisk::HoverR(double hh){
  HoverR_=hh;
}

double Gyoto::Astrobj::SimThickDisk::HoverR() const{
    return HoverR_;
}

double SimThickDisk::operator()(double const coord[4]) {
  // zpos: modulus of altitude above equatorial plane
  // rproj: radius projected in the equatorial plane
  double zpos=0., rproj=0.;
  
  switch (gg_ -> coordKind()) {
  case GYOTO_COORDKIND_SPHERICAL:
    rproj  = coord[1]*sin(coord[2]);
    zpos  = fabs(coord[1]*cos(coord[2]));
    break;
  case GYOTO_COORDKIND_CARTESIAN:
    zpos  = fabs(coord[3]);
    rproj  = sqrt(coord[1]*coord[1]+coord[2]*coord[2]);
    break;
  default:
    GYOTO_ERROR("SimBridge::operator(): unknown COORDKIND");
  }
  double zdisk = HoverR_*rproj; 
  return zpos - zdisk; // >0 outside, <0 inside flared disk 
}

void SimThickDisk::filename(std::string const &f){
  SimBridge::filename(f);

  if (!gg_)
    GYOTO_ERROR("Define the metric in the Astrobj before giving the filename.");
  
  double theta_lim, xmax, ymax, zmax, rproj_max;
  switch (gg_ -> coordKind()) {
  case GYOTO_COORDKIND_SPHERICAL:
    theta_lim = abs(x2_array_[0]);
    HoverR(theta_lim);
    break;
  case GYOTO_COORDKIND_CARTESIAN:
    xmax = max(abs(x1_array_[0]), abs(x1_array_[nx1_-1]));
    ymax = max(abs(x2_array_[0]), abs(x2_array_[nx2_-1]));
    zmax = max(abs(x3_array_[0]), abs(x3_array_[nx3_-1]));
    rproj_max = sqrt(xmax*xmax+ymax*ymax);
    HoverR(zmax/rproj_max);
    break;
  default:
    GYOTO_ERROR("SimBridge::operator(): unknown COORDKIND");
  }
}