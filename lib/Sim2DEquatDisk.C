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
#include "GyotoSim2DEquatDisk.h"

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
GYOTO_PROPERTY_START(Sim2DEquatDisk, "Synchrotron-emitting orbiting blob of plasma")
GYOTO_PROPERTY_END(Sim2DEquatDisk, SimBridge::properties)

Sim2DEquatDisk::Sim2DEquatDisk() : 
  Astrobj::SimBridge(),
  HoverR_(0.)
{
  kind_="Sim2DEquatDisk";
# ifdef GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << "done." << endl;
# endif
}

Sim2DEquatDisk::Sim2DEquatDisk(const Sim2DEquatDisk& orig) : 
  Astrobj::SimBridge(orig),
  HoverR_(orig.HoverR_)
{
}

Sim2DEquatDisk* Sim2DEquatDisk::clone() const { return new Sim2DEquatDisk(*this); }

Sim2DEquatDisk::~Sim2DEquatDisk() {
  if (debug()) cerr << "DEBUG: Sim2DEquatDisk::~Sim2DEquatDisk()\n";
}

string Sim2DEquatDisk::className() const { return  string("Sim2DEquatDisk"); }
string Sim2DEquatDisk::className_l() const { return  string("sim2dequatdisk"); }

void Gyoto::Astrobj::Sim2DEquatDisk::HoverR(double hh){
  HoverR_=hh;
}

double Gyoto::Astrobj::Sim2DEquatDisk::HoverR() const{
    return HoverR_;
}

double Sim2DEquatDisk::operator()(double const coord[4]) {
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