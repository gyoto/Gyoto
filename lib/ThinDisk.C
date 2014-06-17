/*
    Copyright 2011 Thibaut Paumard

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

#include "GyotoPhoton.h"
#include "GyotoThinDisk.h"
#include "GyotoUtils.h"
#include "GyotoFactoryMessenger.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <cmath>
#include <limits>
#include <string>

using namespace std;
using namespace Gyoto;
using namespace Gyoto::Astrobj;

ThinDisk::ThinDisk(std::string kin) :
  Generic(kin), rin_(0.), rout_(DBL_MAX), thickness_(1e-3), dir_(1)
{
  GYOTO_DEBUG << "ThinDisk Construction" << endl;
}

ThinDisk::ThinDisk(const ThinDisk& o) :
  Generic(o), Functor::Double_constDoubleArray(o), rin_(o.rin_), rout_(o.rout_),
  thickness_(o.thickness_), dir_(o.dir_)
{
  GYOTO_DEBUG << "ThinDisk Copy" << endl;
}
ThinDisk* ThinDisk::clone() const
{ return new ThinDisk(*this); }

ThinDisk::~ThinDisk() {
  GYOTO_DEBUG << "ThinDisk Destruction" << endl;
}

double ThinDisk::getInnerRadius() const   { return rin_; }
double ThinDisk::getInnerRadius(string unit) const   {
  return Units::FromGeometrical(getInnerRadius(), unit, gg_);
}
void   ThinDisk::setInnerRadius(double r) { rin_ = r;    }
void   ThinDisk::setInnerRadius(double r, string unit) {
  setInnerRadius(Units::ToGeometrical(r, unit, gg_));
}


double ThinDisk::getOuterRadius() const   { return rout_;}
double ThinDisk::getOuterRadius(string unit) const   {
  return Units::FromGeometrical(getOuterRadius(), unit, gg_);
}
void   ThinDisk::setOuterRadius(double r) { rout_ = r;   }
void   ThinDisk::setOuterRadius(double r, string unit) {
  setOuterRadius(Units::ToGeometrical(r, unit, gg_));
}

double ThinDisk::getThickness() const     { return thickness_;}
double ThinDisk::getThickness(string unit) const   {
  return Units::FromGeometrical(getThickness(), unit, gg_);
}
void   ThinDisk::setThickness(double h)   { thickness_ = h;   }
void   ThinDisk::setThickness(double h, string unit)   {
  setThickness(Units::ToGeometrical(h, unit, gg_));
}

int    ThinDisk::getDir() const           { return dir_; }
void   ThinDisk::setDir(int dir)          { dir_ = dir;  }

double ThinDisk::operator()(double const coord[4])  {
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
    throwError("ThinDisk::Impact(): unknown COORDKIND");
    return 0.;
  }
}

double ThinDisk::projectedRadius(double const coord[4]) const {
  switch (gg_ -> coordKind()) {
  case GYOTO_COORDKIND_SPHERICAL:
    return coord[1];
  case GYOTO_COORDKIND_CARTESIAN:
    return sqrt(coord[1]*coord[1]+coord[2]*coord[2]);
  default:
    throwError("ThinDisk::projectedRadius(): unknown COORDKIND");
    return 0.;
  }
}

double ThinDisk::sphericalPhi(double const coord[4]) const {
  switch (gg_ -> coordKind()) {
  case GYOTO_COORDKIND_SPHERICAL:
    return coord[3];
  case GYOTO_COORDKIND_CARTESIAN:
    return atan2(coord[2], coord[1]);
  default:
    throwError("ThinDisk::sphericalPhi(): unknown COORDKIND");
    return 0.;
  }
}

void ThinDisk::getVelocity(double const pos[4], double vel[4]) {
  gg_ -> circularVelocity(pos, vel, dir_);
}

int ThinDisk::Impact(Photon *ph, size_t index,
			       Astrobj::Properties *data) {
  double coord_ph_hit[8], coord_obj_hit[8];
  double rcross;
  double coord1[8], coord2[8];
  double dt=0.;
  ph->getCoord(index, coord1);
  ph->getCoord(index+1, coord2);

  if (gg_ -> coordKind() == GYOTO_COORDKIND_SPHERICAL){
    //Allows theta and phi to be in the correct range
    ph->checkPhiTheta(coord1);
    ph->checkPhiTheta(coord2);
  }
  

  if (gg_ -> coordKind() == GYOTO_COORDKIND_SPHERICAL &&
      fabs(coord2[2]-coord1[2]) > M_PI)
    throwError ("ThinDisk::Impact: fishy heuristic");

  double h1=operator()(coord1), h2=operator()(coord2);
  double r1=projectedRadius(coord1), r2=projectedRadius(coord2);

  if ( r1 > 2.*rout_ && r2 > 2.*rout_) return 0;
  if ( h1 == h2 && h2 != 0 ) return 0;
  if ( (h1 > 0.) == (h2 > 0.) && h1 != 0. && h2 != 0. ) return 0;
  
  double tlow, thigh;
  if (h1 < h2) {
    tlow = coord1[0]; thigh = coord2[0];
  } else {
    tlow = coord2[0]; thigh = coord1[0];
  }
  ph -> findValue(this, 0., tlow, thigh);
  coord_ph_hit[0]=thigh;

  ph -> getCoord(coord_ph_hit, 1, coord_ph_hit+1, coord_ph_hit+2,
		 coord_ph_hit+3, coord_ph_hit+4, coord_ph_hit+5,
		 coord_ph_hit+6, coord_ph_hit+7);

  if ((rcross=projectedRadius(coord_ph_hit)) < rin_ ||
      rcross > rout_) return 0;

  for (int i=0;i<4;i++) coord_obj_hit[i]=coord_ph_hit[i];
  getVelocity(coord_obj_hit, coord_obj_hit+4);

  if (flag_radtransf_) {
    double vel[3];
    gg_->cartesianVelocity(coord_ph_hit, vel);
    dt = (vel[2]==0.)
      ? (coord2[0] - coord1[0]) 
      : sqrt(1.+(vel[0]*vel[0]+vel[1]*vel[1])/(vel[2]*vel[2]))*thickness_;
  }

  if (data) {
    //Store impact time in user1
    if (data->user1) *data->user1=coord_ph_hit[0];
    if (data->user2) *data->user2=coord_ph_hit[1];
    if (data->user3) *data->user3=coord_ph_hit[3];
  }

  processHitQuantities(ph, coord_ph_hit, coord_obj_hit, dt, data);

  return 1;
}

int ThinDisk::setParameter(std::string name,
			   std::string content,
			   std::string unit) {
    char* tc = const_cast<char*>(content.c_str());
    if      (name=="InnerRadius")     setInnerRadius (atof(tc), unit); 
    else if (name=="OuterRadius")     setOuterRadius (atof(tc), unit); 
    else if (name=="Thickness")       setThickness   (atof(tc), unit); 
    else if (name=="CounterRotating") setDir         (-1);
    else return Generic::setParameter(name, content, unit);
    return 0;
}

#ifdef GYOTO_USE_XERCES
void ThinDisk::fillElement(FactoryMessenger *fmp) const {
  GYOTO_DEBUG <<"InnerRadius" << endl;
  if (rin_!=0.) fmp->setParameter("InnerRadius", rin_);
  GYOTO_DEBUG <<"OuterRadius" << endl;
  if (rout_!=DBL_MAX) fmp->setParameter("OuterRadius", rout_);
  GYOTO_DEBUG <<"Thickness" << endl;
  if (flag_radtransf_) fmp->setParameter("Thickness", thickness_);
  GYOTO_DEBUG <<"Dir" << endl;
  if (dir_==-1) fmp -> setParameter("CounterRotating");
  GYOTO_DEBUG <<"Generic" << endl;
  Generic::fillElement(fmp);
  GYOTO_DEBUG <<"done" << endl;
}
#endif
