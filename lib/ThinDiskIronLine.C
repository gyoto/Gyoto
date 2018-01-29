/*
  Copyright 2013, 2018 Frederic Vincent & Thibaut Paumard
  
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

#include "GyotoThinDiskIronLine.h"
#include "GyotoConverters.h"

#include <iostream>
#include <cmath>
#include <cfloat>
#include <cstdlib>

using namespace Gyoto;
using namespace Gyoto::Astrobj;
using namespace std;

/// Properties
#include "GyotoProperty.h"
GYOTO_PROPERTY_START(ThinDiskIronLine)
GYOTO_PROPERTY_DOUBLE(ThinDiskIronLine, PowerLawIndex, PowerLawIndex)
GYOTO_PROPERTY_DOUBLE_UNIT(ThinDiskIronLine, LineFreq, LineFreq)
GYOTO_PROPERTY_DOUBLE_UNIT(ThinDiskIronLine, CutRadius, CutRadius)
GYOTO_PROPERTY_END(ThinDiskIronLine, ThinDisk::properties)

// ACCESSORS
GYOTO_PROPERTY_ACCESSORS(ThinDiskIronLine, double, plindex_, PowerLawIndex)
GYOTO_PROPERTY_ACCESSORS(ThinDiskIronLine, double, cutradius_, CutRadius)
// Define the accessors with unit manually
void ThinDiskIronLine::CutRadius(double v, std::string const &u) {
  CutRadius(Units::ToGeometrical(v, u, gg_));
}
double ThinDiskIronLine::CutRadius(std::string const &u)const{
  return Units::FromGeometrical(CutRadius(), u, gg_);
}

// The following is not completely standard, let's implement it manually:
#define ___local_f 2.417989579752276e+17 //(1e3*1.60217657e-19/GYOTO_PLANCK);
void ThinDiskIronLine::LineFreq(double v) {linefreq_=v*___local_f;}
double ThinDiskIronLine::LineFreq()const{return linefreq_/___local_f;}
void ThinDiskIronLine::LineFreq(double v, std::string const &u) {
  LineFreq(Units::ToHerz(v, u));
}
double ThinDiskIronLine::LineFreq(std::string const &u)const{
  return Units::FromHerz(LineFreq(), u);
}
#undef ___local_f
///




Gyoto::Astrobj::ThinDiskIronLine::ThinDiskIronLine()
  : ThinDisk("ThinDiskIronLine"), plindex_(0.), linefreq_(0.), 
    cutradius_(-DBL_MAX)
{
  
  GYOTO_DEBUG << "Building ThinDiskIronLine" << endl;
}

Gyoto::Astrobj::ThinDiskIronLine::ThinDiskIronLine(const ThinDiskIronLine &o)
  : ThinDisk(o), plindex_(o.plindex_), linefreq_(o.linefreq_),
    cutradius_(o.cutradius_)
{
  GYOTO_DEBUG << "Copying ThinDiskIronLine" << endl;
}
ThinDiskIronLine * ThinDiskIronLine::clone() const { return new ThinDiskIronLine(*this); }

Gyoto::Astrobj::ThinDiskIronLine::~ThinDiskIronLine()
{
  GYOTO_DEBUG << "Destroying dummy ThinDiskIronLine" << endl;
}

double ThinDiskIronLine::emission(double nu_em, double /* dsem */,
				  state_t const &,
				  double const coord_obj[8]) const{
  double rr=projectedRadius(coord_obj);
  if (rr<cutradius_) return 0.;
  double dfreq=linefreq_/100.;
  /*
    NB: this choice of dfreq is related to the
    spectral resolution of e.g. CHANDRA, which is
    around E / (Delta E) = 100 at 6 keV,
    see e.g. iron line review of Reynolds 2003.
   */
  if (abs(nu_em-linefreq_)>dfreq) return 0.;
  
  double Iem = pow(rr,-plindex_);
  return Iem;
}

void ThinDiskIronLine::getVelocity(double const pos[4], double vel[4]) {
  if (projectedRadius(pos)<cutradius_){
    //any velocity, emission=0 anyway
    for (int ii=1;ii<4;ii++) vel[ii]=0;
    vel[0] = 1.;
    //vel[0] = gg_->SysPrimeToTdot(pos, vel+1); // leads to v>c if a>0...
  }else{
    ThinDisk::getVelocity(pos,vel);
  }
}


