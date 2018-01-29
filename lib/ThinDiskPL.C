/*
    Copyright 2012, 2014, 2016, 2018 Frederic Vincent, Thibaut Paumard

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
#include "GyotoThinDiskPL.h"
#include "GyotoProperty.h"
#include "GyotoUtils.h"
#include "GyotoFactoryMessenger.h"
#include "GyotoKerrBL.h"
#include "GyotoKerrKS.h"


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

GYOTO_PROPERTY_START(ThinDiskPL)
GYOTO_PROPERTY_DOUBLE(ThinDiskPL, Slope, Slope)
GYOTO_PROPERTY_DOUBLE(ThinDiskPL, Tinner, Tinner)
GYOTO_PROPERTY_END(ThinDiskPL, ThinDisk::properties)

// ACCESSORS
void ThinDiskPL::Slope(double alpha) {slope_=alpha;}
double ThinDiskPL::Slope()const{return slope_;}
void ThinDiskPL::Tinner(double TT) {Tinner_=TT;}
double ThinDiskPL::Tinner()const{return Tinner_;}
//

ThinDiskPL::ThinDiskPL() :
  ThinDisk("ThinDiskPL"),
  slope_(0.), Tinner_(1.),
  spectrumBB_(NULL)
{
  if (debug()) cerr << "DEBUG: ThinDiskPL Construction" << endl;
  spectrumBB_ = new Spectrum::BlackBody(); 
}

ThinDiskPL::ThinDiskPL(const ThinDiskPL& o) :
  ThinDisk(o),
  slope_(o.slope_), Tinner_(o.Tinner_), 
  spectrumBB_(NULL)
{
  if (o.gg_()) gg_=o.gg_->clone();
  if (o.spectrumBB_()) spectrumBB_=o.spectrumBB_->clone();
  Generic::gg_=gg_;
}
ThinDiskPL* ThinDiskPL::clone() const
{ return new ThinDiskPL(*this); }

bool ThinDiskPL::isThreadSafe() const {
  return ThinDisk::isThreadSafe()
    && (!spectrumBB_ || spectrumBB_ -> isThreadSafe());
}

ThinDiskPL::~ThinDiskPL() {
  if (debug()) cerr << "DEBUG: ThinDiskPL Destruction" << endl;
}

double ThinDiskPL::emission(double nu, double,
			    state_t const &,
			    double const coord_obj[8]) const{  
  double rcur=projectedRadius(coord_obj);
  double TT = Tinner_*pow(rcur/rin_,slope_);
  //  cout << "In ThinPL rin, slope, Tinner, TT= " << rin_ << " " << slope_ << " " << Tinner_ << " " << TT << endl;
  spectrumBB_->temperature(TT);
  return (*spectrumBB_)(nu);
}
