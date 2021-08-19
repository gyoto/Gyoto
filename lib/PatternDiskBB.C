/*
    Copyright 2012-2014, 2016, 2018 Frederic Vincent, Thibaut Paumard

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
#include "GyotoPatternDiskBB.h"
#include "GyotoProperty.h"
#include "GyotoUtils.h"
#include "GyotoFactoryMessenger.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <fstream>
#include <cstring>
#include <cmath>
#include <limits>

using namespace std;
using namespace Gyoto;
using namespace Gyoto::Astrobj;

GYOTO_PROPERTY_START(PatternDiskBB)
GYOTO_PROPERTY_BOOL(PatternDiskBB,
		    SpectralEmission, BolometricEmission, spectralEmission)
GYOTO_PROPERTY_END(PatternDiskBB, PatternDisk::properties)

bool PatternDiskBB::spectralEmission() const {return SpectralEmission_;}
void PatternDiskBB::spectralEmission(bool t) {SpectralEmission_=t;}

PatternDiskBB::PatternDiskBB() :
  PatternDisk(),
  spectrumBB_(NULL),
  SpectralEmission_(0)
{
  kind_= "PatternDiskBB";
  GYOTO_DEBUG << "PatternDiskBB Construction" << endl;
  spectrumBB_ = new Spectrum::BlackBody(); 
}

PatternDiskBB::PatternDiskBB(const PatternDiskBB& o) :
  PatternDisk(o),
  spectrumBB_(NULL),
  SpectralEmission_(o.SpectralEmission_)
{
  GYOTO_DEBUG << "PatternDiskBB Copy" << endl;
  if (o.spectrumBB_()) spectrumBB_=o.spectrumBB_->clone();
}
PatternDiskBB* PatternDiskBB::clone() const
{ return new PatternDiskBB(*this); }

bool PatternDiskBB::isThreadSafe() const {
  // spectrumBB_ is not a Property
  return PatternDisk::isThreadSafe()
    && (!spectrumBB_ || spectrumBB_->isThreadSafe());
}

PatternDiskBB::~PatternDiskBB() {
  GYOTO_DEBUG << "PatternDiskBB Destruction" << endl;
}

double PatternDiskBB::emission(double nu, double dsem,
			       state_t const &cp,
			       double const co[8]) const{
  // This reimplementation of emission has 2 goals, ensuring that
  // nothing happens below ISCO, and allowing to compute BB emission
  // when the PatternDisk structure contains temperature
  GYOTO_DEBUG << endl;
  
  double Iem=0.;

  if (!SpectralEmission_){
    /*
      Here the intensity is assumed to be given by the emission
      function of PatternDisk
     */
    Iem = PatternDisk::emission(nu,dsem,cp,co);//NB: 3rd argument is not used
  }else{ //Spectral emission    
    /*
      Here the temperature is assumed to be given by the emission
      function of PatternDisk, intensity is computed assuming BB
     */
    double TT;
    TT = PatternDisk::emission(nu,dsem,cp,co);
    if (TT==0.) Iem=0.; // typically: we are outside grid radial range
    else{
      spectrumBB_->temperature(TT);
      Iem=(*spectrumBB_)(nu);
    }
    //cout << "In pattern BB nu, T, Bnu= " << nu << " " << TT << " " << Iem << endl;
  }

  if (!flag_radtransf_) return Iem;
  else GYOTO_ERROR("In PatternDiskBB::emission: should be optically thick!");
  // The PatternDisk::emission function called above will return
  // nonsense for the temperature in case the object is optically thin.

  // Should never reach this
  return 0.;
}
