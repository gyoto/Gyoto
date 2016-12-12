/*
    Copyright 2012-2014, 2016 Frederic Vincent, Thibaut Paumard

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
#include "GyotoKerrBL.h"
#include "GyotoKerrKS.h"

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
GYOTO_PROPERTY_DOUBLE(PatternDiskBB, Risco, risco)
GYOTO_PROPERTY_END(PatternDiskBB, PatternDisk::properties)

bool PatternDiskBB::spectralEmission() const {return SpectralEmission_;}
void PatternDiskBB::spectralEmission(bool t) {SpectralEmission_=t;}

double PatternDiskBB::risco() const {
  if (risco_>0.) return risco_;
  switch (gg_->coordKind()) {
  case GYOTO_COORDKIND_SPHERICAL:
    return static_cast<SmartPointer<Metric::KerrBL> >(gg_) -> getRms();
  default:
    throwError("PatternDiskBB::getVelocity: bad COORDKIND");
  }
  return 0.; // avoid warning, never reached
}
void PatternDiskBB::risco(double r) {risco_=r;}

PatternDiskBB::PatternDiskBB() :
  PatternDisk(),
  spectrumBB_(NULL),
  SpectralEmission_(0), risco_(0.)
{
  GYOTO_DEBUG << "PatternDiskBB Construction" << endl;
  spectrumBB_ = new Spectrum::BlackBody(); 
}

PatternDiskBB::PatternDiskBB(const PatternDiskBB& o) :
  PatternDisk(o),
  spectrumBB_(NULL),
  SpectralEmission_(o.SpectralEmission_), risco_(o.risco_)
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

double const * PatternDiskBB::getVelocity() const { return PatternDisk::getVelocity(); }

void PatternDiskBB::getVelocity(double const pos[4], double vel[4]) {
  // The only use of this reimplementation: ensure nothing happens below ISCO

  double const * const rad=getGridRadius();
  size_t i[3]; // {i_nu, i_phi, i_r}
  getIndices(i, pos, 0.); //NB: last arg should be nu, don't care here
  double rgrid=rad[i[2]-1]; // this is the smallest radius used
                          // when dealing with the current r value

  if (rgrid<risco()){
    //default velocity, emission will be 0 there anyway
    vel[0]=1.;
    for (int ii=1;ii<4;ii++)
      vel[ii]=0.;
  }else{
    PatternDisk::getVelocity(pos, vel);
  }
}

double PatternDiskBB::emission(double nu, double dsem,
			       double *,
			       double co[8]) const{
  // This reimplementation of emission has 2 goals, ensuring that
  // nothing happens below ISCO, and allowing to compute BB emission
  // when the PatternDisk structure contains temperature

  GYOTO_DEBUG << endl;
  
  size_t i[3]; // {i_nu, i_phi, i_r}
  getIndices(i, co, nu);
  double const * const rad=getGridRadius();
  double rgrid=rad[i[2]-1];

  // no emission in any case above rmax_:
  if (rgrid > rmax_ || rgrid < risco()) return 0.; 

  double Iem=0.;
  size_t naxes[3];
  getIntensityNaxes(naxes);
  size_t nnu=naxes[0], nphi=naxes[1];
  if (!SpectralEmission_){
    /*
      Here the intensity is assumed to be given by the emission
      function of PatternDisk
     */
    Iem = PatternDisk::emission(nu,dsem,co,co);//NB: 3rd argument is not used
  }else{ //Spectral emission    
    /*
      Here the temperature is assumed to be given by the emission
      function of PatternDisk, intensity is computed assuming BB
     */
    double TT;
    TT = PatternDisk::emission(nu,dsem,co,co);
    spectrumBB_->temperature(TT);
    Iem=(*spectrumBB_)(nu);
  }

  if (!flag_radtransf_) return Iem;

  double thickness;
  double const * const op = opacity();
  if (op && (thickness=op[i[2]*(nphi*nnu)+i[1]*nnu+i[0]]*dsem))
    return Iem * (1. - exp (-thickness)) ;
  return 0.;
}

void PatternDiskBB::metric(SmartPointer<Metric::Generic> gg) {
  //Metric must be KerrBL or alike
  string kin = gg->kind();
  if ((kin != "KerrBL") && (kin != "ChernSimons"))
    throwError
      ("PatternDiskBB::metric(): metric must be KerrBL or CS");
  ThinDisk::metric(gg);
}
