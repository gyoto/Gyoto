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
#define throwCfitsioError(status) \
    { fits_get_errstatus(status, ermsg); throwError(ermsg); }

#include "GyotoPhoton.h"
#include "GyotoPatternDiskBB.h"
#include "GyotoUtils.h"
#include "GyotoFactoryMessenger.h"
#include "GyotoKerrBL.h"
#include "GyotoKerrKS.h"


#include <fitsio.h>
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

PatternDiskBB::~PatternDiskBB() {
  GYOTO_DEBUG << "PatternDiskBB Destruction" << endl;
}

double const * PatternDiskBB::getVelocity() const { return PatternDisk::getVelocity(); }

void PatternDiskBB::getVelocity(double const pos[4], double vel[4]) {
  // The only use of this reimplementation: ensure nothing happens below ISCO

  double risco;
  if (risco_>0.) risco=risco_;
  else {
    switch (gg_->coordKind()) {
    case GYOTO_COORDKIND_SPHERICAL:
      risco = static_cast<SmartPointer<Metric::KerrBL> >(gg_) -> getRms();
      break;
    default:
      throwError("PatternDiskBB::getVelocity: bad COORDKIND");
      risco=0.;
    }
  }
  
  double const * const radius=getGridRadius();
  size_t i[3]; // {i_nu, i_phi, i_r}
  getIndices(i, pos, 0.); //NB: last arg should be nu, don't care here
  double rgrid=radius[i[2]-1]; // this is the smallest radius used
                          // when dealing with the current r value

  if (rgrid<risco){
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
  
  double risco;
  if (risco_>0.) risco=risco_;
  else {
    switch (gg_->coordKind()) {
    case GYOTO_COORDKIND_SPHERICAL:
      risco = static_cast<SmartPointer<Metric::KerrBL> >(gg_) -> getRms();
      break;
    default:
      throwError("PatternDiskBB::emission: bad COORDKIND");
      risco=0.;
    }
  }

  size_t i[3]; // {i_nu, i_phi, i_r}
  getIndices(i, co, nu);
  double const * const radius=getGridRadius();
  double rgrid=radius[i[2]-1];

  // no emission in any case above rmax_:
  if (rgrid > rmax_ || rgrid < risco) return 0.; 

  double Iem=0.;
  double const * const emiss = getIntensity();
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
    spectrumBB_->setTemperature(TT);
    Iem=(*spectrumBB_)(nu);
  }

  if (!flag_radtransf_) return Iem;

  double thickness;
  double const * const opacity = getOpacity();
  if (opacity && (thickness=opacity[i[2]*(nphi*nnu)+i[1]*nnu+i[0]]*dsem))
    return Iem * (1. - exp (-thickness)) ;
  return 0.;
}

void PatternDiskBB::setMetric(SmartPointer<Metric::Generic> gg) {
  //Metric must be KerrBL or alike
  string kind = gg->kind();
  if ((kind != "KerrBL") && (kind != "ChernSimons"))
    throwError
      ("PatternDiskBB::setMetric(): metric must be KerrBL or CS");
  ThinDisk::setMetric(gg);
}

int PatternDiskBB::setParameter(std::string name,
				std::string content,
				std::string unit) {
  if (name=="SpectralEmission") SpectralEmission_=1;
  else if (name=="Risco") risco_=atof(content.c_str());
  else return PatternDisk::setParameter(name, content, unit);
  return 0;
}

#ifdef GYOTO_USE_XERCES
void PatternDiskBB::fillElement(FactoryMessenger *fmp) const {
  fmp -> setParameter ( SpectralEmission_? "SpectralEmission" : "BolometricEmission");
  PatternDisk::fillElement(fmp);
}

#endif
