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

#include "GyotoSpectrometer.h"
#include "GyotoUniformSpectrometer.h"
#include "GyotoComplexSpectrometer.h"
#include "GyotoUtils.h"
#include "GyotoFactoryMessenger.h"
#include "GyotoConverters.h"
#include "GyotoMetric.h"

#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <float.h> // DBL_MAX
using namespace Gyoto;
using namespace Gyoto::Spectrometer;
using namespace std;

#if defined GYOTO_USE_XERCES
Register::Entry* Gyoto::Spectrometer::Register_ = NULL;
void Spectrometer::initRegister() {
  if (Gyoto::Spectrometer::Register_) delete Gyoto::Spectrometer::Register_;
  Gyoto::Spectrometer::Register_ = NULL;
  // statically fill the register
  Register("wave",    &(Subcontractor<Uniform>));
  Register("wavelog", &(Subcontractor<Uniform>));
  Register("freq",    &(Subcontractor<Uniform>));
  Register("freqlog", &(Subcontractor<Uniform>));
  Register("Complex", &(Subcontractor<Complex>));
}

void Gyoto::Spectrometer::Register(std::string name, Subcontractor_t* scp){
  Register::Entry* ne =
    new Register::Entry(name, (SmartPointee::Subcontractor_t*)scp, Gyoto::Spectrometer::Register_);
  Gyoto::Spectrometer::Register_ = ne;
}

Gyoto::Spectrometer::Subcontractor_t*
Spectrometer::getSubcontractor(std::string name, int errmode) {
  if (!Gyoto::Spectrometer::Register_) throwError("No Spectrometer kind registered!");
  return (Subcontractor_t*)Gyoto::Spectrometer::Register_
    -> getSubcontractor(name, errmode);
}

void Spectrometer::Generic::setParameters(Gyoto::FactoryMessenger *fmp)  {
  string name="", content="", unit="";
  if (fmp)
    while (fmp->getNextParameter(&name, &content, &unit))
      setParameter(name, content, unit);
}
#endif



Generic::Generic() :
  SmartPointee(),
  Teller(),
  kind_(NULL),
  nsamples_(0),
  nboundaries_(0),
  boundaries_(NULL),
  chanind_(NULL),
  midpoints_(NULL),
  widths_(NULL)
{}
Generic::Generic(kind_t kind) :
  SmartPointee(),
  Teller(),
  kind_(kind),
  nsamples_(0),
  nboundaries_(0),
  boundaries_(NULL),
  chanind_(NULL),
  midpoints_(NULL),
  widths_(NULL)
{}
Generic::Generic(const Generic& o) :
  SmartPointee(o),
  Teller(o),
  kind_(o.kind_),
  nsamples_(o.nsamples_),
  nboundaries_(o.nboundaries_),
  boundaries_(NULL),
  chanind_(NULL),
  midpoints_(NULL),
  widths_(NULL)
{
  if (o.boundaries_) boundaries_=new double[nboundaries_];
  memcpy(boundaries_, o.boundaries_, nboundaries_*sizeof(double));
  if (o.widths_) widths_=new double[nsamples_];
  memcpy(widths_, o.widths_, nsamples_*sizeof(double));
  if (o.midpoints_) midpoints_=new double[nsamples_];
  memcpy(midpoints_, o.midpoints_, nsamples_*sizeof(double));
  if (o.chanind_) chanind_=new size_t[2*nsamples_];
  memcpy(chanind_, o.chanind_, 2*nsamples_*sizeof(size_t));
}
Generic::~Generic() {
  if (boundaries_) delete [] boundaries_;
  if (widths_) delete [] widths_;
  if (midpoints_) delete [] midpoints_;
  if (chanind_) delete [] chanind_;
}

char const * Generic::getKind() const {return kind_;}
void Generic::setKind(char const * k) {kind_=k; tellListeners();}

size_t Generic::getNSamples() const { return nsamples_; }
size_t Generic::getNBoundaries() const { return nboundaries_; }

double const * Generic::getMidpoints() const { return midpoints_; }
void Generic::getMidpoints( double data[], std::string unit) {
  for (size_t i=0; i<nsamples_; ++i)
    data[i]=Units::FromHerz(midpoints_[i], unit);
}
double const * Generic::getChannelBoundaries() const { return boundaries_;}
void Generic::getChannelBoundaries( double data[], std::string unit) {
  for (size_t i=0; i<nboundaries_; ++i)
    data[i]=Units::FromHerz(boundaries_[i], unit);
}
size_t const * Generic::getChannelIndices() const { return chanind_; }
double const * Generic::getWidths() const { return widths_; }
void Generic::getWidths( double data[], std::string unit) {
  double cbound[nboundaries_];
  getChannelBoundaries(cbound, unit);
  for(size_t i=0; i<nsamples_; ++i)
    data[i]=fabs(cbound[chanind_[2*i+1]]-cbound[chanind_[2*i]]);
}

void Spectrometer::Generic::setParameter(string name,
					 string content,
					 string unit)
{
  // nothing to do... yet
}

#ifdef GYOTO_USE_XERCES

void Generic::fillElement(FactoryMessenger *fmp) const {
  fmp -> setSelfAttribute( "kind", getKind() );
}
#endif
