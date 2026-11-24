/*
    Copyright 2011-2016 Thibaut Paumard

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

/// Properties

// There is no generic properties for spectrometers. Nevertheless, we
// define this to derived classes can point to
// Spectrometer::Generic::properties rather than Object::properties

#include "GyotoProperty.h"
GYOTO_PROPERTY_START(Spectrometer::Generic,
		     "Spectrometric capabilities of a Gyoto::Screen.")
GYOTO_PROPERTY_END(Spectrometer::Generic, Object::properties)

///

Register::Entry* Gyoto::Spectrometer::Register_ = NULL;
void Spectrometer::initRegister() {
  if (Gyoto::Spectrometer::Register_) delete Gyoto::Spectrometer::Register_;
  Gyoto::Spectrometer::Register_ = NULL;
  // statically fill the register
  Gyoto::Spectrometer::Register("wave",    &(Subcontractor<Uniform>));
  Gyoto::Spectrometer::Register("wavelog", &(Subcontractor<Uniform>));
  Gyoto::Spectrometer::Register("freq",    &(Subcontractor<Uniform>));
  Gyoto::Spectrometer::Register("freqlog", &(Subcontractor<Uniform>));
  Gyoto::Spectrometer::Register("Complex", &(Subcontractor<Complex>));
}

void Gyoto::Spectrometer::Register(std::string name, Subcontractor_t* scp){
  Register::Entry* ne =
    new Register::Entry(name, (SmartPointee::Subcontractor_t*)scp, Gyoto::Spectrometer::Register_);
  Gyoto::Spectrometer::Register_ = ne;
}

GYOTO_GETSUBCONTRACTOR(Spectrometer)
GYOTO_REGISTEREDPLUGINSSLASHKINDS(Spectrometer)

Generic::Generic() :
  SmartPointee(),
  Object(),
  Teller(),
  kindid_(NULL),
  nsamples_(0),
  nboundaries_(0),
  boundaries_(NULL),
  chanind_(NULL),
  midpoints_(NULL),
  widths_(NULL)
{}

Generic::Generic(kind_t kin) :
  SmartPointee(),
  Object(kin),
  Teller(),
  kindid_(kin),
  nsamples_(0),
  nboundaries_(0),
  boundaries_(NULL),
  chanind_(NULL),
  midpoints_(NULL),
  widths_(NULL)
{}

Generic::Generic(const Generic& o) :
  SmartPointee(o),
  Object(o),
  Teller(o),
  kindid_(o.kindid_),
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

char const * Generic::kindid() const {return kindid_;}
void Generic::kindid(char const * k) {kindid_=k; kind_=k; tellListeners();}

size_t Generic::nSamples() const { return nsamples_; }
size_t Generic::getNBoundaries() const { return nboundaries_; }

double const * Generic::getMidpoints() const { return midpoints_; }
void Generic::getMidpoints( double data[], std::string unit) const {
  for (size_t i=0; i<nsamples_; ++i)
    data[i]=Units::FromHerz(midpoints_[i], unit);
}
std::vector<double> Generic::midpoints() const {
  std::vector<double> out(midpoints_, midpoints_+nsamples_);
  return out;
}
std::vector<double> Generic::midpoints(std::string unit) const {
  std::vector<double> out(nsamples_);
  getMidpoints(out.data(), unit);
  return out;
}
double const * Generic::getChannelBoundaries() const { return boundaries_;}
void Generic::getChannelBoundaries( double data[], std::string unit) const {
  for (size_t i=0; i<nboundaries_; ++i)
    data[i]=Units::FromHerz(boundaries_[i], unit);
}
std::vector<double> Generic::channelBoundaries() const {
  std::vector<double> out(boundaries_, boundaries_+nboundaries_);
  return out;
}
std::vector<double> Generic::channelBoundaries(std::string unit) const {
  std::vector<double> out(nboundaries_);
  getChannelBoundaries(out.data(), unit);
  return out;
}
size_t const * Generic::getChannelIndices() const { return chanind_; }
std::vector<unsigned long> Generic::channelIndices() const {
  std::vector<unsigned long> out;
  out.reserve(2*nsamples_);
  for (size_t i = 0; i < 2*nsamples_; ++i)
    out.push_back(chanind_[i]);
  return out;
}
double const * Generic::getWidths() const { return widths_; }
void Generic::getWidths( double data[], std::string unit) const {
  double * cbound = new double[nboundaries_];
  getChannelBoundaries(cbound, unit);
  for(size_t i=0; i<nsamples_; ++i)
    data[i]=fabs(cbound[chanind_[2*i+1]]-cbound[chanind_[2*i]]);
  delete [] cbound;
}
std::vector<double> Generic::widths() const {
  std::vector<double> out(widths_, widths_+nsamples_);
  return out;
}
std::vector<double> Generic::widths(std::string unit) const {
  std::vector<double> out(nsamples_);
  getWidths(out.data(), unit);
  return out;
}
