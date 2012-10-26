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
#include "GyotoUtils.h"
#include "GyotoFactoryMessenger.h"
#include "GyotoConverters.h"
#include "GyotoMetric.h"

#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <float.h> // DBL_MAX
using namespace Gyoto;
using namespace std;

Spectrometer::Spectrometer() :
  kind_(GYOTO_SPECTRO_KIND_NONE),
  nsamples_(0),
  boundaries_(NULL),
  midpoints_(NULL),
  widths_(NULL)
{
  band_[0]=0.; band_[1]=0.;
}
Spectrometer::Spectrometer(size_t nsamples, double band_min, double band_max,
			   SpectroKind_t kind) :
  kind_(kind),
  nsamples_(nsamples),
  boundaries_(NULL),
  midpoints_(NULL),
  widths_(NULL)
{
  band_[0]=band_min; band_[1]=band_max;
  if (nsamples && kind) reset_();
}

Spectrometer::Spectrometer(const Spectrometer& o) :
  SmartPointee(o),
  kind_(o.kind_),
  nsamples_(o.nsamples_),
  boundaries_(NULL),
  midpoints_(NULL),
  widths_(NULL)
{
  band_[0]=o.band_[0]; band_[1]=o.band_[1];
  reset_();
}

Spectrometer* Spectrometer::clone() const { return new Spectrometer(*this); }

Spectrometer::~Spectrometer() {
  if (boundaries_) delete [] boundaries_;
  if (midpoints_)  delete [] midpoints_;
  if (widths_)     delete [] widths_;
}

void Spectrometer::reset_() {
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << endl;
# endif
  if (boundaries_) delete [] boundaries_;
  if (midpoints_)  delete [] midpoints_;
  if (widths_)     delete [] widths_;
  boundaries_ = NULL;
  midpoints_ = NULL;
  widths_ = NULL;
  GYOTO_DEBUG << endl;
  if (!nsamples_ || !kind_) return;

  boundaries_ = new double[nsamples_+1];
  midpoints_  = new double[nsamples_];
  widths_     = new double[nsamples_];

  size_t i=0 ;
# if GYOTO_DEBUG_ENABLED
  GYOTO_IF_DEBUG;
  GYOTO_DEBUG_EXPR(band_);
  if (band_) GYOTO_DEBUG_ARRAY(band_, 2);
  GYOTO_DEBUG_EXPR(kind_);
  GYOTO_ENDIF_DEBUG
# endif
  for (i=0; i<=nsamples_; ++i){
#   if GYOTO_DEBUG_ENABLED
    GYOTO_DEBUG << "boundaries_[" <<i<<"]=";
#   endif
    boundaries_[i]=band_[0]+double(i)*(band_[1]-band_[0])/double(nsamples_);
    if (kind_==GYOTO_SPECTRO_KIND_FREQLOG ||
	kind_==GYOTO_SPECTRO_KIND_WAVELOG)
      boundaries_[i]=pow(10.,boundaries_[i]);
    if (kind_==GYOTO_SPECTRO_KIND_WAVE ||
	kind_==GYOTO_SPECTRO_KIND_WAVELOG)
      boundaries_[i]=boundaries_[i]?GYOTO_C/boundaries_[i]:DBL_MAX;
#   if GYOTO_DEBUG_ENABLED
    if (debug()) cerr << boundaries_[i]<< endl;
#   endif
  }
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG_ARRAY(boundaries_, nsamples_+1);
# endif

  for (i=0; i< nsamples_; ++i){
#   if GYOTO_DEBUG_ENABLED
    GYOTO_DEBUG << "width_[" <<i<<"]=";
#   endif
    widths_[i] = fabs(boundaries_[i+1] - boundaries_[i]);
#   if GYOTO_DEBUG_ENABLED
    GYOTO_IF_DEBUG;
    cerr << widths_[i] << endl;
    GYOTO_DEBUG << "midpoints_[" <<i<<"]=";
    GYOTO_ENDIF_DEBUG
#   endif
      midpoints_[i] = (boundaries_[i+1]*0.5 + boundaries_[i]*0.5);
                             // avoid overflow
#   if GYOTO_DEBUG_ENABLED
    if (debug()) cerr << midpoints_[i] << endl;
#   endif
  }

}

void Spectrometer::setKind(SpectroKind_t k) {
  kind_ = k;
  reset_();
}

void Spectrometer::setKind(std::string str) {
  SpectroKind_t s;

  if      (str == "none"   ) s = GYOTO_SPECTRO_KIND_NONE;
  else if (str == "freq"   ) s = GYOTO_SPECTRO_KIND_FREQ;
  else if (str == "freqlog") s = GYOTO_SPECTRO_KIND_FREQLOG;
  else if (str == "wave"   ) s = GYOTO_SPECTRO_KIND_WAVE;
  else if (str == "wavelog") s = GYOTO_SPECTRO_KIND_WAVELOG;
  else {
    throwError("unknow Spectrometer kind"); s=GYOTO_SPECTRO_KIND_NONE;
  }

  kind_ = s;
  reset_();
}

void Spectrometer::setNSamples(size_t n) { nsamples_ = n; reset_(); }
void Spectrometer::setBand(double nu[2]) {
  band_[0] = nu[0];
  band_[1] = nu[1];
  reset_();
}

void Spectrometer::setBand(double nu[2], string unit, string kind) {
  if (kind != "") setKind(kind);
  double band[2] = {nu[0], nu[1]};

  switch (kind_) {
  case GYOTO_SPECTRO_KIND_FREQ: 
    if (unit != "" && unit != "Hz")
      for (size_t i=0; i<=1; ++i)
	band[i] = Units::ToHerz(nu[i], unit);
    break;
  case GYOTO_SPECTRO_KIND_FREQLOG:
    if (unit != "" && unit != "Hz")
      for (size_t i=0; i<=1; ++i)
	band[i] = log10(Units::ToHerz(pow(10., nu[i]), unit));
    break;
  case GYOTO_SPECTRO_KIND_WAVE:
    if (unit != "" && unit != "m")
      for (size_t i=0; i<=1; ++i)
	band[i] = Units::ToMeters(nu[i], unit);
    break;
  case GYOTO_SPECTRO_KIND_WAVELOG:
    if (unit != "" && unit != "m")
      for (size_t i=0; i<=1; ++i)
	band[i] = log10(Units::ToMeters(pow(10., nu[i]), unit));
    break;
  default:
    throwError("Spectrometer::setBand(double, string, string) at loss: "
	       "please specify Spectrometer kind");
  }

  setBand(band);
}

SpectroKind_t Spectrometer::getKind() const {
  return kind_;
}

std::string Spectrometer::getKindStr() const {
  std::string skind = "";
  stringstream ss;
  switch (kind_) {
  case GYOTO_SPECTRO_KIND_NONE:
    skind = ("none"); break;
  case GYOTO_SPECTRO_KIND_FREQ:
    skind = ("freq"); break;
  case GYOTO_SPECTRO_KIND_FREQLOG:
    skind = ("freqlog"); break;
  case GYOTO_SPECTRO_KIND_WAVE:
    skind = ("wave"); break;
  case GYOTO_SPECTRO_KIND_WAVELOG:
    skind = ("wavelog"); break;
  default:
    ss << "Unknown spectrometer kind: " << kind_;
    throwError( ss.str() );
  }
  return skind;
}

size_t Spectrometer::getNSamples() const { return nsamples_; }

double const * Spectrometer::getBand() const { return band_; }

double const * Spectrometer::getMidpoints() const { return midpoints_; }
double const * Spectrometer::getChannels() const { return boundaries_; }
double const * Spectrometer::getWidths() const { return widths_; }

#ifdef GYOTO_USE_XERCES

void Spectrometer::fillElement(FactoryMessenger *fmp) {
  fmp -> setSelfAttribute( "kind", getKindStr() );
  fmp -> setSelfAttribute( "nsamples", nsamples_ );
  ostringstream ss;
  ss << setprecision(GYOTO_PREC) << setw(GYOTO_WIDTH) << band_[0] << " "
     << setprecision(GYOTO_PREC) << setw(GYOTO_WIDTH) << band_[1];
  fmp -> setFullContent(ss.str()); 
}

SmartPointer<Spectrometer>
Gyoto::SpectrometerSubcontractor(FactoryMessenger* fmp) {

  string skind = fmp -> getSelfAttribute( "kind" );
  size_t nsamples = atol( fmp -> getSelfAttribute( "nsamples" ) . c_str () );
  string unit = fmp -> getSelfAttribute( "unit" );
  string content = fmp -> getFullContent();
  char * tc = const_cast<char*>(content.c_str());
  double band[2];
  band[0]=strtod(tc, &tc);
  band[1]=strtod(tc, &tc);

  Spectrometer* spr =new Spectrometer();
  spr -> setBand(band, unit, skind);
  spr -> setNSamples(nsamples);

  return spr;

}


#endif
