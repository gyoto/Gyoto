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
  Register("wave", &(Uniform::Subcontractor));
  Register("wavelog", &(Uniform::Subcontractor));
  Register("freq", &(Uniform::Subcontractor));
  Register("freqlog", &(Uniform::Subcontractor));
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
#endif



Generic::Generic() :
  SmartPointee(),
  Teller(),
  kind_(GYOTO_SPECTRO_KIND_NONE),
  nsamples_(0),
  nboundaries_(0),
  boundaries_(NULL),
  chanind_(NULL),
  midpoints_(NULL),
  widths_(NULL)
{}
Generic::Generic(SpectroKind_t kind) :
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
void Generic::setKind(char const * k) {kind_=k;}

size_t Generic::getNSamples() const { return nsamples_; }
size_t Generic::getNBoundaries() const { return nboundaries_; }

double const * Generic::getMidpoints() const { return midpoints_; }
double const * Generic::getChannelBoundaries() const { return boundaries_;}
size_t const * Generic::getChannelIndices() const { return chanind_; }
double const * Generic::getWidths() const { return widths_; }

///////////// UNIFORM /////////////////


Uniform::Uniform() :
  Generic()
{
  band_[0]=0.; band_[1]=0.;
}
Uniform::Uniform(size_t nsamples, double band_min, double band_max,
			   SpectroKind_t kind) :
  Generic(kind)
{
  nsamples_=nsamples;
  band_[0]=band_min; band_[1]=band_max;
  if (nsamples && kind) reset_();
}

Uniform::Uniform(const Uniform& o) :
  Generic(o)
{
  band_[0]=o.band_[0]; band_[1]=o.band_[1];
  reset_();
}

Generic* Uniform::clone() const { return new Uniform(*this); }

Uniform::~Uniform() {
}

void Uniform::reset_() {
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << endl;
# endif
  if (boundaries_) delete [] boundaries_;
  if (chanind_)    delete [] chanind_;
  if (midpoints_)  delete [] midpoints_;
  if (widths_)     delete [] widths_;
  boundaries_ = NULL;
  chanind_ = NULL;
  midpoints_ = NULL;
  widths_ = NULL;
  GYOTO_DEBUG << endl;
  if (!nsamples_ || !kind_) return;

  boundaries_ = new double[nsamples_+1];
  chanind_    = new size_t[nsamples_*2];
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
    if (i<nsamples_) {
      chanind_[2*i]=i;
      chanind_[2*i+1]=i+1;
    }
  }
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG_ARRAY(boundaries_, nsamples_+1);
  GYOTO_DEBUG_ARRAY(chanind_, nsamples_*2);
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
  tellListeners();
}

void Uniform::setKind(SpectroKind_t k) {
  kind_ = k;
  reset_();
}

void Uniform::setKind(std::string str) {
  SpectroKind_t s;

  if      (str == "none"   ) s = GYOTO_SPECTRO_KIND_NONE;
  else if (str == "freq"   ) s = GYOTO_SPECTRO_KIND_FREQ;
  else if (str == "freqlog") s = GYOTO_SPECTRO_KIND_FREQLOG;
  else if (str == "wave"   ) s = GYOTO_SPECTRO_KIND_WAVE;
  else if (str == "wavelog") s = GYOTO_SPECTRO_KIND_WAVELOG;
  else {
    throwError("unknown Spectrometer::Uniform kind"); s=GYOTO_SPECTRO_KIND_NONE;
  }

  kind_ = s;
  reset_();
}

void Uniform::setNSamples(size_t n) {
  nsamples_ = n;
  nboundaries_=nsamples_+1;
  reset_();
}
void Uniform::setBand(double nu[2]) {
  band_[0] = nu[0];
  band_[1] = nu[1];
  reset_();
}

void Uniform::setBand(double nu[2], string unit, string kind) {
  if (kind != "") setKind(kind);
  double band[2] = {nu[0], nu[1]};

  if (kind_==GYOTO_SPECTRO_KIND_FREQ) {
    if (unit != "" && unit != "Hz")
      for (size_t i=0; i<=1; ++i)
	band[i] = Units::ToHerz(nu[i], unit);
  } else if (kind_== GYOTO_SPECTRO_KIND_FREQLOG) {
    if (unit != "" && unit != "Hz")
      for (size_t i=0; i<=1; ++i)
	band[i] = log10(Units::ToHerz(pow(10., nu[i]), unit));
  } else if (kind_== GYOTO_SPECTRO_KIND_WAVE) {
    if (unit != "" && unit != "m")
      for (size_t i=0; i<=1; ++i)
	band[i] = Units::ToMeters(nu[i], unit);
  } else if (kind_ == GYOTO_SPECTRO_KIND_WAVELOG) {
    if (unit != "" && unit != "m")
      for (size_t i=0; i<=1; ++i)
	band[i] = log10(Units::ToMeters(pow(10., nu[i]), unit));
  } else {
    throwError("Uniform::setBand(double, string, string) at loss: "
	       "please specify Spectrometer kind");
  }

  setBand(band);
}

double const * Uniform::getBand() const { return band_; }

std::string Gyoto::Spectrometer::Uniform::getKindStr() const { return kind_; }

#ifdef GYOTO_USE_XERCES

void Generic::fillElement(FactoryMessenger *fmp) const {
  fmp -> setSelfAttribute( "kind", getKind() );
}

void Uniform::fillElement(FactoryMessenger *fmp) const {
  fmp -> setSelfAttribute( "nsamples", nsamples_ );
  ostringstream ss;
  ss << setprecision(GYOTO_PREC) << setw(GYOTO_WIDTH) << band_[0] << " "
     << setprecision(GYOTO_PREC) << setw(GYOTO_WIDTH) << band_[1];
  fmp -> setFullContent(ss.str()); 
  Generic::fillElement(fmp);
}

SmartPointer<Generic>
Gyoto::Spectrometer::Uniform::Subcontractor(FactoryMessenger* fmp) {

  string skind = fmp -> getSelfAttribute( "kind" );
  size_t nsamples = atol( fmp -> getSelfAttribute( "nsamples" ) . c_str () );
  string unit = fmp -> getSelfAttribute( "unit" );
  string content = fmp -> getFullContent();
  char * tc = const_cast<char*>(content.c_str());
  double band[2];
  band[0]=strtod(tc, &tc);
  band[1]=strtod(tc, &tc);

  Uniform* spr =new Uniform();
  spr -> setBand(band, unit, skind);
  spr -> setNSamples(nsamples);

  return spr;

}


#endif

char const * const Uniform::WaveKind = "wave";
char const * const Uniform::WaveLogKind = "wavelog";
char const * const Uniform::FreqKind = "freq";
char const * const Uniform::FreqLogKind = "freqlog";
