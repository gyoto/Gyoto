/*
    Copyright 2013 Thibaut Paumard

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

#include "GyotoUniformSpectrometer.h"
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

#include "GyotoProperty.h"
GYOTO_PROPERTY_START(Uniform,
		     "Spectrometer with uniform spacing (log or linear, wave. of freq.).")
GYOTO_PROPERTY_VECTOR_DOUBLE_UNIT(Uniform, Band, band,
				  "Spectral band (default unit: Hz, m, Log(Hz) or Log(m)).")
GYOTO_PROPERTY_STRING(Uniform, Kind, kind,
		      "\"freq\", \"wave\", \"freqlog\" or \"wavelog\".")
GYOTO_PROPERTY_SIZE_T(Uniform, NSamples, nSamples,
		      "Number of spectral channels.")
GYOTO_PROPERTY_END(Uniform, Generic::properties)

void Uniform::fillProperty(Gyoto::FactoryMessenger *fmp,
			   Property const &p) const {
  if (p.type == Property::size_t_t
      && p.name == "NSamples")
    fmp -> setSelfAttribute("nsamples", nsamples_);
  else if (p.type == Property::string_t        && p.name == "Kind")
    ; // do nothing
  else if (p.type == Property::vector_double_t && p.name == "Band") {
    ostringstream ss;
    ss << setprecision(GYOTO_PREC) << setw(GYOTO_WIDTH) << band_[0] << " "
       << setprecision(GYOTO_PREC) << setw(GYOTO_WIDTH) << band_[1];
    fmp -> setFullContent(ss.str());
  } else GYOTO_ERROR("Unsupported Property in Spectrometer::Uniform"); 
}

#ifdef GYOTO_USE_XERCES

void Gyoto::Spectrometer::Uniform::setParameters(FactoryMessenger* fmp) {
  string skind = fmp -> getSelfAttribute( "kind" );
  size_t nsamples = atol( fmp -> getSelfAttribute( "nsamples" ) . c_str () );
  string unit = fmp -> getSelfAttribute( "unit" );

  string content = fmp -> getFullContent();
  double band[2];
  if (FactoryMessenger::parseArray(content, band, 2) != 2)
    GYOTO_ERROR("Spectromecter::Uniform requires exactly 2 tokens");

  Uniform::band(band, unit, skind);
  nSamples(nsamples);
}
#endif

///

Uniform::Uniform() :
  Generic(WaveKind)
{
  band_[0]=0.; band_[1]=0.;
}
Uniform::Uniform(size_t nsamples, double band_min, double band_max, kind_t kin)
: Generic(kin)
{
  nsamples_=nsamples;
  band_[0]=band_min; band_[1]=band_max;
  if (nsamples && kin) reset_();
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
  if (!nsamples_ || !kindid_) return;

  boundaries_ = new double[nsamples_+1];
  chanind_    = new size_t[nsamples_*2];
  midpoints_  = new double[nsamples_];
  widths_     = new double[nsamples_];

  size_t i=0 ;
# if GYOTO_DEBUG_ENABLED
  GYOTO_IF_DEBUG;
  GYOTO_DEBUG_EXPR(band_);
  GYOTO_DEBUG_ARRAY(band_, 2);
  GYOTO_DEBUG_EXPR(kindid_);
  GYOTO_ENDIF_DEBUG
# endif
  for (i=0; i<=nsamples_; ++i){
#   if GYOTO_DEBUG_ENABLED
    GYOTO_DEBUG << "boundaries_[" <<i<<"]=";
#   endif
    boundaries_[i]=band_[0]+double(i)*(band_[1]-band_[0])/double(nsamples_);
    if (kindid_==FreqLogKind ||
	kindid_==WaveLogKind)
      boundaries_[i]=pow(10.,boundaries_[i]);
    if (kindid_==WaveKind ||
	kindid_==WaveLogKind)
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

void Uniform::kindid(kind_t k) {
  kind_=k;
  kindid_ = k;
  reset_();
}

void Uniform::kind(std::string const &str) {
  kind_ = str;
  kind_t s;

  if (str == "freq"   ) s = FreqKind;
  else if (str == "freqlog") s = FreqLogKind;
  else if (str == "wave"   ) s = WaveKind;
  else if (str == "wavelog") s = WaveLogKind;
  else {
    GYOTO_ERROR("unknown Spectrometer::Uniform kind"); s=NULL;
  }

  kindid_ = s;
  reset_();
}

std::string Uniform::kind() const {return Object::kind_;}

void Uniform::nSamples(size_t n) {
  nsamples_ = n;
  nboundaries_=nsamples_+1;
  reset_();
}

void Uniform::band(std::vector<double>const &vnu) {
  if (vnu.size() != 2) GYOTO_ERROR("Band needs exactly two elements");
  double nu[] = {vnu[0], vnu[1]};
  band(nu);
}

void Uniform::band(std::vector<double>const &vnu, std::string const &u) {
  if (vnu.size() != 2) GYOTO_ERROR("Band needs exactly two elements");
  double nu[] = {vnu[0], vnu[1]};
  band(nu, u);
}

std::vector<double> Uniform::band() const {
  std::vector<double> vnu(2, band_[0]);
  vnu[1]=band_[1];
  return vnu;
}

std::vector<double> Uniform::band(std::string const &unit) const {
  std::vector<double> vnu = band();

  if (kindid_== FreqKind) {
    if (unit != "" && unit != "Hz")
      for (size_t i=0; i<=1; ++i)
	vnu[i] = Units::FromHerz(vnu[i], unit);
  } else if (kindid_== FreqLogKind) {
    if (unit != "" && unit != "Hz")
      for (size_t i=0; i<=1; ++i)
	vnu[i] = pow(10., Units::FromHerz(log10(vnu[i]), unit));
  } else if (kindid_== WaveKind) {
    if (unit != "" && unit != "m")
      for (size_t i=0; i<=1; ++i)
	vnu[i] = Units::FromMeters(vnu[i], unit);
  } else if (kindid_ == WaveLogKind) {
    if (unit != "" && unit != "m")
      for (size_t i=0; i<=1; ++i)
	vnu[i] = pow(10., Units::FromMeters(log10(vnu[i]), unit));
  } else {
    GYOTO_ERROR("Uniform::band(string const &unit) at loss: "
	       "please specify Spectrometer kind");
  }

  return vnu;
}

void Uniform::band(double nu[2]) {
  band_[0] = nu[0];
  band_[1] = nu[1];
  reset_();
}

void Uniform::band(double nu[2], string const &unit, string const &skind) {
  if (skind != "") kind(skind);
  band(nu, unit);
}

void Uniform::band(double nu[2], string const &unit) {
  double band[2] = {nu[0], nu[1]};

  if (kindid_== FreqKind) {
    if (unit != "" && unit != "Hz")
      for (size_t i=0; i<=1; ++i)
	band[i] = Units::ToHerz(nu[i], unit);
  } else if (kindid_== FreqLogKind) {
    if (unit != "" && unit != "Hz")
      for (size_t i=0; i<=1; ++i)
	band[i] = log10(Units::ToHerz(pow(10., nu[i]), unit));
  } else if (kindid_== WaveKind) {
    if (unit != "" && unit != "m")
      for (size_t i=0; i<=1; ++i)
	band[i] = Units::ToMeters(nu[i], unit);
  } else if (kindid_ == WaveLogKind) {
    if (unit != "" && unit != "m")
      for (size_t i=0; i<=1; ++i)
	band[i] = log10(Units::ToMeters(pow(10., nu[i]), unit));
  } else {
    GYOTO_ERROR("Uniform::band(double, string) at loss: "
	       "please specify Spectrometer kind");
  }

  Uniform::band(band);
}


double const * Uniform::getBand() const { return band_; }

char const * const Uniform::WaveKind = "wave";
char const * const Uniform::WaveLogKind = "wavelog";
char const * const Uniform::FreqKind = "freq";
char const * const Uniform::FreqLogKind = "freqlog";
