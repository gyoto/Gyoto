/*
    Copyright 2011-2016, 2018-2020 Frederic Vincent, Thibaut Paumard

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

// GYOTO HEADERS
#include "GyotoUtils.h"
#include "GyotoAstrobj.h"
#include "GyotoMetric.h"
#include "GyotoPhoton.h"
#include "GyotoRegister.h"
#include "GyotoFactoryMessenger.h"
#include "GyotoProperty.h"

// SYSTEM HEADERS
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <cfloat>
#include <cmath>
#include <sstream>
#include <limits>

// NAMESPACES
using namespace std;
using namespace Gyoto;
using namespace Gyoto::Astrobj;
using namespace Eigen;

Register::Entry* Gyoto::Astrobj::Register_ = NULL;

GYOTO_PROPERTY_START(Gyoto::Astrobj::Generic,
   "Whatever emits or absorbs light.")
GYOTO_PROPERTY_METRIC(Generic, Metric, metric,
   "The geometry of space-time at this end of the Universe.")
GYOTO_PROPERTY_DOUBLE_UNIT(Generic, RMax, rMax,
   "Maximum distance from the centre of mass (geometrical units).")
GYOTO_PROPERTY_DOUBLE_UNIT(Generic, DeltaMaxInsideRMax, deltaMaxInsideRMax,
   "Maximum step for Photon integration inside RMax (geometrical units).")
GYOTO_PROPERTY_BOOL(Generic, Redshift, NoRedshift, redshift,
    "Whether to take redshift into account.")
GYOTO_PROPERTY_BOOL(Generic, ShowShadow, NoShowShadow, showshadow,
    "Whether to highlight the shadow region on the image.")
GYOTO_PROPERTY_BOOL(Generic, OpticallyThin, OpticallyThick, opticallyThin,
    "Whether the object should be considered optically thin or thick.")
GYOTO_PROPERTY_END(Generic, Object::properties)

Generic::Generic(string kin) :
  SmartPointee(), Object(kin),
  __defaultfeatures(0),
  gg_(NULL), rmax_(DBL_MAX), deltamaxinsidermax_(1.), flag_radtransf_(0),
  noredshift_(0), shadow_(0)
{
#if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << endl;
#endif
}

Generic::Generic() :
  SmartPointee(), Object("Default"),
  __defaultfeatures(0),
  gg_(NULL), rmax_(DBL_MAX), deltamaxinsidermax_(1.), flag_radtransf_(0),
  noredshift_(0), shadow_(0)
{
#if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << endl;
#endif
}

Generic::Generic(double radmax) :
  SmartPointee(), Object("Default"),
  __defaultfeatures(0),
  gg_(NULL), rmax_(radmax), deltamaxinsidermax_(1.), flag_radtransf_(0),
  noredshift_(0), shadow_(0)
{
#if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << endl;
#endif
}

Generic::Generic(const Generic& orig) :
  SmartPointee(orig), Object(orig),
  __defaultfeatures(orig.__defaultfeatures),
  gg_(NULL),
  rmax_(orig.rmax_),
  deltamaxinsidermax_(orig.deltamaxinsidermax_),
  flag_radtransf_(orig.flag_radtransf_),
  noredshift_(orig.noredshift_), shadow_(orig.shadow_)
{
#if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << endl;
#endif
    if (orig.gg_()) {
#if GYOTO_DEBUG_ENABLED
      GYOTO_DEBUG << "orig had a metric, cloning"<< endl;
#endif
      gg_=orig.gg_->clone();
    }
#if GYOTO_DEBUG_ENABLED
    GYOTO_DEBUG << "done" << endl;
#endif
}

Generic::~Generic() {
#if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << endl;
#endif
}

SmartPointer<Metric::Generic> Generic::metric() const { return gg_; }
void Generic::metric(SmartPointer<Metric::Generic> gg) {gg_=gg;}

const string Generic::kind() const { return kind_; }

double Generic::rMax() { return rmax_; }
double Generic::rMax() const { return rmax_; }
double Generic::rMax(string const &unit) {
  return Units::FromGeometrical(rMax(), unit, gg_); }
double Generic::rMax(string const &unit) const {
  return Units::FromGeometrical(rMax(), unit, gg_); }
void Generic::rMax(double val) { rmax_=val; }
void Generic::rMax(double val, string const &unit) {
  rMax(Units::ToGeometrical(val, unit, gg_)); }

double Generic::deltaMaxInsideRMax() const { return deltamaxinsidermax_; }
double Generic::deltaMaxInsideRMax(string const &unit) const {
  return Units::FromGeometrical(deltaMaxInsideRMax(), unit, gg_); }
void Generic::deltaMaxInsideRMax(double val) { deltamaxinsidermax_=val; }
void Generic::deltaMaxInsideRMax(double val, string const &unit) {
  deltaMaxInsideRMax(Units::ToGeometrical(val, unit, gg_)); }

#ifdef GYOTO_USE_XERCES
void Generic::setParameters(FactoryMessenger *fmp) {
  if (fmp) metric(fmp->metric());
  Object::setParameters(fmp);
}
#endif


void Generic::opticallyThin(bool flag) {flag_radtransf_=flag;}
bool Generic::opticallyThin() const {return flag_radtransf_;}

void Generic::showshadow(bool flag) {shadow_=flag;}
bool Generic::showshadow() const {return shadow_;}

void Generic::redshift(bool flag) {noredshift_=!flag;}
bool Generic::redshift() const {return !noredshift_;}

void Generic::processHitQuantities(Photon * ph, state_t const &coord_ph_hit,
				     double const * coord_obj_hit, double dt,
				     Properties* data) const {
#if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << endl;
#endif
  /*
      NB: freqObs is the observer's frequency chosen in
      Screen::getRayCoord for the actual computation of the geodesic ;
      the physical value of nuobs will be used in spectrum
      computations by resorting to the xml specifications of the user
      (see below) ; this freqObs is used to transform the null
      worldline parameter dlambda (see below)
  */
  double freqObs=ph->freqObs(); // this is a useless quantity, always 1
  SmartPointer<Spectrometer::Generic> spr = ph -> spectrometer();
  size_t nbnuobs = spr() ? spr -> nSamples() : 0 ;
  double const * const nuobs = nbnuobs ? spr -> getMidpoints() : NULL;
  double dlambda = dt/coord_ph_hit[4]; //dlambda = dt/tdot
  //  cout << "freqObs="<<freqObs<<endl;
  double ggredm1 = -gg_->ScalarProd(&coord_ph_hit[0],coord_obj_hit+4,
				    &coord_ph_hit[4]);// / 1.;
                                       //this is nu_em/nu_obs
                                       // for nuobs=1. Hz
  if (noredshift_) ggredm1=1.;
  double ggred = 1./ggredm1;           //this is nu_obs/nu_em
  double dsem = dlambda*ggredm1; // * 1Hz ?
  double inc =0.;

  if (data) {
#if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << "data requested. "
	      << "freqObs=" << freqObs << ", ggredm1=" << ggredm1
	      << ", ggred=" << ggred
	      << endl;
#endif

    if (data->redshift) {
      *data->redshift=ggred;
#if GYOTO_DEBUG_ENABLED
      GYOTO_DEBUG_EXPR(*data->redshift);
#endif
    }
    if (data->time) {
      *data->time=coord_ph_hit[0];
#if GYOTO_DEBUG_ENABLED
      GYOTO_DEBUG_EXPR(*data->time);
#endif
    }
    if (data->impactcoords && data->impactcoords[0]==DBL_MAX) {
      if (coord_ph_hit.size() > 8) GYOTO_ERROR("ImpactCoords is incompatible with parallel transport");
      memcpy(data->impactcoords, coord_obj_hit, 8 * sizeof(double));
      memcpy(data->impactcoords+8, &coord_ph_hit[0], 8 * sizeof(double));
    }
#if GYOTO_DEBUG_ENABLED
    GYOTO_DEBUG << "dlambda = (dt="<< dt << ")/(tdot="<< coord_ph_hit[4]
		<< ") = " << dlambda << ", dsem=" << dsem << endl;
#endif
    if (data->intensity) {
      /*
	Comment on intensity integration:
	Eq of radiative transfer used:
	dI_nu_em/ds_em = j_nu_em (assumes no absorption for
	simplicity) where : nu_em is measured by the emitter, ds_em
	is the proper length measured by the emitter corresponding
	to an increase dlambda of the parameter of the null
	worldline of the photon ; We have (see manual for demo) :
	ds_em = dlambda*nu_em

	BUT: dlambda is depending on the choice of freqObs (defined
	above) ; indeed we have: dlambda(freqObs) * freqObs = cst
	(i.e. independent of the choice of freqObs). Thus, for
	spectra computations where the observed frequency varies, we
	have to use the dlambda corresponding to the physical
	frequency chosen by the user, dlambda(nuobs), instead of
	dlambda(freqObs).  The conversion is easy : dlambda(nuobs) =
	dlambda(freqObs)*freqObs/nuobs Thus: ds_em =
	dlambda(nuobs)*nu_em = dlambda(freqObs)*freqObs*ggredm1

	Then, j_nu_em is computed by the emission() function of the
	Astrobj [NB: with rad. transfer emission() computes
	j_nu*dsem, without rad. transfer it computes I_nu, thus
	the result is always homogenous to intensity] Finally:
	I_nu_obs = I_nu_em*(nu_obs/nu_em)^3
      */

      //Intensity increment :
      GYOTO_DEBUG_EXPR(freqObs);
      GYOTO_DEBUG_EXPR(freqObs*ggredm1);
	inc = (emission(freqObs*ggredm1, dsem, coord_ph_hit, coord_obj_hit))
	  * (ph -> getTransmission(size_t(-1)))
	  * ggred*ggred*ggred; // I_nu/nu^3 invariant
#     ifdef HAVE_UDUNITS
      if (data -> intensity_converter_)
	inc = (*data -> intensity_converter_)(inc);
#     endif
      *data->intensity += inc;

#     if GYOTO_DEBUG_ENABLED
	GYOTO_DEBUG
	  << "intensity +=" << *data->intensity
	  << "= emission((dsem=" << dsem << "))="
	  << (emission(freqObs*ggredm1,dsem,coord_ph_hit, coord_obj_hit))
	  << ")*(ggred=" << ggred << ")^3*(transmission="
	  << (ph -> getTransmission(size_t(-1))) << ")"
	  << endl;
#     endif

    }
    if (data->binspectrum) {
      size_t nbounds = spr-> getNBoundaries();
      double const * const channels = spr -> getChannelBoundaries();
      size_t const * const chaninds = spr -> getChannelIndices();
      double * I  = new double[nbnuobs];
      double * boundaries = new double[nbounds];

      for (size_t ii=0; ii<nbounds; ++ii)
	boundaries[ii]=channels[ii]*ggredm1;
      integrateEmission(I, boundaries, chaninds, nbnuobs,
			dsem, coord_ph_hit, coord_obj_hit);
      for (size_t ii=0; ii<nbnuobs; ++ii) {
	inc = I[ii] * ph -> getTransmission(ii) * ggred*ggred*ggred*ggred;
#       ifdef HAVE_UDUNITS
	if (data -> binspectrum_converter_)
	  inc = (*data -> binspectrum_converter_)(inc);
#       endif
	data->binspectrum[ii*data->offset] += inc ;
#       if GYOTO_DEBUG_ENABLED
	GYOTO_DEBUG
	       << "nuobs[" << ii << "]="<< channels[ii]
	       << ", nuem=" << boundaries[ii]
	       << ", binspectrum[" << ii+data->offset << "]="
	       << data->binspectrum[ii*data->offset]<< endl;
#       endif

	if (!data->spectrum) // else it will be done in spectrum
	  ph -> transmit(ii,transmission(nuobs[ii]*ggredm1,dsem,coord_ph_hit, coord_obj_hit));
      }
      delete [] I;
      delete [] boundaries;
    }
    if (data->spectrum||data->stokesQ||data->stokesU||data->stokesV) {
      if (ph -> parallelTransport()) { // Compute polarization
	double * Inu          = new double[nbnuobs];
	double * Qnu          = new double[nbnuobs];
	double * Unu          = new double[nbnuobs];
	double * Vnu          = new double[nbnuobs];
	double * nuem         = new double[nbnuobs];
	Matrix4d * Onu        = new Matrix4d[nbnuobs];

	for (size_t ii=0; ii<nbnuobs; ++ii) {
	  nuem[ii]=nuobs[ii]*ggredm1;
	}
	GYOTO_DEBUG_ARRAY(nuobs, nbnuobs);
	GYOTO_DEBUG_ARRAY(nuem, nbnuobs);
	radiativeQ(Inu, Qnu, Unu, Vnu,
		   Onu, nuem, nbnuobs, dsem,
		   coord_ph_hit, coord_obj_hit);
	ph -> transfer(Inu, Qnu, Unu, Vnu, Onu);
	double ggred3 = ggred*ggred*ggred;
	for (size_t ii=0; ii<nbnuobs; ++ii) {
	  if (data-> spectrum) {
	    inc = Inu[ii] * ggred3;
#           ifdef HAVE_UDUNITS
	    if (data -> spectrum_converter_)
	      inc = (*data -> spectrum_converter_)(inc);
#           endif
	    data->spectrum[ii*data->offset] += inc;
	  }
	  if (data-> stokesQ) {
	    inc = Qnu[ii] * ggred3;
#           ifdef HAVE_UDUNITS
	    if (data -> spectrum_converter_)
	      inc = (*data -> spectrum_converter_)(inc);
#           endif
	    data->stokesQ [ii*data->offset] += inc;
	  }
	  if (data-> stokesU) {
	    inc = Unu[ii] * ggred3;
#           ifdef HAVE_UDUNITS
	    if (data -> spectrum_converter_)
	      inc = (*data -> spectrum_converter_)(inc);
#           endif
	    data->stokesU [ii*data->offset] += inc;
	  }
	  if (data-> stokesV) {
	    inc = Vnu[ii] * ggred3;
#           ifdef HAVE_UDUNITS
	    if (data -> spectrum_converter_)
	      inc = (*data -> spectrum_converter_)(inc);
#           endif
	    data->stokesV [ii*data->offset] += inc;
	  }

#         if GYOTO_DEBUG_ENABLED
	  {
	    double t=ph -> getTransmission(ii);
	    double o = t>0?-log(t):(-std::numeric_limits<double>::infinity());
	  GYOTO_DEBUG
	    << "DEBUG: Generic::processHitQuantities(): "
	    << "nuobs[" << ii << "]="<< nuobs[ii]
	    << ", nuem=" << nuem[ii]
	    << ", dsem=" << dsem
	    << ", Inu * GM/c2="
	    << Inu[ii]
	    << ", spectrum[" << ii*data->offset << "]="
	    << data->spectrum[ii*data->offset]
	    << ", transmission=" << t
	    << ", optical depth=" << o
	    << ", redshift=" << ggred << ")\n" << endl;
	  }
#         endif
	}
	delete [] Inu;
	delete [] Qnu;
	delete [] Unu;
	delete [] Vnu;
	delete [] Onu;
	delete [] nuem;
      } else { // No polarization
	double * Inu          = new double[nbnuobs];
	double * Taunu        = new double[nbnuobs];
	double * nuem         = new double[nbnuobs];

	for (size_t ii=0; ii<nbnuobs; ++ii) {
	  nuem[ii]=nuobs[ii]*ggredm1;
	}
	GYOTO_DEBUG_ARRAY(nuobs, nbnuobs);
	GYOTO_DEBUG_ARRAY(nuem, nbnuobs);
	radiativeQ(Inu, Taunu, nuem, nbnuobs, dsem,
		   coord_ph_hit, coord_obj_hit);
	for (size_t ii=0; ii<nbnuobs; ++ii) {
	  inc = Inu[ii] * ph -> getTransmission(ii) * ggred*ggred*ggred;
#         ifdef HAVE_UDUNITS
	  if (data -> spectrum_converter_)
	    inc = (*data -> spectrum_converter_)(inc);
#         endif
	  data->spectrum[ii*data->offset] += inc;
	  ph -> transmit(ii,Taunu[ii]);

#         if GYOTO_DEBUG_ENABLED
	  {
	    double t = ph -> getTransmission(ii);
	    double o = t>0?-log(t):(-std::numeric_limits<double>::infinity());
	  GYOTO_DEBUG
	    << "DEBUG: Generic::processHitQuantities(): "
	    << "nuobs[" << ii << "]="<< nuobs[ii]
	    << ", nuem=" << nuem[ii]
	    << ", dsem=" << dsem
	    << ", Inu * GM/c2="
	    << Inu[ii]
	    << ", spectrum[" << ii*data->offset << "]="
	    << data->spectrum[ii*data->offset]
	    << ", transmission=" << t
	    << ", optical depth=" << o
	    << ", redshift=" << ggred << ")\n" << endl;
	  }
#         endif
	}
	delete [] Inu;
	delete [] Taunu;
	delete [] nuem;
      }
    }
    /* update photon's transmission */
    ph -> transmit(size_t(-1),
		   transmission(freqObs*ggredm1, dsem,coord_ph_hit, coord_obj_hit));
  } else {
#   if GYOTO_DEBUG_ENABLED
    GYOTO_DEBUG << "NO data requested!" << endl;
#   endif
  }
}

/*
  User code is free to provide any or none of the various versions of
  emission(), transmission() and radiativeQ(). The default
  implementations call one another to try and find user-provided
  code. In order to avoid infinite recursion as well as for
  efficiency, several of those methose set a flag in __defaultfeatures
  if they are called to inform the other methods. This is what each
  method will try:
    - polarized radiativeQ:
      + unpolarized radiativeQ;
    - unpolarized radiativeQ:
      + polarized radiativeQ;
      + emission(double*, ...) and transmission;
    - emission(double*, ...):
      + unpolarized radiativeQ;
      + polarized radiativeQ;
      + emission(double, ...);
    - emission(double, ...):
      + emission(double*, ...);
      + unpolarized radiativeQ;
      + fall back to uniform, unit emission;
    - transmission:
      + unpolarized radiativeQ;
      + fall-back to uniform, unit opacity.
 */

#define __default_radiativeQ_polar 1
#define __default_radiativeQ       2
#define __default_emission_vector  4
double Generic::transmission(double nuem, double dsem, state_t const &coord_ph, double const coord_obj[8]) const {
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG_EXPR(flag_radtransf_);
  GYOTO_DEBUG_EXPR(__defaultfeatures);
# endif
  //cout << (__defaultfeatures & __default_radiativeQ) << endl;
  //cout << (__defaultfeatures & __default_radiativeQ_polar) << endl;
  //cout << __defaultfeatures << "," << __default_radiativeQ << "," << __default_radiativeQ_polar << endl;
  if ((!(__defaultfeatures & __default_radiativeQ)) ||
      (!(__defaultfeatures & __default_radiativeQ_polar))) {
    // We don't know (yet?) whether both radiativeQ are the default
    // implementations. Let's call the unpolarized version, which will
    // call the polarized version if the former is not reimplemented.
    double Inu, Taunu;
    radiativeQ(&Inu, &Taunu, &nuem, 1, dsem, coord_ph, coord_obj);
    return Taunu;
    // If radiativeQ is also the default implementation, it will set
    // __defaultfeatures and recurse back here
  }
  return double(flag_radtransf_);
}

double Generic::emission(double nuem, double dsem, state_t const &cph, double const *co) const
{
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG_EXPR(__defaultfeatures);
# endif

  if (!(__defaultfeatures & __default_emission_vector)) {
    // We don't know (yet?) whether emission(double* ,...) is the
    // default implementation, let's call it
    double Inu;
    emission(&Inu, &nuem , 1, dsem, cph, co);
    return Inu;
  } else if ((!(__defaultfeatures & __default_radiativeQ)) ||
	     (!(__defaultfeatures & __default_radiativeQ_polar))) {
    // emission(double*, ...) is the default implementation, but
    // we don't know (yet?) about radiativeQ(); let's call it
    double Inu, Taunu;
    radiativeQ(&Inu, &Taunu, &nuem, 1, dsem, cph, co);
    return Inu;
  }

  // emission(double*, ...) and radiativeQ are the default
  // implementations: we are on our own. Fall-back to uniform
  // emission.

# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG_EXPR(flag_radtransf_);
# endif
  if (flag_radtransf_) return dsem;
  return 1.;
}

void Generic::emission(double * Inu, double const * nuem , size_t nbnu,
			 double dsem, state_t const &cph, double const *co) const
{
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG_EXPR(__defaultfeatures);
# endif

  // Inform emission(double, ...)  that emission(double *, ...)  is
  // the default implementation and that it should not recurse back
  const_cast<Generic*>(this)->__defaultfeatures |= __default_emission_vector;

  if (!(__defaultfeatures & __default_radiativeQ)) {
    // We don't know (yet?) whether unpolarized radiativeQ is the
    // default implementation, let's call it
    double * Taunu = new double[nbnu];
    radiativeQ(Inu, Taunu, nuem, nbnu, dsem, cph, co);
    delete [] Taunu;
    // If radiativeQ is the default implementation, it will recurse
    // back here and we are going to skip to the next case below
    return;
  } else if (!(__defaultfeatures & __default_radiativeQ_polar)) {
    // We don't know (yet?) whether polarized radiativeQ is the
    // default implementation, let's call it
    double * Taunu = new double[nbnu];
    double * Qnu = new double[nbnu];
    double * Unu = new double[nbnu];
    double * Vnu = new double[nbnu];
    Matrix4d * Onu = new Matrix4d[nbnu];
    double Chi=0;
    radiativeQ(Inu, Qnu, Unu, Vnu,
	       Onu, nuem , nbnu, dsem,
	       cph, co);
    // in all cases, clean up
    delete [] Qnu;
    delete [] Unu;
    delete [] Vnu;
    delete [] Taunu;
    delete [] Onu;
    return;
  }

  // If both radiativeQ() versions are the default implementations,
  // fall back to emission(double, ...)
  for (size_t i=0; i< nbnu; ++i) Inu[i]=emission(nuem[i], dsem, cph, co);
}

void Generic::radiativeQ(double * Inu, double * Taunu,
			 double const * nuem , size_t nbnu,
			 double dsem, state_t const &cph, double const *co) const
{
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG_EXPR(flag_radtransf_);
# endif
  // Inform emission() and transmission() that radiativeQ is the
  // default implementation and that they should not recurse back here
  const_cast<Generic*>(this)->__defaultfeatures |= __default_radiativeQ;

  if (!(__defaultfeatures & __default_radiativeQ_polar)) {
    // We don't know (yet?) whether polarized radiativeQ is the
    // default implementation, let's call it
    double * Qnu = new double[nbnu];
    double * Unu = new double[nbnu];
    double * Vnu = new double[nbnu];
    double * alphaInu = new double[nbnu];
    Matrix4d * Onu = new Matrix4d[nbnu];
    radiativeQ(Inu, Qnu, Unu, Vnu,
	       Onu, nuem , nbnu, dsem,
	       cph, co);
    if (!(__defaultfeatures & __default_radiativeQ_polar)) {
      for (size_t i=0; i<nbnu; ++i) {
				Taunu[i] = Onu[i](0,0);
      }
    } else {
      for (size_t i=0; i<nbnu; ++i) {
				Taunu[i]=transmission(nuem[i], dsem, cph, co);
      }
    }
    // in all cases, clean up
    delete [] Qnu;
    delete [] Unu;
    delete [] Vnu;
    delete [] alphaInu;
    delete [] Onu;
    return;
    // If polarized radiativeQ is not implemented, the default
    // implementation will recurse here.
  }

  // Polarized radiativeQ for sure not implemented, fall-back to
  // emission() and transmission().

  // Call emission()
  emission(Inu, nuem, nbnu, dsem, cph, co);

  // Call transmission()
  for (size_t i=0; i< nbnu; ++i)
    Taunu[i]=transmission(nuem[i], dsem, cph, co);
}

void Generic::radiativeQ(double *Inu, double *Qnu, double *Unu, double *Vnu,
       Matrix4d *Onu,
			 double const *nuem , size_t nbnu, double dsem,
       state_t const &cph, double const *co) const
{
  // cph has 16 elements, 4 elements for each one of
  // X, Xdot, Ephi, Etheta
  //
  // If this method is not reimplemented, the emission is not polarized.

  // Inform non-polarized radiativeQ() that this method is the default
  // implementation and that they should not recurse back here
  const_cast<Generic*>(this)->__defaultfeatures |= __default_radiativeQ_polar;

  // Compute the output from the non-polarized radiativeQ().
  double * Taunu = new double[nbnu];
  double * alphaInu = new double[nbnu];
  Matrix4d identity;
  identity << 1, 0, 0, 0,
              0, 1, 0, 0,
              0, 0, 1, 0,
              0, 0, 0, 1;
  radiativeQ(Inu, Taunu, nuem, nbnu, dsem, cph, co);
  for (size_t i=0; i<nbnu; ++i) {
    // Inu[i] = Inu[i];
    Qnu[i] = 0.;
    Unu[i] = 0.;
    Vnu[i] = 0.;
    GYOTO_DEBUG_EXPR(Taunu[i]);
    if (Taunu[i]<1e-6){
    	//alphaInu[i] = std::numeric_limits<double>::infinity(); // Cause floatting point exception
    	alphaInu[i] = 1.e300; // something very big 
    }
    else alphaInu[i] = -log(Taunu[i])/(dsem*gg_->unitLength());
    GYOTO_DEBUG_EXPR(alphaInu[i]);
    Onu[i]=identity*alphaInu[i]; //Default transmission matrix with all polarisation set to 0, MUST be reimplemented
  }
  delete [] Taunu;
  delete [] alphaInu;
}

Matrix4d Generic::Omatrix(double alphanu[4], double rnu[3], double Chi, double dsem) const{
  
  return Omatrix(alphanu[0], alphanu[1], alphanu[2], alphanu[3], rnu[0], rnu[1], rnu[2], sin(2.*Chi), cos(2.*Chi), dsem);
}

Matrix4d Generic::Omatrix(double alphanu[4], double rnu[3], double sin2Chi, double cos2Chi, double dsem) const{

	return Omatrix(alphanu[0], alphanu[1], alphanu[2], alphanu[3], rnu[0], rnu[1], rnu[2], sin2Chi, cos2Chi, dsem);
}

Matrix4d Generic::Omatrix(double alphaInu, double alphaQnu, double alphaUnu, double alphaVnu,
        double rQnu, double rUnu, double rVnu, double Chi, double dsem) const{

	return Omatrix(alphaInu, alphaQnu, alphaUnu, alphaVnu, rQnu, rUnu, rVnu, sin(2.*Chi), cos(2.*Chi), dsem);
}

Matrix4d Generic::Omatrix(double alphaInu, double alphaQnu, double alphaUnu, double alphaVnu,
        double rQnu, double rUnu, double rVnu, double sin2Chi, double cos2Chi, double dsem) const{
	/** Function which compute the O matrix (see RadiativeTransfertVadeMecum) which represent the exponential
	*		of the Mueller Matrix containing the absorption and Faraday coefficients
	*/
	Matrix4d Onu;
	double alphasqrt, rsqrt, lambda1, lambda2, Theta, sigma;
  
  double aI=alphaInu, aV=alphaVnu;
  double aQ=alphaQnu*cos2Chi+alphaUnu*sin2Chi;
  double aU=alphaUnu*cos2Chi-alphaQnu*sin2Chi;
  double rQ=rQnu*cos2Chi+rUnu*sin2Chi;
  double rU=rUnu*cos2Chi-rQnu*sin2Chi;
  double rV=rVnu;

  alphasqrt = aQ*aQ+aU*aU+aV*aV;
  rsqrt = rQ*rQ+rU*rU+rV*rV;
  lambda1 = pow(pow(pow(alphasqrt-rsqrt,2.)/4.+pow(aQ*rQ+aU*rU+aV*rV,2.),0.5)+pow(alphasqrt-rsqrt,2.)/2.,0.5);
  lambda2 = pow(pow(pow(alphasqrt-rsqrt,2.)/4.+pow(aQ*rQ+aU*rU+aV*rV,2.),0.5)-pow(alphasqrt-rsqrt,2.)/2.,0.5);
  Theta = pow(lambda1,2)+pow(lambda2,2);
  sigma = (aQ*rQ+aU*rU+aV*rV) < 0 ? -1. : 1.;

  double coshlb1=cosh(lambda1*dsem*gg_->unitLength()),
  	sinhlb1=sinh(lambda1*dsem*gg_->unitLength()),
  	coslb2=cos(lambda2*dsem*gg_->unitLength()),
  	sinlb2=sin(lambda2*dsem*gg_->unitLength());

  if (coshlb1==coshlb1+1. || sinhlb1==sinhlb1+1.)
  	GYOTO_ERROR("In Omatrix : the cosh or sinh is infinite, one of the coefficient is to large !");

  Matrix4d zero;
  zero <<  0, 0, 0, 0,
           0, 0, 0, 0,
           0, 0, 0, 0,
           0, 0, 0, 0;
  Matrix4d M1, M2, M3, M4;

  // Fill of M1
  M1 = zero;
  for (int ii=0;ii<4;ii++){
    M1(ii,ii)=1.;
  }
  //cout << "M1 :\n" << M1 << endl;

  // Fill of M2
  M2 = zero;
  M2(0,1)=lambda2*aQ-sigma*lambda1*rQ;
  M2(0,2)=lambda2*aU-sigma*lambda1*rU;
  M2(0,3)=lambda2*aV-sigma*lambda1*rV;
  M2(1,0)=M2(0,1);
  M2(1,2)=sigma*lambda1*aV+lambda2*rV;
  M2(1,3)=-sigma*lambda1*aU-lambda2*rU;
  M2(2,0)=M2(0,2);
  M2(2,1)=-M2(1,2);
  M2(2,3)=sigma*lambda1*aQ+lambda2*rQ;
  M2(3,0)=M2(0,3);
  M2(3,1)=-M2(1,3);
  M2(3,2)=-M2(2,3);
  //cout << "M2 :\n" << M2 << endl;

  // Fill of M3
  M3 = zero;
  M3(0,1)=lambda1*aQ+sigma*lambda2*rQ;
  M3(0,2)=lambda1*aU+sigma*lambda2*rU;
  M3(0,3)=lambda1*aV+sigma*lambda2*rV;
  M3(1,0)=M3(0,1);
  M3(1,2)=-sigma*lambda2*aV+lambda1*rV;
  M3(1,3)=sigma*lambda2*aU-lambda1*rU;
  M3(2,0)=M3(0,2);
  M3(2,1)=-M3(1,2);
  M3(2,3)=-sigma*lambda2*aQ+lambda1*rQ;
  M3(3,0)=M3(0,3);
  M3(3,1)=-M3(1,3);
  M3(3,2)=-M3(2,3);
	//cout << "M3 :\n" << M3 << endl;

  // Fill of M4
  M4 = zero;
  M4(0,0)= (alphasqrt+rsqrt)/2.;
  M4(0,1)=aV*rU-aU*rV;
  M4(0,2)=aQ*rV-aV*rQ;
  M4(0,3)=aU*rQ-aQ*rU;
  M4(1,0)=-M4(0,1);
  M4(1,1)=pow(aQ,2)+pow(rQ,2)-(alphasqrt+rsqrt)/2.;
  M4(1,2)=aQ*aU+rQ*rU;
  M4(1,3)=aV*aQ+rV*rQ;
  M4(2,0)=-M4(0,2);
  M4(2,1)=M4(1,2);
  M4(2,2)=pow(aU,2)+pow(rU,2)-(alphasqrt+rsqrt)/2.;
  M4(2,3)=aU*aV+rU*rV;
  M4(3,0)=-M4(0,3);
  M4(3,1)=M4(1,3);
  M4(3,2)=M4(2,3);
  M4(3,3)=pow(aV,2)+pow(rV,2)-(alphasqrt+rsqrt)/2.;
	//cout << "M4 :\n" << M4 << endl;

  GYOTO_DEBUG
  	<< "alphaS : " << aI << ", " << aQ << ", " << aU << ", " << aV << "\n"
  	<< "rhoS   : " << rQ << ", " << rU << ", " << rV << "\n"
  	<< "alphasqrt : " << alphasqrt << "\n"
  	<< "rsqrt : " << rsqrt << "\n"
  	<< "lambda : " << lambda1 << ", " << lambda2 << "\n"
  	<< "Theta, sigma : " << Theta << ", " << sigma << "\n"
  	<< "dsem*gg_ : " << dsem*gg_->unitLength() << endl;


  Theta = (Theta==0.)?1.:Theta; // Theta equal zero means all coefficients are zero thus the O matrix is Identity and Theta should be 1.  
  // Filling O matrix, output
  Onu=exp(-aI*dsem*gg_->unitLength())*\
  			((coshlb1+coslb2)*M1/2. \
        -sinlb2*M2/Theta \
        -sinhlb1*M3/Theta \
        +(coshlb1-coslb2)*M4/Theta);
  return Onu;
}


Vector4d Generic::rotateJs(double jInu, double jQnu, double jUnu, double jVnu, double Chi) const{
	return rotateJs(jInu, jQnu, jUnu, jVnu, sin(2.*Chi),  cos(2.*Chi));
}

Vector4d Generic::rotateJs(double jInu, double jQnu, double jUnu, double jVnu, double sin2Chi, double cos2Chi) const{
	Matrix4d rot;
    rot << 1.,     0.  ,     0.  , 0.,
           0.,  cos2Chi, -sin2Chi, 0.,
           0.,  sin2Chi,  cos2Chi, 0.,
           0.,     0.  ,     0.  , 1.; // See RadiativeTransfertVadeMecum.pdf
  Vector4d jStokes;
    jStokes(0)=jInu;
    jStokes(1)=jQnu;
    jStokes(2)=jUnu;
    jStokes(3)=jVnu;
    jStokes = rot * jStokes;
    return jStokes;
}

double Generic::getChi(double const Bfourvect[4], state_t const &cph, double const vel[4]) const{
	double Chi=0;
	double sin2Chi, cos2Chi;
	getSinCos2Chi(Bfourvect, cph, vel, &sin2Chi, &cos2Chi);
 	Chi=atan2(sin2Chi,cos2Chi)/2.;

	return Chi;
}

void Generic::getSinCos2Chi(double const Bfourvect[4], state_t const &cph, double const vel[4], double* sin2Chi, double* cos2Chi) const{
	if (cph.size()!=16)
		GYOTO_ERROR("Ephi and Etheta not defined. Enable parrallel transport or implement the non polarised case in polarised RadiativeQ (see exemple in SimplePolarStar.C) ");

	double * Ephi = new double[4];
	double * Etheta = new double[4];
	double * Bproj = new double[4];
	memcpy(Bproj, Bfourvect, 4*sizeof(double));
	double photon_emframe[4]; // photon tgt vector projected in comoving frame; is it k projected ortho to u ?

	for (int ii=0;ii<4;ii++){
    photon_emframe[ii]=cph[ii+4]; // photon wave vector
    Ephi[ii]=cph[ii+8]; // polarization basis vector 1
    Etheta[ii]=cph[ii+12]; // polarization basis vector 2
  }

  // Check that wave vector, Ephi and Etheta are orthogonal:
  double test_tol=1e-2;
  if (fabs(gg_->ScalarProd(&cph[0],Ephi,Etheta))>test_tol
      or fabs(gg_->ScalarProd(&cph[0],Ephi,photon_emframe))>test_tol
      or fabs(gg_->ScalarProd(&cph[0],Etheta,photon_emframe))>test_tol){
      //or fabs(gg_->ScalarProd(&cph[0],Ephi,Ephi)-1.)>test_tol
      //or fabs(gg_->ScalarProd(&cph[0],Etheta,Etheta)-1.)>test_tol
      //or fabs(gg_->ScalarProd(&cph[0],photon_emframe,photon_emframe))>test_tol){
    cerr << "Ephi.Etheta, Ephi.K, Etheta.K)= " << fabs(gg_->ScalarProd(&cph[0],Ephi,Etheta)) << " " << fabs(gg_->ScalarProd(&cph[0],Ephi,photon_emframe)) << " " << fabs(gg_->ScalarProd(&cph[0],Etheta,photon_emframe)) << endl;
    throwError("Polarization basis is not properly parallel transported!");
  }

  // *** Projection into the rest frame of the emitter ***
  
  // NB: projectFourVect(pos, res, u) projects the initial res
  // orthogonally to u; so res is modified in the process.
  
  gg_->projectFourVect(&cph[0], photon_emframe, vel); // so from here on,
  // photon_emframe contains the projection orthogonal to u of the
  // wave vector. Hence the name: photon tangent vector in the rest
  // frame of the emitter, the space orthogonal to its 4vel.
  // Normalizing this projection (which is not by default):
  double norm=sqrt(gg_->ScalarProd(&cph[0],photon_emframe,photon_emframe));
  gg_->multiplyFourVect(photon_emframe,1./norm); // normalize basis vector
  // Check that the tetrad components are well normalized (assuming
  // a rotating emitter)
  double gtt=gg_->gmunu(&cph[0],0,0), grr=gg_->gmunu(&cph[0],1,1),
    gthth=gg_->gmunu(&cph[0],2,2), gpp=gg_->gmunu(&cph[0],3,3),
    Kr_tetrad=sqrt(grr)*photon_emframe[1],
    Kth_tetrad=sqrt(gthth)*photon_emframe[2],
    Kp_tetrad=sqrt(-gtt*gpp)*(photon_emframe[3]*vel[0]-photon_emframe[0]*vel[3]);
  //cout << "K.K in tetrad formalism= " << Kr_tetrad*Kr_tetrad + Kth_tetrad*Kth_tetrad + Kp_tetrad*Kp_tetrad << endl;
  
  gg_->projectFourVect(&cph[0],Bproj,vel); // Here we project the
  // magnetic 4vector orthogonal to u.
  // Normalize it:
  norm=sqrt(gg_->ScalarProd(&cph[0],Bproj,Bproj));
  gg_->multiplyFourVect(Bproj,1./norm);
  //cout << "Bproj BL compo= " << Bproj[0] << " " << Bproj[1] << " " << Bproj[2] << " " << Bproj[3] << endl;
  // Check the tetrd components:
  double Br_tetrad=sqrt(grr)*Bproj[1],
    Bth_tetrad=sqrt(gthth)*Bproj[2],
    Bp_tetrad=sqrt(-gtt*gpp)*(Bproj[3]*vel[0]-Bproj[0]*vel[3]);
  //cout << "Tetrad Bproj compo= " << Br_tetrad << " " << Bth_tetrad << " " << Bp_tetrad << endl;
  //cout << "B.B in tetrad formalism= " << Kr_tetrad*Kr_tetrad + Kth_tetrad*Kth_tetrad + Kp_tetrad*Kp_tetrad << endl;
  
  //gg_->projectFourVect(&cph[0],Ephi,vel); // Same for first polar basis vector

  // Modify the polarization basis vectors by adding to them a multiple
  // of the wavevector, which does not  affect the polarization (gauge freedom,
  // see Eric's notes on optics). This allows to obtain a well-defined
  // orthonormal triad in the emitter's rest frame, see FV rad transfer
  // notes for details.
  double * tmp = new double[4];
  for (int ii=0;ii<4;ii++){
    tmp[ii]=cph[ii+4]; // initialize to wave vector
  }
  gg_->multiplyFourVect(tmp,-gg_->ScalarProd(&cph[0],Ephi,vel)/gg_->ScalarProd(&cph[0],tmp,vel)); // this is well defined because the denominator cannot be zero,
  // it is the scalar prod of a timelike by a null vector
  gg_->addFourVect(Ephi, tmp); // At this point, Ephi lives in the rest frame
  // of the emitter and is orthogonal to photon_emframe.
  norm=sqrt(gg_->ScalarProd(&cph[0],Ephi,Ephi));
  gg_->multiplyFourVect(Ephi,1./norm); // normalize basis vector
  if (fabs(gg_->ScalarProd(&cph[0],Ephi,photon_emframe))>test_tol
      or fabs(gg_->ScalarProd(&cph[0],Ephi,vel))>test_tol
      or fabs(gg_->ScalarProd(&cph[0],Ephi,Ephi)-1.)>test_tol){
    throwError("Bad transformation of the polarization basis Ephi!");
  }

  // Same game for the second polarization basis vector:
  //gg_->projectFourVect(&cph[0],Etheta,vel); // and for second basis vector
  for (int ii=0;ii<4;ii++){
    tmp[ii]=cph[ii+4]; // initialize to wave vector
  }
  gg_->multiplyFourVect(tmp,-gg_->ScalarProd(&cph[0],Etheta,vel)/gg_->ScalarProd(&cph[0],tmp,vel)); 
  gg_->addFourVect(Etheta, tmp);
  norm=sqrt(gg_->ScalarProd(&cph[0],Etheta,Etheta));
  gg_->multiplyFourVect(Etheta,1./norm); // normalize basis vector
  if (fabs(gg_->ScalarProd(&cph[0],Etheta,photon_emframe))>test_tol
      or fabs(gg_->ScalarProd(&cph[0],Etheta,vel))>test_tol
      or fabs(gg_->ScalarProd(&cph[0],Etheta,Ephi))>test_tol
      or fabs(gg_->ScalarProd(&cph[0],Etheta,Etheta)-1.)>test_tol
      or fabs(gg_->ScalarProd(&cph[0],photon_emframe,photon_emframe)-1.)>test_tol){
    //cout << "K.K= " << gg_->ScalarProd(&cph[0],photon_emframe,photon_emframe) << endl;
    throwError("Bad transformation of the polarization basis Etheta!");
  }
  
  // So at this point, all vectors have been projected
  // in the rest frame of the emitter, orthogonal to its 4vel,
  // and we have a well-defined "observer's" orthonormal basis
  // (Ephi, Etheta, photon_emframe). We now need to define
  // the magnetic-field related basis in the rest frame,
  // and compute the rotation angle between the two.
  
  gg_->projectFourVect(&cph[0], Bproj, photon_emframe); // projection
  // orthogonally to photon_emframe of Bproj.
  // This is the second vector of the "emitter's magnetic basis"
  // (polarization vector, Bproj, K), where Bproj is what we have
  // just computed, ie the projection orth to the rest-frame
  // wave vector (K) of the rest-frame magnetic field, and the
  // polarization vector is along K*Bproj (* is a cross product here)
  // and completes the triad.
  // Check that it is still normal to 4vel:
  //cout << "Bperp.u= " << gg_->ScalarProd(&cph[0],Bproj,vel) << endl;
  // Normalize it:
  norm=sqrt(gg_->ScalarProd(&cph[0],Bproj,Bproj));
  gg_->multiplyFourVect(Bproj,1./norm); 
  // Check the tetrd components:
  double Bperp_r_tetrad=sqrt(grr)*Bproj[1],
    Bperp_th_tetrad=sqrt(gthth)*Bproj[2],
    Bperp_p_tetrad=sqrt(-gtt*gpp)*(Bproj[3]*vel[0]-Bproj[0]*vel[3]);
  // For checking, find the same compo from tetrad Euclidian computation:
  double KdotB=Kr_tetrad*Br_tetrad+Kth_tetrad*Bth_tetrad+Kp_tetrad*Bp_tetrad,
    Bperp_r_tetrad_BIS=Br_tetrad - KdotB*Kr_tetrad,
    Bperp_th_tetrad_BIS=Bth_tetrad - KdotB*Kth_tetrad,
    Bperp_p_tetrad_BIS=Bp_tetrad - KdotB*Kp_tetrad,
    Bperp_norm=sqrt(Bperp_r_tetrad_BIS*Bperp_r_tetrad_BIS+Bperp_th_tetrad_BIS*Bperp_th_tetrad_BIS+Bperp_p_tetrad_BIS*Bperp_p_tetrad_BIS);
  //cout << "Comparing tetrad components of B perp= " << endl;
  //cout << "First: " << Bperp_r_tetrad << " " << Bperp_th_tetrad << " " << Bperp_p_tetrad << endl;
  //cout << "Thenn: " << 1./Bperp_norm*Bperp_r_tetrad_BIS << " " << 1./Bperp_norm*Bperp_th_tetrad_BIS << " " << 1./Bperp_norm*Bperp_p_tetrad_BIS << endl;
  
  //cout << "photon_emframe: " << photon_emframe[0] << "," << photon_emframe[1] << "," << photon_emframe[2] << "," << photon_emframe[3] << endl;
  //cout << "Bproj: " << Bproj[0] << "," << Bproj[1] << "," << Bproj[2] << "," << Bproj[3] << endl;
  //cout << "Ephi: " << Ephi[0] << "," << Ephi[1] << "," << Ephi[2] << "," << Ephi[3] << endl;

  // Finally compute the scalar products between Bproj and the
  // observer's polarization basis vectors
  double BperpEtheta=gg_->ScalarProd(&cph[0],Bproj,Etheta),
    BperpEphi=gg_->ScalarProd(&cph[0],Bproj,Ephi);
  // and deduce the angle chi between the polarization basis and
  // the projected magnetic field, which is the final quantity of interest
  // (see FV rad transfer notes for details on the formulas).
  *cos2Chi=(pow(BperpEtheta,2.)-pow(BperpEphi,2.))
    /(pow(BperpEtheta,2.)+pow(BperpEphi,2.));
  *sin2Chi=-2*(BperpEtheta*BperpEphi)/(pow(BperpEtheta,2.)+pow(BperpEphi,2.));

  // For a check let's compare cos(chi) computed in BL and tetrad:
  double Eth_r_tetrad=sqrt(grr)*Etheta[1],
    Eth_th_tetrad=sqrt(gthth)*Etheta[2],
    Eth_p_tetrad=sqrt(-gtt*gpp)*(Etheta[3]*vel[0]-Etheta[0]*vel[3]);
  double Ephi_r_tetrad=sqrt(grr)*Ephi[1],
    Ephi_th_tetrad=sqrt(gthth)*Ephi[2],
    Ephi_p_tetrad=sqrt(-gtt*gpp)*(Ephi[3]*vel[0]-Ephi[0]*vel[3]);
  Bperp_norm=sqrt(Bperp_r_tetrad*Bperp_r_tetrad+Bperp_th_tetrad*Bperp_th_tetrad+Bperp_p_tetrad*Bperp_p_tetrad);
  //cout << "In Astrobj At position: " << cph[0] << " " << cph[1] << " " << cph[2] << " " << cph[3] << endl;
  //cout << "cos chi = " << BperpEtheta << " " << 1./Bperp_norm*(Bperp_r_tetrad*Eth_r_tetrad + Bperp_th_tetrad*Eth_th_tetrad + Bperp_p_tetrad*Eth_p_tetrad) << endl;

  // Define the polarization vector in tetrad formalism:
  double polar_vec[3]={Kth_tetrad*Bperp_p_tetrad - Kp_tetrad*Bperp_th_tetrad,
		       -Kr_tetrad*Bperp_p_tetrad+Kp_tetrad*Bperp_r_tetrad,
		       Kr_tetrad*Bperp_th_tetrad-Kth_tetrad*Bperp_r_tetrad};
  //cout << "Check prod scal: " << polar_vec[0]*Bperp_r_tetrad+polar_vec[1]*Bperp_th_tetrad+polar_vec[2]*Bperp_p_tetrad << " " << polar_vec[0]*Kr_tetrad+polar_vec[1]*Kth_tetrad+polar_vec[2]*Kp_tetrad << " " << polar_vec[0]*polar_vec[0]+polar_vec[1]*polar_vec[1]+polar_vec[2]*polar_vec[2] << endl;

  //cout << "P.Etheta= " << polar_vec[0]*Eth_r_tetrad + polar_vec[1]*Eth_th_tetrad + polar_vec[2]*Eth_p_tetrad << endl;
  //cout << "Bperp.Etheta= " << Bperp_r_tetrad*Eth_r_tetrad + Bperp_th_tetrad*Eth_th_tetrad + Bperp_p_tetrad*Eth_p_tetrad << endl;

  //cout << "P.Ephi= " << polar_vec[0]*Ephi_r_tetrad + polar_vec[1]*Ephi_th_tetrad + polar_vec[2]*Ephi_p_tetrad << endl;
  //cout << "Bperp.Ephi= " << Bperp_r_tetrad*Ephi_r_tetrad + Bperp_th_tetrad*Ephi_th_tetrad + Bperp_p_tetrad*Ephi_p_tetrad << endl;

 	delete [] Bproj;
 	delete [] Ephi;
 	delete [] Etheta;
}

double Generic::getChi(double const basis1[4], double const basis2[4], double const basis3[4], state_t const &cph, double const vel[4]) const{
	double Chi=0;
	double sin2Chi, cos2Chi;
	getSinCos2Chi(basis1, basis2, basis3, cph, vel, &sin2Chi, &cos2Chi);
 	Chi=atan2(sin2Chi,cos2Chi)/2.;

	return Chi;
}

void Generic::getSinCos2Chi(double const basis1[4], double const basis2[4], double const basis3[4], state_t const &cph, double const vel[4], double* sin2Chi, double* cos2Chi) const{
	/**
	 * Compute the angle between the emission basis provided by the user and the parallel transported observer polarisation basis.
	*/

	if (cph.size()!=16)
		GYOTO_ERROR("Ephi and Etheta not defined. Enable parrallel transport or implement the non polarised case in polarised RadiativeQ (see exemple in SimplePolarStar.C) ");

	double * Ephi = new double[4];
	double * Etheta = new double[4];
	double * photon_emframe = new double[4];
	double * e1 = new double[4];
	double * e2 = new double[4];
	double * e3 = new double[4];

	for (int ii=0;ii<4;ii++){
		photon_emframe[ii]=cph[ii+4]; // copy of photon velocity
		Ephi[ii]=cph[ii+8];
		Etheta[ii]=cph[ii+12];
		e1[ii]=basis1[ii];
		e2[ii]=basis2[ii];
		e3[ii]=basis3[ii];
	}
	// Projection into the rest frame of the emitter
	gg_->projectFourVect(&cph[0], Ephi, vel);
	gg_->projectFourVect(&cph[0], Etheta, vel);
  gg_->projectFourVect(&cph[0], photon_emframe, vel);
  gg_->projectFourVect(&cph[0], e1, vel);
  gg_->projectFourVect(&cph[0], e2, vel);
  gg_->projectFourVect(&cph[0], e3, vel);


  // Check that the third vector of the basis is colinear to k projected along vel
  if(gg_->ScalarProd(&cph[0], photon_emframe, e3) != gg_->norm(&cph[0],photon_emframe)*gg_->norm(&cph[0],e3))
  	GYOTO_ERROR("The third vector of the provided basis is not colinear to the photon tangent vector projected in the rest frame of the emitter.");
  

  double e2Etheta=gg_->ScalarProd(&cph[0],e2,Etheta),
  	e2Ephi=gg_->ScalarProd(&cph[0],e2,Ephi);

	*cos2Chi=(pow(e2Etheta,2.)-pow(e2Ephi,2.))/(pow(e2Etheta,2.)+pow(e2Ephi,2.));
	*sin2Chi=-2*(e2Etheta*e2Ephi)/(pow(e2Etheta,2.)+pow(e2Ephi,2.));

  delete [] photon_emframe;
 	delete [] Ephi;
 	delete [] Etheta;
 	delete [] e1;
 	delete [] e2;
 	delete [] e3;
}

void Generic::integrateEmission(double * I, double const * boundaries,
				size_t const * chaninds, size_t nbnu,
				double dsem, state_t const &cph, double const *co) const
{
  for (size_t i=0; i<nbnu; ++i)
    I[i] = integrateEmission(boundaries[chaninds[2*i]],
			     boundaries[chaninds[2*i+1]],
			     dsem, cph, co);
}

double Generic::integrateEmission (double nu1, double nu2, double dsem,
				   state_t const &coord_ph, double const coord_obj[8])
  const {
  double nu;
  if(nu1>nu2) {nu=nu1; nu1=nu2; nu2=nu;}
  double Inu1 = emission(nu1, dsem, coord_ph, coord_obj);
  double Inu2 = emission(nu2, dsem, coord_ph, coord_obj);
  double dnux2 = ((nu2-nu1)*2.);
  double Icur = (Inu2+Inu1)*dnux2*0.25;
  double Iprev;

# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG_EXPR(Icur);
# endif

  do {
    Iprev = Icur;
    dnux2 *= 0.5;
    for (nu = nu1 + 0.5*dnux2; nu < nu2; nu += dnux2) {
      Icur += emission(nu, dsem, coord_ph, coord_obj) * dnux2;
    }
    Icur *= 0.5;
#   if GYOTO_DEBUG_ENABLED
    GYOTO_DEBUG_EXPR(Icur);
#   endif
  } while( fabs(Icur-Iprev) > (1e-2 * Icur) );

# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << "dnu=" << dnux2*0.5
	      << "=(nu2-nu1)/" << (nu2-nu1)/(dnux2*0.5) << endl;
# endif

  return Icur;
}

Quantity_t Generic::getDefaultQuantities() { return GYOTO_QUANTITY_INTENSITY; }

void Astrobj::initRegister() {
  if (Gyoto::Astrobj::Register_) delete Gyoto::Astrobj::Register_;
  Gyoto::Astrobj::Register_ = NULL;
}

double Generic::deltaMax(double coord[8]) {
  double rr=0.,h1max;

  if (!gg_)
    GYOTO_ERROR("Please set metric before calling Astrobj::Generic::deltaMax()");

  switch (gg_ -> coordKind()) {
  case GYOTO_COORDKIND_SPHERICAL:
    rr=coord[1];
    break;
  case GYOTO_COORDKIND_CARTESIAN:
    rr=sqrt(coord[1]*coord[1]+coord[2]*coord[2]+coord[3]*coord[3]);
    break;
  default:
    GYOTO_ERROR("Incompatible coordinate kind in Astrobj.C");
  }

  if (rr<rMax()) h1max=deltamaxinsidermax_; else h1max=rr*0.5;
  return h1max;
}

void Gyoto::Astrobj::Register(std::string name, Subcontractor_t* scp){
  Register::Entry* ne =
    new Register::Entry(name, (SmartPointee::Subcontractor_t*)scp, Gyoto::Astrobj::Register_);
  Gyoto::Astrobj::Register_ = ne;
}

GYOTO_GETSUBCONTRACTOR(Astrobj)

Astrobj::Properties::Properties() :
  intensity(NULL), time(NULL), distance(NULL),
  first_dmin(NULL), first_dmin_found(0),
  redshift(NULL), nbcrosseqplane(NULL),
  spectrum(NULL), stokesQ(NULL), stokesU(NULL), stokesV(NULL),
  binspectrum(NULL), offset(1), impactcoords(NULL),
  user1(NULL), user2(NULL), user3(NULL), user4(NULL), user5(NULL)
# ifdef HAVE_UDUNITS
  , intensity_converter_(NULL), spectrum_converter_(NULL),
  binspectrum_converter_(NULL)
# endif
  , alloc(false)
{}

Astrobj::Properties::Properties( double * I, double * t) :
  intensity(I), time(t), distance(NULL),
  first_dmin(NULL), first_dmin_found(0),
  redshift(NULL), nbcrosseqplane(NULL),
  spectrum(NULL), stokesQ(NULL), stokesU(NULL), stokesV(NULL),
  binspectrum(NULL), offset(1), impactcoords(NULL),
  user1(NULL), user2(NULL), user3(NULL), user4(NULL), user5(NULL)
# ifdef HAVE_UDUNITS
  , intensity_converter_(NULL), spectrum_converter_(NULL),
  binspectrum_converter_(NULL)
# endif
  , alloc(false)
{}

void Astrobj::Properties::init(size_t nbnuobs) {
  if (intensity)      *intensity  = 0.;
  if (time)           *time       = DBL_MAX;
  if (distance)       *distance   = DBL_MAX;
  if (first_dmin){    *first_dmin = DBL_MAX; first_dmin_found=0;}
  if (redshift)       *redshift   = 0.;
  if (nbcrosseqplane) *nbcrosseqplane   = 0.;
  if (spectrum) for (size_t ii=0; ii<nbnuobs; ++ii) spectrum[ii*offset]=0.;
  if (stokesQ)  for (size_t ii=0; ii<nbnuobs; ++ii) stokesQ [ii*offset]=0.;
  if (stokesU)  for (size_t ii=0; ii<nbnuobs; ++ii) stokesU [ii*offset]=0.;
  if (stokesV)  for (size_t ii=0; ii<nbnuobs; ++ii) stokesV [ii*offset]=0.;
  if (binspectrum) for (size_t ii=0; ii<nbnuobs; ++ii)
		     binspectrum[ii*offset]=0.;
  if (impactcoords) for (size_t ii=0; ii<16; ++ii) impactcoords[ii]=DBL_MAX;
  if (user1)          *user1=0.;
  if (user2)          *user2=0.;
  if (user3)          *user3=0.;
  if (user4)          *user4=0.;
  if (user5)          *user5=0.;
}

#ifdef HAVE_UDUNITS
void Astrobj::Properties::intensityConverter(SmartPointer<Units::Converter> conv) {
  intensity_converter_ = conv ;
}

void Astrobj::Properties::intensityConverter(string unit) {
  intensity_converter_ =
    new Units::Converter("J.m-2.s-1.sr-1.Hz-1",
			 unit!=""?unit:"J.m-2.s-1.sr-1.Hz-1");
}

void Astrobj::Properties::spectrumConverter(SmartPointer<Units::Converter> conv) {
  spectrum_converter_ = conv;
}

void Astrobj::Properties::spectrumConverter(string unit) {
  spectrum_converter_ =
    new Units::Converter("J.m-2.s-1.sr-1.Hz-1",
			 unit!=""?unit:"J.m-2.s-1.sr-1.Hz-1");
}

void Astrobj::Properties::binSpectrumConverter(SmartPointer<Units::Converter> conv) {
  binspectrum_converter_ = conv;
}
void Astrobj::Properties::binSpectrumConverter(string unit) {
  binspectrum_converter_ =
    new Units::Converter("J.m-2.s-1.sr-1",
			 unit!=""?unit:"J.m-2.s-1.sr-1");
}

#endif

Astrobj::Properties& Astrobj::Properties::operator+=(ptrdiff_t ofset) {
  if (intensity)      intensity      += ofset;
  if (time)           time           += ofset;
  if (distance)       distance       += ofset;
  if (first_dmin)     first_dmin     += ofset;
  if (redshift)       redshift       += ofset;
  if (nbcrosseqplane) nbcrosseqplane += ofset;
  if (spectrum)       spectrum       += ofset;
  if (stokesQ)        stokesQ        += ofset;
  if (stokesU)        stokesU        += ofset;
  if (stokesV)        stokesV        += ofset;
  if (binspectrum)    binspectrum    += ofset;
  if (impactcoords)   impactcoords   += 16*ofset;
  if (user1)          user1          += ofset;
  if (user2)          user2          += ofset;
  if (user3)          user3          += ofset;
  if (user4)          user4          += ofset;
  if (user5)          user5          += ofset;
  return *this;
}

Astrobj::Properties& Astrobj::Properties::operator++() {
  (*this) += 1;
  return *this;
}

Astrobj::Properties::operator Quantity_t () const {
  Quantity_t res=GYOTO_QUANTITY_NONE;

  if (intensity)      res |= GYOTO_QUANTITY_INTENSITY;
  if (time)           res |= GYOTO_QUANTITY_EMISSIONTIME;
  if (distance)       res |= GYOTO_QUANTITY_MIN_DISTANCE;
  if (first_dmin)     res |= GYOTO_QUANTITY_FIRST_DMIN;
  if (redshift)       res |= GYOTO_QUANTITY_REDSHIFT;
  if (nbcrosseqplane) res |= GYOTO_QUANTITY_NBCROSSEQPLANE;
  if (spectrum)       res |= GYOTO_QUANTITY_SPECTRUM;
  if (stokesQ)        res |= GYOTO_QUANTITY_SPECTRUM_STOKES_Q;
  if (stokesU)        res |= GYOTO_QUANTITY_SPECTRUM_STOKES_U;
  if (stokesV)        res |= GYOTO_QUANTITY_SPECTRUM_STOKES_V;
  if (binspectrum)    res |= GYOTO_QUANTITY_BINSPECTRUM;
  if (impactcoords)   res |= GYOTO_QUANTITY_IMPACTCOORDS;
  if (user1)          res |= GYOTO_QUANTITY_USER1;
  if (user2)          res |= GYOTO_QUANTITY_USER2;
  if (user3)          res |= GYOTO_QUANTITY_USER3;
  if (user4)          res |= GYOTO_QUANTITY_USER4;
  if (user5)          res |= GYOTO_QUANTITY_USER5;

  return res;
}
