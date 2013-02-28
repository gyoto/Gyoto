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

// GYOTO HEADERS
#include "GyotoUtils.h"
#include "GyotoAstrobj.h"
#include "GyotoMetric.h"
#include "GyotoPhoton.h"
#include "GyotoRegister.h"
#include "GyotoFactoryMessenger.h"

// SYSTEM HEADERS
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <cfloat>
#include <cmath>
#include <sstream>

// NAMESPACES
using namespace std;
using namespace Gyoto;
using namespace Gyoto::Astrobj;

#if defined GYOTO_USE_XERCES
Register::Entry* Gyoto::Astrobj::Register_ = NULL;
#endif

Generic::Generic(string kind) :

  gg_(NULL), rmax_(DBL_MAX), rmax_set_(0), kind_(kind), flag_radtransf_(0),
  bolometric_(0)
{
#if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << endl;
#endif
}

Generic::Generic() :

  gg_(NULL), rmax_(DBL_MAX), rmax_set_(0), kind_("Default"), flag_radtransf_(0),
  bolometric_(0)
{
#if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << endl;
#endif
}

Generic::Generic(double radmax) :
  gg_(NULL), rmax_(radmax), rmax_set_(1), kind_("Default"), flag_radtransf_(0),
  bolometric_(0)
{
#if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << endl;
#endif
}

Generic::Generic(const Generic& orig) :
  SmartPointee(orig), gg_(NULL),
  rmax_(orig.rmax_), rmax_set_(orig.rmax_set_), kind_(orig.kind_),
  flag_radtransf_(orig.flag_radtransf_),
  bolometric_(orig.bolometric_)
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

SmartPointer<Metric::Generic> Generic::getMetric() const { return gg_; }
void Generic::setMetric(SmartPointer<Metric::Generic> gg) {gg_=gg;}

double Generic::getRmax() {
  return rmax_;
}

double Generic::getRmax(string unit) {
  return Units::FromGeometrical(getRmax(), unit, gg_);
}

const string Generic::getKind() const {
  return kind_;
}

void Generic::setRmax(double val) {
  rmax_set_=1;
  rmax_=val;
}

void Generic::setRmax(double val, string unit) {
  setRmax(Units::ToGeometrical(val, unit, gg_));
}

void Generic::unsetRmax() {
  rmax_set_=0;
}

#ifdef GYOTO_USE_XERCES
void Generic::fillElement(FactoryMessenger *fmp) const {
  if (rmax_set_) fmp -> setParameter ( "RMax", rmax_ ) ;
  fmp -> setMetric(gg_);
  fmp -> setSelfAttribute("kind", kind_);
  fmp -> setParameter ( flag_radtransf_? "OpticallyThin" : "OpticallyThick");
  fmp -> setParameter ( bolometric_? "Bolometric" : "Specific");
}

void Generic::setParameters(FactoryMessenger *fmp) {
  string name="", content="", unit="";
  if (fmp) {
    setMetric(fmp->getMetric());
    while (fmp->getNextParameter(&name, &content, &unit))
      setParameter(name, content, unit);
  }
}
#endif


void Generic::setFlag_radtransf(int flag) {flag_radtransf_=flag;}
int Generic::getFlag_radtransf() const {return flag_radtransf_;}

void Generic::setBolometric(int bolo) {bolometric_=bolo;}
int Generic::getBolometric() const {return bolometric_;}

int Generic::setParameter(string name, string content, string unit)  {
  char* tc = const_cast<char*>(content.c_str());
  if (name=="Flag_radtransf")  flag_radtransf_= atoi(tc);
  else if (name=="OpticallyThin")  flag_radtransf_= 1;
  else if (name=="OpticallyThick")  flag_radtransf_= 0;
  else if (name=="Specific")  bolometric_= 0;
  else if (name=="Bolometric")  bolometric_= 1;
  else if (name=="RMax") setRmax(atof(tc), unit);
  else return 1;
  return 0;
}

void Generic::processHitQuantities(Photon* ph, double* coord_ph_hit,
				     double* coord_obj_hit, double dt,
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
  double freqObs=ph->getFreqObs(); // this is a useless quantity, always 1
  SmartPointer<Spectrometer> spr = ph -> getSpectrometer();
  size_t nbnuobs = spr() ? spr -> getNSamples() : 0 ;
  double const * const nuobs = nbnuobs ? spr -> getMidpoints() : NULL;
  double dlambda = dt/coord_ph_hit[4]; //dlambda = dt/tdot
  double ggredm1 = -gg_->ScalarProd(coord_ph_hit,coord_obj_hit+4,
				    coord_ph_hit+4) / freqObs; 
                                       //this is nu_em/nu_obs
  double ggred = 1./ggredm1;           //this is nu_obs/nu_em
  double dsem = dlambda*freqObs*ggredm1;
  double inc =0.;
  if (data) {
#if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << "data requested" 
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
    if (data->impactcoords) {
      memcpy(data->impactcoords, coord_obj_hit, 8 * sizeof(double));
      memcpy(data->impactcoords+8, coord_ph_hit, 8 * sizeof(double));
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
      if (!bolometric_) // Specific intensity given by emission (I_nu)
	inc = (emission(freqObs*ggredm1, dsem, coord_ph_hit, coord_obj_hit))
	  * (ph -> getTransmission(size_t(-1)))
	  * ggred*ggred*ggred; // I_nu/nu^3 invariant
      else // Bolometric intensity given by emission (I)
	inc = (emission(freqObs*ggredm1, dsem, coord_ph_hit, coord_obj_hit))
	  * (ph -> getTransmission(size_t(-1)))
	  * ggred*ggred*ggred*ggred; // I/nu^4 invariant
#     ifdef HAVE_UDUNITS
      if (data -> intensity_converter_)
	inc = (*data -> intensity_converter_)(inc);
#     endif
      *data->intensity += inc;

      if (!bolometric_)
#     if GYOTO_DEBUG_ENABLED
	GYOTO_DEBUG
	  << "intensity +=" << *data->intensity
	  << "= emission((dsem=" << dsem << "))="
	  << (emission(freqObs*ggredm1,dsem,coord_ph_hit, coord_obj_hit))
	  << ")*(ggred=" << ggred << ")^3*(transmission="
	  << (ph -> getTransmission(size_t(-1))) << ")"
	  << endl;
#     endif
      else
#     if GYOTO_DEBUG_ENABLED
	GYOTO_DEBUG
	  << "intensity +=" << *data->intensity
	  << "= emission((dsem=" << dsem << "))="
	  << (emission(freqObs*ggredm1,dsem,coord_ph_hit, coord_obj_hit))
	  << ")*(ggred=" << ggred << ")^4*(transmission="
	  << (ph -> getTransmission(size_t(-1))) << ")"
	  << endl;
#     endif
      
    }
    if (data->binspectrum) {
      double const * const channels = spr -> getChannels();
      double * I  = new double[nbnuobs];
      double * boundaries = new double[nbnuobs+1];
      for (size_t ii=0; ii<=nbnuobs; ++ii)
	boundaries[ii]=channels[ii]*ggredm1;
      integrateEmission(I, boundaries, nbnuobs,
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
      }
      delete [] I;
      delete [] boundaries;
    }
    if (data->spectrum) {
      double * Inu  = new double[nbnuobs];
      double * nuem = new double[nbnuobs];
      for (size_t ii=0; ii<nbnuobs; ++ii) {
	nuem[ii]=nuobs[ii]*ggredm1;
      }
      emission(Inu, nuem, nbnuobs, dsem, coord_ph_hit, coord_obj_hit);
      for (size_t ii=0; ii<nbnuobs; ++ii) {
	inc = Inu[ii] * ph -> getTransmission(ii) * ggred*ggred*ggred;
#       ifdef HAVE_UDUNITS
	if (data -> spectrum_converter_)
	  inc = (*data -> spectrum_converter_)(inc);
#       endif
	data->spectrum[ii*data->offset] += inc;
	 
#       if GYOTO_DEBUG_ENABLED
	GYOTO_DEBUG
	       << "DEBUG: Generic::processHitQuantities(): "
	       << "nuobs[" << ii << "]="<< nuobs[ii]
	       << ", nuem=" << nuem 
	       << ", dsem=" << dsem
	       << ", Inu * GM/c2="
	       << Inu[ii]
	       << ", spectrum[" << ii*data->offset << "]="
	       << data->spectrum[ii*data->offset]
	       << ", transmission=" << ph -> getTransmission(ii)
	       << ", redshift=" << ggred << ")\n";
#       endif
      }
      delete [] Inu;
      delete [] nuem;
    }
    /* update photon's transmission */
    ph -> transmit(size_t(-1),
		   transmission(freqObs*ggredm1, dsem,coord_ph_hit));
    for (size_t ii=0; ii<nbnuobs; ++ii)
      ph -> transmit(ii,transmission(nuobs[ii]*ggredm1,dsem,coord_ph_hit));
  } else {
#   if GYOTO_DEBUG_ENABLED
    GYOTO_DEBUG << "NO data requested!" << endl;
#   endif
  }
}

double Generic::transmission(double, double, double*) const {
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG_EXPR(flag_radtransf_);
# endif
  return double(flag_radtransf_);
}

double Generic::emission(double , double dsem, double *, double *) const
{
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG_EXPR(flag_radtransf_);
# endif
  if (flag_radtransf_) return dsem;
  return 1.;
}

void Generic::emission(double * Inu, double * nuem , size_t nbnu,
			 double dsem, double *cph, double *co) const
{
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG_EXPR(flag_radtransf_);
# endif
  for (size_t i=0; i< nbnu; ++i) Inu[i]=emission(nuem[i], dsem, cph, co);
}

void Generic::integrateEmission(double * I, double * boundaries, size_t nbnu,
				double dsem, double *cph, double *co) const
{
  for (size_t i=0; i<nbnu; ++i)
    I[i] = integrateEmission(boundaries[i], boundaries[i+1], dsem, cph, co);
}

double Generic::integrateEmission (double nu1, double nu2, double dsem,
				   double coord_ph[8], double coord_obj[8])
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

#if defined GYOTO_USE_XERCES
void Astrobj::initRegister() {
  if (Gyoto::Astrobj::Register_) delete Gyoto::Astrobj::Register_;
  Gyoto::Astrobj::Register_ = NULL;
}

void Gyoto::Astrobj::Register(std::string name, Subcontractor_t* scp){
  Register::Entry* ne =
    new Register::Entry(name, (SmartPointee::Subcontractor_t*)scp, Gyoto::Astrobj::Register_);
  Gyoto::Astrobj::Register_ = ne;
}

Gyoto::Astrobj::Subcontractor_t*
Astrobj::getSubcontractor(std::string name, int errmode) {
  if (!Gyoto::Astrobj::Register_) throwError("No Astrobj kind registered!");
  return (Subcontractor_t*)Gyoto::Astrobj::Register_
    -> getSubcontractor(name, errmode);
}
#endif


Astrobj::Properties::Properties() :
  intensity(NULL), time(NULL), distance(NULL),
  first_dmin(NULL), first_dmin_found(0),
  redshift(NULL),
  spectrum(NULL), binspectrum(NULL), offset(1), impactcoords(NULL),
  user1(NULL), user2(NULL), user3(NULL), user4(NULL), user5(NULL)
# ifdef HAVE_UDUNITS
  , intensity_converter_(NULL), spectrum_converter_(NULL),
  binspectrum_converter_(NULL)
# endif
{}

Astrobj::Properties::Properties( double * I, double * t) :
  intensity(I), time(t), distance(NULL),
  first_dmin(NULL), first_dmin_found(0),
  redshift(NULL),
  spectrum(NULL), binspectrum(NULL), offset(1), impactcoords(NULL),
  user1(NULL), user2(NULL), user3(NULL), user4(NULL), user5(NULL)
# ifdef HAVE_UDUNITS
  , intensity_converter_(NULL), spectrum_converter_(NULL),
  binspectrum_converter_(NULL)
# endif
{}

void Astrobj::Properties::init(size_t nbnuobs) {
  if (intensity)  *intensity  = 0.;
  if (time)       *time       = DBL_MAX;
  if (distance)   *distance   = DBL_MAX;
  if (first_dmin){*first_dmin = DBL_MAX; first_dmin_found=0;}
  if (redshift)   *redshift   = 0.;
  if (spectrum) for (size_t ii=0; ii<nbnuobs; ++ii) spectrum[ii*offset]=0.; 
  if (binspectrum) for (size_t ii=0; ii<nbnuobs; ++ii)
		     binspectrum[ii*offset]=0.; 
  if (impactcoords) for (size_t ii=0; ii<16; ++ii) impactcoords[ii]=DBL_MAX;
  if (user1)      *user1=0.;
  if (user2)      *user2=0.;
  if (user3)      *user3=0.;
  if (user4)      *user4=0.;
  if (user5)      *user5=0.;
}

#ifdef HAVE_UDUNITS
void Astrobj::Properties::setIntensityConverter(SmartPointer<Units::Converter> conv) {
  intensity_converter_ = conv ;
}

void Astrobj::Properties::setIntensityConverter(string unit) {
  intensity_converter_ =
    new Units::Converter("J.m-2.s-1.sr-1.Hz-1",
			 unit!=""?unit:"J.m-2.s-1.sr-1.Hz-1");
}

void Astrobj::Properties::setSpectrumConverter(SmartPointer<Units::Converter> conv) {
  spectrum_converter_ = conv;
}

void Astrobj::Properties::setSpectrumConverter(string unit) {
  spectrum_converter_ =
    new Units::Converter("J.m-2.s-1.sr-1.Hz-1",
			 unit!=""?unit:"J.m-2.s-1.sr-1.Hz-1");
}

void Astrobj::Properties::setBinSpectrumConverter(SmartPointer<Units::Converter> conv) {
  binspectrum_converter_ = conv;
}
void Astrobj::Properties::setBinSpectrumConverter(string unit) {
  binspectrum_converter_ =
    new Units::Converter("J.m-2.s-1.sr-1",
			 unit!=""?unit:"J.m-2.s-1.sr-1");
}

#endif

Astrobj::Properties Astrobj::Properties::operator++() {
  if (intensity)  ++intensity;
  if (time)       ++time;
  if (distance)   ++distance;
  if (first_dmin) ++first_dmin;
  if (redshift)   ++redshift;
  if (spectrum)   ++spectrum;
  if (binspectrum)++binspectrum;
  if (impactcoords) impactcoords += 16;
  if (user1)      ++user1;
  if (user2)      ++user2;
  if (user3)      ++user3;
  if (user4)      ++user4;
  if (user5)      ++user5;
  return *this;
}
