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

// SYSTEM HEADERS
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <float.h>
#include <cmath>
#include <sstream>

// NAMESPACES
using namespace std;
using namespace Gyoto;

Register::Entry* Gyoto::Astrobj::Register_ = NULL;

Astrobj::Astrobj(string kind) :

  gg_(NULL), rmax_(DBL_MAX), rmax_set_(0), kind_(kind), flag_radtransf_(0),
  critical_value_(DBL_MIN), safety_value_(DBL_MAX)
{
  if (debug()) cerr << "Astrobj Construction" << endl;
  
}

Astrobj::Astrobj() :

  gg_(NULL), rmax_(DBL_MAX), rmax_set_(0), kind_("Default"), flag_radtransf_(0),
  critical_value_(DBL_MIN), safety_value_(DBL_MAX)
{
  if (debug()) cerr << "Astrobj Construction" << endl;
}

Astrobj::Astrobj(double radmax) :
  gg_(NULL), rmax_(radmax), rmax_set_(1), kind_("Default"), flag_radtransf_(0),
  critical_value_(DBL_MIN)
{
  if (debug()) cerr << "Astrobj Construction" << endl;
}
Astrobj::Astrobj(const Astrobj& orig) :
  SmartPointee(orig), gg_(NULL),
  rmax_(orig.rmax_), rmax_set_(orig.rmax_set_), kind_(orig.kind_),
  flag_radtransf_(orig.flag_radtransf_), critical_value_(orig.critical_value_),
  safety_value_(orig.safety_value_)
{
    if (debug()) cerr << "DEBUG: in Astrobj::Astrobj (Copy)" << endl;
    if (orig.gg_()) {
      if (debug())
	cerr << "DEBUG: orig had a metric, cloning" << endl;
      gg_=orig.gg_->clone();
    }
    if (debug()) cerr << "DEBUG: out of Astrobj::Astrobj (Copy)" << endl;
}
Astrobj * Astrobj::clone() const {
  string msg = "Astrobj::clone() called: cloning not supported "
    "for astrobj kind ";
  msg += getKind();
  throwError (msg);
  return const_cast<Astrobj*>(this); // to avoid warning
}

Astrobj::~Astrobj() {
  if (debug()) cerr << "Astrobj Destruction" << endl;
}

SmartPointer<Metric> Astrobj::getMetric() const { return gg_; }
void Astrobj::setMetric(SmartPointer<Metric> gg) {gg_=gg;}

double Astrobj::getRmax() {
  return rmax_;
}

const string Astrobj::getKind() const {
  return kind_;
}

void Astrobj::setRmax(double val) {
  rmax_set_=1;
  rmax_=val;
}

void Astrobj::unsetRmax() {
  rmax_set_=0;
}

double Astrobj::operator()(double const * const) {
  stringstream ss;
  ss << getKind() << "::operator()() is unimplemented";
  throwError(ss.str());
  return 0.;
}

void Astrobj::getVelocity(double const * const, double *) {
  stringstream ss;
  ss << getKind() << "::getVelocity() is unimplemented";
  throwError(ss.str());
}

int Astrobj::Impact(Photon* ph, size_t index, AstrobjProperties *data){
  double p1[8], p2[8];
  ph->getCoord(index, p1);
  ph->getCoord(index+1, p2);
  double tmin, minval;

  double t1 = p1[0], t2=p2[0];
  double val1=(*this)(p1), val2=(*this)(p2);

  if (val1 > critical_value_) {
    if (val2 > critical_value_) {
      if ( val1 > safety_value_ && val2 > safety_value_) {
	if (val1 < val2) {
	  minval = val1; tmin = t1;
	} else {
	  minval = val2; tmin = t2;
	}
      } else
	minval = ph -> findMin(this, p1[0], p2[0], tmin, critical_value_) ;
      if (minval>critical_value_) {
	if (data) {
	  /* EmissionTime */
	  if (data->time) *data->time=tmin;
	  /* MinDistance */
	  if ((data->distance) && (*(data->distance)>minval) )
	    *data->distance=minval;
	  /* FirstMinDist */
	  if (data->first_dmin) { 
	    if (!data->first_dmin_found) {
	      if (*(data->first_dmin)>minval) *(data->first_dmin)=minval;
	      else data->first_dmin_found=1;
	    }
	  }
	}
	return 0;
      }
      ph -> findValue(this, critical_value_, tmin, t2);
    }
    ph -> findValue(this, critical_value_, t2, t1);
  } else if (val2 > critical_value_)
    ph -> findValue(this, critical_value_, t1, t2);

  double cph[8] = { t2 };
  ph -> getCoord(&t2, 1, cph+1, cph+2, cph+3,
		 cph+4, cph+5, cph+6, cph+7);

  double coh[8] = { t2 , cph[1], cph[2], cph[3] };
  getVelocity(coh, coh+4);

  processHitQuantities(ph, cph, coh, t2-t1, data);

  return 1;

}

void Astrobj::fillElement(factoryMessenger *fmp) const {
  fmp -> setSelfAttribute("kind", kind_);
  //  fmp -> setParameter ("Flag_radtransf", flag_radtransf_);
  fmp -> setParameter ( flag_radtransf_? "OpticallyThin" : "OpticallyThick");
}

void Astrobj::setFlag_radtransf(int flag) {flag_radtransf_=flag;}
int Astrobj::getFlag_radtransf() const {return flag_radtransf_;}

void Astrobj::setGenericParameter(string name, string content)  {
  char* tc = const_cast<char*>(content.c_str());
  if (name=="Flag_radtransf")  flag_radtransf_= atoi(tc);
  else if (name=="OpticallyThin")  flag_radtransf_= 1;
  else if (name=="OpticallyThick")  flag_radtransf_= 0;
  else if (name=="RMax")  {
    rmax_ = atof(tc); rmax_set_=1;
  }
}
void Astrobj::processHitQuantities(Photon* ph, double* coord_ph_hit,
				     double* coord_obj_hit, double dt,
				     AstrobjProperties* data) const {
  if (debug())
    cerr << "DEBUG: in Astrobj::processHitQuantities:" << endl;
  /*
      NB: freqObs is the observer's frequency chosen in
      Screen::getRayCoord for the actual computation of the geodesic ;
      the physical value of nuobs will be used in spectrum
      computations by resorting to the xml specifications of the user
      (see below) ; this freqObs is used to transform the null
      worldline parameter dlambda (see below)
  */
  double freqObs=ph->getFreqObs(), nuem; 
  SmartPointer<Spectrometer> spr = ph -> getSpectrometer();
  size_t nbnuobs = spr() ? spr -> getNSamples() : 0 ;
  double const * const nuobs = nbnuobs ? spr -> getMidpoints() : NULL;
  double dlambda = dt/coord_ph_hit[4]; //dlambda = dt/tdot
  double ggredm1 = -gg_->ScalarProd(coord_ph_hit,coord_obj_hit+4,
				    coord_ph_hit+4) / freqObs; 
                                       //this is nu_em/nu_obs
  double ggred = 1./ggredm1;           //this is nu_obs/nu_em
  double dsem = dlambda*freqObs*ggredm1;
  if (data) {
    if (debug())
      cerr << "DEBUG: Astrobj::processHitQuantities: data requested" << endl;

    if (debug())
      cerr << "DEBUG: Astrobj::processHitQuantities: data requested, "
	   << "freqObs=" << freqObs << ", ggredm1=" << ggredm1
	   << ", ggred=" << ggred
	   << endl;


    if (data->redshift) {
      *data->redshift=ggred;
      if (debug())
	cerr << "DEBUG: Astrobj::processHitQuantities(): "
	     << "redshift=" << *data->redshift << endl;
    }
    if (data->time) {
      *data->time=coord_ph_hit[0];
      if (debug())
	cerr << "DEBUG: Astrobj::processHitQuantities(): "
	     << "time=" << *data->time << endl;
    }
    if (debug())
      cerr << "DEBUG: Astrobj::processHitQuantities: "
	   << "dlambda = (dt="<< dt << ")/(tdot="<< coord_ph_hit[4]
	   << ") = " << dlambda << ", dsem=" << dsem << endl;
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
	Astrobj [NB: so far, with rad. transfer emission() computes
	j_nu, without rad. transfer it computes I_nu] Finally:
	I_nu_obs = I_nu_em*(nu_obs/nu_em)^3
      */

      //I_nu_obs increment :
      *data->intensity += 
	(emission(freqObs*ggredm1, dsem, coord_ph_hit, coord_obj_hit))
	* (ph -> getTransmission(size_t(-1)))
	* ggred*ggred*ggred;

      if (debug())
	cerr << "DEBUG: Astrobj::processHitQuantities(): "
	     << "intensity +=" << *data->intensity
	     << "= emission((dsem=" << dsem << "))="
	     << (emission(freqObs*ggredm1,dsem,coord_ph_hit, coord_obj_hit))
	     << ")*(ggred=" << ggred << ")^3*(transmission="
	     << (ph -> getTransmission(size_t(-1))) << ")"
	     << endl;
    }
    if (data->binspectrum) {
      double const * const channels = spr -> getChannels();
      double nuem1, nuem2;
      for (size_t ii=0; ii<nbnuobs; ++ii) {
	nuem1=channels[ii]*ggredm1;
	nuem2=channels[ii+1]*ggredm1;
	data->binspectrum[ii*data->offset] +=
	  integrateEmission(nuem1, nuem2, dsem,
			    coord_ph_hit, coord_obj_hit)
	  * ph -> getTransmission(ii)
	  *ggred*ggred*ggred*ggred;
	if (debug())
	  cerr << "DEBUG: Astrobj::processHitQuantities(): "
	       << "nuobs[" << ii << "]="<< channels[ii]
	       << ", nuem=" << nuem1 
	       << ", binspectrum[" << ii+data->offset << "]="
	       << data->binspectrum[ii*data->offset]<< endl;
      }
    }
    if (data->spectrum) {
      for (size_t ii=0; ii<nbnuobs; ++ii) {
	nuem=nuobs[ii]*ggredm1;
	data->spectrum[ii*data->offset] +=
	  emission(nuem, dsem, coord_ph_hit, coord_obj_hit)
	  * ph -> getTransmission(ii)
	  *ggred*ggred*ggred;
	if (debug())
	  cerr << "DEBUG: Astrobj::processHitQuantities(): "
	       << "nuobs[" << ii << "]="<< nuobs[ii]
	       << ", nuem=" << nuem 
	       << ", dsem=" << dsem
	       << ", Inu * GM/c2="
	       << emission(nuem, dsem, coord_ph_hit, coord_obj_hit)
	       << ", spectrum[" << ii*data->offset << "]="
	       << data->spectrum[ii*data->offset]
	       << ", transmission=" << ph -> getTransmission(ii)
	       << ", redshift=" << ggred << ")\n";
      }
    }
  } else {
    if (debug())
      cerr << "DEBUG: Astrobj::processHitQuantities: NO data requested!\n";
  }
  /* update photon's transmission */
  ph -> transmit(size_t(-1),
		 transmission(freqObs*ggredm1, dsem,coord_ph_hit));
  for (size_t ii=0; ii<nbnuobs; ++ii)
    ph -> transmit(ii,transmission(nuobs[ii]*ggredm1,dsem,coord_ph_hit));
}

double Astrobj::transmission(double, double, double*) const {
  if (debug())
    cerr << "DEBUG: Astrobj::transmission(): flag_radtransf_="
	 << flag_radtransf_ << endl;
  return double(flag_radtransf_);
}

double Astrobj::integrateEmission (double nu1, double nu2, double dsem,
				   double coord_ph[8], double coord_obj[8])
  const {
  double nu;
  if(nu1>nu2) {nu=nu1; nu1=nu2; nu2=nu;}
  double Inu1 = emission(nu1, dsem, coord_ph, coord_obj);
  double Inu2 = emission(nu2, dsem, coord_ph, coord_obj);
  double dnux2 = ((nu2-nu1)*2.);
  double Icur = (Inu2+Inu1)*dnux2*0.25;
  double Iprev;

  if (debug())
      cerr << "DEBUG: Astrobj::integrateEmission(): "
	   << "Icur=" << Icur << endl;

  do {
    Iprev = Icur; 
    dnux2 *= 0.5;
    for (nu = nu1 + 0.5*dnux2; nu < nu2; nu += dnux2) {
      Icur += emission(nu, dsem, coord_ph, coord_obj) * dnux2;
    }
    Icur *= 0.5;
    if (debug())
      cerr << "DEBUG: Astrobj::integrateEmission(): "
	   << "Icur=" << Icur << endl;
  } while( fabs(Icur-Iprev) > (1e-2 * Icur) );

  if (debug())
    cerr << "DEBUG: Astrobj::integrateEmission(): "
	 << "dnu=" << dnux2*0.5
	 << "=(nu2-nu1)/" << (nu2-nu1)/(dnux2*0.5) << endl;

  return Icur;
}

Quantity_t Astrobj::getDefaultQuantities() { return GYOTO_QUANTITY_INTENSITY; }


void Astrobj::initRegister() {
  if (Gyoto::Astrobj::Register_) delete Gyoto::Astrobj::Register_;
  Gyoto::Astrobj::Register_ = NULL;
}

void Gyoto::Astrobj::Register(std::string name, Subcontractor_t* scp){
  Register::Entry* ne =
    new Register::Entry(name, (SmartPointee::Subcontractor_t*)scp, Gyoto::Astrobj::Register_);
  Gyoto::Astrobj::Register_ = ne;
}

Gyoto::Astrobj::Subcontractor_t* Astrobj::getSubcontractor(std::string name) {
  if (!Gyoto::Astrobj::Register_) throwError("No Astrobj kind registered!");
  return (Subcontractor_t*)Gyoto::Astrobj::Register_
    -> getSubcontractor(name);
}


AstrobjProperties::AstrobjProperties() :
  intensity(NULL), time(NULL), distance(NULL),
  first_dmin(NULL), first_dmin_found(0),
  redshift(NULL), rimpact(NULL),
  spectrum(NULL), binspectrum(NULL), offset(1), x(NULL), y(NULL), z(NULL),
  user1(NULL), user2(NULL), user3(NULL), user4(NULL), user5(NULL)
{}

AstrobjProperties::AstrobjProperties( double * I, double * t) :
  intensity(I), time(t), distance(NULL),
  first_dmin(NULL), first_dmin_found(0),
  redshift(NULL), rimpact(NULL),
  spectrum(NULL), binspectrum(NULL), offset(1), x(NULL), y(NULL), z(NULL),
  user1(NULL), user2(NULL), user3(NULL), user4(NULL), user5(NULL)
{}

void AstrobjProperties::init(size_t nbnuobs) {
  if (intensity)  *intensity  = 0.;
  if (time)       *time       = DBL_MAX;
  if (distance)   *distance   = DBL_MAX;
  if (first_dmin){*first_dmin = DBL_MAX; first_dmin_found=0;}
  if (redshift)   *redshift   = 0.;
  if (rimpact)    *rimpact    = 0.;
  if (spectrum) for (size_t ii=0; ii<nbnuobs; ++ii) spectrum[ii*offset]=0.; 
  if (binspectrum) for (size_t ii=0; ii<nbnuobs; ++ii)
		     binspectrum[ii*offset]=0.; 
  if (x)          *x=0.;
  if (y)          *y=0.;
  if (z)          *z=0.;
  if (user1)      *user1=0.;
  if (user2)      *user2=0.;
  if (user3)      *user3=0.;
  if (user4)      *user4=0.;
  if (user5)      *user5=0.;
}

AstrobjProperties AstrobjProperties::operator++() {
  if (intensity)  ++intensity;
  if (time)       ++time;
  if (distance)   ++distance;
  if (first_dmin) ++first_dmin;
  if (redshift)   ++redshift;
  if (rimpact)    ++rimpact;
  if (spectrum)   ++spectrum;
  if (binspectrum)++binspectrum;
  if (x)          ++x;
  if (y)          ++y;
  if (z)          ++z;
  if (user1)      ++user1;
  if (user2)      ++user2;
  if (user3)      ++user3;
  if (user4)      ++user4;
  if (user5)      ++user5;
  return *this;
}
