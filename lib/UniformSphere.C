/*
    Copyright 2011 Thibaut Paumard, Frederic Vincent

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

#include "GyotoUtils.h"
#include "GyotoUniformSphere.h"
#include "GyotoPhoton.h"
#include "GyotoPowerLawSpectrum.h"
#include "GyotoBlackBodySpectrum.h"
#include "GyotoFactoryMessenger.h"
#include "GyotoConverters.h"
#include "GyotoProperty.h"

#include <iostream>
#include <cmath>
#include <string>
#include <cstdlib>
#include <float.h>
#include <sstream>
#include <string.h>

using namespace std;
using namespace Gyoto;
using namespace Gyoto::Astrobj;

GYOTO_PROPERTY_DOUBLE(UniformSphere, Radius, radius, Standard::properties);
GYOTO_PROPERTY_DOUBLE(UniformSphere, Alpha, alpha, &Radius);
GYOTO_PROPERTY_DOUBLE(UniformSphere,
		      DeltaMaxOverRadius, deltaMaxOverRadius, &Alpha);
GYOTO_PROPERTY_DOUBLE(UniformSphere,
		      DeltaMaxOverDistance, deltaMaxOverDistance,
		      &DeltaMaxOverRadius);
GYOTO_PROPERTY_BOOL(UniformSphere,
		    IsotropicEmittedIntensity,
		    TrueEmittedIntensity,
		    isotropic, &DeltaMaxOverDistance);
GYOTO_PROPERTY_SPECTRUM(UniformSphere,
			Opacity, opacity, &IsotropicEmittedIntensity);
GYOTO_PROPERTY_SPECTRUM(UniformSphere,
			Spectrum, spectrum, &Opacity);
GYOTO_PROPERTY_FINALIZE(UniformSphere, &::Spectrum);


#define GYOTO_USPH_DELTAMAX_OVER_RAD 0.1
#define GYOTO_USPH_DELTAMAX_OVER_DST 0.1

UniformSphere::UniformSphere(string kin) :
  Astrobj::Standard(kin),
  // radius_(0.),
  isotropic_(0),
  alpha_(1),
  spectrum_(NULL),
  opacity_(NULL),
  dltmor_(GYOTO_USPH_DELTAMAX_OVER_RAD),
  dltmod_(GYOTO_USPH_DELTAMAX_OVER_DST)
{
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << endl;
# endif

  // also initial safety_value_ etc.
  radius(0.);

  spectrum(new Spectrum::BlackBody()); 
  opacity(new Spectrum::PowerLaw(0., 1.)); 
  opticallyThin(false);
}

UniformSphere::UniformSphere(string kin,
			     SmartPointer<Metric::Generic> met, double rad) :
  Astrobj::Standard(kin),
  //radius_(rad),
  isotropic_(0),
  alpha_(1),
  spectrum_(NULL), opacity_(NULL),
  dltmor_(GYOTO_USPH_DELTAMAX_OVER_RAD),
  dltmod_(GYOTO_USPH_DELTAMAX_OVER_DST)
{
  // also initialize safety_value_ etc.
  radius(rad);

  spectrum(new Spectrum::BlackBody()); 
  opacity(new Spectrum::PowerLaw(0., 1.)); 
  opticallyThin(false);
  gg_=met;

}

UniformSphere::UniformSphere(const UniformSphere& orig) :
  Astrobj::Standard(orig),
  radius_(orig.radius_),
  isotropic_(orig.isotropic_),
  alpha_(orig.alpha_),
  spectrum_(NULL), opacity_(NULL),
  dltmor_(orig.dltmor_),
  dltmod_(orig.dltmod_)

{
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << endl;
# endif
  if (orig.spectrum_()) spectrum_=orig.spectrum_->clone();
  if (orig.opacity_()) opacity_=orig.opacity_->clone();
}

UniformSphere::~UniformSphere() {
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << endl;
# endif
}

string UniformSphere::className() const { return  string("UniformSphere"); }
string UniformSphere::className_l() const { return  string("uniformsphere"); }

SmartPointer<Spectrum::Generic> UniformSphere::spectrum() const { return spectrum_; }
void UniformSphere::spectrum(SmartPointer<Spectrum::Generic> sp) {spectrum_=sp;}

SmartPointer<Spectrum::Generic> UniformSphere::opacity() const { return opacity_; }
void UniformSphere::opacity(SmartPointer<Spectrum::Generic> sp) {
  opticallyThin(sp);
  opacity_=sp;
}


double UniformSphere::operator()(double const coord[4]) {
  double coord_st[4] = {coord[0]};
  double coord_ph[4] = {coord[0]};
  double sintheta;
  getCartesian(coord_st, 1, coord_st+1, coord_st+2, coord_st+3);
  switch (gg_->coordKind()) {
  case GYOTO_COORDKIND_CARTESIAN:
    memcpy(coord_ph+1, coord+1, 3*sizeof(double));
    break;
  case GYOTO_COORDKIND_SPHERICAL:
    coord_ph[1] = coord[1] * (sintheta=sin(coord[2])) * cos(coord[3]);
    coord_ph[2] = coord[1] * sintheta * sin(coord[3]);
    coord_ph[3] = coord[1] * cos(coord[2]) ;
    break;
  default:
    throwError("unsupported coordkind");
  }
  double dx = coord_ph[1]-coord_st[1];
  double dy = coord_ph[2]-coord_st[2];
  double dz = coord_ph[3]-coord_st[3];

  return dx*dx + dy*dy + dz*dz;
}

double UniformSphere::deltaMax(double * coord) {
  return max(dltmod_*sqrt((*this)(coord)), dltmor_*radius_);
}

double UniformSphere::emission(double nu_em, double dsem, double *, double *) const {
  if (isotropic_){
    if (flag_radtransf_){
      return dsem;
    }else{
      return 1.;
    }
  }
  if (flag_radtransf_) return (*spectrum_)(nu_em, (*opacity_)(nu_em), dsem);
  return (*spectrum_)(nu_em);
}

void UniformSphere::processHitQuantities(Photon* ph, double* coord_ph_hit,
					 double* coord_obj_hit, double dt,
					 Properties* data) const {
  if (alpha_==1) {
    // then I_nu \propto nu^0, standard case
    Generic::processHitQuantities(ph,coord_ph_hit,coord_obj_hit,dt,data);
    return;
  }

  // Here nu*I_nu \propto nu^alpha, alpha!=1
  // Emission is assumed to deliver
  // then I_nu integrated over a band is \propto g^(4-alpha_)
  // not simply g^3 as in the standard case 
  double freqObs=ph->freqObs(); // this is a useless quantity, always 1
  SmartPointer<Spectrometer::Generic> spr = ph -> spectrometer();
  double dlambda = dt/coord_ph_hit[4]; //dlambda = dt/tdot
  double ggredm1 = -gg_->ScalarProd(coord_ph_hit,coord_obj_hit+4,
				    coord_ph_hit+4);// / 1.; 
  //this is nu_em/nu_obs
  double ggred = 1./ggredm1;           //this is nu_obs/nu_em
  double dsem = dlambda*ggredm1; // *1.
  double inc =0.;
  if (data) {
    if (data->redshift) throwError("unimplemented");
    if (data->time) throwError("unimplemented");
    if (data->impactcoords) throwError("unimplemented");
    if (data->user4) throwError("unimplemented");
    if (data->binspectrum) throwError("unimplemented");
    if (data->spectrum) throwError("unimplemented");
    if (data->intensity) {
      //Intensity increment :
      inc = (emission(freqObs*ggredm1, dsem, coord_ph_hit, coord_obj_hit))
	* (ph -> getTransmission(size_t(-1)))
	* pow(ggred,4-alpha_);
      *data->intensity += inc;
    }
    
    /* update photon's transmission */
    ph -> transmit(size_t(-1),
		   transmission(freqObs*ggredm1, dsem,coord_ph_hit));
  } else {
#   if GYOTO_DEBUG_ENABLED
    GYOTO_DEBUG << "NO data requested!" << endl;
#   endif
  }
}  
      
double UniformSphere::transmission(double nuem, double dsem, double*) const {
  if (!flag_radtransf_) return 0.;
  double opac = (*opacity_)(nuem);
  
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG <<  "(nuem="    << nuem
	      << ", dsem="    << dsem
	      << "), opacity=" << opac << endl;
# endif

  if (!opac) return 1.;
  return exp(-opac*dsem);
}

double UniformSphere::integrateEmission(double nu1, double nu2, double dsem,
			       double *, double *) const {
  if (flag_radtransf_)
    return spectrum_->integrate(nu1, nu2, opacity_(), dsem);
  return spectrum_->integrate(nu1, nu2);
}


double UniformSphere::radius() const {
  return radius_;
}

void UniformSphere::radius(double r) {
  radius_=r;
  critical_value_ = r*r;
  safety_value_ = critical_value_*1.1+0.1;
}

double UniformSphere::radius(std::string unit) const {
  return Units::FromGeometrical(radius(), unit, gg_);
}

void UniformSphere::radius(double r, std::string unit) {
  radius(Units::ToGeometrical(r, unit, gg_));
}

double UniformSphere::deltaMaxOverRadius() const {return dltmor_;}
void UniformSphere::deltaMaxOverRadius(double f) {dltmor_=f;}

double UniformSphere::deltaMaxOverDistance() const {return dltmod_;}
void UniformSphere::deltaMaxOverDistance(double f) {dltmod_=f;}

double UniformSphere::alpha() const { return alpha_; }
void UniformSphere::alpha(double a) { alpha_ = a; }

bool UniformSphere::isotropic() const { return isotropic_; }
void UniformSphere::isotropic(bool a) { isotropic_ = a; }
