/*
    Copyright 2011-2015, 2018-2019 Thibaut Paumard, Frederic Vincent

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

GYOTO_PROPERTY_START(Gyoto::Astrobj::UniformSphere,
 "Coordinate sphere with uniform emission and absorption.")
GYOTO_PROPERTY_SPECTRUM(UniformSphere, Spectrum, spectrum,
 "Emission law.")
GYOTO_PROPERTY_SPECTRUM(UniformSphere,Opacity, opacity,
 "Absorption law.")
GYOTO_PROPERTY_BOOL(UniformSphere,
		    IsotropicEmittedIntensity, TrueEmittedIntensity,
		    isotropic,
 "If true, then emission just returns 1.")
GYOTO_PROPERTY_DOUBLE(UniformSphere,
		      DeltaMaxOverDistance, deltaMaxOverDistance,
 "Maximum value of step/distance from centre of sphere for photons.")
GYOTO_PROPERTY_DOUBLE(UniformSphere,
		      DeltaMaxOverRadius, deltaMaxOverRadius,
 "Maximum value of step/radius of sphere for photons.")
GYOTO_PROPERTY_DOUBLE(UniformSphere, Alpha, alpha, "Deprecated")
GYOTO_PROPERTY_DOUBLE_UNIT(UniformSphere, Radius, radius, "Sphere radius (geometrical units).")
GYOTO_PROPERTY_END(UniformSphere, Standard::properties)


#define GYOTO_USPH_DELTAMAX_OVER_RAD 0.1
#define GYOTO_USPH_DELTAMAX_OVER_DST 0.1

UniformSphere::UniformSphere(string kin) :
  Astrobj::Standard(kin),
  // radius_(0.),
  isotropic_(0),
  spectrum_(NULL),
  opacity_(NULL),
  dltmor_(GYOTO_USPH_DELTAMAX_OVER_RAD),
  dltmod_(GYOTO_USPH_DELTAMAX_OVER_DST)
{
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << endl;
# endif
  //cout << "in creaor " << Generic::radiativeQ() << endl;
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
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << endl;
# endif
  double coord_st[4] = {coord[0]};
  double coord_ph[4] = {coord[0]};
  double sintheta;
  getCartesian(coord_st, 1, coord_st+1, coord_st+2, coord_st+3);

  // Special treatment for SchwarzschildHarmonic: define star
  // as a sphere of radius r_BL=R_star and not r_harmonic=R_star,
  // in order to ease comparison between coordinate systems.
  if (gg_->kind()=="SchwarzschildHarmonic"){
    double r_st = sqrt(coord_st[1]*coord_st[1]+coord_st[2]*coord_st[2]+coord_st[3]*coord_st[3]);
    double theta = acos(coord_st[3]/r_st), phi = atan(coord_st[2]/coord_st[1]);
    coord_st[1]+= sin(theta)*cos(phi);
    coord_st[2]+= sin(theta)*sin(phi);
    coord_st[3]+= cos(theta);
  }
  switch (gg_->coordKind()) {
  case GYOTO_COORDKIND_CARTESIAN:
    memcpy(coord_ph+1, coord+1, 3*sizeof(double));
    break;
  case GYOTO_COORDKIND_SPHERICAL:
    coord_ph[1] = (coord[1]) * (sintheta=sin(coord[2])) * cos(coord[3]);
    coord_ph[2] = (coord[1]) * sintheta * sin(coord[3]);
    coord_ph[3] = (coord[1]) * cos(coord[2]) ;
    
    // Special treatment for SchwarzschildHarmonic: define star
    // as a sphere of radius r_BL=R_star and not r_harmonic=R_star,
    // in order to ease comparison between coordinate systems.
    if (gg_->kind()=="SchwarzschildHarmonic"){
      coord_ph[1] = (coord[1]+1.) * sintheta * cos(coord[3]);
      coord_ph[2] = (coord[1]+1.) * sintheta * sin(coord[3]);
      coord_ph[3] = (coord[1]+1.) * cos(coord[2]) ;
    }
    break;
  default:
    GYOTO_ERROR("unsupported coordkind");
  }
  //cout << "testcoord: " << coord_ph[1] << " " << coord_st[1] << " " << coord_ph[1] - coord_st[1] <<  endl;
  double dx = coord_ph[1]-coord_st[1];
  double dy = coord_ph[2]-coord_st[2];
  double dz = coord_ph[3]-coord_st[3];
  //cout << "unif= " << dx*dx << " " << dy*dy << " " << dz*dz << endl;

  //double rstar = sqrt(coord_st[1]*coord_st[1] + coord_st[2]*coord_st[2] +coord_st[3]*coord_st[3]);
  //cout << "tph, rph, thph, phph= " << coord[0]<< " " << coord[1] << " " << coord[2] << " " << coord[3] << " " << endl;
  //cout << "tst, rst, thst, phst= " << coord_st[0]<< " " << r_st << " " << theta << " " << phi << endl;
  //cout << "d2 Rstar= " << dx*dx + dy*dy + dz*dz << " " << radius_*radius_<< endl;
  //cout << "trthph ph + st, rsp= " << coord[0] << " " << coord[1] << " " << coord[2] << " " << coord[3] << " ; " << coord_st[0] << " " << rstar << " " << acos(coord_st[3]/rstar) << " " << atan(coord_st[2]/coord_st[1]) << " " << dx*dx + dy*dy + dz*dz << endl;
  return dx*dx + dy*dy + dz*dz;
}

double UniformSphere::deltaMax(double * coord) {
  double r;
  switch (gg_->coordKind()) {
  case GYOTO_COORDKIND_CARTESIAN:
    r=sqrt(coord[1]*coord[1]+coord[2]*coord[2]+coord[3]*coord[3]);
    break;
  case GYOTO_COORDKIND_SPHERICAL:
    r=coord[1];
    break;
  default:
    r=0.;
    GYOTO_ERROR("unsupported coordkind");
  }
  if (rmax_!=DBL_MAX && r>rmax_) return r*0.5; 
  return max(dltmod_*sqrt((*this)(coord)), dltmor_*radius_);
}

double UniformSphere::emission(double nu_em, double dsem, state_t const &, double const *) const {
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << endl;
# endif
  if (isotropic_){
    if (flag_radtransf_){
      return dsem;
    }else{
      //cout << "returning 1 in unif sph" << endl;
      return 1.;
    }
  }
  if (flag_radtransf_) return (*spectrum_)(nu_em, (*opacity_)(nu_em), dsem);
  return (*spectrum_)(nu_em);
}

double UniformSphere::transmission(double nuem, double dsem, state_t const &, double const *) const {
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << endl;
# endif
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
					state_t const &, double const *) const {
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << endl;
# endif
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

double UniformSphere::radius(std::string const &unit) const {
  return Units::FromGeometrical(radius(), unit, gg_);
}

void UniformSphere::radius(double r, std::string const &unit) {
  radius(Units::ToGeometrical(r, unit, gg_));
}

double UniformSphere::deltaMaxOverRadius() const {return dltmor_;}
void UniformSphere::deltaMaxOverRadius(double f) {dltmor_=f;}

double UniformSphere::deltaMaxOverDistance() const {return dltmod_;}
void UniformSphere::deltaMaxOverDistance(double f) {dltmod_=f;}

double UniformSphere::alpha() const { return 1.; }
void UniformSphere::alpha(double a) {
  if (a != 1.) GYOTO_ERROR("property 'Alpha' is deprecated");
}

bool UniformSphere::isotropic() const { return isotropic_; }
void UniformSphere::isotropic(bool a) { isotropic_ = a; }
