/*
    Copyright 2011, 2018 Thibaut Paumard, Frederic Vincent

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
GYOTO_PROPERTY_DOUBLE(UniformSphere, Alpha, alpha)
GYOTO_PROPERTY_DOUBLE_UNIT(UniformSphere, Radius, radius, "Sphere radius (geometrical units).")
GYOTO_PROPERTY_DOUBLE_UNIT(UniformSphere, NumberDensity, numberDensity)
GYOTO_PROPERTY_DOUBLE(UniformSphere, Temperature, temperature)
GYOTO_PROPERTY_DOUBLE_UNIT(UniformSphere, TimeRef, timeRef)
GYOTO_PROPERTY_DOUBLE_UNIT(UniformSphere, TimeSigma, timeSigma)
GYOTO_PROPERTY_DOUBLE(UniformSphere, MagneticParticlesEquipartitionRatio,
		      magneticParticlesEquipartitionRatio)
GYOTO_PROPERTY_DOUBLE(UniformSphere, KappaIndex, kappaIndex)
GYOTO_PROPERTY_END(UniformSphere, Standard::properties)


#define GYOTO_USPH_DELTAMAX_OVER_RAD 0.1
#define GYOTO_USPH_DELTAMAX_OVER_DST 0.1

UniformSphere::UniformSphere(string kin) :
  Astrobj::Standard(kin),
  // radius_(0.),
  isotropic_(0),
  alpha_(1),
  numberDensity_cgs_(1.),
  temperature_(1.),
  timeRef_M_(1.),
  timeSigma_M_(1.),
  magneticParticlesEquipartitionRatio_(1.),
  kappaIndex_(1.),
  spectrum_(NULL),
  opacity_(NULL),
  spectrumKappaSynch_(NULL),
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
  spectrumKappaSynch_ = new Spectrum::KappaDistributionSynchrotron();
  opticallyThin(false);
}

UniformSphere::UniformSphere(string kin,
			     SmartPointer<Metric::Generic> met, double rad) :
  Astrobj::Standard(kin),
  //radius_(rad),
  isotropic_(0),
  alpha_(1),
  numberDensity_cgs_(1.),
  temperature_(1.),
  timeRef_M_(1.),
  timeSigma_M_(1.),
  kappaIndex_(1.),
  magneticParticlesEquipartitionRatio_(1.),
  spectrum_(NULL), opacity_(NULL), spectrumKappaSynch_(NULL),
  dltmor_(GYOTO_USPH_DELTAMAX_OVER_RAD),
  dltmod_(GYOTO_USPH_DELTAMAX_OVER_DST)
{
  // also initialize safety_value_ etc.
  radius(rad);

  spectrum(new Spectrum::BlackBody()); 
  opacity(new Spectrum::PowerLaw(0., 1.));
  spectrumKappaSynch_ = new Spectrum::KappaDistributionSynchrotron();
  opticallyThin(false);
  gg_=met;

}

UniformSphere::UniformSphere(const UniformSphere& orig) :
  Astrobj::Standard(orig),
  radius_(orig.radius_),
  isotropic_(orig.isotropic_),
  alpha_(orig.alpha_),
  numberDensity_cgs_(orig.numberDensity_cgs_),
  temperature_(orig.temperature_),
  timeRef_M_(orig.timeRef_M_),
  timeSigma_M_(orig.timeSigma_M_),
  kappaIndex_(orig.kappaIndex_),
  magneticParticlesEquipartitionRatio_(orig.magneticParticlesEquipartitionRatio_),
  spectrum_(NULL), opacity_(NULL), spectrumKappaSynch_(NULL),
  dltmor_(orig.dltmor_),
  dltmod_(orig.dltmod_)

{
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << endl;
# endif
  if (orig.spectrum_()) spectrum_=orig.spectrum_->clone();
  if (orig.opacity_()) opacity_=orig.opacity_->clone();
  if (orig.spectrumKappaSynch_()) spectrumKappaSynch_=orig.spectrumKappaSynch_->clone();
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
  //cout << "rsp= " << sqrt(coord_st[1]*coord_st[1]+coord_st[2]*coord_st[2]+coord_st[3]*coord_st[3]) << endl;
  double dx = coord_ph[1]-coord_st[1];
  double dy = coord_ph[2]-coord_st[2];
  double dz = coord_ph[3]-coord_st[3];

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
    throwError("unsupported coordkind");
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
      return 1.;
    }
  }
  if (flag_radtransf_) return (*spectrum_)(nu_em, (*opacity_)(nu_em), dsem);
  return (*spectrum_)(nu_em);
}

void UniformSphere::radiativeQ(double Inu[], // output
			       double Taunu[], // output
			       double const nu_ems[], size_t nbnu, // input
			       double dsem,
			       state_t const &coord_ph,
			       double const coord_obj[8]) const {
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << endl;
# endif
  double tcur=coord_ph[0]; //*GMoc3/60.; // in min
  double modulation = exp(-pow((tcur-timeRef_M_)/timeSigma_M_,2));
  double temperature = modulation*temperature_,
    number_density = modulation*numberDensity_cgs_;
  cout << "spot tcur, time_ref, time_sigma, modulation, number_density=" << tcur << " " << timeRef_M_ << " " << timeSigma_M_ << " " << modulation << " " << numberDensity_cgs_ << " " << temperature_ << " " << number_density << " " << temperature << " " << kappaIndex_ << " " << magneticParticlesEquipartitionRatio_ << endl;
  double thetae = GYOTO_BOLTZMANN_CGS*temperature
    /(GYOTO_ELECTRON_MASS_CGS*GYOTO_C2_CGS);
  
  double hypergeom = Gyoto::hypergeom(kappaIndex_, thetae);
  
  double BB = sqrt(8.*M_PI*magneticParticlesEquipartitionRatio_
		   *GYOTO_PROTON_MASS_CGS * GYOTO_C_CGS * GYOTO_C_CGS
		   *number_density);
  double nu0 = GYOTO_ELEMENTARY_CHARGE_CGS*BB
    /(2.*M_PI*GYOTO_ELECTRON_MASS_CGS*GYOTO_C_CGS); // cyclotron freq
  
  // Defining jnus, anus
  double jnu_synch_kappa[nbnu], anu_synch_kappa[nbnu];
  
  for (size_t ii=0; ii<nbnu; ++ii){
    // Initializing to <0 value to create errors if not updated
    // [ exp(-anu*ds) will explose ]
    jnu_synch_kappa[ii]=-1.;
    anu_synch_kappa[ii]=-1.;
  }
  
  // THERMAL SYNCHRO
  spectrumKappaSynch_->kappaindex(kappaIndex_);
  spectrumKappaSynch_->numberdensityCGS(number_density);
  spectrumKappaSynch_->angle_averaged(1);
  spectrumKappaSynch_->angle_B_pem(0.); // avg so we don't care
  spectrumKappaSynch_->cyclotron_freq(nu0);
  spectrumKappaSynch_->thetae(thetae);
  spectrumKappaSynch_->hypergeometric(hypergeom);

  spectrumKappaSynch_->radiativeQ(jnu_synch_kappa,anu_synch_kappa,
				  nu_ems,nbnu);

  // RETURNING TOTAL INTENSITY AND TRANSMISSION
  for (size_t ii=0; ii<nbnu; ++ii){
    double jnu_tot = jnu_synch_kappa[ii],
      anu_tot = anu_synch_kappa[ii];
    //cout << "At r,th= " << coord_ph[1] << " " << coord_ph[2] << endl;
    //cout << "in unif stuff: " << number_density << " " << nu0 << " " << thetae << " " << hypergeom << " " << jnu_tot << " " << anu_tot << " " << dsem << endl;

    // expm1 is a precise implementation of exp(x)-1
    double em1=std::expm1(-anu_tot * dsem * gg_->unitLength());
    Taunu[ii] = em1+1.;
    Inu[ii] = anu_tot == 0. ? jnu_tot * dsem * gg_->unitLength() :
      -jnu_tot / anu_tot * em1;
    
    if (Inu[ii]<0.)
      throwError("In UniformSphere::radiativeQ: Inu<0");
    if (Inu[ii]!=Inu[ii] or Taunu[ii]!=Taunu[ii])
      throwError("In UniformSphere::radiativeQ: Inu or Taunu is nan");
    if (Inu[ii]==Inu[ii]+1. or Taunu[ii]==Taunu[ii]+1.)
      throwError("In UniformSphere::radiativeQ: Inu or Taunu is infinite");
    
  }

}

void UniformSphere::processHitQuantities(Photon* ph, state_t const &coord_ph_hit,
					 double const coord_obj_hit[8], double dt,
					 Properties* data) const {
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << endl;
# endif
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
  double ggredm1 = -gg_->ScalarProd(&coord_ph_hit[0],coord_obj_hit+4,
				    &coord_ph_hit[4]);// / 1.; 
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
      
double UniformSphere::transmission(double nuem, double dsem, state_t const &) const {
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

double UniformSphere::alpha() const { return alpha_; }
void UniformSphere::alpha(double a) { alpha_ = a; }

double UniformSphere::numberDensity() const {
  // Converts internal cgs central enthalpy to SI
  double dens=numberDensity_cgs_;
# ifdef HAVE_UDUNITS
  dens = Units::Converter("cm-3", "m-3")(dens);
# else
  GYOTO_WARNING << "Units ignored, please recompile Gyoto with --with-udunits"
		<< endl ;
# endif
  return dens; }
double UniformSphere::numberDensity(string const &unit) const
{
  double dens = numberDensity();
  if (unit != "") {
# ifdef HAVE_UDUNITS
    dens = Units::Converter("m-3", unit)(dens);
# else
    GYOTO_WARNING << "Units ignored, please recompile Gyoto with --with-udunits"
		  << endl ;
# endif
  }
  return dens;
}
void UniformSphere::numberDensity(double dens) {
# ifdef HAVE_UDUNITS
  dens = Units::Converter("m-3", "cm-3")(dens);
# else
  GYOTO_WARNING << "Units ignored, please recompile Gyoto with --with-udunits"
		<< endl ;
# endif
  numberDensity_cgs_=dens;
}
void UniformSphere::numberDensity(double dens, string const &unit) {
  if (unit != "") {
# ifdef HAVE_UDUNITS
    dens = Units::Converter(unit, "m-3")(dens);
# else
    GYOTO_WARNING << "Units ignored, please recompile Gyoto with --with-udunits"
		  << endl ;
# endif
  }
  numberDensity(dens);
}

double UniformSphere::temperature() const { return temperature_; }
void UniformSphere::temperature(double tt) { temperature_ = tt; }

double UniformSphere::timeRef() const {
  // Converts internal M-unit time to SI
  double tt=timeRef_M_;
# ifdef HAVE_UDUNITS
  tt = Units::ToSeconds(tt,"geometrical_time",gg_);
# else
  GYOTO_WARNING << "Units ignored, please recompile Gyoto with --with-udunits"
		<< endl ;
# endif
  return tt; }
double UniformSphere::timeRef(string const &unit) const
{
  double tt = timeRef();
  if (unit != "") {
# ifdef HAVE_UDUNITS
    tt = Units::Converter("s", unit)(tt);
# else
    GYOTO_WARNING << "Units ignored, please recompile Gyoto with --with-udunits"
		  << endl ;
# endif
  }
  return tt;
}
void UniformSphere::timeRef(double tt) {
# ifdef HAVE_UDUNITS
  tt = Units::ToGeometricalTime(tt, "s", gg_);
# else
  GYOTO_WARNING << "Units ignored, please recompile Gyoto with --with-udunits"
		<< endl ;
# endif
  timeRef_M_ = tt; }
void UniformSphere::timeRef(double tt, string const &unit) {
  if (unit != "") {
# ifdef HAVE_UDUNITS
    tt = Units::ToSeconds(tt,unit,gg_);
# else
    GYOTO_WARNING << "Units ignored, please recompile Gyoto with --with-udunits"
		  << endl ;
# endif
  }
  timeRef(tt);
}

double UniformSphere::timeSigma() const {
  // Converts internal M-unit time to SI
  double tt=timeSigma_M_;
# ifdef HAVE_UDUNITS
  tt = Units::ToSeconds(tt,"geometrical_time",gg_);
# else
  GYOTO_WARNING << "Units ignored, please recompile Gyoto with --with-udunits"
		<< endl ;
# endif
  return tt; }
double UniformSphere::timeSigma(string const &unit) const
{
  double tt = timeSigma();
  if (unit != "") {
# ifdef HAVE_UDUNITS
    tt = Units::Converter("s", unit)(tt);
# else
    GYOTO_WARNING << "Units ignored, please recompile Gyoto with --with-udunits"
		  << endl ;
# endif
  }
  return tt;
}
void UniformSphere::timeSigma(double tt) {
# ifdef HAVE_UDUNITS
  tt = Units::ToGeometricalTime(tt, "s", gg_);
# else
  GYOTO_WARNING << "Units ignored, please recompile Gyoto with --with-udunits"
		<< endl ;
# endif
  timeSigma_M_ = tt; }
void UniformSphere::timeSigma(double tt, string const &unit) {
  if (unit != "") {
# ifdef HAVE_UDUNITS
    tt = Units::ToSeconds(tt,unit,gg_);
# else
    GYOTO_WARNING << "Units ignored, please recompile Gyoto with --with-udunits"
		  << endl ;
# endif
  }
  timeSigma(tt);
}

void UniformSphere::magneticParticlesEquipartitionRatio(double rr) {
  magneticParticlesEquipartitionRatio_=rr;}
double UniformSphere::magneticParticlesEquipartitionRatio()const{
  return magneticParticlesEquipartitionRatio_;}

double UniformSphere::kappaIndex() const { return kappaIndex_; }
void UniformSphere::kappaIndex(double ind) { kappaIndex_ = ind; }

bool UniformSphere::isotropic() const { return isotropic_; }
void UniformSphere::isotropic(bool a) { isotropic_ = a; }
