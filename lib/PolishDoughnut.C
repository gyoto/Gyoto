/*
    Copyright (c) 2012-2016, 2018-2019 Frederic Vincent, Odele Straub,
                                       Thibaut Paumard

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
#include "GyotoPolishDoughnut.h"
#include "GyotoProperty.h"
#include "GyotoPhoton.h"
#include "GyotoFactoryMessenger.h"
#include "GyotoDefs.h"
#include <iostream>
#include <cmath>
#include <string>
#include <cstdlib>
#include <sstream>
using namespace std;
using namespace Gyoto;
using namespace Gyoto::Astrobj;

GYOTO_PROPERTY_START(PolishDoughnut)
GYOTO_PROPERTY_DOUBLE(PolishDoughnut, Lambda, lambda)
GYOTO_PROPERTY_VECTOR_DOUBLE(PolishDoughnut, AngMomRinner, angmomrinner)
GYOTO_PROPERTY_DOUBLE_UNIT(PolishDoughnut, CentralEnthalpyPerUnitVolume, centralEnthalpyPerUnitVolume)
GYOTO_PROPERTY_DOUBLE(PolishDoughnut,
		      CentralTemperature, centralTemp)
GYOTO_PROPERTY_DOUBLE(PolishDoughnut, Beta, beta,
		      "one parametrization of the magnetic to particle "
		      "energy density ratio; this is not the standard "
		      "plasma beta")
GYOTO_PROPERTY_DOUBLE(PolishDoughnut, MagnetizationParameter,
		      magnetizationParameter,
		      "another parametrization of the magnetic to particle "
		      "energy density ratio; this is the standard "
		      "magnetization parameter; this is not the standard "
		      "plasma beta")
GYOTO_PROPERTY_SIZE_T(PolishDoughnut,
		      SpectralOversampling, spectralOversampling)
GYOTO_PROPERTY_BOOL(PolishDoughnut,
		    AngleAveraged, NoAngleAveraged,
		    angleAveraged)
GYOTO_PROPERTY_BOOL(PolishDoughnut,
		    Bremsstrahlung, NoBremsstrahlung,
		    bremsstrahlung)
GYOTO_PROPERTY_VECTOR_DOUBLE(PolishDoughnut,
			     NonThermalDeltaExpo, nonThermalDeltaExpo)
// Since adafparams(vector) sets adaf_ to true, ADAF must come after
// ADAFParameters
GYOTO_PROPERTY_VECTOR_DOUBLE(PolishDoughnut, ADAFParameters, adafparams)
GYOTO_PROPERTY_BOOL(PolishDoughnut, ADAF, NonADAF, adaf)
GYOTO_PROPERTY_BOOL(PolishDoughnut,
		    ChangeCusp, KeepCusp, changeCusp)
GYOTO_PROPERTY_END(PolishDoughnut, Standard::properties)

#ifdef GYOTO_USE_XERCES
/*
  Either lambda_ or rintorus_ is defined. Filter out the other one
  when writing properties to XML.
 */
void PolishDoughnut::fillProperty(Gyoto::FactoryMessenger *fmp, Property const &p) const {
  if ((p.name == "Lambda" && !rochelobefilling_) ||
      (p.name == "AngMomRinner" && !defangmomrinner_))
    return; // do nothing
  else
    Standard::fillProperty(fmp, p);
}
#endif

#define CST_POLY_INDEX 1.5//polytropic index n (gamma=1+1/n=5/3)
#define CST_POLY_INDEX_M1 0.666666666666666666666666666666666666666667
#define CST_HYDRO_FRAC 0.75//hydrogen fraction
#define CST_Z_1 1.//atomic number
#define CST_Z_2 2.
#define CST_MU_ION 1.2307692307692308375521861 //(4./(1. + 3. * CST_HYDRO_FRAC))
#define CST_MU_ELEC 1.1428571428571427937015414 //(2./(1. + CST_HYDRO_FRAC))
#define w_tol 1e-9
#define DEFAULT_L0 10.
#define DEFAULT_RIN 10.
PolishDoughnut::PolishDoughnut() :
  Standard("PolishDoughnut"),
  l0_(DEFAULT_L0),
  lambda_(0.5),
  W_surface_(0.),
  W_centre_(0.),
  r_cusp_(0.),
  r_centre_(0.),
  r_torusouter_(0.),
//DeltaWm1_(),
  central_enthalpy_cgs_(1.),
  central_temperature_(1e10),
  beta_(0.),
  magnetizationParameter_(-1.),
  spectral_oversampling_(10),
  angle_averaged_(0),
  bremsstrahlung_(false),
  deltaPL_(0.),
  adaf_(0),
  ADAFtemperature_(0.),
  ADAFdensity_(0.),
  changecusp_(0),
  rochelobefilling_(0),
  defangmomrinner_(0),
  rintorus_(DEFAULT_RIN),
  intersection(this)
{
#ifdef GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << endl;
#endif
  critical_value_=0.; safety_value_=.1; //rmax_=25.;
  spectrumBrems_ = new Spectrum::ThermalBremsstrahlung();
  spectrumSynch_ = new Spectrum::ThermalSynchrotron();
  spectrumPLSynch_ = new Spectrum::PowerLawSynchrotron();
}
PolishDoughnut::PolishDoughnut(const PolishDoughnut& orig) :
  Standard(orig),
  spectrumBrems_(NULL),
  spectrumSynch_(NULL),
  spectrumPLSynch_(NULL),
  l0_(orig.l0_),
  lambda_(orig.lambda_),
  W_surface_(orig.W_surface_),
  W_centre_(orig.W_centre_),
  r_cusp_(orig.r_cusp_),
  r_centre_(orig.r_centre_),
  r_torusouter_(orig.r_torusouter_),
  DeltaWm1_(orig.DeltaWm1_),
  central_enthalpy_cgs_(orig.central_enthalpy_cgs_),
  central_temperature_(orig.central_temperature_),
  beta_(orig.beta_),
  magnetizationParameter_(orig.magnetizationParameter_),
  aa_(orig.aa_),
  aa2_(orig.aa2_),
  spectral_oversampling_(orig.spectral_oversampling_),
  angle_averaged_(orig.angle_averaged_),
  bremsstrahlung_(orig.bremsstrahlung_),
  deltaPL_(orig.deltaPL_),
  adaf_(orig.adaf_),
  ADAFtemperature_(orig.ADAFtemperature_),
  ADAFdensity_(orig.ADAFdensity_),
  changecusp_(orig.changecusp_),
  rochelobefilling_(orig.rochelobefilling_),
  defangmomrinner_(orig.defangmomrinner_),
  rintorus_(orig.rintorus_),
  intersection(orig.intersection)
{
  intersection.papa=this;
  if (gg_) gg_ -> hook(this);
  if (orig.spectrumBrems_()) spectrumBrems_=orig.spectrumBrems_->clone();
  if (orig.spectrumSynch_()) spectrumSynch_=orig.spectrumSynch_->clone();
  if (orig.spectrumPLSynch_()) spectrumPLSynch_=orig.spectrumPLSynch_->clone();
}
PolishDoughnut* PolishDoughnut::clone() const
{return new PolishDoughnut(*this);}
bool PolishDoughnut::isThreadSafe() const {
  return Standard::isThreadSafe()
    && (!spectrumBrems_ || spectrumBrems_->isThreadSafe())
    && (!spectrumSynch_ || spectrumSynch_->isThreadSafe())
    && (!spectrumPLSynch_ || spectrumPLSynch_->isThreadSafe());
}
double PolishDoughnut::getL0() const { return l0_; }
//void PolishDoughnut::setL0(double l0) { l0_ = l0; }
double PolishDoughnut::getWsurface() const { return W_surface_; }
double PolishDoughnut::getWcentre() const { return W_centre_; }
double PolishDoughnut::getRcusp() const { return r_cusp_; }
double PolishDoughnut::getRcentre() const { return r_centre_; }

double PolishDoughnut::lambda() const {
  if (!rochelobefilling_) {
    if (defangmomrinner_)
      GYOTO_ERROR("Lambda is not set because AngMomRinner is.");
    else
      GYOTO_ERROR("Lambda is not set yet.");
  }
  return lambda_;
}
void PolishDoughnut::lambda(double lam) {
  rochelobefilling_=1; // if here, the torus fills its Roche lobe
  if (defangmomrinner_){
    GYOTO_WARNING << "Setting Lambda overrides AngMomRinner previously set\n";
    defangmomrinner_=0;
  }
  if (!gg_) GYOTO_ERROR("Metric but be set before lambda in PolishDoughnut");
  //Computing marginally stable and marginally bound radii:
  lambda_=lam;
  double rms = gg_->getRms();
  double rmb = gg_->getRmb();
  // marginally stable & marginally bound keplerian angular momentum
  // (Polish doughnut l is rescaled)
  double  l_ms = gg_->getSpecificAngularMomentum(rms);
  double  l_mb = gg_->getSpecificAngularMomentum(rmb);
  l0_ = lambda_*(l_mb-l_ms)+l_ms ;//torus angular momentum
  //Computing the potential at the photon position:
  double r1_min = rmb ;
  double r1_max = rms ;
  double r2_min = rms ;
  double r2_max = 1000. ;
  r_cusp_ = intersection.ridders(r1_min, r1_max) ;
  rintorus_ = r_cusp_;
  r_centre_ = intersection.ridders(r2_min, r2_max) ;
  double poss[4]={0.,r_cusp_,M_PI/2.,0.};
  double posc[4]={0.,r_centre_,M_PI/2.,0.};
  W_surface_ = gg_->getPotential(poss,l0_);
  W_centre_  = gg_->getPotential(posc,l0_);
  DeltaWm1_ = 1./(W_centre_ - W_surface_);

  if (changecusp_){
    /*
      CUSP PROBLEM
      For both large spin and large lambda, the potential line w=0
      can leak out of the tube r=rcusp in the funnel zone (i.e. in
      the zone which is not taken into account for a doughnut, and
      which should be considered as a target for photons). Thus, the
      surface of the torus can be badly defined for big a and lambda.
      Using mma, I checked that this problem only arises for a>0.8
      and lambda>0.3. In this range, rcusp is multiplied by 1.25
      to encompass the w=0 line in the funnel, and to be sure not to
      ray-trace any part of he funnel. This is a crude solution: it
      allows to get rid of the funnel, but it will also remove a small
      part of the proper torus for a>~0.8 and lambda>~0.3. However
      this "forgotten part of the torus" being small, I do not think
      this should be a serious problem. Moreover, this removal of
      part of the torus only affects very slender tori (typically
      few r_g wide).
      But it would be nice to solve this problem in a more elegant way...
    */
    r_cusp_*=1.25;
  }
  // Find torus outer radius
  double r3_min=r_centre_, r3_max=5000.;
  if (lambda_>0.99){
    GYOTO_ERROR("In PolishDoughnut: please use a value of"
	       " lambda < 0.99, or else the outer radius"
	       " finding algorithm may crash");
  }
  outerradius_t outerradius;
  outerradius.papa = this;
  r_torusouter_ = outerradius.ridders(r3_min,r3_max);
  GYOTO_IF_DEBUG;
  GYOTO_DEBUG_EXPR(r_cusp_);
  GYOTO_DEBUG_EXPR(r_torusouter_);
  GYOTO_ENDIF_DEBUG;
  if (r_torusouter_!=r_torusouter_ || r_torusouter_==r_torusouter_+1)
    GYOTO_ERROR("In PolishDoughnut::lambda(): bad r_torusouter_");
  GYOTO_IF_DEBUG
    GYOTO_DEBUG_EXPR(lambda_);
  GYOTO_DEBUG_EXPR(l0_);
  GYOTO_DEBUG_EXPR(r_cusp_);
  GYOTO_DEBUG_EXPR(r_centre_);
  GYOTO_DEBUG_EXPR(W_surface_);
  GYOTO_DEBUG_EXPR(W_centre_);
  GYOTO_ENDIF_DEBUG
    }

void PolishDoughnut::angmomrinner(std::vector<double> const &v) {
  defangmomrinner_=1;
  if (rochelobefilling_){
    GYOTO_WARNING << "Setting AngMomRinner overrides Lambda previously set\n";
    rochelobefilling_=0;
  }
  if (v.size() != 2)
    GYOTO_ERROR("Only 2 arguments to define l0 and rin");
  l0_ = v[0];
  rintorus_ = v[1];
  r_cusp_=rintorus_; // NB: the cusp is most probably not at this radius; however we need to define it to avoid having cases where operator() returns "inside torus" when the photon actually is in the funnel at r<r_cusp_. Defining r_cusp_ that way is okay, we won't miss any part of the real torus.
  //cout << "l0,rin= " << l0_ << " " << rintorus_ << endl;
  double posin[4]={0.,rintorus_,M_PI/2.,0.};
  W_surface_ = gg_->getPotential(posin,l0_);
  double rmin=rintorus_, rmax = 1000.;
  //cout << "rmin max= " << rmin << " " << rmax << endl;
  r_centre_ = intersection.ridders(rmin, rmax) ;
  //cout << "rmin center max= " << rmin << " " << r_centre_ << " " << rmax << endl;
  if (r_centre_ < rmin or r_centre_ > rmax)
    GYOTO_ERROR("In PolishDoughnut::angmomrinner: bad r_centre_");
  double posc[4]={0.,r_centre_,M_PI/2.,0.};
  W_centre_  = gg_->getPotential(posc,l0_);
  DeltaWm1_ = 1./(W_centre_ - W_surface_);
  //cout << "Ws Wc rc= " << W_surface_ << " " << W_centre_ << " " << r_centre_ << endl;
  outerradius_t outerradius;
  outerradius.papa = this;
  rmin=r_centre_;
  r_torusouter_ = outerradius.ridders(rmin,rmax);
  //cout << "Torus rinner, rcen, router= " << rintorus_ << " " << r_centre_ << " " << r_torusouter_ << endl;
  GYOTO_IF_DEBUG;
  GYOTO_DEBUG_EXPR(l0_);
  GYOTO_DEBUG_EXPR(r_centre_);
  GYOTO_DEBUG_EXPR(rintorus_);
  GYOTO_DEBUG_EXPR(W_surface_);
  GYOTO_DEBUG_EXPR(W_centre_);
  GYOTO_ENDIF_DEBUG
}
std::vector<double> PolishDoughnut::angmomrinner() const {
  if (!defangmomrinner_) {
    if (rochelobefilling_)
      GYOTO_ERROR("AngMomRinner is not set because Lambda has been set.");
    else
      GYOTO_ERROR("AngMomRinner is not set yet.");
  }
  std::vector<double> v (2, 0.);
  v[0]=l0_; v[1]=rintorus_;
  return v;
}

double PolishDoughnut::centralEnthalpyPerUnitVolume() const
{
  // Converts internal cgs central enthalpy to SI
  double dens=central_enthalpy_cgs_;
# ifdef HAVE_UDUNITS
  dens = Units::Converter("erg/cm3", "J/m3")(dens);
# else
  GYOTO_WARNING << "Units ignored, please recompile Gyoto with --with-udunits"
		<< endl ;
# endif
  return dens;}

double PolishDoughnut::centralEnthalpyPerUnitVolume(string const &unit) const
{
  double dens = centralEnthalpyPerUnitVolume();
  if (unit != "") {
# ifdef HAVE_UDUNITS
    dens = Units::Converter("J/m3", unit)(dens);
# else
    GYOTO_WARNING << "Units ignored, please recompile Gyoto with --with-udunits"
		  << endl ;
# endif
  }
  return dens;
}

void PolishDoughnut::centralEnthalpyPerUnitVolume(double dens) {
# ifdef HAVE_UDUNITS
  dens = Units::Converter("J/m3", "erg/cm3")(dens);
# else
  GYOTO_WARNING << "Units ignored, please recompile Gyoto with --with-udunits"
		<< endl ;
# endif
  central_enthalpy_cgs_=dens;
}

void PolishDoughnut::centralEnthalpyPerUnitVolume(double dens,
						  string const &unit) {
  if (unit != "") {
# ifdef HAVE_UDUNITS
    dens = Units::Converter(unit, "J/m3")(dens);
# else
    GYOTO_WARNING << "Units ignored, please recompile Gyoto with --with-udunits"
		  << endl ;
# endif
  }
  centralEnthalpyPerUnitVolume(dens);
}

double PolishDoughnut::centralTemp() const
{return central_temperature_;}
void PolishDoughnut::centralTemp(double val)
{central_temperature_=val;}
double PolishDoughnut::beta() const { return beta_; }
void PolishDoughnut::beta(double b) { beta_ = b; }
void PolishDoughnut::magnetizationParameter(double rr) {
  magnetizationParameter_=rr;}
double PolishDoughnut::magnetizationParameter()const{
  return magnetizationParameter_;}
size_t PolishDoughnut::spectralOversampling() const
{ return spectral_oversampling_; }
void PolishDoughnut::spectralOversampling(size_t val)
{ spectral_oversampling_ = val; }
bool PolishDoughnut::changeCusp() const {return changecusp_;}
void PolishDoughnut::changeCusp(bool t) {changecusp_=t;}
bool PolishDoughnut::angleAveraged() const
{return angle_averaged_;}
void PolishDoughnut::angleAveraged(bool ang)
{
  angle_averaged_=ang;
  spectrumSynch_->angle_averaged(ang);
  spectrumPLSynch_->angle_averaged(ang);
}
bool PolishDoughnut::bremsstrahlung() const
{return bremsstrahlung_;}
void PolishDoughnut::bremsstrahlung(bool brems)
{bremsstrahlung_=brems;}
void PolishDoughnut::nonThermalDeltaExpo(std::vector<double> const &v) {
  if (v.size() != 2)
    GYOTO_ERROR("nonThermalDeltaExpo must have exactly 2 elements");
  deltaPL_= v[0];
  double expoPL = v[1];
  spectrumPLSynch_->PLindex(expoPL);
}
std::vector<double> PolishDoughnut::nonThermalDeltaExpo() const {
  std::vector<double> v (2, deltaPL_);
  v[1]=spectrumPLSynch_->PLindex();
  return v;
}

void PolishDoughnut::adafparams(std::vector<double> const &v) {
  if (v.size() != 2)
    GYOTO_ERROR("ADAF must have exactly 2 elements");
  adaf(true);
  ADAFtemperature_ = v[0];
  ADAFdensity_ = v[1];
}
std::vector<double> PolishDoughnut::adafparams() const {
  std::vector<double> v (2, ADAFtemperature_);
  v[1]=ADAFdensity_;
  return v;
}

void PolishDoughnut::adaf(bool t) {adaf_=t;}
bool PolishDoughnut::adaf() const {return adaf_;}

void PolishDoughnut::setParameter(Property const &p,
				  string const & name,
				  string const & content,
				  string const & unit) {
  // Override default behaviour to support obsolete format where
  // ADAFParameters was in ADAF
  if (name=="ADAF") {
    std::vector<double> v=FactoryMessenger::parseArray(content);
    if (v.size()) adafparams(v);
    return ;
  }
  Standard::setParameter(p, name, content, unit);
}

PolishDoughnut::~PolishDoughnut() {
  GYOTO_DEBUG << "PolishDoughnut Destruction" << endl;
  if (gg_) gg_ -> unhook(this);
}

void PolishDoughnut::metric(Gyoto::SmartPointer<Gyoto::Metric::Generic> met)
{
  if (gg_) gg_ -> unhook(this);
  Standard::metric(met);
  if (gg_) gg_ -> hook(this);
  GYOTO_DEBUG << "Metric set, calling lambda\n";
  

  // Initialize other members only if lambda(val) or
  // angmomrinner(vect) has been called already. Mutually exclusive.
  if (defangmomrinner_) angmomrinner(angmomrinner());
  else if (rochelobefilling_) lambda(lambda());
  
  GYOTO_DEBUG << "done\n";
}
void PolishDoughnut::tell(Hook::Teller * met) {
  if (met == gg_) {
    // Initialize other members only if lambda(val) or
    // angmomrinner(vect) has been called already. Mutually exclusive.
    if (defangmomrinner_) angmomrinner(angmomrinner());
    else if (rochelobefilling_) lambda(lambda());
  }
  else GYOTO_ERROR("BUG: PolishDoughnut::tell(Hook::Teller * met) called with"
		  "wrong metric");
}
int PolishDoughnut::Impact(Photon *ph, size_t index,
			   Astrobj::Properties *data) {
  if (beta_==1.) GYOTO_ERROR("Please set beta to != 1.");
  if (adaf_){
    // This is the Impact function for the Yuan+, Broderick+
    // ADAF model, this is actually no longer a Polish doughnut
    // -> only for comparison
    //cout << "ICI1" << endl;
    state_t coord;
    ph->getCoord(index, coord);
    double rr = coord[1], th = coord[2];
    // The outer boundary of the ADAF is simply RMax_ in xml
    // Setting an inner boundary at the ISCO (in projection)
    if (rr*sin(th) < gg_->getRms()) return 0;
    // This allows to reject the points close to the axis
    // such that the cylindrical radius is smaller than Sch ISCO ;
    // there, the Keplerian velocity is not defined
    state_t p1, p2;
    ph->getCoord(index, p1);
    ph->getCoord(index+1, p2);
    double t1 = p1[0], t2=p2[0];
    state_t cph;
    ph -> getCoord(t2, cph);
    double delta=giveDelta(&cph[0]);
    double coh[8];
    while (t2>t1){
      ph -> getCoord(t2, cph);
      for (int ii=0;ii<4;ii++)
	coh[ii] = cph[ii];
      getVelocity(coh, coh+4);
      processHitQuantities(ph, cph, coh, delta, data);
      t2 -= delta;
    }
    return 1;
  }
  return Standard::Impact(ph, index, data);
}
double PolishDoughnut::operator()(double const coord[4]) {
  // w1 = ((potential(r1, theta1, aa) - W_surface_)
  // /(W_centre_ - W_surface_));
  //
  // w1 < 0. outside polishdoughnut, anything inside funnel, 0<w<1
  // inside doughnut.
  //
  // so: operator()() < 0. <=> inside PolishDoughnut.
  double pos[4];
  for (int ii=0;ii<4;ii++) pos[ii]=coord[ii];
  double tmp =  W_surface_ - gg_->getPotential(pos,l0_);
  double rproj = coord[1] * sin(coord[2]);
  if (rproj<r_cusp_) {
    tmp = fabs(tmp)+(r_cusp_-rproj);
  }
  return tmp;
}
void PolishDoughnut::getVelocity(double const pos[4], double vel[4])
{
  if (adaf_) {
    // This will return the circular velocity at the
    // radius projected on the equat plane, or it's Keplerian approximation
    return gg_->circularVelocity(pos,vel,1);
  }
  double gtt=gg_->gmunu(pos,0,0);
  double gtph=gg_->gmunu(pos,0,3);
  double gphph=gg_->gmunu(pos,3,3);
  double Omega=-(l0_*gtt+gtph)/(l0_*gtph+gphph);
  double ut2=-1./(gtt+2.*gtph*Omega+gphph*Omega*Omega);
  if (ut2 < 0.) {
    stringstream ss;
    ss << "PolishDoughnut::getVelocity(pos=[";
    for (int i=0; i<3; ++i) ss << pos[i] << ", ";
    ss << pos[3] << "]): ut^2 is negative.";
    GYOTO_ERROR(ss.str());
  }
  vel[0] = sqrt(ut2);
  vel[1] = vel[2] = 0.;
  vel[3] = Omega*sqrt(ut2);
}
void PolishDoughnut::integrateEmission
(double * I, double const * boundaries,
 size_t const * chaninds, size_t nbnu,
 double dsem, state_t const &cph, double const *co) const
{
  // The original channels may or may not be contiguous. We split
  // each original channels into spectral_oversampling_ subchannels.
  // All we know is that each chunk of spectral_oversampling_
  // subchannels are contiguous. Don't try to recover contiguousness
  // in the original channels, it's too hard for now.
  double som1=1./double(spectral_oversampling_);
  size_t onbnu=nbnu*spectral_oversampling_; // number of subchannels
  size_t onbb = onbnu+nbnu; // number of subchannel boundaries : most
  // are used twice as subchannels are
  // contiguous in each channel.
  double * Inu = new double[onbb];
  double * bo = new double[onbb];
  size_t * ii = new size_t[2*onbnu]; // two indices for each subchannel
  double dnu;
  size_t k=0;
  for (size_t i=0; i<nbnu; ++i) {
    dnu=(boundaries[chaninds[2*i+1]]-boundaries[chaninds[2*i]])*som1;
    for (size_t j=0; j<spectral_oversampling_; ++j) {
      k=i*spectral_oversampling_+j;
      ii[2*k]=k+i;
      ii[2*k+1]=k+i+1;
      bo[ii[2*k]]=boundaries[chaninds[2*i]]+double(j)*dnu;
    }
    bo[ii[2*(i*spectral_oversampling_+spectral_oversampling_-1)+1]]
      =boundaries[chaninds[2*i+1]];
  }
  emission(Inu, bo, onbb, dsem, cph, co);
  for (size_t i=0; i<nbnu; ++i) {
    I[i]=0.;
    for (size_t j=0; j<spectral_oversampling_; ++j) {
      k=i*spectral_oversampling_+j;
      I[i]+=(Inu[ii[2*k+1]]+Inu[ii[2*k]])*0.5*fabs(bo[ii[2*k+1]]-bo[ii[2*k]]);
    }
  }
  delete [] Inu;
  delete [] bo;
  delete [] ii;
}

void PolishDoughnut::radiativeQ(double Inu[], // output
				double Taunu[], // output
				double const nu_ems[], size_t nbnu, // input
				double dsem,
				state_t const &coord_ph,
				double const coord_obj[8]) const {
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << endl;
# endif
  // This function computes the emission and transmission
  // for the Komissarov model, with both thermal and
  // non-thermal electron populations, with proper emission
  // and absorption.
  /* COMPUTING PHYS QUANTITIES */
  double rr = coord_ph[1], theta = coord_ph[2];//NB: rr is units of GM/c^2
  double Msgr = gg_->mass()*1e3; // Gyoto speaks in SI --> here cgs
  double T_electron=0., number_density=0.,
    bnorm = 0., theta_mag=0.;
  if (adaf_){
    if (!angle_averaged_){
      GYOTO_ERROR("In PolishDoughnut: ADAF should be called"
		 " only with angle averaging");
    }
    double zz = rr * fabs(cos(theta)), rcyl = rr * sin(theta);
    // fabs in zz: it is a distance, not an altitude
    if (zz>10.*rcyl) {
      // then exp factor will be
      // vanishingly small, can lead to bad behavior
      for (size_t ii=0; ii<nbnu; ++ii) {Inu[ii]=0.;Taunu[ii]=1.;}
      return;
    }
    double T0 = ADAFtemperature_, nth0 = ADAFdensity_;
    // From Broderick+11:
    number_density = nth0*pow(rr/2.,-1.1)*exp(-zz*zz/(2.*rcyl*rcyl));
    //cout << "ADAF ne= " << number_density << endl;
    T_electron = T0*pow(rr/2.,-0.84);
    double beta = 10., rS = 2.;
    bnorm = sqrt(8.*M_PI*1./beta*number_density
		 *GYOTO_PROTON_MASS_CGS * GYOTO_C_CGS * GYOTO_C_CGS
		 * rS / (12. * rr));
    //cout << "r z ne b= " << rr << " " << zz << " " << nth0*pow(rr/2.,-1.1) << " " << exp(-zz*zz/(2.*rcyl*rcyl)) << " " << number_density << " " << bnorm << endl;
  }else{
    double pos[4]={0.,rr,theta,0.};
    double ww = (gg_->getPotential(pos, l0_) - W_surface_)*DeltaWm1_;
    if (ww<=0.){//Will generate nan in computations w must be strictly positive
      if (fabs(ww)<w_tol) {
	if (ww!=0.) ww=fabs(ww);
	else ww=w_tol;//can be the case if w at entrance in doughnut is exactly 0
      }else{
	GYOTO_ERROR("In PolishDoughnut::emission() w<0!");
      }
    }
    double enthalpy_c=central_enthalpy_cgs_; // Warning: central_density_ is here
    // p+rho*c2 (enthalpy), not rho; model is different from std doughnut
    double g_tt=gg_->gmunu(&coord_ph[0],0,0),
      g_pp=gg_->gmunu(&coord_ph[0],3,3),
      g_tp=gg_->gmunu(&coord_ph[0],0,3),
      LL=g_tp*g_tp-g_tt*g_pp;
    double posc[4]={0.,r_centre_,M_PI/2.,0.};
    double g_ttc=gg_->gmunu(posc,0,0),
      g_ppc=gg_->gmunu(posc,3,3),
      g_tpc=gg_->gmunu(posc,0,3),
      LLc=g_tpc*g_tpc-g_ttc*g_ppc;
    double kappa = 0., kappam=0.;
    if (!angle_averaged_){
      kappa= pow(enthalpy_c,-CST_POLY_INDEX_M1)*(W_centre_-W_surface_)
	/((CST_POLY_INDEX+1)*(1+1./beta_));
      kappam = pow(LLc,-CST_POLY_INDEX_M1)/beta_*kappa;
    }else{
      kappa = pow(enthalpy_c,-CST_POLY_INDEX_M1)*(W_centre_-W_surface_)
	/(CST_POLY_INDEX+1);
    }
    double enthalpy = enthalpy_c*
      pow(
	  ww*
	  (kappa+kappam*pow(LLc,CST_POLY_INDEX_M1))
	  /(kappa+kappam*pow(LL,CST_POLY_INDEX_M1))
	  ,CST_POLY_INDEX
	  );
    number_density = (enthalpy-kappa*pow(enthalpy,1.+CST_POLY_INDEX_M1))
      /(GYOTO_C2_CGS*CST_MU_ELEC*GYOTO_ATOMIC_MASS_UNIT_CGS);
    //cout << "komis ne= " << number_density << endl;
    double number_density_central =
      (enthalpy_c-kappa*pow(enthalpy_c,1.+CST_POLY_INDEX_M1))
      /(GYOTO_C2_CGS*CST_MU_ELEC*GYOTO_ATOMIC_MASS_UNIT_CGS);
    //cout << "central nb density torus= " << number_density_central << endl;
    double magnetic_pressure = 0., fact_b=1.;
    // pm = b^2/fact_b
    if (!angle_averaged_){
      magnetic_pressure = kappam*pow(LL,CST_POLY_INDEX_M1)
	*pow(enthalpy,1.+CST_POLY_INDEX_M1);
      fact_b = 8.*M_PI;
    }else{
      double gas_pressure = kappa*pow(enthalpy,1.+CST_POLY_INDEX_M1);
      magnetic_pressure = gas_pressure/beta_;
      fact_b = 24.*M_PI;
    }
    bnorm = sqrt(fact_b*magnetic_pressure);
    // Redefining bnorm if magnetizationParameter_ is defined;
    // this is for compatibility with Jet.C
    if (magnetizationParameter_!=-1.){
      bnorm = sqrt(4.*M_PI*magnetizationParameter_
		   *GYOTO_PROTON_MASS_CGS * GYOTO_C_CGS * GYOTO_C_CGS
		   *number_density);
    }
    //cout << "ne_c, ne, Bc, B= " << number_density_central << " " << number_density << " " << magnetizationParameter_ << " " << sqrt(4.*M_PI*magnetizationParameter_*GYOTO_PROTON_MASS_CGS * GYOTO_C_CGS * GYOTO_C_CGS * number_density_central)  << " " << bnorm << endl;
    //GYOTO_ERROR("test pol");
    double bphi = bnorm/sqrt(g_pp+2*l0_*g_tp+l0_*l0_*g_tt);
    //NB: in Komissarov it is 2 p_mag in the numerator, but he uses
    // p_mag = B^2/2, and here we use the cgs p_mag = B^2/24pi
    double b4vec[4]={bphi*l0_,0,0,bphi}; // B 4-vector in BL frame
    // this vector is orthogonal to the fluid 4-vel, so it already
    // leaves in the comoving rest space, no need to project
    double vel[4]; // 4-velocity of emitter
    const_cast<PolishDoughnut*>(this)->getVelocity(coord_obj, vel);
    double photon_emframe[4]; // photon tgt vector projected in comoving frame
    for (int ii=0;ii<4;ii++){
      photon_emframe[ii]=coord_ph[ii+4]
	+vel[ii]*gg_->ScalarProd(&coord_ph[0],&coord_ph[4],vel);
    }
    double lnorm = gg_->ScalarProd(&coord_ph[0],photon_emframe,photon_emframe);
    if (lnorm<=0.) GYOTO_ERROR("In PolishDoughnut::radiativeq"
			      " photon_emframe should be spacelike");
    lnorm=sqrt(lnorm);
    double lscalb = gg_->ScalarProd(&coord_ph[0],photon_emframe,b4vec);
    theta_mag = acos(lscalb/(lnorm*bnorm));
    double sth = sin(theta_mag);//, cth = cos(theta_mag);
    if (sth==0.) GYOTO_ERROR("In PolishDoughnut::radiativeq: "
			    "theta_mag is zero leads to undefined emission");

    // doughnut's central temperature
    double T0 = central_temperature_;
    double kappabis = T0*pow(number_density_central,-CST_POLY_INDEX_M1);
    T_electron = kappabis*pow(number_density,CST_POLY_INDEX_M1);
    //cout << "Te= " << T_electron << endl;
  } // End of the switch between doughnut and adaf
  double Theta_elec = GYOTO_BOLTZMANN_CGS*T_electron
    /(GYOTO_ELECTRON_MASS_CGS*GYOTO_C2_CGS);
  
  double coef_ther=0.;
  // coef_ther: see e.g. Ozel+2000, eq. 6
  // here multiplied by Theta_elec coz there would be later
  // a multiplication by Theta_elec anyway
  double besselK3 = bessk(3, 1./Theta_elec),
    besselK2 = bessk(2, 1./Theta_elec),
    besselK1 = bessk1(1./Theta_elec);

  if (Theta_elec > 0.01){
    coef_ther = (3.*besselK3+besselK1)/(4.*besselK2)-1.;
  }else if (Theta_elec > 1e-5){
    // For small Theta_elec, Bessel functions become
    // very small, so I use a linear fit, correct to 1%
    // at theta_e=0.01, and even better for smaller values
    coef_ther=1.5*Theta_elec;
  }else{
    // too low Theta_e leads to Bnu being nan...
    for (size_t ii=0; ii<nbnu; ++ii) {Inu[ii]=0.;Taunu[ii]=1.;}
    return;
  }
  
  double expoPL = spectrumPLSynch_->PLindex();
  //cout << "expopl delta avg in PD= "<< expoPL << " " << deltaPL_ << " " << angle_averaged_ << endl;

  double number_density_PL =
    (expoPL-2.)/(expoPL-1.)*deltaPL_*coef_ther*number_density;
  double nuc = GYOTO_ELEMENTARY_CHARGE_CGS*bnorm
    /(2.*M_PI*GYOTO_ELECTRON_MASS_CGS*GYOTO_C_CGS);
  if (bnorm < 1e-5){
    // too low magnetic field leads to nan in emission
    // synchrotron is anyway vanishingly small
    for (size_t ii=0; ii<nbnu; ++ii) {Inu[ii]=0.;Taunu[ii]=1.;}
    return;
  }
  //cout << "r, ne, npl, nuc= " << rr << " " << number_density << " " << number_density_PL << " " << nuc << endl;
  //cout << "ne, delta, npl= " << number_density << " " << deltaPL_ << " " << number_density_PL << endl;

  // Defining jnus, anus
  double jnu_synch_ther[nbnu], anu_synch_ther[nbnu],
    jnu_synch_PL[nbnu], anu_synch_PL[nbnu],
    jnu_brems[nbnu], anu_brems[nbnu];

  for (size_t ii=0; ii<nbnu; ++ii){
    // Initializing to <0 value to create errors if not updated
    // [ exp(-anu*ds) will explose ]
    jnu_synch_ther[ii]=-1.;
    anu_synch_ther[ii]=-1.;
  } 
  
  // THERMAL SYNCHRO
  //cout << "doughnut stuff= " << T_electron << " " <<  number_density << " " << theta_mag << " " << nuc << " " << bnorm << " " << besselK2 << endl;
  spectrumSynch_->temperature(T_electron);
  spectrumSynch_->numberdensityCGS(number_density);
  spectrumSynch_->angle_B_pem(theta_mag);
  spectrumSynch_->cyclotron_freq(nuc);
  spectrumSynch_->besselK2(besselK2);

  spectrumSynch_->radiativeQ(jnu_synch_ther,anu_synch_ther,
  			     nu_ems,nbnu);

  // NONTHERMAL SYNCHRO
  if (deltaPL_!=0.){
    for (size_t ii=0; ii<nbnu; ++ii){
      // Initializing to <0 value to create errors if not updated
      jnu_synch_PL[ii]=-1.;
      anu_synch_PL[ii]=-1.;
    } 
    spectrumPLSynch_->numberdensityCGS(number_density_PL);
    spectrumPLSynch_->angle_B_pem(theta_mag);
    spectrumPLSynch_->cyclotron_freq(nuc);
    
    spectrumPLSynch_->radiativeQ(jnu_synch_PL,anu_synch_PL,
    				 nu_ems,nbnu);
  }

  // THERMAL BREMSSTRAHLUNG
  if (bremsstrahlung_){
    for (size_t ii=0; ii<nbnu; ++ii){
      // Initializing to <0 value to create errors if not updated
      jnu_brems[ii]=-1.;
      anu_brems[ii]=-1.;
    } 
    spectrumBrems_->temperature(T_electron);
    spectrumBrems_->numberdensityCGS(number_density);

    spectrumBrems_->radiativeQ(jnu_brems,anu_brems,
			       nu_ems,nbnu);
  }      

  // RETURNING TOTAL INTENSITY AND TRANSMISSION
  for (size_t ii=0; ii<nbnu; ++ii){
    double jnu_tot = jnu_synch_ther[ii],
      anu_tot = anu_synch_ther[ii];
    //cout << "at r,th= " << coord_ph[1] << " " << coord_ph[2] << endl;
    //cout << "torus jnu anu synch ther= " << jnu_tot << " " << anu_tot << endl;
    if (deltaPL_>0.){
      jnu_tot += jnu_synch_PL[ii];
      anu_tot += anu_synch_PL[ii];
    }
    if (bremsstrahlung_){
      jnu_tot += jnu_brems[ii];
      anu_tot += anu_brems[ii];
    }

    // expm1 is a precise implementation of exp(x)-1
    double em1=std::expm1(-anu_tot * dsem * gg_->unitLength());
    Taunu[ii] = em1+1.;
    Inu[ii] = anu_tot == 0. ? jnu_tot * dsem * gg_->unitLength() :
      -jnu_tot / anu_tot * em1; 
    
    if (Inu[ii]<0.)
      GYOTO_ERROR("In PolishDoughnut::radiativeQ: Inu<0");
    if (Inu[ii]!=Inu[ii] or Taunu[ii]!=Taunu[ii])
      GYOTO_ERROR("In PolishDoughnut::radiativeQ: Inu or Taunu is nan");
    if (Inu[ii]==Inu[ii]+1. or Taunu[ii]==Taunu[ii]+1.)
      GYOTO_ERROR("In PolishDoughnut::radiativeQ: Inu or Taunu is infinite");
    
  }
  
}
// Intersection of the constant angular momentum l0 with the Keplerian one
//double PolishDoughnut::intersection(double rr) const
PolishDoughnut::intersection_t::intersection_t(PolishDoughnut*parent)
  : papa(parent)
{
}
double PolishDoughnut::intersection_t::operator()(double rr) const
{
  double y = papa->gg_->getSpecificAngularMomentum(rr) - papa->l0_;

  return y ; // y = 0 gives 2 intersections,
  //the cusp and the central radius of the torus
}

double PolishDoughnut::outerradius_t::operator()(double rr) const
{
  double theta = M_PI/2.;
  double pos[4]={0.,rr,theta,0.};
  double ww = (papa->gg_->getPotential(pos,papa->l0_) - papa->W_surface_)*papa->DeltaWm1_;
  return ww;
}
