#include "GyotoDeformedTorus.h"
#include "GyotoFactoryMessenger.h"
#include "GyotoUtils.h"
#include "GyotoPhoton.h"
#include "GyotoDefs.h"
#include "GyotoSpectrum.h"

#include <iostream>
#include <cmath>
#include <string>
#include <cstdlib>
#include <sstream>
using namespace Gyoto;
using namespace Gyoto::Astrobj;
using namespace std;

/// Properties
#include "GyotoProperty.h"
GYOTO_PROPERTY_START(DeformedTorus, "Slender torus subject to simple time-periodic deformations")
GYOTO_PROPERTY_SPECTRUM(DeformedTorus, Spectrum, spectrum)
GYOTO_PROPERTY_DOUBLE(DeformedTorus, LargeRadius, largeRadius)
GYOTO_PROPERTY_DOUBLE(DeformedTorus, Beta, beta)
GYOTO_PROPERTY_DOUBLE(DeformedTorus, BetaSt, betaSt)
GYOTO_PROPERTY_DOUBLE(DeformedTorus, Eta, eta)
GYOTO_PROPERTY_UNSIGNED_LONG(DeformedTorus, Mode, mode)
GYOTO_PROPERTY_STRING(DeformedTorus, PerturbKind, perturbKind)
GYOTO_PROPERTY_END(DeformedTorus, Standard::properties)

// Accessors
GYOTO_PROPERTY_ACCESSORS(DeformedTorus, SmartPointer<Spectrum::Generic>,
			 spectrum_, spectrum)
GYOTO_PROPERTY_ACCESSORS(DeformedTorus, double, c_, largeRadius)
GYOTO_PROPERTY_ACCESSORS_SPECIAL(DeformedTorus, double, param_beta_, beta,
				 if (param_beta_>=1.) GYOTO_ERROR("In DeformedTorus.C: beta should be << 1"); , )
GYOTO_PROPERTY_ACCESSORS(DeformedTorus, double,
			 param_beta_st_, betaSt)
GYOTO_PROPERTY_ACCESSORS(DeformedTorus, double, param_eta_, eta)
GYOTO_PROPERTY_ACCESSORS(DeformedTorus, unsigned long, mode_, mode)

void DeformedTorus::perturbKind(std::string const &k) {
  if      (k == "RadialTranslation")   perturb_kind_ = RadialTranslation;
  else if (k == "VerticalTranslation") perturb_kind_ = VerticalTranslation;
  else if (k == "Rotation")            perturb_kind_ = Rotation;
  else if (k == "Expansion")           perturb_kind_ = Expansion;
  else if (k == "RadialShear")         perturb_kind_ = RadialShear;
  else if (k == "VerticalShear")       perturb_kind_ = VerticalShear;
  else if (k == "PureShear")           perturb_kind_ = PureShear;
  else {
    string errmsg="unknown perturbation kind: '";
    errmsg += k + "'";
    GYOTO_ERROR(errmsg.c_str());
  }
}
std::string DeformedTorus::perturbKind() const {
  switch (perturb_kind_) {
  case RadialTranslation:   return "RadialTranslation";
  case VerticalTranslation: return "VerticalTranslation";
  case Rotation:            return "Rotation";
  case Expansion:           return "Expansion";
  case RadialShear:         return "RadialShear";
  case VerticalShear:       return "VerticalShear";
  case PureShear:           return "PureShear";
  default: GYOTO_ERROR("Unknown perturbation kind");
  }
  return "";
}

///

Gyoto::Astrobj::DeformedTorus::DeformedTorus()
  : Standard("DeformedTorus"),
    gg_(NULL),
    spectrum_(NULL),
    c_(10.8),
    mode_(0),
    param_beta_(0.01),
    param_beta_st_(0.01),
    param_eta_(0.01),
    perturb_kind_(RadialTranslation)
{
  GYOTO_DEBUG << "Building DeformedTorus" << endl;
}

Gyoto::Astrobj::DeformedTorus::DeformedTorus(const DeformedTorus &orig)
  : Standard(orig),
    gg_(NULL),
    spectrum_(NULL),
    c_(orig.c_),
    mode_(orig.mode_),
    param_beta_(orig.param_beta_),
    param_beta_st_(orig.param_beta_st_),
    param_eta_(orig.param_eta_),
    perturb_kind_(orig.perturb_kind_)
{
  if (orig.gg_()) {
    gg_=orig.gg_->clone();
    Standard::gg_ = gg_;
  }
  if (orig.spectrum_()) spectrum_ = orig.spectrum_->clone();
  GYOTO_DEBUG << "Copying DeformedTorus" << endl;
}
DeformedTorus * DeformedTorus::clone() const { return new DeformedTorus(*this); }

Gyoto::Astrobj::DeformedTorus::~DeformedTorus()
{
  GYOTO_DEBUG << "Destroying DeformedTorus" << endl;
}

double DeformedTorus::operator()(double const pos[4]) {

  // needed: operator()() < 0. <=> inside torus
  double posc[4]={0.,c_,M_PI/2.,0.};//don't care about t and phi
  double g_rr=gg_->gmunu(posc,1,1);// covar components
  double g_thth=gg_->gmunu(posc,2,2); 
  double aa=gg_->spin();
  double Omegac=1./(pow(c_,1.5)+aa);
  double omr2=1.-6./c_+8.*aa*pow(c_,-1.5)
    -3.*aa*aa/(c_*c_);
  double omth2=1.-4*aa*pow(c_,-1.5)
    +3.*aa*aa/(c_*c_);
  double x_bar=1./param_beta_
    *sqrt(g_rr)*(pos[1]-c_)/c_;
  double xb2=x_bar*x_bar;
  double y_bar=1./param_beta_
    *sqrt(g_thth)*(M_PI/2.-pos[2])/c_;
  double yb2=y_bar*y_bar;

  double a1=0., a2=0., a3=0.,
    b1=0., b2=0., b3=0.;
  switch (perturb_kind_) {
  case 1: // Radial translation
    a1=1.;a3=param_eta_*sin(Omegac*pos[0]);
    b2=1.;
    break;
  case 2: // Vertical translation
    a1=1.;
    b2=1.;b3=param_eta_*sin(Omegac*pos[0]);
    break;
  case 3: // Rotation
    a1=cos(Omegac*pos[0]);a2=sin(Omegac*pos[0]);
    b1=-sin(Omegac*pos[0]);b2=cos(Omegac*pos[0]);
    break;
  case 4: // Expansion
    a1=1.+param_eta_*sin(Omegac*pos[0]);
    b2=1.+param_eta_*sin(Omegac*pos[0]);
    break;
  case 5: // Simple shear radial
    a1=1.;a2=param_eta_*sin(Omegac*pos[0]);
    b2=1.;
    break;
  case 6: // Simple shear vertical
    a1=1.;
    b1=param_eta_*sin(Omegac*pos[0]);b2=1.;
    break;
  case 7: // Pure shear
    a1=1.+param_eta_*sin(Omegac*pos[0]);
    b2=1./a1;
    break;
  default:
    GYOTO_ERROR("In DeformedTorus.C::operator():"
	       "Unrecognized perturbation kind");
  }
  double deforx = a1*x_bar+a2*y_bar+a3;
  double defory = b1*x_bar+b2*y_bar+b3;
  double ff     = omr2*deforx*deforx
    +omth2*defory*defory
    -1.;

  return ff;
}

void DeformedTorus::getVelocity(double const pos[4], double vel[4]) 
{
  //cout << "pos in getvel: " << pos[2] << " " << pos[3] << endl;
  double aa=gg_->spin();

  // Computations at the torus center for Omegac, lc
  double posc[4]={0.,c_,M_PI/2.,0.};//don't care about t and phi
  double g_tt=gg_->gmunu(posc,0,0);// covar components
  double g_rr=gg_->gmunu(posc,1,1);
  double g_thth=gg_->gmunu(posc,2,2);
  double g_tp=gg_->gmunu(posc,0,3); 
  double g_pp=gg_->gmunu(posc,3,3);
  double Omegac=1./(pow(c_,1.5)+aa); // Kepler rotation vel
  double lc=-(Omegac*g_pp+g_tp)/(Omegac*g_tp+g_tt); // Rescaled ang mom

  // Now computations at the torus surface for gmunu_up coef
  double gtt=gg_->gmunu_up(pos,0,0);
  double gthth=gg_->gmunu_up(pos,2,2);
  double grr=gg_->gmunu_up(pos,1,1);
  double gpp=gg_->gmunu_up(pos,3,3);
  double gtp=gg_->gmunu_up(pos,0,3);  

  // xbar and ybar unperturbed
  double xbar=1./param_beta_
    *sqrt(g_rr)*(pos[1]-c_)/c_;
  double ybar=1./param_beta_
    *sqrt(g_thth)*(M_PI/2.-pos[2])/c_;

  // dr/dt and dtheta/dt depending on transformation

  double drdt=0., dthdt=0.;
  switch (perturb_kind_) {
  case 1: // Radial translation
    drdt=c_/sqrt(g_rr)*param_beta_*param_eta_
      *Omegac*cos(Omegac*pos[0]);
    break;
  case 2: // Vertical translation
    dthdt=-c_/sqrt(g_thth)*param_beta_*param_eta_
      *Omegac*cos(Omegac*pos[0]);
    break;
  case 3: // Rotation
    {
      double x0 = xbar*cos(Omegac*pos[0])-ybar*sin(Omegac*pos[0]),
	y0 = xbar*sin(Omegac*pos[0])+ybar*cos(Omegac*pos[0]);
      drdt = c_/sqrt(g_rr)*param_beta_*Omegac
	*(-sin(Omegac*pos[0])*x0+cos(Omegac*pos[0])*y0);
      dthdt = -c_/sqrt(g_thth)*param_beta_*Omegac
	*(-cos(Omegac*pos[0])*x0-sin(Omegac*pos[0])*y0);
    }
    break;
  case 4: // Expansion
    {
      double x0 = xbar/(1.+param_eta_*sin(Omegac*pos[0])),
	y0 = ybar/(1.+param_eta_*sin(Omegac*pos[0]));
      drdt = c_/sqrt(g_rr)*param_beta_*Omegac*param_eta_
	*cos(Omegac*pos[0])*x0;
      dthdt = -c_/sqrt(g_thth)*param_beta_*Omegac*param_eta_
	*cos(Omegac*pos[0])*y0;
    }
    break;
  case 5: // Simple shear radial
    {
      double y0 = ybar;
      drdt = c_/sqrt(g_rr)*param_beta_*Omegac*param_eta_
	*cos(Omegac*pos[0])*y0;
    }
    break;
  case 6: // Simple shear vertical
    {
      double x0 = xbar;
      dthdt = -c_/sqrt(g_thth)*param_beta_*Omegac*param_eta_
	*cos(Omegac*pos[0])*x0;
    }
    break;
  case 7: // Pure shear
    {
      double x0 = xbar/(1.+param_eta_*sin(Omegac*pos[0])),
	y0 = ybar*(1.+param_eta_*sin(Omegac*pos[0]));
      drdt = c_/sqrt(g_rr)*param_beta_*Omegac*param_eta_
	*cos(Omegac*pos[0])*x0;
      dthdt = -c_/sqrt(g_thth)*param_beta_*Omegac*param_eta_
	*(
	  -cos(Omegac*pos[0])
	  /( (1.+param_eta_*sin(Omegac*pos[0]))*
	     (1.+param_eta_*sin(Omegac*pos[0])) )
	  )*y0;
    }
    break;
  default:
    GYOTO_ERROR("In DeformedTorus.C::operator():"
	       "Unrecognized perturbation kind");
  }

  double u_t2 = -1./
    (gtt+lc*lc*gpp-2.*lc*gtp+drdt*drdt*grr+dthdt*dthdt*gthth);

  if (u_t2 < 0.) {
    stringstream ss;
    ss << "DeformedTorus::getVelocity(pos=[";
    for (int i=0; i<3; ++i) ss << pos[i] << ", ";
    ss << pos[3] << "]): u_t^2 is negative.";
    GYOTO_ERROR(ss.str());
  }

  double u_t=-sqrt(u_t2);
  double u_p=-lc*u_t;
  
  vel[0] = gtt*u_t+gtp*u_p;
  vel[1] = drdt*vel[0];
  vel[2] = dthdt*vel[0];
  vel[3] = gpp*u_p+gtp*u_t;
}

void DeformedTorus::metric(Gyoto::SmartPointer<Gyoto::Metric::Generic> met)
{
  if (met->kind() != "KerrBL")
    GYOTO_ERROR("DeformedTorus::metric(): only KerrBL, please");
  //if (gg_) gg_ -> unhook(this);
  gg_ = SmartPointer<Metric::KerrBL>(met);
  Generic::gg_ = gg_;
  //if (gg_) gg_ -> hook(this);
}




/*int DeformedTorus::Impact(Gyoto::Photon* ph, size_t index,
			 Astrobj::Properties *data) {
  double p1[8], p2[8];
  ph->getCoord(index, p1);
  ph->getCoord(index+1, p2);
  double tmin, minval;

  if (gg_ -> coordKind() == GYOTO_COORDKIND_SPHERICAL){
    //Allows theta and phi to be in the correct range
    checkPhiTheta(p1);
    checkPhiTheta(p2);
  }

  double t1 = p1[0], t2 = p2[0];
  double val1=(*this)(p1), val2=(*this)(p2);
  return 0;
  }*/

double DeformedTorus::emission(double nu_em, double, state_t const &, 
			      double const *) const{
  if (flag_radtransf_)
    GYOTO_ERROR("Radiative transfer not implemented for DeformedTorus.");
  //cout << "in flaring emission" << endl;
  //return (*spectrum_)(nu_em);
  return 1.;
}
