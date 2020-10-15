/*
    Copyright 2016, 2018-2020 Frederic Vincent & Thibaut Paumard

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

#include "GyotoOscilTorus.h"
#include "GyotoFactoryMessenger.h"
#include "GyotoUtils.h"
#include "GyotoPhoton.h"
#include "GyotoDefs.h"
#include "GyotoProperty.h"

#include <iostream>
#include <istream>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <limits> 
using namespace Gyoto;
using namespace Gyoto::Astrobj;
using namespace std;

GYOTO_PROPERTY_START(OscilTorus,
		     "Geometrical Torus with oscillations.")
GYOTO_PROPERTY_DOUBLE(OscilTorus, LargeRadius, largeRadius,
		      "Major radius, distance from centre of tube to centre of torus.")
GYOTO_PROPERTY_UNSIGNED_LONG(OscilTorus,  Mode, mode,
			     "Mode number of oscillations m.")
GYOTO_PROPERTY_DOUBLE(OscilTorus,PolyCst, polyCst,
		      "Polytropic constant kappa.")
GYOTO_PROPERTY_DOUBLE(OscilTorus,PolyIndex, polyIndex,
		      "Polytropic index n.")
GYOTO_PROPERTY_DOUBLE(OscilTorus,CentralDensity, centralDensity,
		      "Central density.")
GYOTO_PROPERTY_STRING(OscilTorus,PerturbKind, perturbKind,
		      "One of: Radial Vertical X Plus Breathing")
/*
  Nomenclature note: an old nomencalture is sometimes used in my
  notes, refering to "minus" and "plus" modes. This is due to the
  expression of the eigenfreq for these modes that differ by a sign.
  However, this is not the proper nomenclature and one should use
  Blaes et al. (2006) naming, so:
  old minus mode --> plus mode
  old plus mode  --> breathing mode
 */
GYOTO_PROPERTY_DOUBLE(OscilTorus, PerturbIntens, perturbIntens,
		      "Perturbations intensity.")
GYOTO_PROPERTY_FILENAME(OscilTorus, EmittingArea, emittingArea,
		      "Only for mode=0, file containing time series of cross section area")
GYOTO_PROPERTY_END(OscilTorus, Standard::properties)

GYOTO_PROPERTY_ACCESSORS_GEOMETRICAL_SPECIAL(OscilTorus, c_, largeRadius, gg_,
				 updateCachedValues(); , )
GYOTO_PROPERTY_ACCESSORS(OscilTorus, unsigned long, mode_, mode)
GYOTO_PROPERTY_ACCESSORS(OscilTorus, double, polycst_, polyCst)
GYOTO_PROPERTY_ACCESSORS_SPECIAL(OscilTorus, double, polyindex_, polyIndex,
				 updateCachedValues(); , )
GYOTO_PROPERTY_ACCESSORS(OscilTorus, double, central_density_, centralDensity)
GYOTO_PROPERTY_ACCESSORS(OscilTorus, double, perturb_intens_, perturbIntens)

void OscilTorus::perturbKind(std::string const &k) {
  if      (k == "Radial")    perturb_kind_ = Radial;
  else if (k == "Vertical")  perturb_kind_ = Vertical;
  else if (k == "X")         perturb_kind_ = X;
  else if (k == "Plus")      perturb_kind_ = Plus;
  else if (k == "Breathing") perturb_kind_ = Breathing;
  else {
    string errmsg="unknown perturbation kind: '";
    errmsg += k + "'";
    GYOTO_ERROR(errmsg.c_str());
  }
  updateCachedValues();
}
std::string OscilTorus::perturbKind() const {
  switch (perturb_kind_) {
  case Radial:    return "Radial";
  case Vertical:  return "Vertical";
  case X:         return "X";
  case Plus:      return "Plus";
  case Breathing: return "Breathing";
  default: GYOTO_ERROR("Unknown kind");
  }
  return "Should not reach this";
}

std::string OscilTorus::emittingArea() const {return emitting_area_;}
void OscilTorus::emittingArea(std::string const &f)  {
  if (f=="" || f.substr(f.size()-1) == "/") {
    emitting_area_ = "";
    with_cross_=0;
    tt_.clear();
    area_.clear();
    return;
  }
  ifstream file(f, ios::in);
  if (file) {
    with_cross_=1;
    double tt, area;
    tt_.clear();
    area_.clear();
    while (!file.eof()) {
      file >> tt >> area;
      if (area) {
	tt_.push_back(tt);
	area_.push_back(area);
      } else {
	//this means that last line of file was blank, area is
	//never 0, so just forget last line of data
	break;
      }
      file.ignore(numeric_limits<streamsize>::max(), '\n');//next line
    }
    nbt_=tt_.size();
    emitting_area_ = f;
  } else GYOTO_ERROR("Unable to read " + f);
}

Gyoto::Astrobj::OscilTorus::OscilTorus()
  : Standard("OscilTorus"),
    c_(10.8),
    mode_(0),
    polycst_(0.01),
    polyindex_(0.01),
    central_density_(0.01),
    perturb_kind_(Radial),
    perturb_intens_(0.1),
    kerrbl_(NULL),
    tt_(),
    area_(),
    nbt_(0),
    with_cross_(0),
    sigma_(0.),
    alpha_(0.),
    w1_(0.),
    w2_(0.),
    omr2_(0.),
    omth2_(0.),
    Omegac_(0.),
    lc_(0.),
    g_rr_(0.),
    g_thth_(0.),
    hold_(false)
{
  GYOTO_DEBUG << "Building OscilTorus" << endl;
}

Gyoto::Astrobj::OscilTorus::OscilTorus(const OscilTorus &orig)
  : Standard(orig),
    c_(orig.c_),
    mode_(orig.mode_),
    polycst_(orig.polycst_),
    polyindex_(orig.polyindex_),
    central_density_(orig.central_density_),
    perturb_kind_(orig.perturb_kind_),
    perturb_intens_(orig.perturb_intens_),
    kerrbl_(NULL),
    tt_(orig.tt_),
    area_(orig.area_),
    nbt_(orig.nbt_),
    with_cross_(orig.with_cross_),
    sigma_(orig.sigma_),
    alpha_(orig.alpha_),
    w1_(orig.w1_),
    w2_(orig.w2_),
    omr2_(orig.omr2_),
    omth2_(orig.omth2_),
    Omegac_(orig.Omegac_),
    lc_(orig.lc_),
    g_rr_(orig.g_rr_),
    g_thth_(orig.g_thth_),
    hold_(orig.hold_)
{
  GYOTO_DEBUG << "Copying OscilTorus" << endl;
  if (gg_) {
    kerrbl_=SmartPointer<Metric::KerrBL>(gg_);
    gg_->hook(this);
  }
}
OscilTorus * OscilTorus::clone() const { return new OscilTorus(*this); }

Gyoto::Astrobj::OscilTorus::~OscilTorus()
{
  GYOTO_DEBUG << "Destroying OscilTorus" << endl;
  if (gg_) gg_->unhook(this);
}

double OscilTorus::operator()(double const pos[4]) {
  //  cout << "pos= " << pos[1] << " " << pos[2] << " " << pos[3] << endl;
  // if far from the torus center, return any >0 value,
  // no hit
  //  if (fabs(pos[1]-c_)/c_*100. > 50.) return 100.;
  // USE RMax IN XML INSTEAD

  double x_bar=0., y_bar=0.;
  computeXbYb(pos,x_bar,y_bar);
  //  cout << "xb,yb= " << x_bar << " " << y_bar << endl;
  double uu=0.; // perturbation-dependent factor of surface function
  switch (perturb_kind_) {
  case Radial:
    uu = x_bar;
    break;
  case Vertical:
    uu = y_bar;
    break;
  case X:
    uu = x_bar*y_bar;
    break;
  case Plus:
  case Breathing:
    uu = 1+w1_*x_bar*x_bar+w2_*y_bar*y_bar;
    break;
  default:
    GYOTO_ERROR("In OscilTorus.C::operator():"
	       "Unrecognized perturbation kind");
  }
  // non-perturbed torus f
  double fnoperturb = omr2_*x_bar*x_bar + omth2_*y_bar*y_bar - 1.; 
  // correction
  double correc=perturb_intens_*sigma_*alpha_*uu
    *cos(mode_*pos[3] - Omegac_*(sigma_+mode_)*pos[0]);

  double ff     = fnoperturb + correc; // perturbed torus f
  //cout << "ff= " << ff << endl;

  return ff;
}

void OscilTorus::getVelocity(double const pos[4], double vel[4]) 
{

  //cout << "in getveloc" << endl;

  double gtt=kerrbl_->gmunu_up(pos,0,0);
  double gthth=kerrbl_->gmunu_up(pos,2,2);
  double grr=kerrbl_->gmunu_up(pos,1,1);
  double gpp=kerrbl_->gmunu_up(pos,3,3);
  double gtp=kerrbl_->gmunu_up(pos,0,3);  

  // Beta parameter
  double poly=(polyindex_+1.)/polyindex_;
  double centralp = polycst_*pow(central_density_,poly);

  double x_bar=0., y_bar=0.;
  computeXbYb(pos,x_bar,y_bar);

  double vr=0., vth=0.;
  switch (perturb_kind_) {
  case Radial:
    vr = 1.;
    break;
  case Vertical:
    vth = 1.;
    break;
  case X:
    vr = y_bar;
    vth = x_bar;
    break;
  case Plus:
  case Breathing:
    vr = 2.*w1_*x_bar;
    vth = 2.*w2_*y_bar;
    break;
  default:
    GYOTO_ERROR("In OscilTorus.C::operator():"
	       "Unrecognized perturbation kind");
  }

  double u_r = -perturb_intens_*sqrt(centralp/central_density_)*sqrt(g_rr_)
    *alpha_*vr
    *sin(mode_*pos[3] - Omegac_*(sigma_+mode_)*pos[0]);

  double u_th = perturb_intens_*sqrt(centralp/central_density_)*sqrt(g_thth_)
    *alpha_*vth
    *sin(mode_*pos[3] - Omegac_*(sigma_+mode_)*pos[0]);

  double u_t2 = (-1.-grr*u_r*u_r-gthth*u_th*u_th)/
    (gtt+lc_*lc_*gpp-2.*lc_*gtp);

  //cout << "num,deno= " << -1.-grr*u_r*u_r-gthth*u_th*u_th << " " << gtt+lc_*lc_*gpp-2.*lc_*gtp << endl;

  if (u_t2 < 0.) {
    stringstream ss;
    ss << "OscilTorus::getVelocity(pos=[";
    for (int i=0; i<3; ++i) ss << pos[i] << ", ";
    ss << pos[3] << "]): u_t^2 is negative.";
    GYOTO_ERROR(ss.str());
  }

  double u_t=-sqrt(u_t2);
  double u_p=-lc_*u_t;
  
  vel[0] = gtt*u_t+gtp*u_p;
  vel[1] = grr*u_r;
  vel[2] = gthth*u_th;
  vel[3] = gpp*u_p+gtp*u_t;
}

void OscilTorus::computeXbYb(const double * pos, double & xb, double & yb){
  double aa=kerrbl_->spin();
  
  // Computations at the torus center for Omegac_, lc_
  double posc[4]={0.,c_,M_PI/2.,0.};//don't care about t and phi
  double g_tt=gg_->gmunu(posc,0,0);// covar components
  double g_rr_=gg_->gmunu(posc,1,1);
  double g_thth_=gg_->gmunu(posc,2,2);
  double g_tp=gg_->gmunu(posc,0,3); 
  double g_pp=gg_->gmunu(posc,3,3);
  double Omegac_=1./(pow(c_,1.5)+aa); // Kepler rotation vel
  double lc_=-(Omegac_*g_pp+g_tp)/(Omegac_*g_tp+g_tt); // Rescaled ang mom

  // Now computations at the torus surface for gmunu_up coef
  double gtt=kerrbl_->gmunu_up(pos,0,0);
  double gthth=kerrbl_->gmunu_up(pos,2,2);
  double grr=kerrbl_->gmunu_up(pos,1,1);
  double gpp=kerrbl_->gmunu_up(pos,3,3);
  double gtp=kerrbl_->gmunu_up(pos,0,3);  

  // Beta parameter
  double poly=(polyindex_+1.)/polyindex_;
  double centralp = polycst_*pow(central_density_,poly);
  double cs2 = poly*centralp/central_density_;
  double ut_central2 = -1/(g_tt+g_pp*Omegac_*Omegac_+2.*g_tp*Omegac_);
  double beta2 = 2.*polyindex_*cs2
    /(c_*c_*ut_central2*Omegac_*Omegac_);
  if (beta2<=0.) {
    GYOTO_ERROR("In OscilTorus::computeXbYb(): "
	       "bad beta parameter");
  }
  double beta=sqrt(beta2);

  xb=1./beta
    *sqrt(g_rr_)*(pos[1]-c_)/c_;
  yb=1./beta
    *sqrt(g_thth_)*(M_PI/2.-pos[2])/c_;
}

void OscilTorus::metric(Gyoto::SmartPointer<Gyoto::Metric::Generic> met)
{
  if (!met) {
    if (gg_) gg_->unhook(this);
    kerrbl_=NULL;
    gg_=NULL;
    return;
  }
  kerrbl_ = Gyoto::SmartPointer<Gyoto::Metric::KerrBL>(met);
  if (!kerrbl_) GYOTO_ERROR("OscilTorus::metric(): only KerrBL, please");
  if (gg_) gg_->unhook(this);
  Standard::metric(met);
  gg_->hook(this);
  updateCachedValues();
}

void OscilTorus::tell(Hook::Teller* msg) {
  if (msg==gg_) updateCachedValues();
}

void OscilTorus::updateCachedValues() {
  if (hold_ || !gg_ || !c_) return;
  double aa=kerrbl_->spin();
  double posc[4]={0.,c_,M_PI/2.,0.};//don't care about t and phi
  double g_tt=gg_->gmunu(posc,0,0);// covar components
  g_rr_=gg_->gmunu(posc,1,1);
  g_thth_=gg_->gmunu(posc,2,2);
  double g_tp=gg_->gmunu(posc,0,3); 
  double g_pp=gg_->gmunu(posc,3,3);
  Omegac_=1./(pow(c_,1.5)+aa); // Kepler rotation vel
  lc_=-(Omegac_*g_pp+g_tp)/(Omegac_*g_tp+g_tt); // Rescaled ang mom

  // Epicyclic pulsation
  omr2_=1.-6./c_+8.*aa*pow(c_,-1.5)
    -3.*aa*aa/(c_*c_);
  omth2_=1.-4*aa*pow(c_,-1.5)
    +3.*aa*aa/(c_*c_);

  if (omr2_<=0. || omth2_<=0.) {
    GYOTO_ERROR("In OscilTorus::updateCachedValues(): "
	       "bad epicyclic freq");
  }
  double alpha0 = sqrt(polyindex_*sqrt(omr2_)*sqrt(omth2_)/M_PI);
  switch (perturb_kind_) {
  case Radial: // Radial oscillation
    {
      sigma_ = sqrt(omr2_);
      alpha_ = alpha0*sqrt(2*(polyindex_+1)*omr2_);
      break;
    }
  case Vertical: // Vertical oscillation
    {
      sigma_ = sqrt(omth2_);
      alpha_ = alpha0*sqrt(2*(polyindex_+1)*omth2_);
      break;
    }
  case X: // X mode
    {
      sigma_ = sqrt(omr2_+omth2_);
      alpha_ = alpha0*sqrt(4*(polyindex_+1)*(polyindex_+2)*omr2_*omth2_);
      break;
    }
  case Plus: // + mode
    {
      double sigmaminus2 = (
			    (2.*polyindex_+1)*(omr2_+omth2_) 
			    - sqrt(
				   4.*polyindex_*(polyindex_+1)
				   *(omr2_-omth2_)*(omr2_-omth2_)
				   +(omr2_+omth2_)*(omr2_+omth2_)
				   )
			    ) / (2.*polyindex_);
      sigma_ = sqrt(sigmaminus2);
      w1_ = - (omr2_
	       *(2.*omth2_+2.*polyindex_*omth2_-polyindex_*sigmaminus2)
	       ) / (omth2_-omr2_);
      w2_ = (omth2_
	     *(2.*omr2_+2.*polyindex_*omr2_-polyindex_*sigmaminus2)
	     ) / (omth2_-omr2_);
      alpha_ = alpha0*
	sqrt(
	     (polyindex_+2)*(sigmaminus2-(omr2_+omth2_))
	     /
	     (2.*polyindex_*sigmaminus2-(2.*polyindex_+1)*(omr2_+omth2_))
	     );
      break;
    }
  case Breathing: // breathing mode
    {
      double sigmaplus2 = (
			   (2.*polyindex_+1)*(omr2_+omth2_) 
			   + sqrt(
				  4.*polyindex_*(polyindex_+1)
				  *(omr2_-omth2_)*(omr2_-omth2_)
				  +(omr2_+omth2_)*(omr2_+omth2_)
				  )
			   ) / (2.*polyindex_);
      sigma_ = sqrt(sigmaplus2);
      w1_ = - (omr2_
	       *(2.*omth2_+2.*polyindex_*omth2_-polyindex_*sigmaplus2)
	       ) / (omth2_-omr2_);
      w2_ = (omth2_
	     *(2.*omr2_+2.*polyindex_*omr2_-polyindex_*sigmaplus2)
	     ) / (omth2_-omr2_);
      alpha_ = alpha0*
	sqrt(
	     (polyindex_+2)*(sigmaplus2-(omr2_+omth2_))
	     /
	     (2.*polyindex_*sigmaplus2-(2.*polyindex_+1)*(omr2_+omth2_))
	     );
      break;
    }
  default:
    GYOTO_ERROR("In OscilTorus.C::setParameter():"
	       "Unrecognized perturbation kind");
  }
}

double OscilTorus::emission(double nu_em, double, state_t const &cp, 
			      double const *) const{
  //cout << "r,theta,rcosth= " << cp[1] << " " << cp[2] << " " << cp[1]*cos(cp[2]) << endl;
  if (flag_radtransf_)
    GYOTO_ERROR("Radiative transfer not implemented for OscilTorus.");
  if (with_cross_){
    /*
      If the torus mode under consideration (breathing - and +)
      leads to change of volume, surface emitted intensity must be
      modulated accordingly.
      Simple assumption: I_nu \propto 1/volume
      ie: \propto 1/(cross-section area) if mode number m=0
      So here the area of the cross-section is determined
      at impact time to modulate the emitted intensity.
     */
    if (mode_!=0) {
      GYOTO_ERROR("In OscilTorus.C::emission:"
		 "mode=0 is required for area determination");
    }
    if (perturb_kind_==Vertical || perturb_kind_==X)
      GYOTO_ERROR("In OscilTorus::setParamter: bad perturbation kind");

    // Rescaled time and area determination
    double AA = Omegac_*sigma_; // cos modulation is 2pi/AA periodic
    double tt=cp[0],tmax=tt_[0], area=-1.;
    double myt = tt;
    while (myt>2.*M_PI/AA) myt-=2.*M_PI/AA; // myt is in ]0,2pi/AA]
    int ii=0;
    while (myt>tmax && ii < nbt_-1){
      ii++;
      tmax=tt_[ii];
    }  
    if (ii==0 || ii==nbt_-1){
      area=area_[ii];
    }else{
      area=
	area_[ii-1]+(myt-tt_[ii-1])*(area_[ii]-area_[ii-1])/(tt_[ii-1]-tt_[ii]);
    }
    if (area<=0. || area!=area) GYOTO_ERROR("In OscilTorus::emission:"
					  "bad area value");
    return 1./area;
  }else{
    return 1.;
  }
}

#ifdef GYOTO_USE_XERCES
void OscilTorus::setParameters(Gyoto::FactoryMessenger *fmp) {
  hold_=true;
  Standard::setParameters(fmp);
  hold_=false;
  updateCachedValues();
}
#endif
