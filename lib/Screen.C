/*
    Copyright 2011-2020 Thibaut Paumard, Frederic Vincent

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
#include "GyotoFactoryMessenger.h"
#include "GyotoScreen.h"
#include "GyotoConverters.h"
#include "GyotoKerrBL.h"

#include <iostream>
#include <cmath>
#include <cstring>
#include <string>
#include <sstream>

#ifdef GYOTO_USE_CFITSIO
#include <fitsio.h>
#define throwCfitsioError(status) \
    { fits_get_errstatus(status, ermsg); GYOTO_ERROR(ermsg); }
#endif

#ifdef HAVE_BOOST_MULTIPRECISION_CPP_DEC_FLOAT_HPP
#include <boost/multiprecision/cpp_dec_float.hpp>
#endif

using namespace std ; 
using namespace Gyoto;

/// Properties
#include "GyotoProperty.h"
GYOTO_PROPERTY_START(Screen,
		     "Pin-hole camera with spectrometer.")
GYOTO_PROPERTY_METRIC(Screen, Metric, metric,
		      "It surrounds us and penetrates us; it binds the galaxy together.")
GYOTO_PROPERTY_SPECTROMETER(Screen, Spectrometer, spectrometer,
			    "Spectrometric capabilities.")
GYOTO_PROPERTY_DOUBLE_UNIT(Screen, Time, time,
			   "Observing date (seconds).")
GYOTO_PROPERTY_DOUBLE_UNIT(Screen, FieldOfView, fieldOfView,
			   "Field-of-view of the camera, in angle or length (radians).")
GYOTO_PROPERTY_DOUBLE(Screen, AzimuthalFieldOfView, azimuthalFieldOfView,
		      "Azimuthal field-of-view of the camera, for Spherical Angles images.")
GYOTO_PROPERTY_DOUBLE_UNIT(Screen, PALN, PALN,
			   "Position angle of the line of nodes of the equatorial plane (radians).")
GYOTO_PROPERTY_DOUBLE_UNIT(Screen, Inclination, inclination,
			   "Angle between the equatorial and sky planes (radians).")
GYOTO_PROPERTY_DOUBLE_UNIT(Screen, Argument, argument,
			   "Angle between the line of nodes and the Ox axis (radians).")
GYOTO_PROPERTY_DOUBLE_UNIT(Screen, Dangle1, dangle1, "(radians)")
GYOTO_PROPERTY_DOUBLE_UNIT(Screen, Dangle2, dangle2, "(radians)")
GYOTO_PROPERTY_SIZE_T(Screen, Resolution, resolution,
		      "Number of rows and columns.")
GYOTO_PROPERTY_DOUBLE_UNIT(Screen, Distance, distance,
			   "Distance from observer to centre of metric (meters).")
GYOTO_PROPERTY_DOUBLE(Screen, DMax, dMax,
		      "If DMax<Distance, use Distance=DMax and rescale (geometrical units).")
GYOTO_PROPERTY_DOUBLE_UNIT(Screen, FreqObs, freqObs,
			   "Observing frequency, wavelength or energy (Hz).")
GYOTO_PROPERTY_STRING(Screen, AngleKind, anglekind,
		      "\"EquatorialAngles\" (default), \"SphericalAngles\" or \"Rectilinear\".")
GYOTO_PROPERTY_STRING(Screen, ObserverKind, observerKind,
		      "\"ObserverAtInfinity\" (default), \"KeplerianObserver\", \"ZAMO\", \"VelocitySpecified\" or \"FullySpecified\".")
GYOTO_PROPERTY_FILENAME(Screen, Mask, maskFile,
			"FITS file to use as mask.")
GYOTO_PROPERTY_VECTOR_DOUBLE(Screen, FourVelocity, fourVel, "4-velocity of observer.")
GYOTO_PROPERTY_VECTOR_DOUBLE(Screen, ScreenVector1, screenVector1, "Screen e1 4-vector.")
GYOTO_PROPERTY_VECTOR_DOUBLE(Screen, ScreenVector2, screenVector2, "Screen e2 4-vector.")
GYOTO_PROPERTY_VECTOR_DOUBLE(Screen, ScreenVector3, screenVector3, "Screen e3 4-vector.")
GYOTO_PROPERTY_END(Screen, Object::properties)
///

// Default constructor
Screen::Screen() : 
  //tobs_(0.), fov_(M_PI*0.1), tmin_(0.), npix_(1001),
tobs_(0.), fov_(M_PI*0.1), azimuthal_fov_(2.*M_PI),
  npix_(1001), mask_(NULL), mask_filename_(""),
  distance_(1.), dmax_(GYOTO_SCREEN_DMAX), anglekind_(equatorial_angles),
  dangle1_(0.), dangle2_(0.),
  gg_(NULL), spectro_(NULL),
  freq_obs_(1.),
  observerkind_(GYOTO_OBSKIND_ATINFINITY)
{
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG_EXPR(dmax_);
# endif
  euler_[0]=euler_[1]=euler_[2]=0.;
  setProjection(0., 0., 0.);
  for (int ii=0;ii<4;ii++) {
    fourvel_[ii]=0.;
    screen1_[ii]=0.;
    screen2_[ii]=0.;
    screen3_[ii]=0.;
  }
}

Screen::Screen(const Screen& o) :
  SmartPointee(o),
  tobs_(o.tobs_), fov_(o.fov_), azimuthal_fov_(o.azimuthal_fov_),
  npix_(o.npix_), mask_(NULL),
  mask_filename_(o.mask_filename_),
  distance_(o.distance_),
  dmax_(o.dmax_),
  anglekind_(o.anglekind_),
  dangle1_(o.dangle1_), dangle2_(o.dangle2_),
  gg_(NULL), spectro_(NULL), freq_obs_(o.freq_obs_),
  observerkind_(o.observerkind_)
{
  if (o.gg_()) gg_=o.gg_->clone();
  if (o.spectro_()) spectro_ = o.spectro_ -> clone();
  int i;
  for (i=0; i<3; ++i) {
    euler_[i] = o.euler_[i];
    ex_[i]=o.ex_[i];
    ey_[i]=o.ey_[i];
    ez_[i]=o.ez_[i];
  }
  for (int ii=0;ii<4;ii++) {
    fourvel_[ii]=o.fourvel_[ii];
    screen1_[ii]=o.screen1_[ii];
    screen2_[ii]=o.screen2_[ii];
    screen3_[ii]=o.screen3_[ii];
  }
  if (o.mask_) {
    mask_=new double[npix_*npix_];
    memcpy(mask_, o.mask_, npix_*npix_*sizeof(double));
  } 

}
Screen * Screen::clone() const { return new Screen(*this); }

bool Screen::isThreadSafe() const {
  return Object::isThreadSafe()
    && (!gg_ || gg_ -> isThreadSafe())
    && (!spectro_ || spectro_ -> isThreadSafe());
}

Screen::~Screen(){if (mask_) delete [] mask_;}

std::ostream& Screen::print( std::ostream& o) const {
  o << "distance="    << distance_ << ", " ;
  o << "paln="        << euler_[0] << ", " ;
  o << "inclination=" << euler_[1] << ", " ;
  o << "argument="    << euler_[2] << "." << endl;
  return o;
}

/***************Definition of the physical scene**************/

void Screen::setProjection(const double paln,
			   const double incl,
			   const double arg) {
  euler_[0]=paln;
  euler_[1]=incl;
  euler_[2]=arg;
  computeBaseVectors();
}

void Screen::setProjection(const double dist,
			   const double paln,
			   const double incl,
			   const double arg) {
  distance_=dist;
  setProjection(paln, incl, arg);
}

void Screen::distance(double dist, const string &unit)    {
  distance(Units::ToMeters(dist, unit, gg_));
}
void Screen::distance(double dist)    {
  distance_=dist;
  computeBaseVectors();
}

void Screen::dMax(double dist) {
  dmax_ = dist;
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG_EXPR(dmax_);
# endif
}
double Screen::dMax() const { return dmax_; }


void Screen::PALN(double paln, const string &unit)        {
  if ((unit=="") || (unit=="rad"));
# ifdef HAVE_UDUNITS
  else paln = Units::Converter(unit, "rad")(paln);
# else
  else if (unit=="degree" || unit=="°") paln *= GYOTO_DEGRAD;
# endif
  PALN(paln);
}
void Screen::PALN(double paln)        {
  euler_[0]=paln; computeBaseVectors();
}

void Screen::inclination(double incl, const string &unit) { 
  if ((unit=="") || (unit=="rad"));
# ifdef HAVE_UDUNITS
  else incl = Units::Converter(unit, "rad")(incl);
# else
  else if (unit=="degree" || unit=="°") incl *= GYOTO_DEGRAD;
# endif
  inclination(incl);
}
void Screen::inclination(double incl) {
  euler_[1]=incl; computeBaseVectors();
}

void Screen::argument(double arg, const string &unit) { 
  if (unit=="" || unit=="rad");
# ifdef HAVE_UDUNITS
  else arg = Units::Converter(unit, "rad")(arg);
# else
  else if (unit=="degree" || unit=="°") arg *= GYOTO_DEGRAD;
# endif
  argument(arg);
}
void Screen::argument(double arg)     {
  euler_[2]=arg;  computeBaseVectors();
}

void Screen::freqObs(double fo, const string &unit) { 
  freqObs(Units::ToHerz(fo, unit));
}
void Screen::freqObs(double fo)     {
  freq_obs_=fo;
  GYOTO_DEBUG_EXPR(freq_obs_);
}
double Screen::freqObs() const {
  return freq_obs_;
}
double Screen::freqObs(const string &unit) const {
  return Units::FromHerz(freq_obs_, unit);
}


void Screen::metric(SmartPointer<Metric::Generic> gg) { gg_ = gg; computeBaseVectors(); }

int Screen::coordKind()      const { return gg_ -> coordKind(); }
double Screen::distance()    const { return distance_; }
double Screen::distance(const std::string& unit)    const {
  return Units::FromMeters(distance(), unit, gg_); 
}
double Screen::PALN()        const { return euler_[0]; }
double Screen::PALN(const string &unit) const {
  double paln = PALN();
  if (unit != "" && unit != "rad" ) {
# ifdef HAVE_UDUNITS
    paln = Units::Converter("rad", unit)(paln);
#else
    GYOTO_WARNING << "unit ignored, please recompile with --with-udunits\n";
#endif
  }
  return paln;
}

double Screen::inclination() const { return euler_[1]; }
double Screen::inclination(const string &unit) const {
  double incl = inclination();
  if (unit != "" && unit != "rad" ) {
# ifdef HAVE_UDUNITS
    incl = Units::Converter("rad", unit)(incl);
#else
    GYOTO_WARNING << "unit ignored, please recompile with --with-udunits\n";
#endif
  }
  return incl;
}

double Screen::argument()    const { return euler_[2]; }
double Screen::argument(const string &unit) const {
  double arg = argument();
  if (unit != "" && unit != "rad" ) {
# ifdef HAVE_UDUNITS
    arg = Units::Converter("rad", unit)(arg);
#else
    GYOTO_WARNING << "unit ignored, please recompile with --with-udunits\n";
#endif
  }
  return arg;
}

SmartPointer<Metric::Generic> Screen::metric() const { return gg_; }

void Screen::setObserverPos(const double coord[4]) {
  tobs_ = coord[0] * gg_ -> unitLength() / GYOTO_C;
  euler_[0]=M_PI;//Par défaut, A CHANGER SI BESOIN
  //dans le cas standard ex=-ephi et ephi dirige la ligne des noeuds d'où le pi
  //NB: c'est -pi dans mes notes, donc OK [2pi]
  //NB : ne décrit que la rotation de la caméra dans le plan x,y

  int coordkind=gg_ -> coordKind();
  switch (coordkind) {
  case GYOTO_COORDKIND_SPHERICAL:
    distance_=coord[1]* gg_ -> unitLength();
    //distance_=coord[1]; DEBUG
    euler_[1]=M_PI-coord[2];
    euler_[2]=-M_PI/2-coord[3];
    break;
  case GYOTO_COORDKIND_CARTESIAN:
    {
      double rks=sqrt(coord[1]*coord[1]+coord[2]*coord[2]+coord[3]*coord[3]);
      double thetaks=acos(coord[3]/rks);
      double phiks=atan2(coord[2],coord[1]);
      distance_=rks * gg_ -> unitLength();
      euler_[1]=M_PI-thetaks;
      euler_[2]=-M_PI/2-phiks;
    }
    break;
  default:
    GYOTO_ERROR("Incompatible coordinate kind in Screen::setObserverCoord");
  }
  computeBaseVectors();
}

void Screen::observerKind(const string &kind) {
  if(kind == "ObserverAtInfinity")
    observerkind_ = GYOTO_OBSKIND_ATINFINITY;
  else if(kind == "KeplerianObserver")
    observerkind_ = GYOTO_OBSKIND_KEPLERIAN;
  else if(kind == "ZAMO")
    observerkind_ = GYOTO_OBSKIND_ZAMO;
  else if(kind == "VelocitySpecified")
    observerkind_ = GYOTO_OBSKIND_VELOCITYSPECIFIED;
  else if(kind == "FullySpecified")
    observerkind_ = GYOTO_OBSKIND_FULLYSPECIFIED;
  else
    throwError("unknown observer kind");
}
string Screen::observerKind() const {
  switch (observerkind_) {
  case GYOTO_OBSKIND_ATINFINITY:
    return "ObserverAtInfinity";
  case GYOTO_OBSKIND_KEPLERIAN:
    return "KeplerianObserver";
  case GYOTO_OBSKIND_ZAMO:
    return "ZAMO";
  case GYOTO_OBSKIND_VELOCITYSPECIFIED:
    return "VelocitySpecified";
  case GYOTO_OBSKIND_FULLYSPECIFIED:
    return "FullySpecified";
  default:
    throwError("unknown observer kind tag");
  }
  return "will not reach here, this line to avoid compiler warning";
}

void Screen::setFourVel(const double coord[4]) {
  for (int ii=0;ii<4;ii++)
    fourvel_[ii]=coord[ii];
}

void Screen::fourVel(std::vector<double> const &coord) {
  if (coord.size() != 4)
    GYOTO_ERROR("base screen vectors require 4 elements");
  for (int ii=0;ii<4;ii++)
    fourvel_[ii]=coord[ii];
}
std::vector<double> Screen::fourVel() const {
  std::vector<double> output(4, 0.);
  for (int ii=0;ii<4;ii++) output[ii]=fourvel_[ii];
  return output;
}

void Screen::screenVector1(std::vector<double> const &coord) {
  if (coord.size() != 4)
    GYOTO_ERROR("base screen vectors require 4 elements");
  for (int ii=0;ii<4;ii++)
    screen1_[ii]=coord[ii];
}
std::vector<double> Screen::screenVector1() const {
  std::vector<double> output(4, 0.);
  for (int ii=0;ii<4;ii++) output[ii]=screen1_[ii];
  return output;
}

void Screen::screenVector2(std::vector<double> const &coord) {
  if (coord.size() != 4)
    GYOTO_ERROR("base screen vectors require 4 elements");
  for (int ii=0;ii<4;ii++)
    screen2_[ii]=coord[ii];
}
std::vector<double> Screen::screenVector2() const {
  std::vector<double> output(4, 0.);
  for (int ii=0;ii<4;ii++) output[ii]=screen2_[ii];
  return output;
}

void Screen::screenVector3(std::vector<double> const &coord) {
  if (coord.size() != 4)
    GYOTO_ERROR("base screen vectors require 4 elements");
  for (int ii=0;ii<4;ii++)
    screen3_[ii]=coord[ii];
}
std::vector<double> Screen::screenVector3() const {
  std::vector<double> output(4, 0.);
  for (int ii=0;ii<4;ii++) output[ii]=screen3_[ii];
  return output;
}

void Screen::setScreen1(const double coord[4]) {
  for (int ii=0;ii<4;ii++)
    screen1_[ii]=coord[ii];
}

void Screen::setScreen2(const double coord[4]) {
  for (int ii=0;ii<4;ii++)
    screen2_[ii]=coord[ii];
}

void Screen::setScreen3(const double coord[4]) {
  for (int ii=0;ii<4;ii++)
    screen3_[ii]=coord[ii];
}

void Screen::getObserverPos(double coord[]) const
{
  double r0     = distance_ / gg_ -> unitLength();
  //remark : Pb here if mass=0 (unitLength=0...) -> pb for flat metric

  //if (debug()) cout << "distance_ in screen= " << distance_ << endl;
  //double r0     = distance_ ;DEBUG
  double theta0 = M_PI-euler_[1];
  double phi0 = -M_PI/2-euler_[2];
    
  int coordkind = gg_ -> coordKind();

  stringstream ss;

  switch (coordkind) {
  case GYOTO_COORDKIND_CARTESIAN:
    coord[0] = tobs_ * GYOTO_C / gg_ -> unitLength();
    coord[1] = r0*cos(phi0)*sin(theta0);
    coord[2] = r0*sin(phi0)*sin(theta0);
    coord[3] = r0*cos(theta0);
    break;
  case GYOTO_COORDKIND_SPHERICAL:
    coord[0] = tobs_ * GYOTO_C / gg_ -> unitLength();
    coord[1] = r0;
    coord[2] = theta0;
    coord[3] = phi0;
    break;
  default:
    ss << "Incompatible coordinate kind in Screen::getObserverPos: "
       << coordkind;
    GYOTO_ERROR(ss.str());
  }
}

void Screen::getFourVel(double fourvel[]) const{
  for (int ii=0;ii<4;ii++) fourvel[ii]=fourvel_[ii];
}

void Screen::getScreen1(double output[]) const{
  for (int ii=0;ii<4;ii++) output[ii]=screen1_[ii];
}

void Screen::getScreen2(double output[]) const{
  for (int ii=0;ii<4;ii++) output[ii]=screen2_[ii];
}

void Screen::getScreen3(double output[]) const{
  for (int ii=0;ii<4;ii++) output[ii]=screen3_[ii];
}

/* SPECTROMETER */

void Screen::spectrometer(SmartPointer<Spectrometer::Generic> spr) { spectro_=spr; }
SmartPointer<Spectrometer::Generic> Screen::spectrometer() const { return spectro_; }

void Screen::getRayTriad(const size_t i, const size_t j,
			 double coord[],
			 bool compute_polar_basis,
			 double Ephi[4], double Etheta[4]) const {
  double xscr, yscr;
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << "(i=" << i << ", j=" << j << ", coord)" << endl;
# endif
  switch (anglekind_) {
  case Screen::spherical_angles:{
    /*
      GYOTO screen labelled by spherical
      angles a and b (see Fig. in user guide)
     */
    xscr = double(i-1)*fov_/(2.*double(npix_-1));
    yscr = M_PI-(double(j-1)*azimuthal_fov_/double(npix_-1));
    
    // NB: here xscr and yscr are the spherical angles
    // a and b ; the b->pi-b transformation boils down
    // to performing X->-X, just as below for equat angles.
    break;
  }
  case Screen::equatorial_angles:{
    /*
      GYOTO screen labelled by equatorial
      angles alpha and delta (see Fig. in user guide)
    */
    const double delta= fov_/double(npix_);
    yscr=delta*(double(j)-double(npix_+1)/2.);
    xscr=-delta*(double(i)-double(npix_+1)/2.);
    break;
  }
    // transforming X->-X (X being coord along e_1 observer vector)
    // this is due to the orientation convention of the screen 
    // (cf InitialConditions.pdf)
  case Screen::rectilinear: {
    if (fov_ >= M_PI)
      GYOTO_ERROR("Rectilinear projection requires fov_ < M_PI");
    const double xfov=2.*tan(fov_*0.5);
    const double delta= xfov/double(npix_);
    yscr=delta*(double(j)-double(npix_+1)/2.);
    xscr=-delta*(double(i)-double(npix_+1)/2.);
    break;
  }
  default:
    xscr=yscr=0.;
    GYOTO_ERROR("Unrecognized anglekind_");
  }
  getRayTriad(xscr, yscr, coord,
	      compute_polar_basis, Ephi, Etheta); 
}

void Screen::getRayTriad(double angle1, double angle2,			 
			 double coord[],
			 bool compute_polar_basis,
			 double Ephi[4], double Etheta[4]) const

{
  double normtol=1e-10;
  int i; // dimension : 0, 1, 2
  double pos[4];
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << "(angle1="<<angle1<<",angle2="<<angle2<<",coord)" << endl;
# endif
  getObserverPos(coord);

  double robs = distance_ / gg_->unitLength();

  if (robs > dmax_) {
    // scale
    coord[0] -= robs - dmax_;
    double scale = robs / dmax_;
    switch(gg_->coordKind()) {
    case GYOTO_COORDKIND_SPHERICAL:
      coord[1] = dmax_;
      break;
    case GYOTO_COORDKIND_CARTESIAN:
      coord[1] /= scale;
      coord[2] /= scale;
      coord[3] /= scale;
      break;
    default:
      throwError("Unimplemented coordkind");
    }
    angle1 *= scale;
    angle2 *= scale;
  }

  coord[4]=coord[5]=coord[6]=coord[7]=0.;//initializing 4-velocity
  double spherical_angle_a,
    spherical_angle_b;
  
  switch (anglekind_) {
  case spherical_angles:
    /*
      GYOTO screen labelled by spherical
      angles a and b (see Fig. in user guide)
     */
    spherical_angle_a = angle1+dangle1_;
    spherical_angle_b = angle2+dangle2_;
    break;
  case rectilinear:
    spherical_angle_a = atan(sqrt(angle1*angle1+angle2*angle2));
    spherical_angle_b = atan2(angle2, angle1);
    break;
  case equatorial_angles:
    /*
      GYOTO screen labelled by equatorial
      angles alpha and delta (see Fig. in user guide)
      Must compute spherical angles a and b
      from these.
    */
    
    /*
      NB: spherical_angles = spherical angles associated with the
      orthonormal 3-frame of the observer's local rest space
    */
    
    /*
      NBB: for spherical_angle_a, the relation comes from the spherical
      law of cosines of spherical trigonometry, the arcs
      spherical_angle_a, alpha and delta forming a spherical triangle,
      with alpha othogonal to delta.
      NBBB: for spherical_angle_b, the relation also comes from the spherical
      law of cosine. 
      --> Following transformations are OK even for non-small alpha, delta
    */

#ifdef HAVE_BOOST_MULTIPRECISION_CPP_DEC_FLOAT_HPP
    // using boost multiprecision to avoid information loss in trigonometry
    {
      boost::multiprecision::cpp_dec_float_100
	angle1100=angle1, angle2100=angle2, a, b;
      a=acos(cos(angle1100)*cos(angle2100));
      b=atan2(tan(angle2100),sin(angle1100));
      spherical_angle_a=a.convert_to<double>();
      spherical_angle_b=b.convert_to<double>();
    }
#else
    if (abs(angle1)<1e-6 || abs(angle2) < 1e-6) {
      spherical_angle_a = sqrt(angle1*angle1+angle2*angle2);
    } else {
      spherical_angle_a = acos(cos(angle1)*cos(angle2));
    }
    spherical_angle_b =
		   (angle1==0. && angle2==0.) ? 0. : atan2(tan(angle2),sin(angle1));
#endif
    break;
  default:
    spherical_angle_a=spherical_angle_b=0.;
    GYOTO_ERROR("Unknown angle type");
  }

  // Move these two angles to [0,pi], [0,2pi]
  double s1tmp=spherical_angle_a, s2tmp=spherical_angle_b;
  while (s1tmp>M_PI) s1tmp-=2.*M_PI;
  while (s1tmp<-M_PI) s1tmp+=2.*M_PI;//then s1 in [-pi,pi]
  if (s1tmp<0.) {
    s1tmp=-s1tmp;//then s1 in [0,pi]
    s2tmp+=M_PI;//thus, same direction
  }
  while (s2tmp>=2.*M_PI) s2tmp-=2.*M_PI;
  while (s2tmp<0.) s2tmp+=2.*M_PI;//then s2 in [0,2pi[
  spherical_angle_a=s1tmp;
  spherical_angle_b=s2tmp;

  /* 
     Tangent vector of incident photon in observer's local frame
     vel is thus orthogonal to observer's 4-velocity
     vel = vel[0]*screen1_+vel[1]*screen2_+vel[2]*screen3_
     where screen1,2,3_ is a triad in observer's rest space
  */

  //NB: following minus signs because the photon doesn't leave the screen, but
  //is heading towards it!
  double vel[3]={-sin(spherical_angle_a)*cos(spherical_angle_b),
		 -sin(spherical_angle_a)*sin(spherical_angle_b),
		 -cos(spherical_angle_a)};
  if (anglekind_ != spherical_angles) {
    // Rotate around Y by dangle1, then around new X by dangle2
    double c, s; sincos(dangle1_, &s, &c);
    double vel_rot[3]={c*vel[0]+s*vel[2],
		       vel[1],
		       -s*vel[0]+c*vel[2]};
    sincos(dangle2_, &s, &c);
    vel[0]=vel_rot[0];
    vel[1]=c*vel_rot[1]-s*vel_rot[2];
    vel[2]=s*vel_rot[1]+c*vel_rot[2];
    if (observerkind_!=GYOTO_OBSKIND_ATINFINITY){
      // Apply PALN
      sincos(M_PI-euler_[0], &s, &c);
      vel_rot[1]=vel[1];
      vel_rot[2]=vel[2];
      vel[0]=(c*vel_rot[0]-s*vel_rot[1]);
      vel[1]=(s*vel_rot[0]+c*vel_rot[1]);
    }
  }
  // 4-vector tangent to photon geodesic
  
  if (observerkind_==GYOTO_OBSKIND_ATINFINITY){
    /* 
       ---> Observer local frame not given in XML <---
       Assume observer static at infinity ("standard Gyoto")
       Treatment depending on coordinate system
    */
    switch (gg_ -> coordKind()) {
    case GYOTO_COORDKIND_CARTESIAN:
      {
	double rr=coord[1]*coord[1]+
	  coord[2]*coord[2]+
	  coord[3]*coord[3],
	  rinf=20.; // this rinf is just a very crude test
               	    // I take rinf=10*rhor in Sch metric
	if (rr<rinf)
	  GYOTO_ERROR("In Screen::getRayTriad: "
		     "observer is not at spatial infinity "
		     "and should be");

	//Transforming to KS:
	for (i=0;i<3;++i) {
	  coord[5]+=vel[i]*ex_[i];//xdot0
	  coord[6]+=vel[i]*ey_[i];//ydot0
	  coord[7]+=vel[i]*ez_[i];//zdot0
	}
	
      }
      break;
    case GYOTO_COORDKIND_SPHERICAL:
      {
	if (coord[2]==0. || coord[2]==M_PI)
	  GYOTO_ERROR("Please move Screen away from z-axis");

	double rr=coord[1],
	  rinf=20.; // this rinf is just a very crude test
               	    // I take rinf=10*rhor in Sch metric
	if (rr<rinf)
	  GYOTO_ERROR("In Screen::getRayTriad: "
		     "observer is not at spatial infinity "
		     "and should be");

	pos[0]=coord[0];pos[1]=coord[1];pos[2]=coord[2];pos[3]=coord[3];
	double grr=gg_->gmunu(pos,1,1), 
	  gthth=gg_->gmunu(pos,2,2), gphph=gg_->gmunu(pos,3,3);
	coord[5]=-vel[2]/sqrt(grr);
	double sp=sin(euler_[0]);
	double cp=cos(euler_[0]);
	coord[6]=(-sp*vel[0]+cp*vel[1])/sqrt(gthth);
	coord[7]=( cp*vel[0]+sp*vel[1])/sqrt(gphph);
      }
      break;
    default:
      GYOTO_ERROR("Incompatible coordinate kind in Screen::getRayTriad()");
      break;
    }

    // 0-component of photon tangent 4-vector found by normalizing
    gg_ -> nullifyCoord(coord);
    
  } else {
    /* 
       ---> Observer local frame given in XML <---
       Express photon tangent 4-vector in the observer basis
       Treatment is coordinate independent 
       (except for z-axis check right below)
    */

    if (gg_ -> coordKind() == GYOTO_COORDKIND_SPHERICAL){
      if (coord[2]==0. || coord[2]==M_PI)
	GYOTO_ERROR("Please move Screen away from z-axis");
    }

    // Implemented observer specifid in XML, local tetrad computed by Metric
    const double fourpos[4]={coord[0],coord[1],coord[2],coord[3]};
    double fourvel[4], screen1[4], screen2[4], screen3[4];
    memcpy(fourvel, fourvel_, 4*sizeof(double));
    memcpy(screen1, screen1_, 4*sizeof(double));
    memcpy(screen2, screen2_, 4*sizeof(double));
    memcpy(screen3, screen3_, 4*sizeof(double));
    gg_ -> observerTetrad(observerkind_,fourpos,fourvel,screen1,
			  screen2,screen3);

    /* Photon tagent 4-vector l defined by: l = p + fourvel_
       where p gives the direction of the photon
       in the observer's rest space (orthogonal to fourvel).
       Here we choose the particular tangent 4-vector l
       that satisfies l.fourvel_=-1
       Then l = fourvel_ + (orthogonal proj of l onto rest space)
       = fourvel_ + p
       and p = vel[0]*screen1_ + vel[1]*screen2_ + vel[2]*screen3_
       with p.p = 1, thus l.l = 0 as it should
    */

    coord[4]=vel[0]*screen1[0]
      +vel[1]*screen2[0]
      +vel[2]*screen3[0]
      +fourvel[0];
    coord[5]=vel[0]*screen1[1]
      +vel[1]*screen2[1]
      +vel[2]*screen3[1]
      +fourvel[1];
    coord[6]=vel[0]*screen1[2]
      +vel[1]*screen2[2]
      +vel[2]*screen3[2]
      +fourvel[2];
    coord[7]=vel[0]*screen1[3]
      +vel[1]*screen2[3]
      +vel[2]*screen3[3]
      +fourvel[3];
  }
  if (fabs(gg_->ScalarProd(coord,coord+4,coord+4))>normtol){
    GYOTO_SEVERE << "In Screen::getRayTriad: "
		 << "tangent 4-vector to photon not properly normalized: "
		 << "norm = "
		 << gg_->ScalarProd(coord,coord+4,coord+4)
		 << endl;
  }

  if (compute_polar_basis==true){
    switch (gg_ -> coordKind()) {
    case GYOTO_COORDKIND_SPHERICAL:
      {
	double ca, sa; sincos(spherical_angle_a, &sa, &ca);
	double cb, sb; sincos(spherical_angle_b, &sb, &cb);
	double Ephi_screenBasis[3] = {-sb*sb*(1-ca)-ca,
				      sb*cb*(1-ca),
				      cb*sa};
	double Etheta_screenBasis[3] = {-sb*cb*(1-ca),
					(cb*cb*(1-ca)+ca),
					-sb*sa};
	
	double cp, sp; sincos(euler_[0], &sp, &cp);
	
	if (observerkind_==GYOTO_OBSKIND_ATINFINITY){
	  double grr=gg_->gmunu(coord,1,1), 
	    gthth=gg_->gmunu(coord,2,2), gphph=gg_->gmunu(coord,3,3);
	  // Ephi
	  Ephi[0]=0.;
	  Ephi[1]=-Ephi_screenBasis[2]/sqrt(grr);      
	  Ephi[2]=(-sp*Ephi_screenBasis[0]
		   +cp*Ephi_screenBasis[1])/sqrt(gthth);
	  Ephi[3]=( cp*Ephi_screenBasis[0]
		    +sp*Ephi_screenBasis[1])/sqrt(gphph);
	  // Etheta
	  Etheta[0]=0.;
	  Etheta[1]=-Etheta_screenBasis[2]/sqrt(grr);
	  Etheta[2]=(-sp*Etheta_screenBasis[0]
		     +cp*Etheta_screenBasis[1])/sqrt(gthth);
	  Etheta[3]=( cp*Etheta_screenBasis[0]
		      +sp*Etheta_screenBasis[1])/sqrt(gphph);
	  //cout << "In Screen init Ephi= " << Ephi[0] << " " << Ephi[1] << " " << coord[1]*Ephi[2] << " " << coord[1]*abs(sin(coord[2]))*Ephi[3] << endl;
	  //cout << "In Screen init Etheta= " << Etheta[0] << " " << Etheta[1] << " " << coord[1]*Etheta[2] << " " << coord[1]*abs(sin(coord[2]))*Etheta[3] << endl;
	  //throwError("test Eth");
	}else{
	  throwError("Observer should be at infinity");
	}
      
	// double k_phi = gg_->gmunu(coord,3,3)*coord[7]
	// 	+ gg_->gmunu(coord,0,3)*coord[4]; // phi covariant compo
	//                      // of tangent vector to null geodesic
	// double ktheta = coord[6];
	// double rr = coord[1];
	// double sth=sin(coord[2]);
	// if (sth==0.)
	// 	GYOTO_ERROR("Please move Screen away from z-axis");
	// double rsm1 = 1./(rr*sth);
	
	// double sp=sin(euler_[0]);
	// double cp=cos(euler_[0]);
	// // Ephi
	// Ephi[0]=0.;
	// Ephi[1]=-k_phi*rsm1;      
	// Ephi[2]=sp/rr;
	// Ephi[3]=-cp*rsm1;
	// // Etheta
	// Etheta[0]=0.;
	// Etheta[1]=rr*ktheta;
	// Etheta[2]=cp/rr;
	// Etheta[3]=sp*rsm1;
	break;
      }
      
    default:
      GYOTO_ERROR("Non implemented coord kind for polarization");
    }
  }
}

/************** MASK ******************/

void Screen::mask(double const * const mm, size_t res) {
  if (res) npix_=res;
  if (mask_) { delete[] mask_; mask_=NULL; mask_filename_=""; }
  if (mm) {
    mask_ = new double[npix_*npix_];
    memcpy(mask_, mm, npix_*npix_*sizeof(double));
  }
}

double const * Screen::mask() const { return mask_; }

void Screen::maskFile(std::string const &fname) {
# ifdef GYOTO_USE_CFITSIO
  if (fname != "") fitsReadMask(fname);
# else
  GYOTO_WARNING << "No FITS i/o, Screen mask file ignored" << endl;
# endif 
}

std::string Screen::maskFile() const {return mask_filename_;}

#ifdef GYOTO_USE_CFITSIO
void Screen::fitsReadMask(std::string const &filename) {
  GYOTO_DEBUG << "Screen::fitsReadMask(\"" << filename << "\")"<<endl;
  char*     pixfile   = const_cast<char*>(filename.c_str());
  fitsfile* fptr      = NULL;
  int       status    = 0;
  int       anynul    = 0;
  long      naxes []  = {1, 1, 1};
  long      fpixel[]  = {1, 1, 1};
  long      inc   []  = {1, 1, 1};
  char      ermsg[31] = ""; // ermsg is used in throwCfitsioError()
  if (fits_open_file(&fptr, pixfile, 0, &status)) {
    GYOTO_WARNING << "Unable to read Screen mask file '"
		  << filename << "', ignoring." << endl;
    return;
  }
  if (fits_get_img_size(fptr, 3, naxes, &status)) throwCfitsioError(status) ;
  if (naxes[0] != naxes[1])
    GYOTO_ERROR("Screen::fitsReadMask(): mask must be square");
  //  if (naxes[2] > 1)
  //    GYOTO_ERROR("Screen::fitsReadMask(): mask must have only one plane");
  npix_=naxes[0];
  if (mask_) { delete[] mask_; mask_=NULL; mask_filename_="";}
  mask_ = new double[npix_*npix_];
  inc[2]=naxes[2]; // read only first plane!
  if (fits_read_subset(fptr, TDOUBLE, fpixel, naxes, inc,
		       0, mask_,&anynul,&status)) {
    GYOTO_DEBUG << " error, trying to free pointer" << endl;
    delete [] mask_; mask_=NULL;
    throwCfitsioError(status) ;
  }
  if (fits_close_file(fptr, &status)) throwCfitsioError(status) ;
  mask_filename_=filename;
  fptr = NULL;
}

void Screen::fitsWriteMask(string const &fname) {
  std::string filename = fname==""?mask_filename_:fname;
  if (filename=="") filename=mask_filename_;
  if (filename=="") GYOTO_ERROR("no filename specified");
  if (!mask_) GYOTO_ERROR("Screen::fitsWriteMask(filename): nothing to save!");
  char*     pixfile   = const_cast<char*>(filename.c_str());
  fitsfile* fptr      = NULL;
  int       status    = 0;
  long      naxes []  = {long(npix_), long(npix_)};
  long      fpixel[]  = {1,1};

  char      ermsg[31] = ""; // ermsg is used in throwCfitsioError()

  ////// CREATE FILE
  GYOTO_DEBUG << "creating file \"" << pixfile << "\"... ";
  fits_create_file(&fptr, pixfile, &status);
  if (debug()) cerr << "done." << endl;
  fits_create_img(fptr, DOUBLE_IMG, 2, naxes, &status);
  if (status) throwCfitsioError(status) ;

  ////// SAVE EMISSION IN PRIMARY HDU ///////
  GYOTO_DEBUG << "saving emission_\n";
  fits_write_pix(fptr, TDOUBLE, fpixel, npix_*npix_, mask_, &status);
  if (status) throwCfitsioError(status) ;

  ////// CLOSING FILE ///////
  GYOTO_DEBUG << "close FITS file\n";
  if (fits_close_file(fptr, &status)) throwCfitsioError(status) ;
  mask_filename_=filename;
  fptr = NULL;
}
#endif

bool Screen::operator()(size_t i, size_t j) {
  GYOTO_DEBUG << "i=" << i << ", j=" << j << endl;
  if ( i<=0 || i> npix_ || j>npix_ || j<=0) GYOTO_ERROR("wrong index");
  if (!mask_) return true;
  return mask_[i-1+npix_*(j-1)];
}

void Screen::computeBaseVectors() {
  // See http://en.wikipedia.org/wiki/Euler_angles

  //  if (debug()) cout << "In compute base vectore" << endl;
  double ca, sa; sincos(euler_[0], &sa, &ca); // paln
  double cb, sb; sincos(euler_[1], &sb, &cb); // inclination
  double cc, sc; sincos(euler_[2], &sc, &cc); // argument
  double unit  = 1.;//1./distance_;

  //NB: here x,y,z are KS coordinates
  ex_[0] = unit*( ca*cc - sa*cb*sc);//component of KS vector eX along observer East direction ex
  ex_[1] = unit*( sa*cc + ca*cb*sc);//component of KS vector eX along observer North direction ey
  ex_[2] = unit*( sb*sc);//etc...
  ey_[0] = unit*(-ca*sc - sa*cb*cc);
  ey_[1] = unit*(-sa*sc + ca*cb*cc);
  ey_[2] = unit*( sb*cc);
  ez_[0] = unit*( sb*sa);
  ez_[1] = unit*(-sb*ca);
  ez_[2] = unit*( cb);

}

void Screen::coordToSky(const double pos[4], double skypos[3],
			bool geometrical) const {
  double xyz[3];
  coordToXYZ(pos, xyz);
  double ul = geometrical?1.:gg_ -> unitLength();

  skypos[0]=(xyz[0]*ex_[0]+xyz[1]*ey_[0]+xyz[2]*ez_[0]) * ul;
  skypos[1]=(xyz[0]*ex_[1]+xyz[1]*ey_[1]+xyz[2]*ez_[1]) * ul;
  skypos[2]=(xyz[0]*ex_[2]+xyz[1]*ey_[2]+xyz[2]*ez_[2]) * ul;
}

void Screen::skyToCoord(const double skypos[3], double pos[4],
			bool geometrical) const {
  // convert from m to ug
  double ulm1 = geometrical?1.:1./gg_ -> unitLength();
  double alpha=skypos[0]*ulm1, delta=skypos[1]*ulm1, zobs=skypos[2]*ulm1;

  // rotate from sky frame to BH frame in Cartesian coordinates
  double ca, sa; sincos(euler_[0], &sa, &ca);
  double cb, sb; sincos(euler_[1], &sb, &cb);
  double cc, sc; sincos(euler_[2], &sc, &cc);
  double xyz[3];
  xyz[0]=(-cb*sa*sc+ca*cc)*alpha+(ca*cb*sc+cc*sa)*delta+sb*sc*zobs;
  xyz[1]=(-cb*cc*sa-ca*sc)*alpha+(ca*cb*cc-sa*sc)*delta+cc*sb*zobs;
  xyz[2]=sa*sb*alpha-ca*sb*delta+cb*zobs;

  // possibly convert to spherical
  switch (gg_->coordKind()) {
  case GYOTO_COORDKIND_CARTESIAN:
    pos[1]=xyz[0];
    pos[2]=xyz[1];
    pos[3]=xyz[2];
    break;
  case GYOTO_COORDKIND_SPHERICAL:
    cartesianToSpherical(xyz, pos+1);
    break;
  default:
    GYOTO_ERROR("Unimplemented coordinate kind in Screen::skyToCoord");
  }
}

std::ostream & Screen::printBaseVectors(std::ostream &o) const {
 
  o << endl;
  o << setprecision(3) << setw(8) << ex_[0] << ", "
    << setprecision(3) << setw(8) << ey_[0] << ", "
    << setprecision(3) << setw(8) << ez_[0] << endl;
  o << setprecision(3) << setw(8) << ex_[1] << ", "
    << setprecision(3) << setw(8) << ey_[1] << ", "
    << setprecision(3) << setw(8) << ez_[1] << endl;
  o << setprecision(3) << setw(8) << ex_[2] << ", "
    << setprecision(3) << setw(8) << ey_[2] << ", "
    << setprecision(3) << setw(8) << ez_[2] << endl;
  return o ;
} 


/***************Useful functions**************/


void Screen::coordToXYZ(const double pos[4], double xyz[3]) const {
  int coordkind = gg_ -> coordKind();
  switch(coordkind) {
  case GYOTO_COORDKIND_SPHERICAL:
    xyz[0] = pos[1]*sin(pos[2])*cos(pos[3]);
    xyz[1] = pos[1]*sin(pos[2])*sin(pos[3]);
    xyz[2] = pos[1]*cos(pos[2]);
    break;
  case GYOTO_COORDKIND_CARTESIAN:
    xyz[0] = pos[1];
    xyz[1] = pos[2];
    xyz[2] = pos[3];
    break;
  default:
    GYOTO_ERROR("Incompatible coordinate kind in Screen::coordToXYZ");
    break;
  }
}

double Screen::time() const { return tobs_ ; }
double Screen::time(const string &unit) const {
  return Units::FromSeconds(time(), unit, gg_) ;
}

void Screen::time(double tobs, const string &unit) {
# ifdef GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG_EXPR(tobs);
  GYOTO_DEBUG_EXPR(unit);
# endif
  time(Units::ToSeconds(tobs, unit, gg_));
}
void Screen::time(double tobs) {
  tobs_ = tobs;
# ifdef GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG_EXPR(tobs_);
# endif
}

//double Screen::getMinimumTime() { return tmin_; }
//void Screen::setMinimumTime(double tmin) { tmin_ = tmin; }
double Screen::fieldOfView() const { return fov_; }

double Screen::fieldOfView(string const &unit) const {
  double fov = fieldOfView();
  if (unit=="" || unit=="rad") ;
  else if (unit=="geometrical") fov *= distance_ / gg_ -> unitLength();
# ifdef HAVE_UDUNITS
  else if (Units::areConvertible(unit, "m"))
    fov = Units::FromMeters(fov*distance_, unit) ;
  else fov = Units::Converter("rad", unit)(fov);
# else
  else if (unit=="degree") fov /= GYOTO_DEGRAD;
  else if (unit=="arcmin") fov /= GYOTO_MINRAD;
  else if (unit=="arcsec") fov /= GYOTO_SECRAD;
  else if (unit=="milliarcsec") fov /= GYOTO_MASRAD;
  else if (unit=="microarcsec") fov /= GYOTO_MUASRAD;
  else {
    stringstream ss;
    ss << "Screen::fieldOfView(): unknown unit: \"" << unit << "\""
       << " (you may have more chance compiling gyoto with --with-udunits)";
    GYOTO_ERROR(ss.str());
  }
# endif
  return fov;
}

void Screen::fieldOfView(double fov, const string &unit) {
  if (unit=="" || unit=="rad") ;
  else if (unit=="geometrical") fov *= gg_ -> unitLength() / distance_ ;
# ifdef HAVE_UDUNITS
  else {
    Units::Unit from (unit);
    if (Units::areConvertible(from, "m"))
      fov = Units::ToMeters(fov, from) / distance_;
    else fov = Units::Converter(from, "rad")(fov);
  }
# else
  else if (unit=="degree" || unit=="°")        fov *= GYOTO_DEGRAD;
  else if (unit=="arcmin")                     fov *= GYOTO_MINRAD;
  else if (unit=="arcsec" || unit=="as")       fov *= GYOTO_SECRAD;
  else if (unit=="milliarcsec" || unit=="mas") fov *= GYOTO_MASRAD;
  else if (unit=="microarcsec" || unit=="µas" || unit=="uas")
                                               fov *= GYOTO_MUASRAD;
  else {
    stringstream ss;
    ss << "Screen::fieldOfView(): unknown unit: \"" << unit << "\""
       << " (you may have more chance compiling gyoto with --with-udunits)";
    GYOTO_ERROR(ss.str());
  }
# endif
  fieldOfView(fov);
}
void Screen::fieldOfView(double fov) { fov_ = fov; }

double Screen::azimuthalFieldOfView() const {return azimuthal_fov_;}
void Screen::azimuthalFieldOfView(double fov) {azimuthal_fov_ = fov;}

void Screen::dangle1(double aa) { dangle1_ = aa; }
double Screen::dangle1() const { return dangle1_; }
void Screen::dangle2(double bb) { dangle2_ = bb; }
double Screen::dangle2() const { return dangle2_; }

double Screen::dangle1(string const &unit) const {
  double fov = dangle1();
  if (unit=="" || unit=="rad") ;
  else if (unit=="geometrical") fov *= distance_ / gg_ -> unitLength();
# ifdef HAVE_UDUNITS
  else if (Units::areConvertible(unit, "m"))
    fov = Units::FromMeters(fov*distance_, unit) ;
  else fov = Units::Converter("rad", unit)(fov);
# else
  else if (unit=="degree") fov /= GYOTO_DEGRAD;
  else if (unit=="arcmin") fov /= GYOTO_MINRAD;
  else if (unit=="arcsec") fov /= GYOTO_SECRAD;
  else if (unit=="milliarcsec") fov /= GYOTO_MASRAD;
  else if (unit=="microarcsec") fov /= GYOTO_MUASRAD;
  else {
    stringstream ss;
    ss << "Screen::dangle1(): unknown unit: \"" << unit << "\""
       << " (you may have more chance compiling gyoto with --with-udunits)";
    GYOTO_ERROR(ss.str());
  }
# endif
  return fov;
}

double Screen::dangle2(string const &unit) const{
  double fov = dangle2();
  if (unit=="" || unit=="rad") ;
  else if (unit=="geometrical") fov *= distance_ / gg_ -> unitLength();
# ifdef HAVE_UDUNITS
  else if (Units::areConvertible(unit, "m"))
    fov = Units::FromMeters(fov*distance_, unit) ;
  else fov = Units::Converter("rad", unit)(fov);
# else
  else if (unit=="degree") fov /= GYOTO_DEGRAD;
  else if (unit=="arcmin") fov /= GYOTO_MINRAD;
  else if (unit=="arcsec") fov /= GYOTO_SECRAD;
  else if (unit=="milliarcsec") fov /= GYOTO_MASRAD;
  else if (unit=="microarcsec") fov /= GYOTO_MUASRAD;
  else {
    stringstream ss;
    ss << "Screen::dangle2(): unknown unit: \"" << unit << "\""
       << " (you may have more chance compiling gyoto with --with-udunits)";
    GYOTO_ERROR(ss.str());
  }
# endif
  return fov;
}

void Screen::dangle1(double fov, const string &unit) {
  if (unit=="" || unit=="rad") ;
  else if (unit=="geometrical") fov *= gg_ -> unitLength() / distance_ ;
# ifdef HAVE_UDUNITS
  else {
    Units::Unit from (unit);
    if (Units::areConvertible(from, "m"))
      fov = Units::ToMeters(fov, from) / distance_;
    else fov = Units::Converter(from, "rad")(fov);
  }
# else
  else if (unit=="degree" || unit=="°")        fov *= GYOTO_DEGRAD;
  else if (unit=="arcmin")                     fov *= GYOTO_MINRAD;
  else if (unit=="arcsec" || unit=="as")       fov *= GYOTO_SECRAD;
  else if (unit=="milliarcsec" || unit=="mas") fov *= GYOTO_MASRAD;
  else if (unit=="microarcsec" || unit=="µas" || unit=="uas")
                                               fov *= GYOTO_MUASRAD;
  else {
    stringstream ss;
    ss << "Screen::fieldOfView(): unknown unit: \"" << unit << "\""
       << " (you may have more chance compiling gyoto with --with-udunits)";
    GYOTO_ERROR(ss.str());
  }
# endif
  dangle1(fov);
}

void Screen::dangle2(double fov, const string &unit) {
  if (unit=="" || unit=="rad") ;
  else if (unit=="geometrical") fov *= gg_ -> unitLength() / distance_ ;
# ifdef HAVE_UDUNITS
  else {
    Units::Unit from (unit);
    if (Units::areConvertible(from, "m"))
      fov = Units::ToMeters(fov, from) / distance_;
    else fov = Units::Converter(from, "rad")(fov);
  }
# else
  else if (unit=="degree" || unit=="°")        fov *= GYOTO_DEGRAD;
  else if (unit=="arcmin")                     fov *= GYOTO_MINRAD;
  else if (unit=="arcsec" || unit=="as")       fov *= GYOTO_SECRAD;
  else if (unit=="milliarcsec" || unit=="mas") fov *= GYOTO_MASRAD;
  else if (unit=="microarcsec" || unit=="µas" || unit=="uas")
                                               fov *= GYOTO_MUASRAD;
  else {
    stringstream ss;
    ss << "Screen::fieldOfView(): unknown unit: \"" << unit << "\""
       << " (you may have more chance compiling gyoto with --with-udunits)";
    GYOTO_ERROR(ss.str());
  }
# endif
  dangle2(fov);
}

void Screen::anglekind(int kind) { anglekind_ = kind; }
void Screen::anglekind(std::string const &skind) {
  if      (skind=="EquatorialAngles") anglekind_=equatorial_angles;
  else if (skind=="SphericalAngles")  anglekind_=spherical_angles;
  else if (skind=="Rectilinear")      anglekind_=rectilinear;
  else GYOTO_ERROR("Invalid string value for anglekind_");
}

std::string Screen::anglekind() const {
  switch (anglekind_) {
  case equatorial_angles:
    return "EquatorialAngles";
    break;
  case spherical_angles:
    return "SphericalAngles";
    break;
  case rectilinear:
    return "Rectilinear";
    break;
  default:
    GYOTO_ERROR("Invalid integer value for Screen::anglekind_");
  }
  return ""; // silence warning
}

size_t Screen::resolution() const { return npix_; }
void Screen::resolution(size_t n) {
  if (mask_ && n != npix_) {
    GYOTO_INFO << "Changing resolution: deleting mask" << endl;
    delete[] mask_; mask_=NULL;
  }
  npix_ = n;
}

#ifdef HAVE_UDUNITS
void Gyoto::Screen::mapPixUnit() {
  Units::Unit radian ("radian");
  double delta = fov_/double(npix_);
  ut_unit * pix = ut_scale(delta, radian);
  ut_status status = UT_BAD_ARG;
  status = ut_map_name_to_unit("pixel", UT_ASCII, pix);
  switch (status) {
  case UT_SUCCESS:
    break;
  case UT_EXISTS:
    GYOTO_ERROR("name \"pixel\" already registered");
    break;
  default:
    GYOTO_ERROR("error initializing \"pixel\" unit");
  }
  status = ut_map_symbol_to_unit("pix", UT_ASCII, pix);
  switch (status) {
  case UT_SUCCESS:
    break;
  case UT_EXISTS:
    GYOTO_ERROR("symbol \"pix\" already registered");
    break;
  default:
    GYOTO_ERROR("error initializing \"pixel\" unit");
  }
  ut_free(pix);
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << "\"pix\" unit mapped\n";
# endif
}
void Gyoto::Screen::unmapPixUnit() {
  ut_system * sys = Units::getSystem();
  if ((ut_unmap_name_to_unit(sys, "pixel", UT_ASCII) != UT_SUCCESS) ||
      (ut_unmap_symbol_to_unit(sys, "pix", UT_ASCII) != UT_SUCCESS))
    GYOTO_ERROR("Error unmapping \"pixel\" unit");
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << "\"pix\" unit unmapped\n";
# endif
}
#endif

#ifdef GYOTO_USE_XERCES
void Screen::fillProperty(Gyoto::FactoryMessenger *fmp,
			  Property const &p) const {
  if (p.name=="Distance") {
    FactoryMessenger* child = NULL;
    double d = distance();
    if (gg_() && (gg_->mass() == 1.)) {
      d /=  gg_->unitLength();
      fmp -> setParameter ("Distance", &d, 1, &child);
      child -> setSelfAttribute("unit", "geometrical");
    } else {
      string unit = "m";
      if (     d >= GYOTO_KPC*1e3   ) { unit = "Mpc"; d = Units::FromMeters(d, unit, gg_);}
      else if (d >= GYOTO_KPC       ) { unit = "kpc"; d = Units::FromMeters(d, unit, gg_);}
      else if (d >= GYOTO_KPC*1e-3  ) { unit = "pc";  d = Units::FromMeters(d, unit, gg_);}
      else if (d >= GYOTO_LIGHT_YEAR) { unit = "ly";  d = Units::FromMeters(d, unit, gg_);}
      else if (d >= GYOTO_ASTRONOMICAL_UNIT)
	{ unit = "AU"; d = Units::FromMeters(d, unit, gg_);}
      else if (d >= GYOTO_SUN_RADIUS)
	{ unit = "sunradius"; d = Units::FromMeters(d, unit, gg_);}
      else if (d >= 1e3)              { unit = "km";  d = Units::FromMeters(d, unit, gg_);}
      else if (d >= 1.) ;//           { d *= 1.;                unit = "m";}
      else if (d >= 1e-2)             { unit = "cm";  d = Units::FromMeters(d, unit, gg_);}
      //else ;                        { d *= 1.;                unit = "m";}
      fmp -> setParameter ("Distance", &d, 1, &child);
      child -> setSelfAttribute("unit", unit);
    }
    child -> setSelfAttribute("dmax", dmax_);
    delete child; child = NULL;
  }
  else if (p.name=="DMax") ; // do nothing
  else if (p.name=="Mask" && mask_filename_!="")
    fmp->setParameter("Mask",
		      (mask_filename_.compare(0,1,"!") ?
		       mask_filename_ :
		       mask_filename_.substr(1)));
  else Object::fillProperty(fmp, p);
}


SmartPointer<Screen> Screen::Subcontractor(FactoryMessenger* fmp) {
  string name="", content="", unit="", tunit="", aunit="", dunit="";
  vector<string> plugin;
  SmartPointer<Screen> scr = new Screen();
  if (!fmp) return scr;
  scr -> metric(fmp->metric());
  int tobs_found=0;
  double tobs_tmp, pos[4] ;
  char * tc;

  // Deal with fov later as we need Inclination
  double fov; string fov_unit; int fov_found=0;
  double dangle1=0.; int dangle1_found=0; 
  double dangle2=0.; int dangle2_found=0;

  while (fmp->getNextParameter(&name, &content, &unit)) {
    tc = const_cast<char*>(content.c_str());
#   ifdef GYOTO_DEBUG_ENABLED
    GYOTO_IF_DEBUG
    GYOTO_DEBUG_EXPR(name);
    GYOTO_DEBUG_EXPR(content);
    GYOTO_DEBUG_EXPR(unit);
    GYOTO_ENDIF_DEBUG
#   endif
    if      (name=="Time")
      {tobs_tmp = Gyoto::atof(tc); tunit=unit; tobs_found=1;}
    else if (name=="Position") {
      if (FactoryMessenger::parseArray(content, pos, 4) != 4)
	GYOTO_ERROR("Screen \"Position\" requires exactly 4 tokens");
      scr -> setObserverPos (pos); 
    }
    else if (name=="Distance")    
      {
	scr -> distance    ( Gyoto::atof(tc), unit );
	string dmax = fmp -> getAttribute("dmax");
	if (dmax != "") scr -> dMax(Gyoto::atof(dmax.c_str()));
      }
    else if (name=="FieldOfView") {
      fov = Gyoto::atof(tc); fov_unit=unit; fov_found=1;
    }
    else if (name=="Spectrometer") {
      scr ->
	spectrometer((Spectrometer::getSubcontractor(fmp->getAttribute("kind"),
						     plugin))
		     (fmp->getChild(), plugin));
    }
    else if (name=="Dangle1"){
      dangle1 = Gyoto::atof(tc); dangle1_found=1; aunit=unit;
    }
    else if (name=="Dangle2"){
      dangle2 = Gyoto::atof(tc); dangle2_found=1; dunit=unit;
    }
    else if (name=="SphericalAngles" ||
	     name=="EquatorialAngles" ||
	     name=="Rectilinear")  scr -> anglekind(name);
    else if (name=="Mask")
      scr -> maskFile(content==""?"":fmp->fullPath(content));
    else if (scr->setParameter(name, content, unit))
      GYOTO_ERROR("no such Screen Property");
  }

  // Must be processed after Position
  if (tobs_found) scr -> time ( tobs_tmp, tunit );

  // Must be processed after Position and Distance so pix_unit is
  // defined
  if (fov_found) scr -> fieldOfView ( fov, fov_unit );
  if (dangle1_found) scr -> dangle1(dangle1, aunit);
  if (dangle2_found) scr -> dangle2(dangle2, dunit);

# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG_EXPR(scr->dMax());
# endif

  return scr;
}
#endif

//////// COORD2DSET

Screen::Coord2dSet::Coord2dSet(CoordType_e k) : kind(k) {}
Screen::Coord2dSet::~Coord2dSet() {}

GYOTO_ARRAY<size_t, 2> Screen::Coord2dSet::operator* () const {
  if (kind==pixel)
    GYOTO_ERROR("BUG: Coord2dSet of kind pixel should implement operator*");
  else
    GYOTO_ERROR("Coord2dSet of kind angle should not be dereferenced");
  // avoid warning
  GYOTO_ARRAY<size_t, 2> a;
  return a;
}

GYOTO_ARRAY<double, 2> Screen::Coord2dSet::angles () const {
  if (kind==Screen::angle)
    GYOTO_ERROR("BUG: Coord2dSet of kind angle should implement angles()");
  else
    GYOTO_ERROR("angles() should not be called on Coord2dSet of kind pixel");
  // avoid warning
  GYOTO_ARRAY<double, 2> a;
  return a;
}

//////// GRID

Screen::Grid::Grid(Coord1dSet &iset, Coord1dSet &jset,
		   const char * const p)
  : Coord2dSet(pixel),
    prefix_(NULL),
    iset_(iset), jset_(jset)
{
  if (p) {
    size_t sz=strlen(p)+1;
    prefix_ = new char[sz];
    memcpy(prefix_, p, sz*sizeof(char));
  }
}
Screen::Grid::~Grid(){
  if (prefix_) delete[] prefix_;
}

GYOTO_ARRAY<size_t, 2> Screen::Grid::operator* () const {
#if defined HAVE_BOOST_ARRAY_HPP
  GYOTO_ARRAY<size_t, 2> ij = {*iset_, *jset_};
#else
  GYOTO_ARRAY<size_t, 2> ij;
  ij[0]=*iset_;
  ij[1]=*jset_;
#endif
  return ij;
}
void Screen::Grid::begin() {iset_.begin(); jset_.begin();}
bool Screen::Grid::valid() {return iset_.valid() && jset_.valid();}
size_t Screen::Grid::size() {return (iset_.size())*(jset_.size());}

Screen::Coord2dSet& Screen::Grid::operator++() {
  if (!valid()) return *this;
  ++iset_;
  if (!iset_.valid()) {
    iset_.begin();
    ++jset_;
    if (prefix_ && verbose() >= GYOTO_QUIET_VERBOSITY && jset_.valid())
      cout << prefix_ << jset_.index()+1 << "/" << jset_.size() << flush;
  }

  return *this;
}

//////// COORD1DSET

Screen::Coord1dSet::Coord1dSet(CoordType_e k) : kind(k) {}
Screen::Coord1dSet::~Coord1dSet() {}

size_t Screen::Coord1dSet::operator* () const {
  if (kind==pixel)
    GYOTO_ERROR("BUG: Coord1dSet of kind pixel should implement operator*");
  else
    GYOTO_ERROR("Coord1dSet of kind angle should not be dereferenced");
  // avoid warning
  return 0;
}

double Screen::Coord1dSet::angle () const {
  if (kind==Screen::angle)
    GYOTO_ERROR("BUG: Coord1dSet of kind angle should implement angle()");
  else
    GYOTO_ERROR("angle() should not be called on Coord1dSet of kind pixel");
  // avoid warning
  return 0.;
}

///////

Screen::Range::Range(size_t mi, size_t ma, size_t d)
  : Coord1dSet(pixel), mi_(mi), ma_(ma), d_(d), sz_((ma-mi)/d+1), cur_(mi)
{}

void Screen::Range::begin() {cur_=mi_;}
Screen::Coord1dSet& Screen::Range::operator++() {
  cur_ += d_; return *this;
}
bool Screen::Range::valid() {return cur_ <= ma_;}
size_t Screen::Range::size() {return sz_;}
size_t Screen::Range::operator*() const {return cur_;}
size_t Screen::Range::index() const {return (cur_-mi_) / d_;}

//////

Screen::Indices::Indices (size_t const*const buf, size_t sz)
  : Coord1dSet(pixel), indices_(NULL), sz_(sz), i_(0)
{
  if (!buf) return;
  indices_ = new size_t[sz];
  memcpy(indices_, buf, sz*sizeof(size_t));
}
Screen::Indices::~Indices(){
  if (!indices_) return;
  delete[] indices_;
}
void Screen::Indices::begin() {i_=0;}
bool Screen::Indices::valid() {return i_ < sz_;}
size_t Screen::Indices::size(){return sz_;}
Screen::Coord1dSet& Screen::Indices::operator++() {++i_; return *this;}
size_t Screen::Indices::operator*() const {return indices_[i_];}
size_t Screen::Indices::index() const {return i_;}


/////

Screen::Angles::Angles (double const*const buf, size_t sz)
  : Coord1dSet(Screen::angle), buf_(NULL), sz_(sz), i_(0)
{
  if (!buf) return;
  buf_ = new double[sz];
  memcpy(buf_, buf, sz*sizeof(double));
}
Screen::Angles::~Angles(){
  if (!buf_) return;
  delete[] buf_;
}
void Screen::Angles::begin() {i_=0;}
bool Screen::Angles::valid() {return i_<sz_;}
size_t Screen::Angles::size(){return sz_;}
Screen::Coord1dSet& Screen::Angles::operator++(){++i_; return *this;}
double Screen::Angles::angle() const {return buf_[i_];}
size_t Screen::Angles::index() const {return i_;}

Screen::RepeatAngle::RepeatAngle (double val, size_t sz)
  : Coord1dSet(Screen::angle), val_(val), sz_(sz), i_(0)
{}
void Screen::RepeatAngle::begin() {i_=0;}
bool Screen::RepeatAngle::valid() {return i_<sz_;}
size_t Screen::RepeatAngle::size(){return sz_;}
Screen::Coord1dSet& Screen::RepeatAngle::operator++(){++i_; return *this;}
double Screen::RepeatAngle::angle() const {return val_;}
size_t Screen::RepeatAngle::index() const {return i_;}

Screen::Bucket::Bucket (Coord1dSet &alp, Coord1dSet &del)
  : Coord2dSet(alp.kind), alpha_(alp), delta_(del)
{
  if (alp.kind != del.kind)
    GYOTO_ERROR("both specifiers must be of same kind");
  if (alp.size() != del.size())
    GYOTO_ERROR("alpha and delta should be of same size"); 
}
void Screen::Bucket::begin() {alpha_.begin(); delta_.begin();}
bool Screen::Bucket::valid() {return alpha_.valid() && delta_.valid();}
size_t Screen::Bucket::size(){return alpha_.size();}
Screen::Coord2dSet& Screen::Bucket::operator++(){
  ++alpha_; ++delta_; return *this;}
GYOTO_ARRAY<double, 2> Screen::Bucket::angles() const {
#if defined HAVE_BOOST_ARRAY_HPP
  GYOTO_ARRAY<double, 2> out {alpha_.angle(), delta_.angle()};
#else
  GYOTO_ARRAY<double, 2> out;
  out[0]=alpha_.angle();
  out[1]=delta_.angle();
#endif
  return out;
}
GYOTO_ARRAY<size_t, 2> Screen::Bucket::operator* () const {
#if defined HAVE_BOOST_ARRAY_HPP
  GYOTO_ARRAY<size_t, 2> ij = {*alpha_, *delta_};
#else
  GYOTO_ARRAY<size_t, 2> ij;
  ij[0]=*alpha_;
  ij[1]=*delta_;
#endif
  return ij;
}

///// Empty

Screen::Empty::Empty () : Coord2dSet(pixel) {}
Screen::Coord2dSet& Screen::Empty::operator++() {return *this;}
void Screen::Empty::begin() {}
bool Screen::Empty::valid() {return false;}
size_t Screen::Empty::size() {return 0;}
