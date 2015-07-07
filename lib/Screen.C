/*
    Copyright 2011-2015 Thibaut Paumard, Frederic Vincent

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
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <string>
#include <sstream>

#ifdef GYOTO_USE_CFITSIO
#include <fitsio.h>
#define throwCfitsioError(status) \
    { fits_get_errstatus(status, ermsg); throwError(ermsg); }
#endif

#ifdef HAVE_BOOST_MULTIPRECISION_CPP_DEC_FLOAT_HPP
#include <boost/multiprecision/cpp_dec_float.hpp>
#endif

using namespace std ; 
using namespace Gyoto;

/// Properties
#include "GyotoProperty.h"
GYOTO_PROPERTY_START(Screen)
GYOTO_PROPERTY_METRIC(Screen, Metric, metric)
GYOTO_PROPERTY_SPECTROMETER(Screen, Spectrometer, spectrometer)
GYOTO_PROPERTY_DOUBLE_UNIT(Screen, Time, time)
GYOTO_PROPERTY_DOUBLE_UNIT(Screen, FieldOfView, fieldOfView)
GYOTO_PROPERTY_DOUBLE_UNIT(Screen, PALN, PALN)
GYOTO_PROPERTY_DOUBLE_UNIT(Screen, Inclination, inclination)
GYOTO_PROPERTY_DOUBLE_UNIT(Screen, Argument, argument)
GYOTO_PROPERTY_DOUBLE_UNIT(Screen, Alpha0, alpha0)
GYOTO_PROPERTY_DOUBLE_UNIT(Screen, Delta0, delta0)
GYOTO_PROPERTY_SIZE_T(Screen, Resolution, resolution)
GYOTO_PROPERTY_DOUBLE_UNIT(Screen, Distance, distance)
GYOTO_PROPERTY_DOUBLE(Screen, DMax, dMax)
GYOTO_PROPERTY_DOUBLE_UNIT(Screen, FreqObs, freqObs)
GYOTO_PROPERTY_STRING(Screen, AngleKind, anglekind)
GYOTO_PROPERTY_STRING(Screen, ObserverKind, observerKind)
GYOTO_PROPERTY_FILENAME(Screen, Mask, maskFile)
GYOTO_PROPERTY_VECTOR_DOUBLE(Screen, FourVelocity, fourVel)
GYOTO_PROPERTY_VECTOR_DOUBLE(Screen, ScreenVector1, screenVector1)
GYOTO_PROPERTY_VECTOR_DOUBLE(Screen, ScreenVector2, screenVector2)
GYOTO_PROPERTY_VECTOR_DOUBLE(Screen, ScreenVector3, screenVector3)
GYOTO_PROPERTY_END(Screen, Object::properties)
///

// Default constructor
Screen::Screen() : 
  //tobs_(0.), fov_(M_PI*0.1), tmin_(0.), npix_(1001),
  tobs_(0.), fov_(M_PI*0.1), npix_(1001), mask_(NULL), mask_filename_(""),
  distance_(1.), dmax_(GYOTO_SCREEN_DMAX), anglekind_(equatorial_angles),
  alpha0_(0.), delta0_(0.),
  gg_(NULL), spectro_(NULL),
  freq_obs_(1.),
  observerkind_("ObserverAtInfinity")
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
  tobs_(o.tobs_), fov_(o.fov_), npix_(o.npix_), mask_(NULL),
  mask_filename_(o.mask_filename_),
  distance_(o.distance_),
  dmax_(o.dmax_),
  anglekind_(o.anglekind_),
  alpha0_(o.alpha0_), delta0_(o.delta0_),
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
    throwError("Incompatible coordinate kind in Screen::setObserverCoord");
  }
  computeBaseVectors();
}

void Screen::observerKind(const string &kind) {
  observerkind_=kind;
}
string Screen::observerKind() const {
  return observerkind_;
}

void Screen::setFourVel(const double coord[4]) {
  for (int ii=0;ii<4;ii++)
    fourvel_[ii]=coord[ii];
}

void Screen::fourVel(std::vector<double> const &coord) {
  if (coord.size() != 4)
    throwError("base screen vectors require 4 elements");
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
    throwError("base screen vectors require 4 elements");
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
    throwError("base screen vectors require 4 elements");
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
    throwError("base screen vectors require 4 elements");
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
    throwError(ss.str());
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

void Screen::getRayCoord(const size_t i, const size_t j, double coord[]) const {
  double xscr, yscr;
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << "(i=" << i << ", j=" << j << ", coord)" << endl;
# endif
  switch (anglekind_) {
  case Screen::spherical_angles:
    /*
      GYOTO screen labelled by spherical
      angles a and b (see Fig. in user guide)
     */
    xscr = double(i-1)*fov_/(2.*double(npix_-1));
    yscr = M_PI-(double(j-1)*2.*M_PI/double(npix_-1));
    // NB: here xscr and yscr are the spherical angles
    // a and b ; the b->pi-b transformation boils down
    // to performing X->-X, just as below for equat angles.
    break;
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
      throwError("Rectilinear projection requires fov_ < M_PI");
    const double xfov=2.*tan(fov_*0.5);
    const double delta= xfov/double(npix_);
    yscr=delta*(double(j)-double(npix_+1)/2.);
    xscr=-delta*(double(i)-double(npix_+1)/2.);
    break;
  }
  default:
    xscr=yscr=0.;
    throwError("Unrecognized anglekind_");
  }
  getRayCoord(xscr, yscr, coord); 
}

void Screen::getRayCoord(double alpha, double delta,
		      double coord[]) const

{
  alpha+=alpha0_; delta+=delta0_; // Screen orientation
  double normtol=1e-10;
  int i; // dimension : 0, 1, 2
  double pos[4];
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << "(alpha="<<alpha<<",delta="<<delta<<",coord)" << endl;
# endif
  getObserverPos(coord);

  if (coord[1] > dmax_) {
    // scale
    coord[0] -= coord[1] - dmax_;
    double scale = coord[1] / dmax_;
    coord[1] = dmax_;
    alpha *= scale;
    delta *= scale;
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
    spherical_angle_a = alpha;
    spherical_angle_b = delta;
    break;
  case rectilinear:
    spherical_angle_a = atan(sqrt(alpha*alpha+delta*delta));
    spherical_angle_b = atan2(delta, alpha);
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
	alpha100=alpha, delta100=delta, a, b;
      a=acos(cos(alpha100)*cos(delta100));
      b=atan2(tan(delta100),sin(alpha100));
      spherical_angle_a=a.convert_to<double>();
      spherical_angle_b=b.convert_to<double>();
    }
#else
    if (abs(alpha)<1e-6 || abs(delta) < 1e-6) {
      spherical_angle_a = sqrt(alpha*alpha+delta*delta);
    } else {
      spherical_angle_a = acos(cos(alpha)*cos(delta));
    }
    spherical_angle_b =
		   (alpha==0. && delta==0.) ? 0. : atan2(tan(delta),sin(alpha));
#endif
    break;
  default:
    spherical_angle_a=spherical_angle_b=0.;
    throwError("Unknown angle type");
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
  
  // 4-vector tangent to photon geodesic
  
  if (fourvel_[0]==0. && observerkind_=="ObserverAtInfinity"){
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
	  throwError("In Screen::getRayCoord: "
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
	  throwError("Please move Screen away from z-axis");

	double rr=coord[1],
	  rinf=20.; // this rinf is just a very crude test
               	    // I take rinf=10*rhor in Sch metric
	if (rr<rinf)
	  throwError("In Screen::getRayCoord: "
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
      throwError("Incompatible coordinate kind in Screen::getRayCoord()");
      break;
    }

    // 0-component of photon tangent 4-vector found by normalizing
    gg_ -> nullifyCoord(coord);
    
  }else{    
    /* 
       ---> Observer local frame given in XML <---
       Express photon tangent 4-vector in the observer basis
       Treatment is coordinate independent 
       (except for z-axis check right below)
    */

    if (gg_ -> coordKind() == GYOTO_COORDKIND_SPHERICAL){
      if (coord[2]==0. || coord[2]==M_PI)
	throwError("Please move Screen away from z-axis");
    }

    if (fourvel_[0]!=0. && observerkind_!="ObserverAtInfinity"){
      throwError("In Screen:getRayCoord: "
		 " choose an implemented observer kind OR"
		 " explicitly give the local tetrad in the XML");
    }

    if (fourvel_[0]==0){
      // Implemented observer specifid in XML, local tetrad computed by Metric
      const double fourpos[4]={coord[0],coord[1],coord[2],coord[3]};
      double fourvel[4], screen1[4], screen2[4], screen3[4];
      gg_ -> observerTetrad(observerkind_,fourpos,fourvel,screen1,
			    screen2,screen3);

      /*cout << "Vectors in Screen: " << setprecision(17) << endl;
      for (int ii=0;ii<4;ii++) cout << fourvel[ii] << " ";
      cout << endl;
      for (int ii=0;ii<4;ii++) cout << screen1[ii] << " ";
      cout << endl;
      for (int ii=0;ii<4;ii++) cout << screen2[ii] << " ";
      cout << endl;
      for (int ii=0;ii<4;ii++) cout << screen3[ii] << " ";
      cout << endl;*/

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
	+fourvel_[2];
      coord[7]=vel[0]*screen1[3]
	+vel[1]*screen2[3]
	+vel[2]*screen3[3]
	+fourvel[3];
    }else{
      // Local tetrad given by the user in the XML file. Check it.
      if (fabs(gg_->ScalarProd(coord,fourvel_,fourvel_)+1.)>normtol ||
	  fabs(gg_->ScalarProd(coord,screen1_,screen1_)-1.)>normtol ||
	  fabs(gg_->ScalarProd(coord,screen2_,screen2_)-1.)>normtol ||
	  fabs(gg_->ScalarProd(coord,screen3_,screen3_)-1.)>normtol){
	cout << "norm= " << gg_->ScalarProd(coord,fourvel_,fourvel_) << " " << gg_->ScalarProd(coord,screen1_,screen1_) << " " << gg_->ScalarProd(coord,screen2_,screen2_) << " " << gg_->ScalarProd(coord,screen3_,screen3_) << endl;
	throwError("In Screen:getRayCoord: observer's local"
		   " basis is not properly normalized");
      }
      
      if (fabs(gg_->ScalarProd(coord,fourvel_,screen1_))>normtol ||
	  fabs(gg_->ScalarProd(coord,fourvel_,screen2_))>normtol ||
	  fabs(gg_->ScalarProd(coord,fourvel_,screen3_))>normtol ||
	  fabs(gg_->ScalarProd(coord,screen1_,screen2_))>normtol ||
	  fabs(gg_->ScalarProd(coord,screen1_,screen3_))>normtol ||
	  fabs(gg_->ScalarProd(coord,screen2_,screen3_))>normtol)
	throwError("In Screen:getRayCoord: observer's local"
		   " basis is not orthogonal");
     
      coord[4]=vel[0]*screen1_[0]
	+vel[1]*screen2_[0]
	+vel[2]*screen3_[0]
	+fourvel_[0];
      coord[5]=vel[0]*screen1_[1]
	+vel[1]*screen2_[1]
	+vel[2]*screen3_[1]
	+fourvel_[1];
      coord[6]=vel[0]*screen1_[2]
	+vel[1]*screen2_[2]
	+vel[2]*screen3_[2]
	+fourvel_[2];
      coord[7]=vel[0]*screen1_[3]
	+vel[1]*screen2_[3]
	+vel[2]*screen3_[3]
	+fourvel_[3];
    }
    if (fabs(gg_->ScalarProd(coord,coord+4,coord+4))>normtol){
      throwError("In Screen::getRayCoord: "
		 " tangent 4-vector to photon not properly normalized");
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
    throwError("Screen::fitsReadMask(): mask must be square");
  //  if (naxes[2] > 1)
  //    throwError("Screen::fitsReadMask(): mask must have only one plane");
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
  if (filename=="") throwError("no filename specified");
  if (!mask_) throwError("Screen::fitsWriteMask(filename): nothing to save!");
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
  if ( i<=0 || i> npix_ || j>npix_ || j<=0) throwError("wrong index");
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

void Screen::coordToSky(const double pos[4], double skypos[3]) const {
  double xyz[3];
  coordToXYZ(pos, xyz);
  double ul = gg_ -> unitLength();

  skypos[0]=(xyz[0]*ex_[0]+xyz[1]*ey_[0]+xyz[2]*ez_[0]) * ul;
  skypos[1]=(xyz[0]*ex_[1]+xyz[1]*ey_[1]+xyz[2]*ez_[1]) * ul;
  skypos[2]=(xyz[0]*ex_[2]+xyz[1]*ey_[2]+xyz[2]*ez_[2]) * ul;
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
    throwError("Incompatible coordinate kind in Screen::coordToXYZ");
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
    throwError(ss.str());
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
    throwError(ss.str());
  }
# endif
  fieldOfView(fov);
}
void Screen::fieldOfView(double fov) { fov_ = fov; }

void Screen::alpha0(double alpha) { alpha0_ = alpha; }
double Screen::alpha0() const { return alpha0_; }
void Screen::delta0(double delta) { delta0_ = delta; }
double Screen::delta0() const { return delta0_; }

double Screen::alpha0(string const &unit) const {
  double fov = alpha0();
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
    ss << "Screen::alpha0(): unknown unit: \"" << unit << "\""
       << " (you may have more chance compiling gyoto with --with-udunits)";
    throwError(ss.str());
  }
# endif
  return fov;
}

double Screen::delta0(string const &unit) const{
  double fov = delta0();
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
    ss << "Screen::delta0(): unknown unit: \"" << unit << "\""
       << " (you may have more chance compiling gyoto with --with-udunits)";
    throwError(ss.str());
  }
# endif
  return fov;
}

void Screen::alpha0(double fov, const string &unit) {
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
    throwError(ss.str());
  }
# endif
  alpha0(fov);
}

void Screen::delta0(double fov, const string &unit) {
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
    throwError(ss.str());
  }
# endif
  delta0(fov);
}

void Screen::anglekind(int kind) { anglekind_ = kind; }
void Screen::anglekind(std::string const &skind) {
  if      (skind=="EquatorialAngles") anglekind_=equatorial_angles;
  else if (skind=="SphericalAngles")  anglekind_=spherical_angles;
  else if (skind=="Rectilinear")      anglekind_=rectilinear;
  else throwError("Invalid string value for anglekind_");
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
    throwError("Invalid integer value for Screen::anglekind_");
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
    throwError("name \"pixel\" already registered");
    break;
  default:
    throwError("error initializing \"pixel\" unit");
  }
  status = ut_map_symbol_to_unit("pix", UT_ASCII, pix);
  switch (status) {
  case UT_SUCCESS:
    break;
  case UT_EXISTS:
    throwError("symbol \"pix\" already registered");
    break;
  default:
    throwError("error initializing \"pixel\" unit");
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
    throwError("Error unmapping \"pixel\" unit");
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
      if (     d >= GYOTO_KPC*1e3   ) { d /= (1000.*GYOTO_KPC); unit = "Mpc"; }
      else if (d >= GYOTO_KPC       ) { d /= GYOTO_KPC;         unit = "kpc"; }
      else if (d >= GYOTO_KPC*1e-3  ) { d /= (GYOTO_KPC*1e-3);  unit = "pc"; }
      else if (d >= GYOTO_LIGHT_YEAR) { d /= GYOTO_LIGHT_YEAR;  unit = "ly"; }
      else if (d >= GYOTO_ASTRONOMICAL_UNIT)
	{ d /= GYOTO_ASTRONOMICAL_UNIT; unit = "AU"; }
      else if (d >= GYOTO_SUN_RADIUS)
	{ d /= GYOTO_SUN_RADIUS;  unit = "sunradius"; }
      else if (d >= 1e3)              { d *= 1e-3;              unit = "km";}
      else if (d >= 1.) ;//           { d *= 1.;                unit = "m";}
      else if (d >= 1e-2)             { d *= 1e2;               unit = "cm";}
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
  SmartPointer<Screen> scr = new Screen();
  scr -> metric(fmp->metric());
  int tobs_found=0;
  double tobs_tmp, pos[4] ;
  char * tc;

  // Deal with fov later as we need Inclination
  double fov; string fov_unit; int fov_found=0;
  double alpha0=0.; int alpha0_found=0; 
  double delta0=0.; int delta0_found=0;

  while (fmp->getNextParameter(&name, &content, &unit)) {
    tc = const_cast<char*>(content.c_str());
#   ifdef GYOTO_DEBUG_ENABLED
    GYOTO_IF_DEBUG
    GYOTO_DEBUG_EXPR(name);
    GYOTO_DEBUG_EXPR(content);
    GYOTO_DEBUG_EXPR(unit);
    GYOTO_ENDIF_DEBUG
#   endif
    if      (name=="Time")     {tobs_tmp = atof(tc); tunit=unit; tobs_found=1;}
    else if (name=="Position") {
      if (FactoryMessenger::parseArray(content, pos, 4) != 4)
	throwError("Screen \"Position\" requires exactly 4 tokens");
      scr -> setObserverPos (pos); 
    }
    else if (name=="KeplerianObserver" ||
	     name=="ZAMO" ||
	     name=="ObserverAtInfinity") {
      scr -> observerKind(name);
    }
    else if (name=="Distance")    
      {
	scr -> distance    ( atof(tc), unit );
	string dmax = fmp -> getAttribute("dmax");
	if (dmax != "") scr -> dMax(atof(dmax.c_str()));
      }
    else if (name=="FieldOfView") {
      fov = atof(tc); fov_unit=unit; fov_found=1;
    }
    else if (name=="Spectrometer") {
      scr -> spectrometer ((Spectrometer::getSubcontractor(fmp->getAttribute("kind")))(fmp->getChild()));
    }
    else if (name=="Alpha0"){
      alpha0 = atof(tc); alpha0_found=1; aunit=unit;
    }
    else if (name=="Delta0"){
      delta0 = atof(tc); delta0_found=1; dunit=unit;
    }
    else if (name=="SphericalAngles" ||
	     name=="EquatorialAngles" ||
	     name=="Rectilinear")  scr -> anglekind(name);
    else if (name=="Mask")
      scr -> maskFile(content==""?"":fmp->fullPath(content));
    else if (scr->setParameter(name, content, unit))
      throwError("no such Screen Property");
  }

  // Must be processed after Position
  if (tobs_found) scr -> time ( tobs_tmp, tunit );

  // Must be processed after Position and Distance so pix_unit is
  // defined
  if (fov_found) scr -> fieldOfView ( fov, fov_unit );
  if (alpha0_found) scr -> alpha0(alpha0, aunit);
  if (delta0_found) scr -> delta0(delta0, dunit);

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
    throwError("BUG: Coord2dSet of kind pixel should implement operator*");
  else
    throwError("Coord2dSet of kind angle should not be dereferenced");
  // avoid warning
  GYOTO_ARRAY<size_t, 2> a;
  return a;
}

GYOTO_ARRAY<double, 2> Screen::Coord2dSet::angles () const {
  if (kind==Screen::angle)
    throwError("BUG: Coord2dSet of kind angle should implement angles()");
  else
    throwError("angles() should not be called on Coord2dSet of kind pixel");
  // avoid warning
  GYOTO_ARRAY<double, 2> a;
  return a;
}

//////// GRID

Screen::Grid::Grid(Coord1dSet &iset, Coord1dSet &jset,
		   const char * const p)
  : Coord2dSet(pixel),
    prefix_(p),
    iset_(iset), jset_(jset)
{}

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
      cout << prefix_ << *jset_ << "/" << jset_.size() << flush;
  }

  return *this;
}

//////// COORD1DSET

Screen::Coord1dSet::Coord1dSet(CoordType_e k) : kind(k) {}
Screen::Coord1dSet::~Coord1dSet() {}

size_t Screen::Coord1dSet::operator* () const {
  if (kind==pixel)
    throwError("BUG: Coord1dSet of kind pixel should implement operator*");
  else
    throwError("Coord1dSet of kind angle should not be dereferenced");
  // avoid warning
  return 0;
}

double Screen::Coord1dSet::angle () const {
  if (kind==Screen::angle)
    throwError("BUG: Coord1dSet of kind angle should implement angle()");
  else
    throwError("angle() should not be called on Coord1dSet of kind pixel");
  // avoid warning
  return 0.;
}

///////

Screen::Range::Range(size_t mi, size_t ma, size_t d)
  : Coord1dSet(pixel), mi_(mi), ma_(ma), d_(d), sz_((ma-mi+1)/d), cur_(mi)
{}

void Screen::Range::begin() {cur_=mi_;}
Screen::Coord1dSet& Screen::Range::operator++() {
  cur_ += d_; return *this;
}
bool Screen::Range::valid() {return cur_ <= ma_;}
size_t Screen::Range::size() {return sz_;}
size_t Screen::Range::operator*() const {return cur_;}

//////

Screen::Indices::Indices (size_t const*const buf, size_t sz)
  : Coord1dSet(pixel), indices_(buf), sz_(sz), i_(0)
{}
void Screen::Indices::begin() {i_=0;}
bool Screen::Indices::valid() {return i_ < sz_;}
size_t Screen::Indices::size(){return sz_;}
Screen::Coord1dSet& Screen::Indices::operator++() {++i_; return *this;}
size_t Screen::Indices::operator*() const {return indices_[i_];}


/////

Screen::Angles::Angles (double const*const buf, size_t sz)
  : Coord1dSet(Screen::angle), buf_(buf), sz_(sz), i_(0)
{}
void Screen::Angles::begin() {i_=0;}
bool Screen::Angles::valid() {return i_<sz_;}
size_t Screen::Angles::size(){return sz_;}
Screen::Coord1dSet& Screen::Angles::operator++(){++i_; return *this;}
double Screen::Angles::angle() const {return buf_[i_];}

Screen::RepeatAngle::RepeatAngle (double val, size_t sz)
  : Coord1dSet(Screen::angle), val_(val), sz_(sz), i_(0)
{}
void Screen::RepeatAngle::begin() {i_=0;}
bool Screen::RepeatAngle::valid() {return i_<sz_;}
size_t Screen::RepeatAngle::size(){return sz_;}
Screen::Coord1dSet& Screen::RepeatAngle::operator++(){++i_; return *this;}
double Screen::RepeatAngle::angle() const {return val_;}

Screen::Bucket::Bucket (Coord1dSet &alp, Coord1dSet &del)
  : Coord2dSet(alp.kind), alpha_(alp), delta_(del)
{
  if (alp.kind != del.kind)
    throwError("both specifiers must be of same kind");
  if (alp.size() != del.size())
    throwError("alpha and delta should be of same size"); 
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
