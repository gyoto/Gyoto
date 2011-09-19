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
#include "GyotoFactoryMessenger.h"
#include "GyotoScreen.h"

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <string>
#include <sstream>

using namespace std ; 
using namespace Gyoto;

// Default constructor
Screen::Screen() : 
  //tobs_(0.), fov_(M_PI*0.1), tmin_(0.), npix_(1001),
  tobs_(0.), fov_(M_PI*0.1), npix_(1001),
  distance_(1.), dmax_(GYOTO_SCREEN_DMAX), gg_(NULL), spectro_(NULL)
{
  if (debug()) cerr << "DEBUG: Screen::Screen()" << endl ;
  if (debug()) cerr << "DEBUG: Screen::dmax_ = " << dmax_ << endl ;
  euler_[0]=euler_[1]=euler_[2]=0.;
  setProjection(0., 0., 0.);
}

Screen::Screen(const Screen& o) :
  SmartPointee(o),
  tobs_(o.tobs_), fov_(o.fov_), npix_(o.npix_), distance_(o.distance_),
  dmax_(o.dmax_), gg_(NULL), spectro_(NULL)
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
}
Screen * Screen::clone() const { return new Screen(*this); }

Screen::~Screen(){}

std::ostream& Screen::print( std::ostream& o) const {
  o << "distance="    << distance_ << ", " ;
  o << "paln="        << euler_[0] << ", " ;
  o << "inclination=" << euler_[1] << ", " ;
  o << "argument="    << euler_[2] ;
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

void Screen::setDistance(double dist, const string unit)    {
  if ((unit=="") || (unit=="m"))     ; // default is SI
  else if (unit=="geometrical") dist *= gg_ -> unitLength();
  else if (unit=="cm")          dist *= 1e-2;
  else if (unit=="km")          dist *= 1e3;
  else if (unit=="sunradius")   dist *= GYOTO_SUN_RADIUS;
  else if ((unit=="astronomicalunit") ||
	   (unit=="AU") ||
	   (unit=="au") ||
	   (unit=="ua"))
                                dist *= GYOTO_ASTRONOMICAL_UNIT;
  else if (unit=="ly")          dist *= GYOTO_LIGHT_YEAR;
  else if (unit=="pc")          dist *= GYOTO_KPC * 1e-3;
  else if (unit=="kpc")         dist *= GYOTO_KPC;
  else if (unit=="Mpc")         dist *= GYOTO_KPC * 1e3;
  setDistance(dist);
}
void Screen::setDistance(double dist)    {
  distance_=dist;
  computeBaseVectors();
}

void Screen::setDmax(double dist) {
  if (debug()) cerr << "DEBUG: Screen::setDmax("<<dist<<")\n";
  dmax_ = dist;
}
double Screen::getDmax() const { return dmax_; }


void Screen::setPALN(double paln, const string &unit)        {
  if ((unit=="") || (unit=="rad"));
  else if (unit=="deg") paln *= GYOTO_DEGRAD;
  setPALN(paln);
}
void Screen::setPALN(double paln)        {
  euler_[0]=paln; computeBaseVectors();
}

void Screen::setInclination(double incl, const string &unit) { 
  if ((unit=="") || (unit=="rad"));
  else if (unit=="deg") incl *= GYOTO_DEGRAD;
  setInclination(incl);
}
void Screen::setInclination(double incl) {
  euler_[1]=incl; computeBaseVectors();
}

void Screen::setArgument(double arg, const string &unit) { 
  if (unit=="" || unit=="rad");
  else if (unit=="deg") arg *= GYOTO_DEGRAD;
  setArgument(arg);
}
void Screen::setArgument(double arg)     {
  euler_[2]=arg;  computeBaseVectors();
}

void Screen::setMetric(SmartPointer<Metric> gg) { gg_ = gg; computeBaseVectors(); }

int Screen::getCoordKind()      const { return gg_ -> getCoordKind(); }
double Screen::getDistance()    const { return distance_; }
double Screen::getPALN()        const { return euler_[0]; }
double Screen::getInclination() const { return euler_[1]; }
double Screen::getArgument()    const { return euler_[2]; }
SmartPointer<Metric> Screen::getMetric() const { return gg_; }

void Screen::setObserverPos(const double coord[4]) {
  tobs_ = coord[0] * gg_ -> unitLength() / GYOTO_C;
  euler_[0]=M_PI;//Par défaut, A CHANGER SI BESOIN
  //dans le cas standard ex=-ephi et ephi dirige la ligne des noeuds d'où le pi
  //NB: c'est -pi dans mes notes, donc OK [2pi]
  //NB : ne décrit que la rotation de la caméra dans le plan x,y

  int coordkind=gg_ -> getCoordKind();
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

void Screen::getObserverPos(double coord[]) const
{
  double r0     = distance_ / gg_ -> unitLength();
  //remark : Pb here if mass=0 (unitLength=0...) -> pb for flat metric

  //if (debug()) cout << "distance_ in screen= " << distance_ << endl;
  //double r0     = distance_ ;DEBUG
  double theta0 = M_PI-euler_[1];
  double phi0 = -M_PI/2-euler_[2];
    
  int coordkind = gg_ -> getCoordKind();

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

/* SPECTROMETER */

void Screen::setSpectrometer(SmartPointer<Spectrometer> spr) { spectro_=spr; }
SmartPointer<Spectrometer> Screen::getSpectrometer() const { return spectro_; }

void Screen::getRayCoord(const size_t i, const size_t j, double coord[]) const {
  const double delta= fov_/npix_;
  double xscr, yscr;
  if (debug()) cerr << "Gyoto::Screen::getRayCoord(i=" << i
		    << ", j=" << j << ", coord)" << endl;
  yscr=delta*(j-(npix_+1)/2.);
  xscr=delta*(i-(npix_+1)/2.);
  getRayCoord(-xscr, yscr, coord);
}

void Screen::getRayCoord(double alpha, double delta,
		      double coord[]) const

{
  //if (debug()) cout << "alpha,delta= " << alpha << " " << delta << endl;
  int i; // dimension : 0, 1, 2
  double pos[4];
  if (debug()) cerr << "getRayCoord(alpha="<<alpha<<",delta="<<delta<<",coord)" << endl;
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

  //3-vector of norm unity directing the photon trajectory
  
  //NB: spherical_angles = spherical angles associated to the
  //cartesian frame centered on the observer's screen,
  //x=East,y=North,z=line-of-sight ; number 1 is from z-axis to
  //vector, number 2 from x axis to projection on x-y plane. See notes
  //for details.

  //NBB: for spherical_angle_1, the relation comes from the spherical
  //law of cosines of spherical trigonometry, the arcs
  //spherical_angle_1, alpha and delta forming a spherical triangle,
  //with alpha othogonal to delta.
  double spherical_angle_1 = acos(cos(alpha)*cos(delta));
  double spherical_angle_2 = (alpha==0. && delta==0.) ? 0. : atan2(delta,alpha);
                 // C++ atan2 seems to behave badly if the two arguments are 0

  //NB: minus signs because the photon doesn't leave the screen, but
  //is heading towards it!
  double vel[3]={-sin(spherical_angle_1)*cos(spherical_angle_2),
		 -sin(spherical_angle_1)*sin(spherical_angle_2),
		 -cos(spherical_angle_1)};

  switch (gg_ -> getCoordKind()) {
  case GYOTO_COORDKIND_CARTESIAN:
    {
      //3 velocity normalization (see manual for details)
      pos[0]=coord[0];pos[1]=coord[1];pos[2]=coord[2];pos[3]=coord[3];
      double gxx=gg_->gmunu(pos,1,1), gyy=gg_->gmunu(pos,2,2), gzz=gg_->gmunu(pos,3,3), gxy=gg_->gmunu(pos,1,2), gxz=gg_->gmunu(pos,1,3), gyz=gg_->gmunu(pos,2,3);

      //Transforming to KS:
      for (i=0;i<3;++i) {
	coord[5]+=vel[i]*ex_[i];//xdot0
	coord[6]+=vel[i]*ey_[i];//ydot0
	coord[7]+=vel[i]*ez_[i];//zdot0
      }

      double vx=coord[5], vy=coord[6], vz=coord[7];
      double denom=gxx*vx*vx+gyy*vy*vy+gzz*vz*vz+2.*gxy*vx*vy+2.*gxz*vx*vz+2.*gyz*vy*vz;
      double lambda;
      double nuobs=1.;
      if (denom < 0.){
	throwError("In Screen.C: impossible to normalize in KS case!");
      }else{
	lambda=nuobs/sqrt(denom);
      }
      coord[5]*=lambda;coord[6]*=lambda;coord[7]*=lambda;
    }
    break;
  case GYOTO_COORDKIND_SPHERICAL:
    {
      pos[0]=coord[0];pos[1]=coord[1];pos[2]=coord[2];pos[3]=coord[3];
      double grr=gg_->gmunu(pos,1,1), gthth=gg_->gmunu(pos,2,2), gphph=gg_->gmunu(pos,3,3);
      coord[5]=-vel[2];
      //coord[5]=-vel[2]/sqrt(grr);
      double sp=sin(euler_[0]);
      double cp=cos(euler_[0]);
      coord[6]=(-sp*vel[0]+cp*vel[1])/coord[1];
      //coord[6]=(-sp*vel[0]+cp*vel[1])/sqrt(gthth);
      coord[7]=( cp*vel[0]+sp*vel[1])/(coord[1]*sin(coord[2]));
      //coord[7]=( cp*vel[0]+sp*vel[1])/sqrt(gphph);

      //Normalization of Boyer-Lindquist projected velocity to 1 (see manual for details). This is a choice of normalization for the photon tangent vector (i.e. 4-momentum) that boils down to taking observed frequency = 1.

      double lambda;
      double nuobs=1.;

      if (grr*coord[5]*coord[5] + gthth*coord[6]*coord[6] + gphph*coord[7]*coord[7] < 0.){
	throwError("In Screen.C: impossible normalization!");
      }else{
	lambda=nuobs/sqrt(grr*coord[5]*coord[5] + gthth*coord[6]*coord[6] + gphph*coord[7]*coord[7]);
      }
      coord[5]*=lambda;coord[6]*=lambda;coord[7]*=lambda;
    }
    break;
  default:
    throwError("Incompatible coordinate kind in Screen::getRayCoord()");
    break;
  }

  // 3-vel [dx1/dlambda,dx2/dlambda,dx3/dlambda] -> 4-vel
  gg_ -> nullifyCoord(coord);
}

void Screen::computeBaseVectors() {
  // See http://en.wikipedia.org/wiki/Euler_angles

  //  if (debug()) cout << "In compute base vectore" << endl;
  double ca = cos(euler_[0]); // paln
  double sa = sin(euler_[0]);
  double cb = cos(euler_[1]); // inclination
  double sb = sin(euler_[1]);
  double cc = cos(euler_[2]); // argument
  double sc = sin(euler_[2]);
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
  int coordkind = gg_ -> getCoordKind();
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

double Screen::getTime() { return tobs_ ; }

void Screen::setTime(double tobs, const string &unit) {
  if (unit=="" || unit=="s") ;
  else if (unit=="geometrical") tobs *= gg_ -> unitLength() / GYOTO_C ;
  else if (unit=="min") tobs *= 60. ;
  else if (unit=="h") tobs *= 3600. ;
  else if (unit=="d") tobs *= 86400. ;
  else if (unit=="y") tobs *= 3.15576e+07;
  else {
    stringstream ss;
    ss << "Screen::setTime(): unkwon unit \"" << unit << "\". Accepted units: "
       << "[s] geometrical min h d y";
    throwError (ss.str());
  }
  setTime(tobs);
}
void Screen::setTime(double tobs) { tobs_ = tobs; }

//double Screen::getMinimumTime() { return tmin_; }
//void Screen::setMinimumTime(double tmin) { tmin_ = tmin; }
double Screen::getFieldOfView() { return fov_; }

void Screen::setFieldOfView(double fov, const string &unit) {
  if (unit=="" || unit=="rad") ;
  else if (unit=="deg") fov *= GYOTO_DEGRAD;
  else if (unit=="arcmin") fov *= GYOTO_MINRAD;
  else if (unit=="arcsec") fov *= GYOTO_SECRAD;
  else if (unit=="mas") fov *= GYOTO_MASRAD;
  else if (unit=="�as" || unit=="uas") fov *= GYOTO_MUASRAD;
  else {
    stringstream ss;
    ss << "Screen::setFieldOfView(): unknown unit: \"" << unit << "\"";
    throwError(ss.str());
  }
  setFieldOfView(fov);
}
void Screen::setFieldOfView(double fov) { fov_ = fov; }

size_t Screen::getResolution() { return npix_; }
void Screen::setResolution(size_t n) { npix_ = n; }

#ifdef GYOTO_USE_XERCES
void Screen::fillElement(FactoryMessenger *fmp) {
  FactoryMessenger* child = NULL;
  if (gg_) fmp -> setMetric (gg_) ;
  fmp -> setParameter ("Time", tobs_);
  fmp -> setParameter ("FieldOfView", fov_);
  fmp -> setParameter ("Resolution", npix_);
  double d = getDistance();
  if (gg_() && (gg_->getMass() == 1.)) {
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
  if (dmax_ != GYOTO_SCREEN_DMAX) child -> setSelfAttribute("dmax", dmax_);
  delete child; child = NULL;
  fmp -> setParameter ("PALN", getPALN());
  fmp -> setParameter ("Inclination", getInclination());
  fmp -> setParameter ("Argument", getArgument());
  if (spectro_ && spectro_ -> getKind() != GYOTO_SPECTRO_KIND_NONE) {
    child = fmp -> makeChild("Spectrometer");
    spectro_ -> fillElement(child) ;
    delete child; child = NULL;
  }
}

SmartPointer<Screen> Screen::Subcontractor(FactoryMessenger* fmp) {

  string name="", content="", unit="", tunit="";
  SmartPointer<Screen> scr = new Screen();
  scr -> setMetric(fmp->getMetric());
  int tobs_found=0;
  double tobs_tmp, pos[4] ;
  char * tc;

  while (fmp->getNextParameter(&name, &content)) {
    tc = const_cast<char*>(content.c_str());
    unit = fmp -> getAttribute("unit");
    if      (name=="Time")     {tobs_tmp = atof(tc); tunit=unit; tobs_found=1;}
    else if (name=="Position") {for (int i=0;i<4;++i) pos[i] = strtod(tc, &tc);
                                  scr -> setObserverPos (pos); }
    else if (name=="Distance")    
      {
	scr -> setDistance    ( atof(tc), unit );
	string dmax = fmp -> getAttribute("dmax");
	if (dmax != "") scr -> setDmax(atof(dmax.c_str()));
      }
    else if (name=="PALN")        scr -> setPALN        ( atof(tc), unit );
    else if (name=="Inclination") scr -> setInclination ( atof(tc), unit );
    else if (name=="Argument")    scr -> setArgument    ( atof(tc), unit );
    else if (name=="FieldOfView") scr -> setFieldOfView ( atof(tc), unit );
    else if (name=="Resolution")  scr -> setResolution  ( atoi(tc) );
    else if (name=="Spectrometer")
      scr -> setSpectrometer (SpectrometerSubcontractor(fmp->getChild()));
  }

  if (tobs_found) scr -> setTime ( tobs_tmp, tunit );

  if (debug())
    cerr << "DEBUG: Screen::Subcontractor(): dmax_ = " << scr->getDmax() <<endl;

  return scr;
}
#endif
