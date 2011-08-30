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

#include "GyotoUtils.h"
#include <GyotoPhoton.h>
#include <GyotoFixedStar.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <cmath>
#include <limits>
#include <cstring>
#include <float.h>

using namespace std;
using namespace Gyoto;

FixedStar::FixedStar() : Astrobj("FixedStar"), radius_(0),
			 use_generic_impact_(0)
{
  if (debug())
    cerr << "DEBUG: in FixedStar::FixedStar(void)" << endl;
  for (int i=0;i<3;++i) pos_[i]=0.;
  setRadius(0.);
  if (debug())
    cerr << "DEBUG: out FixedStar::FixedStar(void)" << endl;
}

FixedStar::FixedStar(SmartPointer<Gyoto::Metric> gg, double StPsn[3],
		     double rad) :
  Astrobj("FixedStar"), radius_(rad), use_generic_impact_(0)
{
  if (debug())
    cerr << "DEBUG: in FixedStar::FixedStar(metric, pos, rad)" << endl;
  gg_ = gg;
  for (int i=0;i<3;++i) pos_[i] = StPsn[i]; 
  setRadius(rad);
  if (debug())
    cerr << "DEBUG: out FixedStar::FixedStar(metric, pos, rad)" << endl;
}

FixedStar::FixedStar(const FixedStar& orig) :
  Astrobj(orig), radius_(orig.radius_),
  use_generic_impact_(orig.use_generic_impact_)
{
  for (int i=0; i<3; ++i) pos_[i] = orig.pos_[i];
}
FixedStar* FixedStar::clone() const { return new FixedStar(*this); }

FixedStar::~FixedStar() {

  if (debug()) cout << "FixedStar Destruction" << endl;

}

void FixedStar::getVelocity(double const pos[4], double vel[4]) {
  for (size_t i=0; i<4; ++i) vel[i]=0.;
  vel[0]=gg_->SysPrimeToTdot(pos, vel+1);
}

double FixedStar::operator()(double const coord[4]) {
  double x, y, z, xs, ys, zs;
  switch (gg_->getCoordKind()) {
  case GYOTO_COORDKIND_CARTESIAN:
    xs= pos_[0];
    ys= pos_[1];
    zs= pos_[2];
    x = coord[1];
    y = coord[2];
    z = coord[3];
    break;
  case GYOTO_COORDKIND_SPHERICAL:
    {
      double rs=pos_[0];
      double ths=pos_[1];
      double phs=pos_[2];
      xs= rs*sin(ths)*cos(phs);
      ys= rs*sin(ths)*sin(phs);
      zs= rs*cos(ths);
    }
    x= coord[1]*sin(coord[2])*cos(coord[3]);
    y= coord[1]*sin(coord[2])*sin(coord[3]);
    z= coord[1]*cos(coord[2]);
    break;
  default:
    throwError("unsupported coordkind");
  }
  double dx = x-xs, dy=y-ys, dz=z-zs;

  return dx*dx + dy*dy + dz*dz;
}
int FixedStar::Impact(Photon *ph, size_t index, AstrobjProperties *data) {
  if (debug())
    cerr << "DEBUG: FixedStar::Impact(): use_generic_impact_="
	 << use_generic_impact_ << endl;
  if (use_generic_impact_) return Astrobj::Impact(ph, index, data);
  Impact_(ph, index, data);
}

int FixedStar::Impact_(Photon *ph, size_t index, AstrobjProperties *data) {
  // coord2 is the coordinate of the photon at index, coord1 the previous location.
  double coord2[8], coord1[8], coord_ph_hit[8], coord_obj_hit[8];
  ph->getCoord(index, coord1);
  ph->getCoord(index+1, coord2);
  int hit=0;

  //TO IMPROVE!
  //for (int ii=0;ii<8;ii++) coord_ph_hit[ii]=coord1[ii];
  
  int coordkind = gg_ -> getCoordKind();
  
  double x1, y1, z1, x2, y2, z2, xs, ys, zs;

  switch (coordkind) {
  case GYOTO_COORDKIND_SPHERICAL:
    {
      double rs=pos_[0];
      double ths=pos_[1];
      double phs=pos_[2];
      xs= rs*sin(ths)*cos(phs);
      ys= rs*sin(ths)*sin(phs);
      zs= rs*cos(ths);
    }
    x2= coord2[1]*sin(coord2[2])*cos(coord2[3]);
    y2= coord2[1]*sin(coord2[2])*sin(coord2[3]);
    z2= coord2[1]*cos(coord2[2]);
    x1= coord1[1]*sin(coord1[2])*cos(coord1[3]);
    y1= coord1[1]*sin(coord1[2])*sin(coord1[3]);
    z1= coord1[1]*cos(coord1[2]);
    break;
  case GYOTO_COORDKIND_CARTESIAN:
    xs= pos_[0];
    ys= pos_[1];
    zs= pos_[2];
    x2=coord2[1];
    y2=coord2[2];
    z2=coord2[3];
    x1=coord1[1];
    y1=coord1[2];
    z1=coord1[3];
    break;
  default:
    throwError("Incompatible coordinate kind in FixedStar::Impact");
  }

  /*cout << "In FixedStar Impact! coords= " << endl;
  cout << xs << " " << ys << " " << zs << endl;
  cout << x1 << " " << y1 << " " << z1 << endl;
  cout << "***" << endl;*/

  if (coord2[0]==coord1[0]) throwError("FixedStar::Impact: t1==t2");

  double dtm1=1./(coord2[0]-coord1[0]); // (t2-t1)^(-1)
  double ax=(x2-x1)*dtm1; // x(t)-xs=ax*t+bx
  double ay=(y2-y1)*dtm1;
  double az=(z2-z1)*dtm1;

  double bx=x1-ax*coord1[0]-xs;
  double by=y1-ay*coord1[0]-ys;
  double bz=z1-az*coord1[0]-zs;

  // dist^2=a*t^2+b*t+c ;
  double a=    ax*ax + ay*ay + az*az;
  double b=2.*(ax*bx + ay*by + az*bz);
  double c=    bx*bx + by*by + bz*bz - radius_*radius_;

  if (a == 0 ) throwError("FixedStar::Impact: a ==0");
  double tmin = -b/(2.*a);
  if (tmin > coord2[0]) tmin=coord2[0];
  if (tmin < coord1[0]) tmin=coord1[0];

  double d2min=a*tmin*tmin+b*tmin+c;

  if (d2min<=0) {
    double discr = b*b-4.*a*c , tentry=DBL_MIN;
    if (debug()) cerr << "DEBUG: FixedStar::Impact(): hit!" << endl;
    hit=1;
    if (data) {
      tmin=(-b+sqrt(discr))/(2.*a);
      tentry=(-b-sqrt(discr))/(2.*a);
      if (tmin > coord2[0]) tmin=coord2[0];
      if (tmin < coord1[0]) throwError("FixedStar::Impact(): fishy heuristic");
      if (tentry > coord2[0])
	throwError("FixedStar::Impact(): fishy heuristic");
      if (tentry < coord1[0]) tentry=coord1[0];
      d2min=0;
      //if (data->intensity) *data->intensity=1.;
    }

    // Update photon coordinate to impact point coordinates:
    coord_ph_hit[0]=tmin;
    ph -> getCoord(coord_ph_hit, 1, coord_ph_hit+1, coord_ph_hit+2,
		   coord_ph_hit+3, coord_ph_hit+4, coord_ph_hit+5,
		   coord_ph_hit+6, coord_ph_hit+7);

    if (data) {
      for (int ii=0; ii<4; ++ii) {
	coord_obj_hit[ii]=coord_ph_hit[ii];
	coord_obj_hit[ii+4]=0.;
      }
      coord_obj_hit[4]=gg_->SysPrimeToTdot(coord_obj_hit, coord_obj_hit+5);
      processHitQuantities(ph, coord_ph_hit, coord_obj_hit, tmin-tentry, data);
    }
  }
  if (data) {
    if (data->time) *data->time=tmin;
    if (data->distance && *(data->distance)>d2min) {
      *data->distance=d2min;
      /*if (data->x) *data->x= ax*(tmin-coord1[0]) + coord1[1];
      if (data->y) *data->y= ay*(tmin-coord1[0]) + coord1[2];
      if (data->z) *data->z= az*(tmin-coord1[0]) + coord1[3];*/ //obsolete version (useful?)
    }
    /* FirstMinDist */
    if (data->first_dmin) { 
      if (!data->first_dmin_found) {
	if (*(data->first_dmin)>d2min) *(data->first_dmin)=d2min;
	else data->first_dmin_found=1;
      }
    }
  }
 
  return hit;
}

double FixedStar::emission(double, double dsem, double*, double*) const{
  if (flag_radtransf_) return dsem;
  return 1.;
}

double FixedStar::getRadius() const { return radius_; }

const double* const FixedStar::getPos() const { return pos_; }

void FixedStar::getPos(double dst[3]) const
{ for (int i=0; i<3;++i) dst[i]=pos_[i]; }

void FixedStar::setMetric(SmartPointer<Metric> gg) {
 if (debug())
   cerr << "DEBUG: in FixedStar::setMetric(gg)\n";
 Astrobj::setMetric(gg);
 setRadius(radius_);
}

void FixedStar::setRadius(double r) {
  radius_ = r;
  critical_value_=r*r;
  safety_value_=1.1*critical_value_;
  if (!gg_()) {
    if (debug())
      cerr << "DEBUG: FixedStar::setRadius(radius): metric is not set yet\n";
    return;
  }
  switch (gg_ -> getCoordKind()) {
  case GYOTO_COORDKIND_SPHERICAL:
    rmax_=3.*(pos_[0]+radius_);
    break;
  case GYOTO_COORDKIND_CARTESIAN:
    rmax_=3.*(sqrt(pos_[0]*pos_[0]+pos_[1]*pos_[1]+pos_[2]*pos_[2])+radius_);
    break;
  default:
    throwError("unimplemented coordinate system in FixedStar");
  } 
}

void FixedStar::setPos(const double p[3])
{ for (int i=0; i<3; ++i) pos_[i]=p[i]; setRadius(radius_);}

void FixedStar::useGenericImpact(int val) {
  if (debug())
    cerr << "DEBUG: FixedStar::useGenericImpact(" << val << ")\n";
  use_generic_impact_=val;
}

#ifdef GYOTO_USE_XERCES
void FixedStar::fillElement(factoryMessenger *fmp) const {
  fmp -> setParameter ("Radius", getRadius());
  fmp -> setParameter ("Position", const_cast<double*>(pos_), 3);
  if (use_generic_impact_) fmp -> setParameter ("UseGenericImpact");
  Astrobj::fillElement(fmp);
}

SmartPointer<Astrobj> Gyoto::FixedStar::Subcontractor(factoryMessenger* fmp) {

  string name="", content="";

  SmartPointer<FixedStar> ao = new FixedStar();
  ao -> setMetric(fmp->getMetric());

  while (fmp->getNextParameter(&name, &content)) {
    char* tc = const_cast<char*>(content.c_str());
    if(name=="Radius") ao -> setRadius(atof(tc));
    else if(name=="Position") {
      double pos[3];
      for (int i=0;i<3;++i) pos[i] = strtod(tc, &tc);
      ao -> setPos(pos);
    }
    else if (name=="UseGenericImpact")
      ao -> useGenericImpact(1);
    else
      ao -> setGenericParameter(name, content);
  }

  return ao;

}

void Gyoto::FixedStar::Init() {
  Gyoto::Astrobj::Register("FixedStar", &Gyoto::FixedStar::Subcontractor);
}
#endif
