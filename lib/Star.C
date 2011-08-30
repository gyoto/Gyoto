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
#include "GyotoStar.h"
#include "GyotoPhoton.h"
#include "GyotoPowerLawSpectrum.h"
#include "GyotoBlackBodySpectrum.h"
#include <iostream>
#include <cmath>
#include <string>
#include <cstdlib>
#include <float.h>
#include <sstream>
#include <string.h>

using namespace std;
using namespace Gyoto;

Star::Star() :
  Astrobj("Star"),
  spectrum_(NULL),
  opacity_(NULL)
{
  if (debug())
    cerr << "DEBUG: in Star::Star()" << endl;
  setRadius(0.);

  metric_=NULL; gg_=NULL;
  spectrum_ = new Spectrum::BlackBody(); 
  opacity_ = new Spectrum::PowerLaw(0., 1.); 
}

Star::Star(SmartPointer<Metric> met, double rad,
	   double pos[4],
	   double v[3]) :
  Astrobj("Star"),
  radius_(rad),
  spectrum_(NULL), opacity_(NULL)
{
  critical_value_ = radius_*radius_;
  safety_value_ = critical_value_*1.1 + 0.1;
  spectrum_ = new Spectrum::BlackBody(); 
  opacity_ = new Spectrum::PowerLaw(0., 1.); 
  if (debug()) {
    cerr << "DEBUG: Star Construction " << endl
	 << "       POS=[" << pos[0];
    for (int i=1; i<4; ++i) cerr << ", " << pos[i];
    cerr << "]\n       VEL=[" << v[0] ;
    for (int i=1; i<3; ++i) cerr << ", " << v[i];
    cerr << "]\n       RADIUS=" << rad << endl;

  }

  metric_=met;
  gg_=met;

  double tdot0=metric_->SysPrimeToTdot(pos, v);

  if (debug()) cerr << "       TDOT0=" << tdot0 << endl;

  double coord[8]={pos[0], pos[1], pos[2], pos[3],
		   tdot0, v[0]*tdot0, v[1]*tdot0, v[2]*tdot0};

  Worldline::setInitialCondition(metric_, coord, 1);
    //last number : direction of integration + or -1
}

Star::Star(const Star& orig) :
  Astrobj(orig), Worldline(orig),
  radius_(orig.radius_)
{
  if (debug()) cerr << "Star copy" << endl;
  gg_ = metric_; // we have two distinct clones of the metric, not good...
  if (orig.spectrum_()) spectrum_=orig.spectrum_->clone();
  if (orig.opacity_()) opacity_=orig.opacity_->clone();
}

Star* Star::clone() const { return new Star(*this); }

Star::~Star() {
  if (debug()) cerr << "DEBUG: Star::~Star()\n";
}

string Star::className() const { return  string("Star"); }
string Star::className_l() const { return  string("star"); }

SmartPointer<Metric> Star::getMetric() const { return gg_; }
void Star::setMetric(SmartPointer<Metric> gg) {gg_=gg; metric_=gg;}

SmartPointer<Spectrum::Generic> Star::getSpectrum() const { return spectrum_; }
void Star::setSpectrum(SmartPointer<Spectrum::Generic> sp) {spectrum_=sp;}

SmartPointer<Spectrum::Generic> Star::getOpacity() const { return opacity_; }
void Star::setOpacity(SmartPointer<Spectrum::Generic> sp) {opacity_=sp;}


void Star::setInitialCondition(double coord[8]) {
  if (!metric_) throwError("Please set metric before calling Star::setInitialCondition(double*)");
  Worldline::setInitialCondition(metric_, coord, 1);
}

double Star::getMass() const {return 1. ;}

void Star::getVelocity(double const pos[4], double vel[4]) {
  getCoord(pos, 1, NULL, NULL, NULL, vel, vel+1, vel+2, vel+3);
}

double Star::operator()(double const coord[4]) {
  double coord_st[4] = {coord[0]};
  double coord_ph[4] = {coord[0]};
  double sintheta;
  getCartesian(coord_st, 1, coord_st+1, coord_st+2, coord_st+3);
  switch (gg_->getCoordKind()) {
  case GYOTO_COORDKIND_CARTESIAN:
    memcpy(coord_ph+1, coord+1, 3*sizeof(double));
    break;
  case GYOTO_COORDKIND_SPHERICAL:
    coord_ph[1] = coord[1] * (sintheta=sin(coord[2])) * cos(coord[3]);
    coord_ph[2] = coord[1] * sintheta * sin(coord[3]);
    coord_ph[3] = coord[1] * cos(coord[2]) * sin(coord[3]);
    break;
  default:
    throwError("unsupported coordkind");
  }
  double dx = coord_ph[1]-coord_st[1];
  double dy = coord_ph[2]-coord_st[2];
  double dz = coord_ph[3]-coord_st[3];

  return dx*dx + dy*dy + dz*dz;
}

int Star::Impact_(Photon *ph, size_t index, AstrobjProperties *data) {
  //xyz coord from BL
  double coord_ph_hit[8], coord_obj_hit[8];
  double p1[8], p2[8];
  ph->getCartesianPos(index, p1);
  ph->getCartesianPos(index+1, p2);
  double t1= p1[0];
  double t2= p2[0];

  if (debug())
    cerr << "DEBUG: in Star::Impact(): t1="<<t1<<", t2="<<t2<<endl;

  // Assume the Photon's trajectory is linear over t1 < t2
  double dtm1=1./(t2-t1); // (t2-t1)^(-1)
  double ax=(p2[1]-p1[1])*dtm1; // x(t)=ax*t+bx
  double ay=(p2[2]-p1[2])*dtm1;
  double az=(p2[3]-p1[3])*dtm1;
  double bx=p1[1]-ax*t1;
  double by=p1[2]-ay*t1;
  double bz=p1[3]-az*t1;

  // Integrate Star's geodesic over same period
  xFill(t1); xFill(t2);// has been called at the previous step

  int hit=0;

  // t1 -> t2 may span several integration steps for the Star
  // integrate emission & absorption backwards between t_exit and t_entry.
  // Assume we cross the star only once.
  size_t i2=imin_+1;
  while (x0_[i2]<t2){
    if (i2 == imax_) {
      if (debug())
	cerr << "DEBUG: Star::Impact(): Object doesn't exist at this time\n";
      return 0;
    }
    ++i2;
  }

  size_t i1 = i2-1;

  // d2_1, d2_2: Star <-> Photon distance at i1, i2
  double d2_1, d2_2, o1[4], o2[4];

  getCartesianPos(i1, o1);
  d2_1=(ax*o1[0]+bx-o1[1])*(ax*o1[0]+bx-o1[1])
    +(ay*o1[0]+by-o1[2])*(ay*o1[0]+by-o1[2])
    +(az*o1[0]+bz-o1[3])*(az*o1[0]+bz-o1[3]);

  getCartesianPos(i2, o2);
  d2_2=(ax*o2[0]+bx-o2[1])*(ax*o2[0]+bx-o2[1])
    +(ay*o2[0]+by-o2[2])*(ay*o2[0]+by-o2[2])
    +(az*o2[0]+bz-o2[3])*(az*o2[0]+bz-o2[3]);

  double r2 = radius_*radius_;

  while ( d2_1 <= d2_2 && o1[0] > t1 && d2_1 > r2) {

    i2=i1--;
    for (int i=0; i<4; ++i) o2[i] = o1[i];
    d2_2=d2_1;

    getCartesianPos(i1, o1);
    d2_1=(ax*o1[0]+bx-o1[1])*(ax*o1[0]+bx-o1[1])
      +(ay*o1[0]+by-o1[2])*(ay*o1[0]+by-o1[2])
      +(az*o1[0]+bz-o1[3])*(az*o1[0]+bz-o1[3]);

  }

  // Assume object's trajectory is linear in cartesian 3 coord
  double dtom1=1./(o2[0]-o1[0]); // (t2-t1)^(-1)
  double aox=(o2[1]-o1[1])*dtom1; // xo(t)=aox*t+box
  double aoy=(o2[2]-o1[2])*dtom1;
  double aoz=(o2[3]-o1[3])*dtom1;

  double box=o1[1]-aox*o1[0];
  double boy=o1[2]-aoy*o1[0];
  double boz=o1[3]-aoz*o1[0];

  // Solve for closest approach to star's _center_
  double dax = (ax-aox), day=ay-aoy, daz=az-aoz;
  double dbx = (bx-box), dby=by-boy, dbz=bz-boz;
  double a =   (dax*dax + day*day + daz*daz);
  double b =2.*(dax*dbx + day*dby + daz*dbz);
  double c =   (dbx*dbx + dby*dby + dbz*dbz - r2); 

  double tmin = -b/(2.*a);

  if (tmin > t2) tmin=t2;
  if (tmin < t1) tmin=t1;
  if (tmin > o2[0]) tmin=o2[0];
  if (tmin < o1[0]) tmin=o1[0];

  if (tmin == o2[0] && d2_2 > r2) {
    // in this case the minimum may be in the next sample, do it again
    i1=i2++;
    for (int i=0; i<4; ++i) o1[i] = o2[i];
    d2_1=d2_2;

    getCartesianPos(i2, o2);
    d2_2=(ax*o2[0]+bx-o2[1])*(ax*o2[0]+bx-o2[1])
      +(ay*o2[0]+by-o2[2])*(ay*o2[0]+by-o2[2])
      +(az*o2[0]+bz-o2[3])*(az*o2[0]+bz-o2[3]);

    dtom1=1./(o2[0]-o1[0]); // (t2-t1)^(-1)
    aox=(o2[1]-o1[1])*dtom1; // xo(t)=aox*t+box
    aoy=(o2[2]-o1[2])*dtom1;
    aoz=(o2[3]-o1[3])*dtom1;
    
    box=o1[1]-aox*o1[0];
    boy=o1[2]-aoy*o1[0];
    boz=o1[3]-aoz*o1[0];
    
    dax = (ax-aox); day=ay-aoy; daz=az-aoz;
    dbx = (bx-box); dby=by-boy; dbz=bz-boz;
    a =   (dax*dax + day*day + daz*daz);
    b =2.*(dax*dbx + day*dby + daz*dbz);
    c =   (dbx*dbx + dby*dby + dbz*dbz - r2); 
    
    if (a == 0 ) throwError("Star::Impact: a == 0");
    tmin = -b/(2.*a);
    if (tmin > t2) tmin=t2;
    if (tmin < t1) tmin=t1;
    if (tmin > o2[0]) tmin=o2[0];
    if (tmin < o1[0]) tmin=o1[0];

    if ( tmin == o2[0] ) throwError("Star::Impact: fishy heuristic");

  }

  double d2min=a*tmin*tmin+b*tmin+c;
  if (debug())
    cerr << "DEBUG: Star::Impact(): d2min="<<d2min<<", tmin=" <<tmin<<endl;

  if (d2min<=0) {
    // now solve for point crossing the surface
    double tlast, tfirst;
    tmin=(-b+sqrt(b*b-4.*a*c))/(2.*a);

    if (tmin > t2) tmin=t2;
    if (tmin < t1) tmin=t1;
    if (tmin > o2[0]) tmin=o2[0];
    if (tmin < o1[0]) tmin=o1[0];
    tfirst=tmin;

    double d2last; d2min=a*tmin*tmin+b*tmin+c;
    o1[0]=tmin;

    getCartesian(o1, 1, o1+1, o1+2, o1+3);
    if (debug()) {
      cerr << "DEBUG: Star::Impact: o1=["<< o1[0] << ", " << o1[1]
	   << ", " << o1[2] << ", " << o1[3] << "]\n";
    }
    d2last=(ax*o1[0]+bx-o1[1])*(ax*o1[0]+bx-o1[1])
      +(ay*o1[0]+by-o1[2])*(ay*o1[0]+by-o1[2])
      +(az*o1[0]+bz-o1[3])*(az*o1[0]+bz-o1[3]);

    if (d2last>r2*1.1) {
      stringstream ss;
      ss << "Star::Impact(): d2 min = 0 but d2last="<<d2last;
      throwError(ss.str());
    }
    d2last=r2;
    if (debug())
      cerr << "DEBUG: Star::Impact: hit! at t_photon= " << tmin << endl;
    hit=1;

    // integrate backwards until entry point
    if (flag_radtransf_) {
      while (d2last<=r2 && (tlast=tfirst)>t1) {
 
	for (int i=0; i<4; ++i) o2[i] = o1[i];

	getCartesianPos(i1, o1);
	if(debug())
	  cerr << "DEBUG: Star::Impact(): t1=" << t1
	       << ", tfirst=" << o1[0]
	       << ", tlast=" << tlast
	       << ", t2=" << t2
	       << endl;

	if (o1[0]<t1) {
	  o1[0]=t1;
	  getCartesian(o1, 1, o1+1, o1+2, o1+3);
	  if (debug())
	    cerr << "DEBUG: Star::Impact(): reaching t1. tfirst=t1="
		 << o1[0] << endl;
	}

	d2last=(ax*o1[0]+bx-o1[1])*(ax*o1[0]+bx-o1[1])
	  +(ay*o1[0]+by-o1[2])*(ay*o1[0]+by-o1[2])
	  +(az*o1[0]+bz-o1[3])*(az*o1[0]+bz-o1[3]);

	if (d2last > r2) {
	  // solve to find entry point. There _is_ one.
	  if (tlast==tmin) {
	    o2[0]=tmin;
	    getCartesian(o2, 1, o2+1, o2+2, o2+3);
	  } else getCartesianPos(i2, o2);

	  if (debug())
	    cerr << "DEBUG: Star::Impact(): r2=" << r2
		 << ", d2(tfirst)=" << d2last
		 << ", d2(tlast)="
		 << ((ax*o2[0]+bx-o2[1])*(ax*o2[0]+bx-o2[1])+
		     (ay*o2[0]+by-o2[2])*(ay*o2[0]+by-o2[2])+
		     (az*o2[0]+bz-o2[3])*(az*o2[0]+bz-o2[3]))
		 << endl;

	  dtom1=1./(o2[0]-o1[0]); // (t2-t1)^(-1)
	  aox=(o2[1]-o1[1])*dtom1; // xo(t)=aox*t+box
	  aoy=(o2[2]-o1[2])*dtom1;
	  aoz=(o2[3]-o1[3])*dtom1;
	  box=o1[1]-aox*o1[0];
	  boy=o1[2]-aoy*o1[0];
	  boz=o1[3]-aoz*o1[0];
	  if (debug())
	    cerr << "DEBUG: Star::Impact(): solving for entry point. "
		 << "dt=" << 1./dtom1 
		 <<", x="<<aox<<"*dt+"<<box
		 <<", y="<<aoy<<"*dt+"<<boy
		 <<", z="<<aoz<<"*dt+"<<boz<<endl;
    
	  dax = (ax-aox); day=ay-aoy; daz=az-aoz;
	  dbx = (bx-box); dby=by-boy; dbz=bz-boz;
	  a =   (dax*dax + day*day + daz*daz);
	  b =2.*(dax*dbx + day*dby + daz*dbz);
	  c =   (dbx*dbx + dby*dby + dbz*dbz - r2); 
	  double Delta = b*b-4.*a*c;
	  if (Delta<0) break;
	  o1[0]=(-b-sqrt(b*b-4.*a*c))/(2.*a);
	  if (debug())
	    cerr << "DEBUG: Star::Impact(): solving for entry point. "
		 << "a=" << a << ", b=" << b << ", c=" << c
		 << ", Delta=" << b*b-4.*a*c 
		 << ", tfirst=" << o1[0] << endl;
	}
	tfirst=o1[0];
	if (tfirst > tlast) break;
	coord_ph_hit[0]=coord_obj_hit[0]=tlast;
	getCoord(coord_obj_hit, 1,
		 coord_obj_hit+1, coord_obj_hit+2, coord_obj_hit+3,
		 coord_obj_hit+4, coord_obj_hit+5, coord_obj_hit+6,
		 coord_obj_hit+7);
	ph->getCoord(coord_ph_hit, 1,
		     coord_ph_hit+1, coord_ph_hit+2, coord_ph_hit+3,
		     coord_ph_hit+4, coord_ph_hit+5, coord_ph_hit+6,
		     coord_ph_hit+7);
	if (debug())
	  cerr << "DEBUG: Star::Impact: tlast=" << tlast
	       << ", tfirst=" << tfirst << ", dt=" << tlast-tfirst << endl;

	processHitQuantities(ph, coord_ph_hit, coord_obj_hit,
			     tlast-tfirst, data);
	i2=i1--;
      }
    } else {
      coord_ph_hit[0]=coord_obj_hit[0]=tmin;
      getCoord(coord_obj_hit, 1,
	       coord_obj_hit+1, coord_obj_hit+2, coord_obj_hit+3,
	       coord_obj_hit+4, coord_obj_hit+5, coord_obj_hit+6,
	       coord_obj_hit+7);
      ph->getCoord(coord_ph_hit, 1,
		   coord_ph_hit+1, coord_ph_hit+2, coord_ph_hit+3,
		   coord_ph_hit+4, coord_ph_hit+5, coord_ph_hit+6,
		   coord_ph_hit+7);
      
      processHitQuantities(ph, coord_ph_hit, coord_obj_hit,
			   1., data);
    }
  }

  if (data) {
    /* EmissionTime */
    if (data->time) *data->time=tmin;
    /* MinDistance */
    if ((data->distance) && (*(data->distance)>d2min) )
      *data->distance=hit?0:d2min;
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

double Star::emission(double nu_em, double dsem, double *, double *) const {
  if (flag_radtransf_) return (*spectrum_)(nu_em, (*opacity_)(nu_em), dsem);
  return (*spectrum_)(nu_em);
}

double Star::transmission(double nuem, double dsem, double*) const {
  if (!flag_radtransf_) return 0.;
  double opacity = (*opacity_)(nuem);
  if (debug())
    cerr << "DEBUG: Star::transmission(nuem="<<nuem<<", dsem="<<dsem<<"), "
	 << "opacity=" << opacity << "\n";
  if (!opacity) return 1.;
  return exp(-opacity*dsem);
}

double Star::integrateEmission(double nu1, double nu2, double dsem,
			       double *, double *) const {
  if (flag_radtransf_)
    return spectrum_->integrate(nu1, nu2, opacity_(), dsem);
  return spectrum_->integrate(nu1, nu2);
}

double Star::getRmax() {
  // rmax may not be up to date. Assume it's OK if it's non-0 warning:
  // if line is extended, need to either update rmax or reset it to 0
  // if (debug()) cerr << "DEBUG: Star::getRmax(): rmax_set_==" 
  //                   << rmax_set_ << endl;
  if (!rmax_set_ && !rmax_) {
    size_t i;
    for (i=imin_;i<=imax_;++i) if (x1_[i]>rmax_) rmax_=x1_[i];
    rmax_*=3.;
  }
  return rmax_;
}

void Star::unsetRmax() {
  rmax_set_=0;
  rmax_=DBL_MAX;
}

double Star::getRadius() const {
  return radius_;
}

void Star::setRadius(double r) {
  radius_=r;
  critical_value_ = r*r;
  safety_value_ = critical_value_*1.1+0.1;
}

#ifdef GYOTO_USE_XERCES
void Star::fillElement(factoryMessenger *fmp) const {
  factoryMessenger * childfmp=NULL;

  fmp -> setMetric (getMetric()) ;
  fmp -> setParameter ("Radius", getRadius());
  double coord[8];
  getInitialCoord(coord);
  fmp -> setParameter ("Position", coord, 4);
  double vel[3] = {coord[5]/coord[4], coord[6]/coord[4], coord[7]/coord[4]};
  fmp -> setParameter ("Velocity", vel, 3);

  childfmp = fmp -> makeChild ( "Spectrum" );
  spectrum_ -> fillElement(childfmp);
  delete childfmp;

  childfmp = fmp -> makeChild ( "Opacity" );
  opacity_ -> fillElement(childfmp);
  delete childfmp;

  Astrobj::fillElement(fmp);
}

SmartPointer<Astrobj> Gyoto::Star::Subcontractor(factoryMessenger* fmp) {

  string name="", content="";
  int pos_found=0, vel_found=0;
  double pos[4], v[3], radius=0.;
  SmartPointer<Metric> gg = NULL;
  SmartPointer<Spectrum::Generic> sp = NULL, op = NULL;
  factoryMessenger * child = NULL;

  gg = fmp->getMetric();

  while (fmp->getNextParameter(&name, &content)) {
    char* tc = const_cast<char*>(content.c_str());
    if      (name=="Radius") radius=atof(tc);
    else if (name=="Position") {
      pos_found=1;
      for (int i=0;i<4;++i) pos[i] = strtod(tc, &tc);
    }
    else if (name=="Velocity") {
      vel_found=1;
      for (int i=0;i<3;++i) v[i] = strtod(tc, &tc);
    }
    else if (name=="Spectrum") {
      content = fmp -> getAttribute("kind");
      child = fmp -> getChild();
      sp = (*Spectrum::getSubcontractor(content))(child);
      delete child;
    }
    else if (name=="Opacity") {
      content = fmp -> getAttribute("kind");
      child = fmp -> getChild();
      op = (*Spectrum::getSubcontractor(content))(child);
      delete child;
    }
  }
  if (!pos_found) throwError("Position MUST be set in Star definition");
  if (!vel_found) throwError("Velocity MUST be set in Star definition");

  SmartPointer<Star> st = new Star(gg, radius, pos, v);
  if (sp) st -> setSpectrum(sp);
  if (op) {
    st -> setOpacity(op);
    st -> setFlag_radtransf(1);
  } else
    st -> setFlag_radtransf(0);

  fmp->reset();
  while (fmp->getNextParameter(&name, &content)) {
    st -> setGenericParameter(name, content);
  }

  return st;
}

void Gyoto::Star::Init() {
  Gyoto::Astrobj::Register("Star", &Gyoto::Star::Subcontractor);
}
#endif
