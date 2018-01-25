/*
    Copyright 2011-2018 Frederic Vincent, Thibaut Paumard

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
#include "GyotoPhoton.h"
#include "GyotoScreen.h"
#include "GyotoDefs.h"
#include "GyotoError.h"
#include "GyotoProperty.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <cstring>
#include <cstdlib>
#include <cfloat>

#define GYOTO_LIMIT_TRANSMISSION 1e-6 //0.36788
// NB: 0.36788=exp(-1), 1e-6 transmission is well thicker
// than the standard optical depth = 1

using namespace std;
using namespace Gyoto;

GYOTO_PROPERTY_START(Photon)
GYOTO_PROPERTY_ASTROBJ(Photon, Astrobj, astrobj)
GYOTO_WORLDLINE_PROPERTY_END(Photon, Object::properties)

Photon::Photon() :
  Worldline(),
  Object("Photon"),
  object_(NULL),
  freq_obs_(1.), transmission_freqobs_(1.),
  spectro_(NULL), transmission_(NULL)
 {}

Photon::Photon(const Photon& o) :
  Worldline(o), SmartPointee(o),
  object_(NULL),
  freq_obs_(o.freq_obs_), transmission_freqobs_(o.transmission_freqobs_),
  spectro_(NULL), transmission_(NULL)
{
  if (o.object_()) {
    object_  = o.object_  -> clone();
    object_ -> metric(metric_);
  }
  if (o.spectro_()) {
    spectro_ = o.spectro_ -> clone();
    _allocateTransmission();
    if (size_t nsamples = spectro_->nSamples())
      memcpy(transmission_, o.getTransmission(), nsamples*sizeof(double));
  }
}

Photon * Photon::clone() const { return new Photon(*this); }

bool Photon::isThreadSafe() const {
  // spectro_ is not a Property
  return
    Object::isThreadSafe()
    && (!spectro_ || spectro_ -> isThreadSafe());
}

Photon::Photon(Photon* orig, size_t i0, int dir, double step_max) :
  Worldline(orig, i0, dir, step_max), SmartPointee(),
  object_(orig->object_),
  freq_obs_(orig->freq_obs_),
  transmission_freqobs_(orig->transmission_freqobs_),
  spectro_(orig->spectro_), transmission_(orig->transmission_)
{
}

Photon::Refined::Refined(Photon* orig, size_t i0, int dir, double step_max) :
  Photon(orig, i0, dir, step_max),
  parent_(orig)
{
  freqObs(orig->freqObs());
}

Photon::Photon(SmartPointer<Metric::Generic> met,
	       SmartPointer<Astrobj::Generic> obj,
	       double* coord):
  Worldline(), freq_obs_(1.), transmission_freqobs_(1.), spectro_(NULL), transmission_(NULL)
{
  setInitialCondition(met, obj, coord);
}

Photon::Photon(SmartPointer<Metric::Generic> met, 
	       SmartPointer<Astrobj::Generic> obj, 
	       SmartPointer<Screen> screen, 
	       double d_alpha, double d_delta):
  Worldline(), object_(obj), freq_obs_(screen->freqObs()),
  transmission_freqobs_(1.),
  spectro_(NULL), transmission_(NULL)
{
  double coord[8], Ephi[4], Etheta[4];
  screen -> getRayCoord(d_alpha, d_delta, coord);
  if (parallel_transport_)
    screen -> getRayTriad(coord, Ephi, Etheta);
  Worldline::setInitialCondition(met, coord, -1, Ephi, Etheta);
  spectrometer(screen);
}

Photon::~Photon() {}

/* TRANSMISSION STUFF */
void Photon::_allocateTransmission() {
  if (transmission_) {
    delete [] transmission_;
    transmission_ = NULL;
  }
  if (spectro_()) {
    size_t nsamples = spectro_->nSamples();
    if (nsamples) {
      transmission_ = new double[nsamples];
      resetTransmission();
    }
  }
}

void Photon::resetTransmission() {
  transmission_freqobs_ = 1.;
  if (spectro_() && transmission_) {
    size_t nsamples = spectro_->nSamples();
    for (size_t i = 0; i < nsamples; ++i) transmission_[i] = 1.;
  }
}


double Photon::getMass() const { return 0.; }

void Photon::astrobj(SmartPointer<Astrobj::Generic> ao) {
  if (object_!=ao) {
    if (imin_<=imax_) imin_=imax_=i0_;
    object_=ao;
    if (metric_) object_->metric(metric_);
  }
}

void Photon::metric(SmartPointer<Metric::Generic> met) {
  Worldline::metric(met);
  if (object_) object_->metric(met);
}

void Photon::spectrometer(SmartPointer<Spectrometer::Generic> spr) {
  spectro_=spr;
  _allocateTransmission();
}
SmartPointer<Spectrometer::Generic> Photon::spectrometer() const { return spectro_; }


string Photon::className() const { return  string("Photon"); }
string Photon::className_l() const { return  string("photon"); }

SmartPointer<Astrobj::Generic> Photon::astrobj() const { return object_; }



void Photon::setInitialCondition(SmartPointer<Metric::Generic> met,
				 SmartPointer<Astrobj::Generic> obj,
				 SmartPointer<Screen> screen,
				 const double d_alpha,
				 const double d_delta)
{
  double coord[8], Ephi[4], Etheta[4];
  screen -> getRayCoord(d_alpha, d_delta, coord);
  if (parallel_transport_)
    screen -> getRayTriad(coord, Ephi, Etheta);
  Worldline::setInitialCondition(met, coord, -1, Ephi, Etheta);
  if (obj) object_=obj;

}

void Photon::setInitialCondition(SmartPointer<Metric::Generic> met,
				 SmartPointer<Astrobj::Generic> obj,
				 const double coord[8])
{
  if (!met) met = metric_;
  Worldline::setInitialCondition(met, coord, -1);
  if (obj) object_=obj;
}

void Photon::setInitialCondition(SmartPointer<Metric::Generic> met,
				 SmartPointer<Astrobj::Generic> obj,
				 const double coord[8],
				 const double Ephi[4],
				 const double Etheta[4])
{
  if (!met) met = metric_;
  Worldline::setInitialCondition(met, coord, -1, Ephi, Etheta);
  if (obj) object_=obj;
}

int Photon::hit(Astrobj::Properties *data) {

  /*
    Ray-tracing of the photon until the object_ is hit. Radiative
    transfer inside the object_ may then be performed depending on
    flag_radtransf. Final result (observed flux for instance,
    depending on object_'s Astrobj::Properties) will be stored in data.
   */

  //tmin_=-1000.;//DEBUG //NB: integration stops when t < Worldline::tmin_

  transmission_freqobs_=1.;
  size_t nsamples;
  if (spectro_() && (nsamples = spectro_->nSamples()))
    for (size_t ii=0; ii<nsamples; ++ii) transmission_[ii]=1.;

  double rmax=object_ -> rMax();
  int coordkind = metric_ -> coordKind();

  int hitt=0;
  //hitted=1 if object is hitted at least one time (hitt can be 0 even
  //if the object was hit if cross_max>0) ; hitt_crude=1 if the object
  //is not yet hit with adaptive integration step. A second integration
  //with small fixed step will be performed to determine more precisely
  //the surface point.
  Worldline::IntegState::state_type coord(parallel_transport_?16:8);
  int dir=(tmin_>x0_[i0_])?1:-1;
  size_t ind=i0_;
  stopcond=0;
  double rr=DBL_MAX, rr_prev=DBL_MAX;

  //-------------------------------------------------
  /*
    1-
    Call object_->Impact on the already computed part 
    of the geodesic to check whether the object_ was hit.
   */
  for (ind=i0_+dir;
       ((dir==1)?
	(ind<=imax_ && x0_[ind]<=tmin_): // conditions if dir== 1
	(ind>=imin_ && x0_[ind]>=tmin_)) // conditions if dir==-1
	 && !hitt;                    // condition in all cases
       ind+=dir) {
    switch (coordkind) {
    case GYOTO_COORDKIND_SPHERICAL:
      rr=x1_[ind];
      break;
    case GYOTO_COORDKIND_CARTESIAN:
      rr=sqrt(x1_[ind]*x1_[ind]+x2_[ind]*x2_[ind]+x3_[ind]*x3_[ind]);
      break;
    default:
      throwError("Incompatible coordinate kind in Photon.C");
    }
    if (rr<rmax)
      hitt = object_ -> Impact(this, ind, data);
  }
  if (hitt) {
#   if GYOTO_DEBUG_ENABLED
    GYOTO_DEBUG << "DEBUG: Photon.C: Hit for already computed position; "
		<< "Warning: radiative transfer not implemented "
		<< "for that case" << endl;
#   endif
    return hitt;
  } else if (((dir==1)?
	(ind==imax_ && x0_[ind]>=tmin_): // conditions if dir== 1
	(ind>=imin_ && x0_[ind]<=tmin_)) // conditions if dir==-1
	     && !hitt)
    return hitt;
  if (ind!=i0_) ind-=dir;
  //-------------------------------------------------

  //-------------------------------------------------
  /*
    2-
    Need to compute geodesic further.
    Set up integration.
  */

  // Check whether the arrays need to be expanded first
  if (dir==1 && ind==x_size_) ind=xExpand(1);
  else if (dir==-1 && ind==0) ind=xExpand(-1);

  // Set up integration
  getCoord(ind, coord);
  switch (coordkind) {
  case GYOTO_COORDKIND_SPHERICAL:
    rr=coord[1];
    break;
  case GYOTO_COORDKIND_CARTESIAN:
    rr=sqrt(coord[1]*coord[1]+coord[2]*coord[2]+coord[3]*coord[3]);
    break;
  default:
    throwError("Incompatible coordinate kind in Photon.C");
  }

  state_->init(this, coord, delta_* dir);
  //delta_ = initial integration step (defaults to 0.01)

  size_t count=0;// Must remain below count_max (prevents infinite integration)

  //-------------------------------------------------

  //-------------------------------------------------
  /*
    3- Integration loop: integrate the geodesic until stopcond is 1.
    Possible stopping conditions: 
    - transmission_freqobs_ low [see transmission() function 
       in astrobjs, which defaults to 0 (optically thick) 
       or 1 (optically thin) in Astrobj.C]
    - t < tmin_ (if dir=-1), [NB: tmin_ defaults to -DBL_MAX in Worldline.C]
    - photon is at r>rmax (defined for each object) and goes even further
    - metric tells it's time to stop (eg horizon crossing)
    - t does not evolve [to investigate, metric should have stopped
       integration before, see above]
    - count>count_max [should never be used, just to prevent infinite
    integration in case of a bug]
   */

  double h1max=DBL_MAX;
  while (!stopcond) {
    // Next step along photon's worldline
    h1max=object_ -> deltaMax(&coord[0]);
    stopcond  = state_ -> nextStep(coord, h1max);
    if (!secondary_){ // to compute only primary image
      double sign = x1_[i0_]*cos(x2_[i0_]);
      if (coord[1]*cos(coord[2])*sign<0. && x1_[ind]*cos(x2_[ind])*sign<0.)
	return 0;
    }
    if (stopcond) {
#     if GYOTO_DEBUG_ENABLED
      GYOTO_DEBUG << "stopcond set by integrator\n";
#     endif
      int shadow=object_->showshadow();
      if (shadow && data && data->spectrum){
	SmartPointer<Spectrometer::Generic> spr = spectrometer();
	size_t nbnuobs = spr() ? spr -> nSamples() : 0 ;
	for (size_t ii=0; ii<nbnuobs; ++ii) {
	  data->spectrum[ii*data->offset] = 1e10; // something very big
	}
      }
      break;
    }
    if (coord[0] == x0_[ind]) { // here, ind denotes previous step
      stopcond=1;
#     if GYOTO_DEBUG_ENABLED
      GYOTO_DEBUG << "time did not evolve, break." << endl;
#     endif
      break;
    }
    if((stopcond=metric_->isStopCondition(&coord[0]))) {
#     if GYOTO_DEBUG_ENABLED
      GYOTO_DEBUG << "stopcond step by metric"<<endl;
#     endif
      break;
    }

    if ( ++count > maxiter_ ) {
      GYOTO_SEVERE << "Photon::hit: too many iterations ("<<count<<" vs. "
		   << maxiter_<<"), break" << endl;
      stopcond = 1;
      break;
    }
     
    ind +=dir;
    // store photon's trajectory for later use
    x0_[ind] = coord[0];
    x1_[ind] = coord[1];
    x2_[ind] = coord[2];
    x3_[ind] = coord[3];
    x0dot_[ind] = coord[4];
    x1dot_[ind] = coord[5];
    x2dot_[ind] = coord[6];
    x3dot_[ind] = coord[7];
    if (parallel_transport_) {
      ep0_[ind] = coord[8];
      ep1_[ind] = coord[9];
      ep2_[ind] = coord[10];
      ep3_[ind] = coord[11];
      et0_[ind] = coord[12];
      et1_[ind] = coord[13];
      et2_[ind] = coord[14];
      et3_[ind] = coord[15];
    }


    if (dir==1) ++imax_; else --imin_;

    if (imin_!=ind && imax_!=ind) {
#     if GYOTO_DEBUG_ENABLED
      GYOTO_DEBUG << "imin_=" << imin_ << ", imax_=" << imax_
		  << ", ind=" << ind << endl;
#     endif
      throwError("BUG: Photon.C: bad index evolution, "
		 "ind should be equal to imin or imax");
    }
    //************************************
    /* 
       3-a
       Call to object_ -> Impact 
    */
    // Check if we can reach the object_
    switch (coordkind) {
    case GYOTO_COORDKIND_SPHERICAL:
      rr = x1_[ind];
      break;
    case GYOTO_COORDKIND_CARTESIAN:
      rr=sqrt(x1_[ind]*x1_[ind]+x2_[ind]*x2_[ind]+x3_[ind]*x3_[ind]);
      break;
    default:
      throwError("Incompatible coordinate kind in Photon.C");
    }

#   if GYOTO_DEBUG_ENABLED
    GYOTO_IF_DEBUG
      GYOTO_DEBUG_EXPR(rmax);
      GYOTO_DEBUG_EXPR(rr);
    GYOTO_ENDIF_DEBUG
#   endif

    if (rr<rmax) {

#     if GYOTO_DEBUG_ENABLED
      GYOTO_DEBUG << "calling Astrobj::Impact\n";
#     endif

      hitt |= object_ -> Impact(this, ind, data);
      if (hitt && !data) stopcond=1;

#     if GYOTO_DEBUG_ENABLED
      GYOTO_DEBUG_EXPR(transmission_freqobs_);
#     endif

      if ( getTransmissionMax() < GYOTO_LIMIT_TRANSMISSION ) {
	stopcond=1;

#       if GYOTO_DEBUG_ENABLED
	GYOTO_DEBUG << "stopping because we are optically thick\n";
#       endif

      }
    } else {
      if ( rr > rr_prev ) {

#       if GYOTO_DEBUG_ENABLED
	GYOTO_DEBUG << "Stopping because "
		    << "1) we are far from this object and "
		    << "2) we are flying away" << endl;
#       endif

	// store coordinates of outgoing photon in impactcoords
	if (data && data->impactcoords && data->impactcoords[0]==DBL_MAX) {
	  for (size_t i=0; i<8; ++i) data->impactcoords[i]=DBL_MAX;
	  memcpy(data->impactcoords+8, &coord[0], 8 * sizeof(double));
	}

	stopcond=1;
	break;
      }
    }
    rr_prev=rr;


    //************************************

    //************************************
    /* 
       3-c Checks whether t < tmin_ (with dir=-1) and expands arrays
       if necessary to be able to store next step's results
    */
    switch (dir) {
    case 1:
      if (coord[0]>tmin_) {
#       if GYOTO_DEBUG_ENABLED
	GYOTO_DEBUG << "stopping because time goes beyond time limit\n";
#       endif
	stopcond=1;
      }
      if ((!stopcond) && (ind==x_size_)) {
	imax_=x_size_-1;
	ind=xExpand(1);
      }
      break;
    default:
      if (coord[0]<tmin_) {
#       if GYOTO_DEBUG_ENABLED
	GYOTO_DEBUG << "stopping because time goes beyond time limit\n";
#       endif
	stopcond=1;
      }
      if ((!stopcond) && (imin_==0)) {
	ind=xExpand(-1);
      }
    }
    //************************************

  }
  // End of stopcond loop
  //-------------------------------------------------

  return hitt;
  
}

double Photon::findMin(Functor::Double_constDoubleArray* object,
		       double t1, double t2, double &tmin,
		       double threshold) {
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << endl;
# endif
  double p1[8] = {t1}, p2[8] = {t2};
  getCoord(p1, 1, p1+1, p1+2, p1+3, p1+4, p1+5, p1+6, p1+7);
  getCoord(p2, 1, p2+1, p2+2, p2+3, p2+4, p2+5, p2+6, p2+7);
  double curval = DBL_MAX, pcur[8], val1, val2;

  pcur[0]=t1;
  getCoord(pcur, 1, pcur+1, pcur+2, pcur+3, pcur+4, pcur+5, pcur+6, pcur+7);
  val1=(*object)(pcur);

  pcur[0]=t2;
  getCoord(pcur, 1, pcur+1, pcur+2, pcur+3, pcur+4, pcur+5, pcur+6, pcur+7);
  val2=(*object)(pcur);

  while ( (fabs(t2-t1)>GYOTO_T_TOL) && (curval>threshold) ) {
    pcur[0] = 0.5*(t1+t2);
    if ((pcur[0]==t1) ||pcur[0]==t2) {
      GYOTO_SEVERE << "Photon::findMin(): dt still above GYOTO_T_TOL (t2-t1="
		   << t2-t1 << ")";
      break;
    }
    getCoord(pcur, 1, pcur+1, pcur+2, pcur+3, pcur+4, pcur+5, pcur+6, pcur+7);
    curval=(*object)(pcur);
    if (val1<val2) {
      t2=pcur[0];
      val2=curval;
    } else {
      t1=pcur[0];
      val1=curval;
    }
  }

  if (val1<val2) {
    tmin=t1;
    return val1;
  }
    
  tmin=t2;
  return val2;

}

void Photon::findValue(Functor::Double_constDoubleArray* object,
		       double value,
		       double tinside, double &toutside) {
  double pcur[8];
  while (fabs(toutside-tinside) > GYOTO_T_TOL) {
    pcur[0] = 0.5*(tinside+toutside);
    getCoord(pcur, 1, pcur+1, pcur+2, pcur+3, pcur+4, pcur+5, pcur+6, pcur+7);
    if ( (*object)(pcur) < value ) tinside = pcur[0];
    else toutside = pcur[0];
  }
  toutside = tinside;
}

void Photon::freqObs(double fo) {
  freq_obs_=fo; 
  GYOTO_DEBUG_EXPR(freq_obs_);
}
double Photon::freqObs() const {
  GYOTO_DEBUG_EXPR(freq_obs_);
  return freq_obs_;
}

double Photon::getTransmission(size_t i) const {
  if (i==size_t(-1)) return transmission_freqobs_;
  if (!spectro_() || i>=spectro_->nSamples())
    throwError("Photon::getTransmission(): i > nsamples");
  return transmission_[i];
}

double Photon::getTransmissionMax() const {
  double transmax=transmission_freqobs_;
  if (spectro_()) {
    transmax=0.;
    size_t i=0, imax= spectro_->nSamples();
    for (i=0; i < imax; ++i)
      if (transmission_[i] > transmax)
	transmax = transmission_[i];
  }
# ifdef GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG_EXPR(transmax);
# endif
  return transmax;
}

double const * Photon::getTransmission() const { return transmission_; }
void Photon::transmit(size_t i, double t) {
  if (i==size_t(-1)) { transmission_freqobs_ *= t; return; }
  if (!spectro_() || i>=spectro_->nSamples())
    throwError("Photon::getTransmission(): i > nsamples");
  transmission_[i] *= t;
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << "(i="<<i<< ", transmission="<<t<<"):"
	      << "transmission_[i]="<< transmission_[i]<< "\n";
# endif
}
void Photon::Refined::transmit(size_t i, double t) {
  parent_->transmit(i, t);
  if (i==size_t(-1)) transmission_freqobs_ = parent_->transmission_freqobs_;
}


#ifdef GYOTO_USE_XERCES
void Photon::setParameters(FactoryMessenger* fmp) {
  wait_pos_ = 1;
  metric(fmp->metric());
  Object::setParameters(fmp);
  wait_pos_ = 0;
  if (init_vel_) {
    delete[] init_vel_; init_vel_=NULL;
    throwError("Worldline::setParameters(): "
	       "Velocity was found but not Position");
  }
}

SmartPointer<Photon> Gyoto::Photon::Subcontractor(FactoryMessenger* fmp) {
  SmartPointer<Photon> ph = new Photon();
  ph -> setParameters(fmp);
  return ph;
}
#endif
