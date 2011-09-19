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
#include "GyotoFactoryMessenger.h"
#include "GyotoPhoton.h"
#include "GyotoScreen.h"
#include "GyotoWorldlineIntegState.h"
#include "GyotoDefs.h"
#include "GyotoError.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <cstring>
#include <cstdlib>
#include <cfloat>

#define COUNT_MAX 100000


using namespace std;
using namespace Gyoto;

Photon::Photon() :
  Worldline(), freq_obs_(1.),
  transmission_freqobs_(1.), spectro_(NULL), transmission_(NULL)
 {}

Photon::Photon(const Photon& o) :
  Worldline(o), SmartPointee(o),
  object_(NULL),
  freq_obs_(o.freq_obs_), transmission_freqobs_(o.transmission_freqobs_),
  spectro_(NULL), transmission_(NULL)
{
  if (o.object_())  object_  = o.object_  -> clone();
  if (o.spectro_()) {
    spectro_ = o.spectro_ -> clone();
    _allocateTransmission();
    if (size_t nsamples = spectro_->getNSamples())
      memcpy(transmission_, o.getTransmission(), nsamples*sizeof(double));
  }
}

Photon * Photon::clone() const { return new Photon(*this); }

Photon::Photon(SmartPointer<Metric::Generic> met, SmartPointer<Astrobj::Generic> obj,
	       double* coord):
  Worldline(), transmission_freqobs_(1.), spectro_(NULL), transmission_(NULL)
{
  setInitialCondition(met, obj, coord);
}

Photon::Photon(SmartPointer<Metric::Generic> met, SmartPointer<Astrobj::Generic> obj, 
	       SmartPointer<Screen> screen, double d_alpha, double d_delta):
  Worldline(), object_(obj), transmission_freqobs_(1.),
  spectro_(NULL), transmission_(NULL)
{
  double coord[8];
  screen -> getRayCoord(d_alpha, d_delta, coord);
  Worldline::setInitialCondition(met, coord, -1);
  setSpectrometer(screen);
}

Photon::~Photon() {}

/* TRANSMISSION STUFF */
void Photon::_allocateTransmission() {
  if (transmission_) {
    delete [] transmission_;
    transmission_ = NULL;
  }
  if (spectro_()) {
    size_t nsamples = spectro_->getNSamples();
    if (nsamples) {
      transmission_ = new double[nsamples];
      resetTransmission();
    }
  }
}

void Photon::resetTransmission() {
  transmission_freqobs_ = 1.;
  if (spectro_() && transmission_) {
    size_t nsamples = spectro_->getNSamples();
    for (size_t i = 0; i < nsamples; ++i) transmission_[i] = 1.;
  }
}


double Photon::getMass() const { return 0.; }

void Photon::setAstrobj(SmartPointer<Astrobj::Generic> ao) {
  imin_=imax_=i0_;
  object_=ao;
}

void Photon::setSpectrometer(SmartPointer<Spectrometer> spr) {
  spectro_=spr;
  _allocateTransmission();
}
SmartPointer<Spectrometer> Photon::getSpectrometer() const { return spectro_; }


string Photon::className() const { return  string("Photon"); }
string Photon::className_l() const { return  string("photon"); }

SmartPointer<Astrobj::Generic> Photon::getAstrobj() const { return object_; }



void Photon::setInitialCondition(SmartPointer<Metric::Generic> met,
				 SmartPointer<Astrobj::Generic> obj,
				 SmartPointer<Screen> screen,
				 const double d_alpha,
				 const double d_delta)
{
  double coord[8];
  screen -> getRayCoord(d_alpha, d_delta, coord);
  Worldline::setInitialCondition(met, coord, -1);
  object_=obj;

}

void Photon::setInitialCondition(SmartPointer<Metric::Generic> met,
				 SmartPointer<Astrobj::Generic> obj,
				 const double coord[8])
{
  
  /*if(debug()) cout << coord[0] << " "
       << coord[1] << " "
       << coord[2] << " "
       << coord[3] << " "
       << coord[4] << " "
       << coord[5] << " "
       << coord[6] << " "
       << coord[7] << endl;*/
  double gtt0=met->gmunu(coord,0,0);
  double ObsVel[4]={sqrt(-1./gtt0),0.,0.,0.};
  double sp_rec=met->ScalarProd(coord,coord+4,ObsVel);
  freq_obs_ = -sp_rec;
  //  cout << "ut= " << sqrt(-1./gtt0) << endl;
  //  cout << "In Photon.C: freq obs= " << freq_obs_ << endl;
  Worldline::setInitialCondition(met, coord, -1);
  object_=obj;
}

int Photon::hit(Astrobj::Properties *data) {

  /*
    Ray-tracing of the photon until the object_ is hit. Radiative
    transfer inside the object_ may then be performed depending on
    flag_radtransf. Final result (observed flux for instance,
    depending on object_'s Astrobj::Properties) will be stored in data.
   */

  //tlim_=-1000.;//DEBUG //NB: integration stops when t < Worldline::tlim_

  transmission_freqobs_=1.;
  size_t nsamples;
  if (spectro_() && (nsamples = spectro_->getNSamples()))
    for (size_t ii=0; ii<nsamples; ++ii) transmission_[ii]=1.;

  double rmax=object_ -> getRmax();
  int coordkind = metric_ -> getCoordKind();

  int hitt=0;
  //hitted=1 if object is hitted at least one time (hitt can be 0 even
  //if the object was hit if cross_max>0) ; hitt_crude=1 if the object
  //is not yet hit with adaptive integration step. A second integration
  //with small fixed step will be performed to determine more precisely
  //the surface point.
  double coord[8];
  int dir=(tlim_>x0_[i0_])?1:-1;
  size_t ind=i0_;
  int stopcond=0;
  double rr=DBL_MAX, rr_prev=DBL_MAX;

  //-------------------------------------------------
  /*
    1-
    Call object_->Impact on the already computed part 
    of the geodesic to check whether the object_ was hit.
   */
  for (ind=i0_+dir;
       ((dir==1)?
	(ind<=imax_ && x0_[ind]<=tlim_): // conditions if dir== 1
	(ind>=imin_ && x0_[ind]>=tlim_)) // conditions if dir==-1
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
  }
  if (rr<rmax)
    hitt = object_ -> Impact(this, ind, data);
  if (hitt) {
    if (debug()) cout << "Photon.C: Hit for already computed position; "
		      << "Warning: radiative transfer not implemented "
		      << "for that case" << endl;
    return hitt;
  } else if (((dir==1)?
	(ind==imax_ && x0_[ind]>=tlim_): // conditions if dir== 1
	(ind>=imin_ && x0_[ind]<=tlim_)) // conditions if dir==-1
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

  SmartPointer<WorldlineIntegState> state
    = new WorldlineIntegState(metric_, coord, delta_* dir);
  //delta_ = initial integration step (defaults to 0.01)

  size_t count=0;// Must remain below count_max (prevents infinite integration)

  //-------------------------------------------------

  //-------------------------------------------------
  /*
    3- Integration loop: integrate the geodesic until stopcond is 1.
    Possible stopping conditions: 
    - transmission_freqobs_ low
    - t < tlim_ (if dir=-1) 
    - t does not evolve
    - metric tells it's time to stop
    - count>count_max (should never be used, just to prevent infinite
    integration in case of a bug)
   */

  while (!stopcond) {
    // Next step along photon's worldline
    stopcond  = state -> nextStep(this, coord);
    if (stopcond) {
      if (debug()) cerr<<"DEBUG: Photon::hit(): stopcond set by integrator\n";
      break;
    }
    if (coord[0] == x0_[ind]) { // here, ind denotes previous step
      stopcond=1;
      if (verbose() >= GYOTO_SEVERE_VERBOSITY)
	cerr << "SEVERE: Photon::hit(): time did not evolve, break." << endl;
      break;
    }
    if((stopcond=metric_->isStopCondition(coord))) {
      if (debug()) cerr << "DEBUG: Photo::hit(): stopcond step by metric"<<endl;
      break;
    }

    if ( ++count > COUNT_MAX ) {
      if (verbose()>= GYOTO_SEVERE_VERBOSITY)
	cerr << "***WARNING (severe): Photon::hit: too many iterations, break"
	     << endl;
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


    if (dir==1) ++imax_; else --imin_;

    if (imin_!=ind && imax_!=ind) {
      if (debug()) cerr << "\nimin_=" << imin_
			<< ", imax_=" << imax_
			<< ", ind=" << ind << endl;
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

    if (debug())
      cerr << "DEBUG: Photon::hit(): rmax="<< rmax <<", rr="<<rr<<endl;
    if (rr<rmax) {
      if (debug()) cerr << "DEBUG: Photon::hit() calling Astrobj::Impact\n";
      hitt = object_ -> Impact(this, ind, data);
      if (debug()) cerr << "DEBUG: Photon::hit(): transmission_freqobs_="
			<< transmission_freqobs_ << endl;
      if ( transmission_freqobs_ < 1e-6 ) {
	stopcond=1;
	if (debug()) cerr << "DEBUG: Photon::hit(): stopping because we "
			  << "are optically thick\n";
      }
    } else {
      if ( rr > rr_prev ) {
	if (debug())
	  cerr << "DEBUG: Photon::hit(): Stopping because "
	       << "1) we are far from this object and "
	       << "2) we are flying away" << endl;
	stopcond=1;
	break;
      }
    }
    rr_prev=rr;


    //************************************

    //************************************
    /* 
       3-c Checks whether t < tlim_ (with dir=-1) and expands arrays
       if necessary to be able to store next step's results
    */
    switch (dir) {
    case 1:
      if (coord[0]>tlim_) {
	stopcond=1;
      }
      if ((!stopcond) && (ind==x_size_)) {
	imax_=x_size_-1;
	ind=xExpand(1);
      }
      break;
    default:
      if (coord[0]<tlim_) {
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

double Photon::findMin(Astrobj::Generic* object,
		       double t1, double t2, double &tmin,
		       double threshold) {
  if (debug())
    cerr << "DEBUG: in Photon::findMind()\n";
  double p1[4] = {t1}, p2[4] = {t2};
  getCoord(p1, 1, p1+1, p1+2, p1+3);
  getCoord(p2, 1, p2+1, p2+2, p2+3);
  double curval = DBL_MAX, pcur[4], val1, val2;

  pcur[0]=t1;
  getCoord(pcur, 1, pcur+1, pcur+2, pcur+3);
  val1=(*object)(pcur);

  pcur[0]=t2;
  getCoord(pcur, 1, pcur+1, pcur+2, pcur+3);
  val2=(*object)(pcur);

  while ( (fabs(t2-t1)>GYOTO_T_TOL) && (curval>threshold) ) {
    pcur[0] = 0.5*(t1+t2);
    getCoord(pcur, 1, pcur+1, pcur+2, pcur+3);
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

void Photon::findValue(Astrobj::Generic* object, double value,
		       double tinside, double &toutside) {
  double pcur[4];
  while (fabs(toutside-tinside) > GYOTO_T_TOL) {
    pcur[0] = 0.5*(tinside+toutside);
    getCoord(pcur, 1, pcur+1, pcur+2, pcur+3);
    if ( (*object)(pcur) < value ) tinside = pcur[0];
    else toutside = pcur[0];
  }
  toutside = tinside;
}

double Photon::getFreqObs() const { return freq_obs_; }

double Photon::getTransmission(size_t i) const {
  if (i==size_t(-1)) return transmission_freqobs_;
  if (!spectro_() || i>=spectro_->getNSamples())
    throwError("Photon::getTransmission(): i > nsamples");
  return transmission_[i];
}
double const * Photon::getTransmission() const { return transmission_; }
void Photon::transmit(size_t i, double t) {
  if (i==size_t(-1)) { transmission_freqobs_ *= t; return; }
  if (!spectro_() || i>=spectro_->getNSamples())
    throwError("Photon::getTransmission(): i > nsamples");
  transmission_[i] *= t;
  if (debug())
    cerr << "DEBUG: Photon::transmit(i="<<i<< ", transmission="<<t<<"):"
	 << "transmission_[i]="<< transmission_[i]<< "\n";
}

#ifdef GYOTO_USE_XERCES
void Photon::fillElement(FactoryMessenger *fmp) {
  if (metric_)     fmp -> setMetric (metric_) ;
  if (object_)    fmp -> setAstrobj (object_) ;

  double coord[8];
  getInitialCoord(coord);
  fmp -> setParameter("InitCoord", coord, 8);

  if (delta_ != GYOTO_DEFAULT_DELTA)
    fmp -> setParameter ("Delta", delta_);
}

SmartPointer<Photon> Gyoto::PhotonSubcontractor(FactoryMessenger* fmp) {

  string name="", content="";
  SmartPointer<Metric::Generic> gg = NULL;
  SmartPointer<Astrobj::Generic> ao = NULL;

  SmartPointer<Photon> ph = new Photon();
  ph -> setMetric(  fmp->getMetric() );
  ph -> setAstrobj( fmp->getAstrobj() );


  while (fmp->getNextParameter(&name, &content)) {
    char* tc = const_cast<char*>(content.c_str());
    if(name=="Delta") ph -> setDelta( atof(tc) );
    if(name=="InitCoord") {
      double coord[8];
      for (int i=0;i<8;++i) coord[i] = strtod(tc, &tc);
      ph -> setInitCoord(coord, -1);
    }
  }

  return ph;
}
#endif
