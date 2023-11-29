/*
    Copyright 2011-2016, 2018-2020, 2022 Frederic Vincent, Thibaut Paumard

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
#include "GyotoKerrBL.h"
#include "GyotoValue.h"

#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/jacobi_elliptic.hpp>

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
using namespace Eigen;

// don't add e.g. "using namespace boost::math"
// it is not compatible with "using namespace std" in the general case

GYOTO_PROPERTY_START(Photon)
GYOTO_PROPERTY_ASTROBJ(Photon, Astrobj, astrobj)
GYOTO_WORLDLINE_PROPERTY_END(Photon, Object::properties)

Photon::Photon() :
  Worldline(),
  Object("Photon"),
  object_(NULL),
  freq_obs_(1.), transmission_freqobs_(1.),
  spectro_(NULL), transmission_(NULL), nb_cross_eqplane_(0),
  transmissionMatrix_(NULL), transmissionMatrix_freqobs_()
 {}

Photon::Photon(const Photon& o) :
  Worldline(o), SmartPointee(o),
  object_(NULL),
  freq_obs_(o.freq_obs_), transmission_freqobs_(o.transmission_freqobs_),
  spectro_(NULL), transmission_(NULL), nb_cross_eqplane_(o.nb_cross_eqplane_),
  transmissionMatrix_(NULL), transmissionMatrix_freqobs_(o.transmissionMatrix_freqobs_)
{
  if (o.object_()) {
    object_  = o.object_  -> clone();
    object_ -> metric(metric_);
  }
  if (o.spectro_()) {
    spectro_ = o.spectro_ -> clone();
    _allocateTransmission();
    _allocateTransmissionMatrix();
    if (size_t nsamples = spectro_->nSamples()){
      memcpy(transmission_, o.getTransmission(), nsamples*sizeof(double));
      memcpy(transmissionMatrix_,o.getTransmissionMatrix(), nsamples*sizeof(Matrix4d));
    }
  }
}

Photon * Photon::clone() const {
  return new Photon(*this); }

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
  spectro_(orig->spectro_), transmission_(orig->transmission_),
  nb_cross_eqplane_(orig->nb_cross_eqplane_),
  transmissionMatrix_(orig->transmissionMatrix_),
  transmissionMatrix_freqobs_(orig->transmissionMatrix_freqobs_)
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
  Worldline(), freq_obs_(1.), transmission_freqobs_(1.), spectro_(NULL), transmission_(NULL), nb_cross_eqplane_(0),
  transmissionMatrix_(NULL), transmissionMatrix_freqobs_()
{
  setInitialCondition(met, obj, coord);
}

Photon::Photon(SmartPointer<Metric::Generic> met, 
	       SmartPointer<Astrobj::Generic> obj, 
	       SmartPointer<Screen> screen, 
	       double d_alpha, double d_delta):
  Worldline(), object_(obj), freq_obs_(screen->freqObs()),
  transmission_freqobs_(1.),
  spectro_(NULL), transmission_(NULL),
  nb_cross_eqplane_(0),
  transmissionMatrix_(NULL),
  transmissionMatrix_freqobs_()
{
  double coord[8], Ephi[4], Etheta[4];
  bool compute_polar_basis=false;
  if (parallel_transport_) compute_polar_basis=true;
  screen -> getRayTriad(d_alpha, d_delta,
			coord,
			compute_polar_basis,
			Ephi, Etheta);
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

void Photon::_allocateTransmissionMatrix() {
  if (transmissionMatrix_){
    delete [] transmissionMatrix_;
    transmissionMatrix_ = NULL;
  }
  if (spectro_()){
    size_t nsamples = spectro_->nSamples();
    if (nsamples) {
      transmissionMatrix_ = new Matrix4d[nsamples]();
    }
    resetTransmissionMatrix();
  }
}

void Photon::resetTransmission() {
  transmission_freqobs_ = 1.;
  if (spectro_() && transmission_) {
    size_t nsamples = spectro_->nSamples();
    for (size_t i = 0; i < nsamples; ++i) transmission_[i] = 1.;
  }  
}

void Photon::resetTransmissionMatrix() {
  if (spectro_() && transmissionMatrix_ && &transmissionMatrix_freqobs_) {
    size_t nsamples = spectro_->nSamples();
    Matrix4d identity;
    identity << 1, 0, 0, 0,
                0, 1, 0, 0,
                0, 0, 1, 0,
                0, 0, 0, 1;
    transmissionMatrix_freqobs_ = identity;
    for (size_t i = 0; i < nsamples; ++i){
      transmissionMatrix_[i]=identity;
    }
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
  _allocateTransmissionMatrix();
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
  bool compute_polar_basis=false;
  if (parallel_transport_) compute_polar_basis=true;
  screen -> getRayTriad(d_alpha, d_delta,
			coord,
			compute_polar_basis,
			Ephi, Etheta);
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

  /*transmission_freqobs_=1.;
  size_t nsamples;
  if (spectro_() && (nsamples = spectro_->nSamples()))
    for (size_t ii=0; ii<nsamples; ++ii) transmission_[ii]=1.;
  */
  resetTransmission();
  resetTransmissionMatrix();

  double rmax=object_ -> rMax();
  int coordkind = metric_ -> coordKind();

  int hitt=0;
  //hitted=1 if object is hitted at least one time (hitt can be 0 even
  //if the object was hit if cross_max>0) ; hitt_crude=1 if the object
  //is not yet hit with adaptive integration step. A second integration
  //with small fixed step will be performed to determine more precisely
  //the surface point.
  state_t coord(parallel_transport_?16:8);
  double tau;
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
      GYOTO_ERROR("Incompatible coordinate kind in Photon.C");
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
  tau=tau_[ind];
  switch (coordkind) {
  case GYOTO_COORDKIND_SPHERICAL:
    rr=coord[1];
    break;
  case GYOTO_COORDKIND_CARTESIAN:
    rr=sqrt(coord[1]*coord[1]+coord[2]*coord[2]+coord[3]*coord[3]);
    break;
  default:
    GYOTO_ERROR("Incompatible coordinate kind in Photon.C");
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

  /*
    Image order computation:
    
    An image can be decomposed into various orders depending on the
    number of turns that a photon executes around a black hole before
    reaching the observer (this is broadly speaking, details below).
    The 0th-order (aka direct) image is the map of specific intensity 
    obtained by integrating the radiative transfer along the part of 
    the geodesic that makes between 0 and half a turn around the BH,
    when integrating from the observer towards the BH. The 1st-order
    image is obtained by taking into account the intensity accumulated
    between half a turn and one full turn, etc. Again this is broadly
    speaking and just to have a first naive view. In particular, note that
    one given geodesic is not "of order n", but it rather has some part of it
    which is of order 0, some part of order 1, etc.

    In more details now, there are at least 2 ways that one can encode this:
    - we can keep track of the theta turning points along the geodesic,
    and change the geodesic order when we encounter a theta turning point
    (corresponding to a half orbit around the BH, where an orbit is defined here
    as the geodesic path between one theta turning point back to this
    same turning point).
    This is reasonably easy to implement (but there are few tricks, which
    makes it less obviously correct than the second method) and is 
    spacetime agnostic.
    - a more precise way of dealing with that, which is specific to Kerr,
    is to keep track of the Mino time as defined eg in Gralla & Lupsasca,
    Phys Rev D, 101, 044031 (2020, hereafter GL20a), Eq. 34-36. 
    The Mino time is directly linked to the number of orbits around the BH. 
    So here there is no trick, we simply integrate the
    Mino time, and use it to define the number of orbits, and define the
    image order as above as the number of half orbits. However, this is
    not generalizable to a non-Kerr spacetime.

    I (FV) have compared the two approaches and found 
    excellent agreement in Kerr.

    The image-order code below is thus decomposed into two parts, one 
    which takes care of keeping track of the theta turning points, and
    the other that takes care of integrating the Mino time. One can switch
    between the two by choosing compute_Mino=0 or 1 below (hardcoded).

    Careful: all this is well tested at low inclination but might
    not be well adapted for high inclination, TBC.
   */

  int geodesic_order=0; // will be changed to 1 when the geodesic becomes
  // 1st-order, then to 2 when it becomes second-order.

  // Theta-turning-point image-order-tracking parameters:
  bool theta_has_changed=false; // flag checking the change of theta coordinate
  // needed for tracking various orders of the image if max_cross_eqplane_
  // is set
  bool theta_is_increasing=true; // 1: if theta is increasing along integration,
  // ie from Screen towards Object; 0: means theta is decreasing. Flag
  // needed for tracking various orders of the image if max_cross_eqplane_
  // is set. Default value here, will be updated later.
  int nb_theta_turningpoints=0; // keeping track of number of theta turning points
  // Mino-time image-order-tracking parameters:
  int compute_Mino=0; // IMPORTANT! if put to 1, the theta turning points
  // tracking is ignored and the image-order tracking is done by means
  // of the Mino time
  double G_theta=0.; // angular integral along geodesic = elapsed Mino time
  double I_radial=0.; // radial integral along geodesic = elapsed Mino time also
  double G_theta_1libration=0., nn_Mino_1st_turning=0., KK=0., Fo=0., sign_o=0.,
    uplus=0., uminus=0., uratio=0., eta=0., lambda=0., spin2=0., nn_Mino=0., Einit=0., Linit=0., spin=0.;
  double mytol=1e-6; // it must be possible to decrease this tolerance
  // by using a finer AbsTol and RelTol in xml (checked for one case).
  if ((maxCrossEqplane_<DBL_MAX || (data && data->nbcrosseqplane)) && compute_Mino==1){ // The first two conditions ensure that we are interesting in computing image orders; then we compute here everything that is constant along the null geodesic
    string kin = metric_->kind(); // check that we are in Kerr for Mino computation
    if (kin != "KerrBL")
      GYOTO_ERROR("Photon::hit: KerrBL needed for Mino time computation!");
    double pos[4]={coord[0],coord[1],coord[2],coord[3]};
    double g_tt=metric_->gmunu(pos,0,0), g_thth=metric_->gmunu(pos,2,2),
      g_pp=metric_->gmunu(pos,3,3), g_tp=metric_->gmunu(pos,0,3); 
    double p_t = g_tt*coord[4] + g_tp*coord[7],
      p_ph = g_tp*coord[4] + g_pp*coord[7],
      p_th = g_thth*coord[6];
    //cout << "E,L= " << -p_t << " " << p_ph << endl;
    double theta = coord[2];

    spin= metric_ -> get("Spin");
    // Note: static_cast<SmartPointer<Metric::KerrBL> >(metric_) -> spin()
    // does not work because Photon belongs to libgyoto and KerrBL to
    // libgyoto-stdplug so there is no possible cast between them.
    // The get function is more time consuming, keep that in mind...
    
    spin2 = spin*spin;
    //cout << "pth= " << p_th << endl;
    double Carter = p_th*p_th - cos(theta)*cos(theta)*(spin2*p_t*p_t
						       - p_ph*p_ph/(sin(theta)*sin(theta)));
    Einit = -p_t; Linit = p_ph;
    //cout << setprecision(16) << "E,L init= " << Einit <<  " " << Linit << endl;
    if (Carter<0.) {
      //GYOTO_WARNING << "Carter<0" << endl;
      compute_Mino=0.; // the Mino time expressions of GL20a
      // are restricted to positive Carter cst geodesics (most of them are)
    }
    lambda = -p_ph/p_t;
    eta = Carter/(p_t*p_t);
    //cout << "eta,lambda= " << eta << " " << lambda << endl;

    if (spin==0.) throwError("Photon::hit: Mino-time formalism derived for spin>0!");

    double delta_th = 0.5*(1.-(eta+lambda*lambda)/spin2);
    //cout<< "delta_th= "<< delta_th << endl;
    uplus = delta_th + sqrt(delta_th*delta_th + eta/spin2);
    uminus = delta_th - sqrt(delta_th*delta_th + eta/spin2);
    if (uminus==0.) throwError("Photon::hit: uminus is zero and we need "
			       "to divide by that number!");
    uratio = uplus/uminus; // note that uratio is always <0 for the
    // cases of interest where Carter>0 (then, uminus<0, uplus>0).

    /*
      IMPORTANT NOTE ON ELLIPTIC INTEG:
      F_GrallaLupsasca(phi,k) = ellint_1(sqrt(k),phi) in C++
      careful to the sqrt(k); 
      also note that the C++ implementation only deals with k>0. 
      For k<0 we use Eq. 17.4.17 of Abramowitz & Stegun which expresses
      F(phi,-k) as a function of F(phi,k/(1+k)),
      hence the sqrt(-uratio/(1-uratio)) appearing regularly below,
      because uratio<0.
     */
    
    if (compute_Mino==1){ // compute_Mino may have been changed to zero if
                          // Carter cst is <0
      KK = 1./sqrt(1-uratio)*(boost::math::ellint_1(sqrt(-uratio/(1-uratio)),M_PI/2.)-boost::math::ellint_1(sqrt(-uratio/(1-uratio)),0.)); // using Abramowitz&Stegun formula, this is
      // KK = F(pi/2 | uratio)
      G_theta_1libration = 4.*KK/(spin*sqrt(-uminus)); // Eq. 35 of GL20a
      sign_o = (x2dot_[i0_] > 0) - (x2dot_[i0_] < 0); // this is the sign of the theta component of the photon momentum at the observer

      // Computing F(arcsin(cos(theta_obs)/sqrt(uplus) | uratio):
      double argasin=cos(x2_[i0_])/sqrt(uplus); // uplus is the maximum allowed value of cos^2 theta, so argasin should be <1; it might not be, due to limited precision, cos^2 theta being slightly above uplus. So we take care of this:
      if (argasin-1>0. && argasin-1<mytol)
	argasin=1.;
      if (-1-argasin>0. && -1-argasin<mytol)
	argasin=-1.;
      double Fobs = 1./sqrt(1-uratio)*(boost::math::ellint_1(sqrt(-uratio/(1-uratio)),M_PI/2.)-boost::math::ellint_1(sqrt(-uratio/(1-uratio)),M_PI/2.-asin(argasin))); // using Abramowitz&Stegun formula, this is F(arcsin(cos(theta_obs)/sqrt(uplus) | uratio)

      // Computing F(arcsin(cos(theta_turning)/sqrt(uplus) | uratio):
      double theta_turn=acos(sign_o*sqrt(uplus)); // first turning point encountered along backwards integration
      argasin=cos(theta_turn)/sqrt(uplus); // same stuff as above
      if (argasin-1>0. && argasin-1<mytol)
	argasin=1.;
      if (-1-argasin>0. && -1-argasin<mytol)
	argasin=-1.;
      double Fturn = 1./sqrt(1-uratio)*(boost::math::ellint_1(sqrt(-uratio/(1-uratio)),M_PI/2.)-boost::math::ellint_1(sqrt(-uratio/(1-uratio)),M_PI/2.-asin(argasin)));

      // From these two elliptic integrals, get the value of G_theta
      // at the first turning point:
      double G_theta_1st_turning = sign_o*(Fobs-Fturn)*(-1./(spin*sqrt(-uminus)));    // using Eq. 29 and 36 of Gralla&Lupsasca, PRD, 101, 044032 (2020),
      // hereafter GL20b

      nn_Mino_1st_turning = G_theta_1st_turning/G_theta_1libration;
      // fractional orbit number at first crossing, see GL20a Eq. 34
      
      if (nn_Mino_1st_turning<0.1) nn_Mino_1st_turning+=0.5; // I impose the
      // fractional orbit number to be of order 0.5 at first crossing; if it
      // is not, it means that there is a theta crossing far away from
      // the BH and I don't want to take it into account. See FV notes on
      // image-order tracking for details on these "useless turning points".
      
    }

  }
  double h1max=DBL_MAX;
  while (!stopcond) {
    // Next step along photon's worldline
    h1max=object_ -> deltaMax(&coord[0]);
    stopcond  = state_ -> nextStep(coord, tau, h1max);
    //cout << "IN ph r, z= " << coord[1] << " " << coord[1]*cos(coord[2]) << endl;
    double pos[4]={coord[0],coord[1],coord[2],coord[3]};
    double g_tt=metric_->gmunu(pos,0,0), g_thth=metric_->gmunu(pos,2,2),
      g_pp=metric_->gmunu(pos,3,3), g_tp=metric_->gmunu(pos,0,3); 
    double p_t = g_tt*coord[4] + g_tp*coord[7],
      p_ph = g_tp*coord[4] + g_pp*coord[7],
      p_th = g_thth*coord[6];
    double Ecur=-p_t, Lcur=p_ph;
    //cout << setprecision(16) << "E,L cur= " << -p_t <<  " " << fabs(Einit+p_t) << " " << p_ph << " " << fabs(Linit-p_ph) << endl;
    
    if (!secondary_){ // to compute only primary image (outdated, use MaxCrossEqplane above instead)
      // Thin disk case
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
    xStore(ind, coord, tau);

    if (dir==1) ++imax_; else --imin_;

    if (imin_!=ind && imax_!=ind) {
#     if GYOTO_DEBUG_ENABLED
      GYOTO_DEBUG << "imin_=" << imin_ << ", imax_=" << imax_
		  << ", ind=" << ind << endl;
#     endif
      GYOTO_ERROR("BUG: Photon.C: bad index evolution, "
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
      GYOTO_ERROR("Incompatible coordinate kind in Photon.C");
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
    
    //*****************************
    // *** Image order tracking ***
    //*****************************
    if (maxCrossEqplane_<DBL_MAX || (data && data->nbcrosseqplane)){

      if (compute_Mino==1){
	/*
	  DETERMINE GEODESIC ORDER BY MEANS OF MINO TIME
	  (specific to Kerr)
	 */
	
	double theta_step = -coord[2] + x2_[ind-dir]; // it is indeed like that because we are integrating backwards, the order matters a lot; Gth is the integral from source to observer.
	double costh_cur, sinth_cur; sincos(coord[2], &sinth_cur, &costh_cur);
	double costh_prev, sinth_prev; sincos(x2_[ind-dir], &sinth_prev,
					      &costh_prev);
	double c2_cur=costh_cur*costh_cur, s2_cur=sinth_cur*sinth_cur,	
	  c2_prev=costh_prev*costh_prev, s2_prev=sinth_prev*sinth_prev;
	double Theta_cur = eta + spin2*c2_cur - lambda*lambda*c2_cur/s2_cur;
	double Theta_prev = eta + spin2*c2_prev - lambda*lambda*c2_prev/s2_prev;
	double pth_cur = coord[6], rr_cur= coord[1], mytheta_cur = coord[2];
	double Sigma = rr_cur*rr_cur + spin2*cos(mytheta_cur)*cos(mytheta_cur);
	double Delta = rr_cur*rr_cur + spin2 - 2.*rr_cur;
	double RR_cur = (rr_cur*rr_cur + spin2 - spin*lambda)*(rr_cur*rr_cur + spin2 - spin*lambda) - Delta*(eta + (lambda-spin)*(lambda-spin));
	double pphi_rhs = spin/Delta*(rr_cur*rr_cur + spin2 - spin*lambda) + lambda/(sin(mytheta_cur)*sin(mytheta_cur)) - spin;
	double pt_rhs = (rr_cur*rr_cur + spin2)/Delta*(rr_cur*rr_cur + spin2 - spin*lambda) + spin*(lambda-spin*sin(mytheta_cur)*sin(mytheta_cur));
	double pth_prev = Ecur/Sigma*sqrt(Theta_cur),
	  pr_prev = Ecur/Sigma*sqrt(RR_cur),
	  pt_prev = Ecur/Sigma*pt_rhs,
	  pphi_prev = Ecur/Sigma*pphi_rhs;

	//cout << "Potentials= " << RR_cur <<" " << Theta_cur << endl;

	// Next if loop is for computing the I_r integral in addition
	// to the G_theta integral (see GL20). Both should agree.
	// Only useful for testing.
	int compute_Irad=0;
	if (compute_Irad==1){
	  double Arond = spin2 - eta - lambda*lambda,
	    Brond = 2.*(eta + (lambda - spin)*(lambda - spin)),
	    Crond = -spin2*eta,
	    Prond = -Arond*Arond/12. - Crond,
	    Qrond = -Arond/3.*(Arond*Arond/36. - Crond) - Brond*Brond/8.,
	    delta_cubic = Prond*Prond*Prond/27. + Qrond*Qrond/4.;
	  
	  double xi0;
	  if (delta_cubic<0.){
	    complex<double> omegaplus =
	      pow(-Qrond/2.
		  + complex<double>(0, sqrt(-delta_cubic)),1./3.),
	      omegaminus = pow(-Qrond/2. -
			       complex<double>(0, sqrt(-delta_cubic)),1./3.);
	    if (imag(omegaplus)+imag(omegaminus)!=0.)
	      throwError("Bad cubic roots!");
	    xi0 = real(omegaplus) + real(omegaminus) - Arond/3.;
	  }else{
	    double omegaplus = pow(-Qrond/2. + sqrt(delta_cubic),1./3.),
	      omegaminus = pow(-Qrond/2. - sqrt(delta_cubic),1./3.);
	    xi0 = omegaplus + omegaminus - Arond/3.;
	  }
	  if (xi0<0.) throwError("Bad cubic roots!");
	  double zquartic = sqrt(xi0/2.),
	    r1quartic =
	    -zquartic
	    - sqrt(-Arond/2. - zquartic*zquartic + Brond/(4.*zquartic)),
	    r2quartic =
	    -zquartic
	    + sqrt(-Arond/2. - zquartic*zquartic + Brond/(4.*zquartic)),
	    r3quartic =
	    zquartic
	    - sqrt(-Arond/2. - zquartic*zquartic - Brond/(4.*zquartic)),
	    r4quartic =
	    zquartic
	    + sqrt(-Arond/2. - zquartic*zquartic - Brond/(4.*zquartic));
	  if (fabs(r1quartic+r2quartic+r3quartic+r4quartic)>1e-10)
	    throwError("Bad quartic roots!"); // sum should be zero
	  
	  // Here I assume a type-1 radial evolution, ie r1<r2<r<r3<r4
	  double r32 = r3quartic - r2quartic,
	    r41 = r4quartic - r1quartic,
	    r31 = r3quartic - r1quartic,
	    r42 = r4quartic - r2quartic,
	    kradial = r32*r41/(r31*r42),
	    x2radial_cur = sqrt((rr_cur - r4quartic)/(rr_cur - r3quartic)
				* r31/r41),
	    x2radial_prev = sqrt((x1_[ind-dir] - r4quartic)
				 /(x1_[ind-dir] - r3quartic) * r31/r41);
	  
	  if (kradial<0. || kradial >1. || x2radial_cur<0.
	      || x2radial_cur>1. ||  x2radial_prev<0. || x2radial_prev>1.)
	    throwError("Bad radial potential params");
	  
	  double F2radial_cur = 2./sqrt(r31*r42)
	    *boost::math::ellint_1(sqrt(kradial),asin(x2radial_cur)),
	    F2radial_prev = 2./sqrt(r31*r42)
	    *boost::math::ellint_1(sqrt(kradial),asin(x2radial_prev));
	  double sign_radial_s_cur = (coord[5] > 0) - (coord[5] < 0);
	  double sign_radial_s_prev = (x1dot_[ind-dir] > 0)
	    - (x1dot_[ind-dir] < 0);
	  double sign_radial_s=0., dI_radial=0.;
	  if (sign_radial_s_cur == sign_radial_s_prev) {
	    sign_radial_s = sign_radial_s_cur;
	    dI_radial = sign_radial_s*(F2radial_prev - F2radial_cur);
	  }else{
	    double F2radial_turn = 0.; //2./sqrt(r31*r42)*boost::math::ellint_1(sqrt(kradial),0.); // at r=r4quartic, x2radial is zero, and asin(0)=0, hence the 0 in second slot, which leads to ellint1=0
	    dI_radial = 2.*sign_radial_s_cur*F2radial_turn - sign_radial_s_cur*(F2radial_prev + F2radial_cur);
	  }
	  if (dI_radial<0.) throwError("radial integral should increase");
	  I_radial += dI_radial;

	  // Equatorial plane crossing
	  double argasin=cos(x2_[i0_])/sqrt(uplus);
	  if (argasin-1>0. && argasin-1<mytol)
	    argasin=1.;
	  if (-1-argasin>0. && -1-argasin<mytol)
	    argasin=-1.;
	  double F_obs = 1./sqrt(1-uratio)
	    *(boost::math::ellint_1(sqrt(-uratio/(1-uratio)),M_PI/2.)
	      -boost::math::ellint_1(sqrt(-uratio/(1-uratio)),M_PI/2.-asin(argasin)));
	  double G_theta_eqplane_0turnings=(-sign_o*F_obs)
	    *(1./(spin*sqrt(-uminus))),
	    G_theta_eqplane_1turnings=(2.*KK-sign_o*F_obs)
	    *(1./(spin*sqrt(-uminus))),
	    G_theta_eqplane_2turnings=(2.*2.*KK-sign_o*F_obs)
	    *(1./(spin*sqrt(-uminus)));

	  // Computing all equat crossing when after the other,
	  // erasing the jacobi integrals from one to the other,
	  // only the req_* matter.
	  // Crossing after zero turning point:
	  double jacobi_elliptic_sine =
	    boost::math::jacobi_sn(sqrt(kradial),
		      0.5*sqrt(r31*r42)*G_theta_eqplane_0turnings
		      - boost::math::ellint_1(sqrt(kradial),asin(sqrt(r31/r41)))),
	    jacobi_elliptic_sine_squared =
	    jacobi_elliptic_sine*jacobi_elliptic_sine,
	    req_0turnings = (r4quartic*r31
			     - r3quartic*r41*jacobi_elliptic_sine_squared)
	    /(r31 - r41*jacobi_elliptic_sine_squared);

	  // crossing after 1 turning
	  jacobi_elliptic_sine =
	    boost::math::jacobi_sn(sqrt(kradial),
		      0.5*sqrt(r31*r42)*G_theta_eqplane_1turnings
		      - boost::math::ellint_1(sqrt(kradial),asin(sqrt(r31/r41))));
	  jacobi_elliptic_sine_squared =
	    jacobi_elliptic_sine*jacobi_elliptic_sine;
	  double  req_1turnings =
	    (r4quartic*r31 - r3quartic*r41*jacobi_elliptic_sine_squared)
	    /(r31 - r41*jacobi_elliptic_sine_squared);

	  // crossing after 2 turnings
	  jacobi_elliptic_sine =
	    boost::math::jacobi_sn(sqrt(kradial),
		      0.5*sqrt(r31*r42)*G_theta_eqplane_2turnings
		      - boost::math::ellint_1(sqrt(kradial),asin(sqrt(r31/r41))));
	  jacobi_elliptic_sine_squared =
	    jacobi_elliptic_sine*jacobi_elliptic_sine;
	  double req_2turnings =
	    (r4quartic*r31 - r3quartic*r41*jacobi_elliptic_sine_squared)
	    /(r31 - r41*jacobi_elliptic_sine_squared);
	  
	  //cout << "equat turnings= " << req_0turnings << " " << req_1turnings << " " << req_2turnings << endl;
	} // end of Irad if loop
	
	// Back to main loop, computing G_theta:

	// Same stuff as before essentially
	// Computing elliptic integral F at current point
	double argasin=cos(coord[2])/sqrt(uplus); 
	if (argasin-1>0. && argasin-1<mytol)
	  argasin=1.;
	if (-1-argasin>0. && -1-argasin<mytol)
	  argasin=-1.;
	double Fs_cur = 1./sqrt(1-uratio)*(boost::math::ellint_1(sqrt(-uratio/(1-uratio)),M_PI/2.)-boost::math::ellint_1(sqrt(-uratio/(1-uratio)),M_PI/2.-asin(argasin)));
	// Computing elliptic integral F at previous point
	argasin=cos(x2_[ind-dir])/sqrt(uplus);
	if (argasin-1>0. && argasin-1<mytol)
	  argasin=1.;
	if (-1-argasin>0. && -1-argasin<mytol)
	  argasin=-1.;
	double Fs_prev = 1./sqrt(1-uratio)*(boost::math::ellint_1(sqrt(-uratio/(1-uratio)),M_PI/2.)-boost::math::ellint_1(sqrt(-uratio/(1-uratio)),M_PI/2.-asin(argasin)));//boost::math::ellint_1(uratio,asin(argasin));

	// From the two Fs compute the increment of G_theta
	double sign_s_cur = (coord[6] > 0) - (coord[6] < 0);
	double sign_s_prev = (x2dot_[ind-dir] > 0) - (x2dot_[ind-dir] < 0);
	double sign_s=0., dG_theta=0.;
	if (sign_s_cur == sign_s_prev) {
	  sign_s = sign_s_cur;
	  dG_theta = sign_s*(Fs_prev - Fs_cur)*(-1./(spin*sqrt(-uminus)));
	}else{
	  double Fs_turn = G_theta_1libration/4.;
	  if (coord[2]<M_PI/2.) Fs_turn*=-1.;
	  dG_theta = 2*sign_s_cur*Fs_turn
	    - sign_s_cur*(Fs_prev + Fs_cur)*(-1./(spin*sqrt(-uminus)));
	}
	
	G_theta += dG_theta;
	// and deduce fractional orbit number:
	nn_Mino = G_theta/G_theta_1libration;

	// Now simply update the geodesic order from nn:
	if (geodesic_order==0 && nn_Mino>nn_Mino_1st_turning) {
	  //cout << setprecision(16)<< "***MINO: geod becomes n=1 at rc, z= " << coord[1]*sin(coord[2]) << " " << coord[1]*cos(coord[2]) << ", and nn_Mino= "<< nn_Mino << endl;
	  geodesic_order=1;
	}
	if (geodesic_order==1 && nn_Mino>nn_Mino_1st_turning+0.5) {
	  //cout << "***MINO: geod becomes n=2 at rc, z= " << coord[1]*sin(coord[2]) << " " << coord[1]*cos(coord[2]) << ", and nn_Mino= "<< nn_Mino << endl;
	  geodesic_order=2;
	}
	if (geodesic_order==2 && nn_Mino>nn_Mino_1st_turning+1.) {
	  //cout << "***MINO: geod becomes n=3 at rc, z= " << coord[1]*sin(coord[2]) << " " << coord[1]*cos(coord[2]) << ", and nn_Mino= "<< nn_Mino << endl;
	  geodesic_order=3;
	}

	//cout <<"geodesic order= " << geodesic_order << endl;
	
      }else{

	/*
	  DETERMINE GEODESIC ORDER BY MEANS OF THETA TURNING POINTS
	  (not specific to Kerr)

	  For details on the tricks, see FV notes on image order tracking.
	 */
	
	double zsign=0.;
	switch (coordkind) {
	case GYOTO_COORDKIND_SPHERICAL:
	  // *** z sign change tracking
	  zsign = x1_[i0_]*cos(x2_[i0_]); // sign of first z position
	  if (nb_cross_eqplane_>0) zsign *= pow(-1,nb_cross_eqplane_); // update it when crossing equatorial plane
	  if (coord[1]*cos(coord[2])*zsign<0.){ 
	    nb_cross_eqplane_+=1; // equatorial plane has been just crossed
	    //cout << "***updating nbcross to " << nb_cross_eqplane_ << endl;
	    //cout << "at rc,z= " << coord[1]*sin(coord[2]) << " " << coord[1]*cos(coord[2]) << endl;
	  }
	  
	  // *** theta turning points tracking
	  if (!theta_has_changed){
	    // theta has not yet changed
	    // just keep track of it starting to change
	    // and update theta_is_increasing accordingly
	    if (x2_[i0_]!=coord[2]){
	      theta_has_changed=true;
	      if (coord[2]<x2_[i0_]) theta_is_increasing=false;
	    }
	  }else{
	    // theta is now changing
	    // keep track of turning points
	    if ((theta_is_increasing && coord[2]<x2_[ind-dir])
		|| (!theta_is_increasing && coord[2]>x2_[ind-dir])){
	      // so here theta was increasing in the past of the integration
	      // and it now starts to decrease, or the other way round:
	      // we have a new turning point
	      theta_is_increasing = !theta_is_increasing;
	      if (nb_theta_turningpoints==nb_cross_eqplane_-1){
		// eqplane should be crossed before the theta turning point
		nb_theta_turningpoints+=1;
		if (nb_theta_turningpoints==1) {
		  //cout << "***TURNING: geod becomes n=1 at rc, z= " << coord[1]*sin(coord[2]) << " " << coord[1]*cos(coord[2]) << endl;
		  geodesic_order=1;
		}
		if (nb_theta_turningpoints==2){
		  //cout << "***TURNING: geod becomes n=2 at rc, z= " << coord[1]*sin(coord[2]) << " " << coord[1]*cos(coord[2]) << endl;
		  geodesic_order=2;
		}
		if (nb_theta_turningpoints==3){
		  //cout << "***TURNING: geod becomes n=3 at rc, z= " << coord[1]*sin(coord[2]) << " " << coord[1]*cos(coord[2]) << endl;
		  geodesic_order=3;
		}
	      }
	    }
	    
	  }
	  
	  break;
	case GYOTO_COORDKIND_CARTESIAN:
	  {
	    throwError("to be implemented");

	    break;
	  }
	default:
	  GYOTO_ERROR("Incompatible coordinate kind in Photon.C");
	}
      } // END OF IMAGE ORDER STUFF
      
      GYOTO_DEBUG_EXPR(nb_cross_eqplane_);

      // Store if necessary the image order stuff in the relevant data:
      if (data->nbcrosseqplane && compute_Mino==0 &&
	  nb_cross_eqplane_==nb_theta_turningpoints){
      	*data->nbcrosseqplane=nb_cross_eqplane_;
      }

      if (data->nbcrosseqplane && compute_Mino==1)
	*data->nbcrosseqplane=geodesic_order;

      if (data->user1 &&
	  geodesic_order==0) {
	if (data->spectrum)
	  *data->user1=data->spectrum[0]; // User1 keep tracks of Spectrum until the geodesic becomes of order 1
	if (data->user4)
	  *data->user1=data->user4[0]; 
	//cout << "saving order 0 spectrum of " << data->spectrum[0] << " " << *data->user1 << " at rc,z=" << coord[1]*sin(coord[2]) << " " << coord[1]*cos(coord[2]) << endl;
      }
      if (data->user2 &&
	  geodesic_order==1) {
	if (data->spectrum)
	  *data->user2=data->spectrum[0]-*data->user1;
	if (data->user4)
	  *data->user2=data->user4[0]-*data->user1;
	//cout << "saving order 1 spectrum of " << data->spectrum[0] << " " << *data->user1 << " " << *data->user2 << " at rc,z=" << coord[1]*sin(coord[2]) << " " << coord[1]*cos(coord[2]) << endl;
      }
      if (data->user3 &&
	  geodesic_order==2) {
	if (data->spectrum)
	  *data->user3=data->spectrum[0]-*data->user1-*data->user2;
	if (data->user4)
	  *data->user3=data->user4[0]-*data->user1-*data->user2;
	//cout << "saving order 2 spectrum of " << *data->user3 << " at rc,z=" << coord[1]*sin(coord[2]) << " " << coord[1]*cos(coord[2]) << endl;
      }
      if (geodesic_order == maxCrossEqplane_){
	// stop integration if maxcross reached
	return 0;
      }
    }    

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

void Photon::nb_cross_eqplane(int nb) {
  nb_cross_eqplane_=nb; 
  GYOTO_DEBUG_EXPR(nb_cross_eqplane_);
}
int Photon::nb_cross_eqplane() const {
  GYOTO_DEBUG_EXPR(nb_cross_eqplane_);
  return nb_cross_eqplane_;
}

double Photon::getTransmission(size_t i) const {
  if (i==size_t(-1)) return transmission_freqobs_;
  if (!spectro_() || i>=spectro_->nSamples())
    GYOTO_ERROR("Photon::getTransmission(): i > nsamples");
  return transmission_[i];
}

Matrix4d Photon::getTransmissionMatrix(size_t i) const {
  if (i==size_t(-1)) return transmissionMatrix_freqobs_;
  if (!spectro_() || i>=spectro_->nSamples())
    GYOTO_ERROR("Photon::getTransmission(): i > nsamples");
  return transmissionMatrix_[i];
}

double Photon::getTransmissionMax() const {
  double transmax = 0.;
  if (parallel_transport_ and spectro_()){
    transmax=0.;
    size_t i=0, imax= spectro_->nSamples();
    for (i=0; i < imax; ++i){
      Matrix4d mat=transmissionMatrix_[i];
      if (mat(0,0)>transmax)
        transmax=mat(0,0);
    }
  }else{
    transmax=transmission_freqobs_;
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
  }
  return transmax;
}

double const * Photon::getTransmission() const { return transmission_; }
Matrix4d const * Photon::getTransmissionMatrix() const { return transmissionMatrix_; }

void Photon::transmit(size_t i, double t) {
  if (i==size_t(-1)) { transmission_freqobs_ *= t; return; }
  if (!spectro_() || i>=spectro_->nSamples())
    GYOTO_ERROR("Photon::getTransmission(): i > nsamples");
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

void Photon::transmit(size_t i, Matrix4d mat){
  if (i==size_t(-1)) { transmissionMatrix_freqobs_ *= mat; return; }
  if (!spectro_() || i>=spectro_->nSamples())
    GYOTO_ERROR("Photon::getTransmissionMatrix(): i > nsamples");
  transmissionMatrix_[i] *= mat;
}
void Photon::Refined::transmit(size_t i, Matrix4d mat) {
  parent_ -> transmit(i, mat);
}

void Photon::transfer(double * Inu, double * Qnu, double * Unu, double * Vnu, Matrix4d * Onu){
  // Apply transfer function to I, Q, U and V, then update the transfer function.
  // For the prototype,
  //   * just apply the transmission to Inu;
  //   * only update transmission.
  size_t nbnuobs = spectro_() ? spectro_->nSamples() : 0;
  Matrix4d Tau;
  for (size_t ii=0; ii<nbnuobs; ++ii) {
    Tau=transmissionMatrix_[ii];
    Vector4d Stokes(Inu[ii],Qnu[ii],Unu[ii],Vnu[ii]);
    Vector4d StokesOut = Tau*Stokes;
    Inu[ii]=StokesOut(0);
    Qnu[ii]=StokesOut(1);
    Unu[ii]=StokesOut(2);
    Vnu[ii]=StokesOut(3);
    transmit(ii, Onu[ii]);
  }
}

void Photon::Refined::transfer(double * Inu, double * Qnu, double * Unu, double * Vnu, Matrix4d * Onu){
  parent_ -> transfer(Inu, Qnu, Unu, Vnu, Onu);
}

#ifdef GYOTO_USE_XERCES
void Photon::setParameters(FactoryMessenger* fmp) {
  wait_pos_ = 1;
  metric(fmp->metric());
  Object::setParameters(fmp);
  wait_pos_ = 0;
  if (init_vel_) {
    delete[] init_vel_; init_vel_=NULL;
    GYOTO_ERROR("Worldline::setParameters(): "
	       "Velocity was found but not Position");
  }
}

SmartPointer<Photon> Gyoto::Photon::Subcontractor(FactoryMessenger* fmp) {
  SmartPointer<Photon> ph = new Photon();
  ph -> setParameters(fmp);
  return ph;
}
#endif
