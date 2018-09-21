/*
    Copyright 2011, 2012, 2014, 2015, 2017, 2018 Thibaut Paumard & Frederic Vincent

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

// GYOTO HEADERS
#include "GyotoUtils.h"
#include "GyotoStandardAstrobj.h"
#include "GyotoMetric.h"
#include "GyotoPhoton.h"
#include "GyotoRegister.h"
#include "GyotoFactoryMessenger.h"
#include "GyotoProperty.h"

// SYSTEM HEADERS
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <cstring>
#include <float.h>
#include <cmath>
#include <sstream>

// NAMESPACES
using namespace std;
using namespace Gyoto;
using namespace Gyoto::Astrobj;

GYOTO_PROPERTY_START(Gyoto::Astrobj::Standard,
  "Gyoto::Astrobj whose shape is defined by a scalar function.")
GYOTO_PROPERTY_DOUBLE(Standard, SafetyValue, safetyValue,
  "Value of the function below which to look more carefully.")
GYOTO_PROPERTY_DOUBLE(Standard, DeltaInObj, deltaInObj,
		      "Value of the constant integration step "
		      "inside the astrobj (geometrical units)")
GYOTO_PROPERTY_END(Standard, Generic::properties)

Standard::Standard(string kin) :
  Generic(kin),
  critical_value_(DBL_MIN), safety_value_(DBL_MAX),
  delta_inobj_(0.05)
{
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << endl;
# endif
}

Standard::Standard() :
  Generic(),
  critical_value_(DBL_MIN), safety_value_(DBL_MAX),
  delta_inobj_(0.05)
{
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << endl;
# endif
}

Standard::Standard(double radmax) :
  Generic(radmax),
  critical_value_(DBL_MIN), safety_value_(DBL_MAX)
{
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << endl;
# endif
}

Standard::Standard(const Standard& orig) :
  Generic(orig), Functor::Double_constDoubleArray(orig),
  critical_value_(orig.critical_value_), safety_value_(orig.safety_value_),
  delta_inobj_(orig.delta_inobj_)
{
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << endl;
# endif
}

Standard::~Standard() {
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << endl;
# endif
}

int Standard::Impact(Photon* ph, size_t index, Properties *data){
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG_EXPR(kind());
# endif
  size_t sz = ph -> parallelTransport()?16:8;
  state_t p1(sz), p2(sz);
  ph->getCoord(index, p1);
  ph->getCoord(index+1, p2);
  double tmin, minval;

  if (gg_ -> coordKind() == GYOTO_COORDKIND_SPHERICAL){
    //Allows theta and phi to be in the correct range
    ph->checkPhiTheta(&p1[0]);
    ph->checkPhiTheta(&p2[0]);
  }

  double t1 = p1[0], t2=p2[0];
  double val1=(*this)(&p1[0]), val2=(*this)(&p2[0]);

  if (val1 > critical_value_) {
    if (val2 > critical_value_) {
      if ( val1 > safety_value_ && val2 > safety_value_) {
	if (val1 < val2) {
	  minval = val1; tmin = t1;
	} else {
	  minval = val2; tmin = t2;
	}
      } else
	minval = ph -> findMin(this, p1[0], p2[0], tmin, critical_value_) ;
      if (minval>critical_value_) {
	if (data) {
	  /* EmissionTime */
	  if (data->time) *data->time=tmin;
	  /* MinDistance */
	  if ((data->distance) && (*(data->distance)>minval) )
	    *data->distance=minval;
	  /* FirstMinDist */
	  if (data->first_dmin) { 
	    if (!data->first_dmin_found) {
	      if (*(data->first_dmin)>minval) *(data->first_dmin)=minval;
	      else data->first_dmin_found=1;
	    }
	  }
	}
	return 0;
      }
      ph -> findValue(this, critical_value_, tmin, t2);
    } else tmin=t2;
    ph -> findValue(this, critical_value_, tmin, t1);
  } else if (val2 > critical_value_)
    ph -> findValue(this, critical_value_, t1, t2);

  state_t cph(sz);
  ph -> getCoord(t2, cph);

  double coh[8] = {cph[0], cph[1], cph[2], cph[3]};
  getVelocity(coh, coh+4);
  bool current_is_inside = true; // by construction, t2 is always inside

  double delta=giveDelta(&cph[0]);
  double dt;
  double coh_next[8];
  state_t cph_next(sz);
  bool next_is_inside;

  while (cph[0]>t1){
    // Warning: Impact must not extend the Worldline!
    // never call get Coord with anything outside [t1, t2].
    ph -> getCoord(max(cph[0] - delta, t1), cph_next);

    memcpy(coh_next, &cph_next[0], 4*sizeof(cph_next[0]));
    getVelocity(coh_next, coh_next+4);

    next_is_inside = ((*this)(coh_next) <= critical_value_);

    if (current_is_inside) {
      if (next_is_inside) {
	// Both points are inside
	dt = cph[0]-cph_next[0];
      } else {
	// Late point in object, early outside
	// Find date of surface crossing and update dt.
	double t_out = cph_next[0];
	ph -> findValue(this, critical_value_, cph[0], t_out);
	dt = cph[0] - t_out;
      }
    } else {
      if (next_is_inside) {
	// Early point in object, late point outside
	// Place cph and coh inside, near surface; update dt
	double t_out = cph[0];
	ph -> findValue(this, critical_value_, cph_next[0], t_out);
	ph -> getCoord(t_out, cph);
	memcpy(coh, &cph[0], 4*sizeof(cph_next[0]));
	getVelocity(coh, coh+4);
	dt = cph[0]-cph_next[0];
      } else {
	// Two points outside
	dt = 0.;
      }
    }

    // dt == 0 means the two points are outside
    if (dt != 0.)
      processHitQuantities(ph, cph, coh, dt, data);

    // Copy next to current
    cph = cph_next;
    memcpy(coh, coh_next, 8*sizeof(coh[0]));
    current_is_inside = next_is_inside;

  }

  return 1;

}

void Standard::safetyValue(double val) {safety_value_ = val; }
double Standard::safetyValue() const { return safety_value_; }

double Standard::deltaInObj() const { return delta_inobj_; }
void   Standard::deltaInObj(double val) { delta_inobj_ = val; }

double Standard::giveDelta(double *) { return deltaInObj(); }
