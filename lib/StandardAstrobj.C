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

// GYOTO HEADERS
#include "GyotoUtils.h"
#include "GyotoStandardAstrobj.h"
#include "GyotoMetric.h"
#include "GyotoPhoton.h"
#include "GyotoRegister.h"
#include "GyotoFactoryMessenger.h"

// SYSTEM HEADERS
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <float.h>
#include <cmath>
#include <sstream>

// NAMESPACES
using namespace std;
using namespace Gyoto;
using namespace Gyoto::Astrobj;

Standard::Standard(string kind) :
  Generic(kind),
  critical_value_(DBL_MIN), safety_value_(DBL_MAX)
{
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << endl;
# endif
}

Standard::Standard() :
  Generic(),
  critical_value_(DBL_MIN), safety_value_(DBL_MAX)
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
  critical_value_(orig.critical_value_), safety_value_(orig.safety_value_)
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
  GYOTO_DEBUG_EXPR(getKind());
# endif
  double p1[8], p2[8];
  ph->getCoord(index, p1);
  ph->getCoord(index+1, p2);
  double tmin, minval;

  if (gg_ -> getCoordKind() == GYOTO_COORDKIND_SPHERICAL){
    //Allows theta and phi to be in the correct range
    ph->checkPhiTheta(p1);
    ph->checkPhiTheta(p2);
  }

  double t1 = p1[0], t2=p2[0];
  double val1=(*this)(p1), val2=(*this)(p2);

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
    }
    ph -> findValue(this, critical_value_, t2, t1);
  } else if (val2 > critical_value_)
    ph -> findValue(this, critical_value_, t1, t2);

  double cph[8] = { t2 };
  ph -> getCoord(&t2, 1, cph+1, cph+2, cph+3,
		 cph+4, cph+5, cph+6, cph+7);

  double delta=giveDelta(cph);
  double coh[8];
  while (cph[0]>t1){
    ph -> getCoord(cph, 1, cph+1, cph+2, cph+3,
		   cph+4, cph+5, cph+6, cph+7);
    for (int ii=0;ii<4;ii++) 
      coh[ii] = cph[ii];
    
    getVelocity(coh, coh+4);
    //Next test to insure every point given to process
    //is inside objetc. Not obvious as the worldline between
    //t1 and t2 is not necessarily straight (at small r in particular)
    if ((*this)(coh)<critical_value_)
      processHitQuantities(ph, cph, coh, delta, data);
    cph[0]-=delta;
  }

  return 1;

}

void Standard::setSafetyValue(double val) {safety_value_ = val; }
double Standard::getSafetyValue() const { return safety_value_; }

double Standard::giveDelta(double *) {return 0.05;}

int Standard::setParameter(string name, string content, string unit)  {
  if (name == "SafetyValue") safety_value_ = atof(content.c_str());
  else return Generic::setParameter(name, content, unit);
  return 0;
}

#ifdef GYOTO_USE_XERCES
void Standard::fillElement(FactoryMessenger* fmp) const {
  fmp -> setParameter("SafetyValue", safety_value_);
  Generic::fillElement(fmp);
}
#endif
