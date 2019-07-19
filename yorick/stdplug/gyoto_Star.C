/*
    Copyright 2011, 2013, 2018 Thibaut Paumard

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

#include "../ygyoto.h"
#include "yapi.h"
#include "pstdlib.h"
#include <GyotoStar.h>
#include <GyotoStarTrace.h>
#ifdef GYOTO_USE_XERCES
#include <GyotoFactory.h>
#endif

#include <iostream>
#include <cstring>
using namespace std;

using namespace Gyoto;
using namespace Gyoto::Astrobj;

#define OBJ ao

// on_eval worker
void ygyoto_Star_eval(SmartPointer<Astrobj::Generic>* ao_, int argc) {
  GYOTO_DEBUG << endl;

  // Define keywords
  static char const * knames[]={
    "unit",
    "radius", "metric", "initcoord", "spectrum", "opacity", "delta", "adaptive",
    "deltamaxoverradius", "deltamaxoverdistance",
    "maxiter",
    "integrator", "deltamin", "deltamax", "deltamaxoverr", "abstol", "reltol",
    "reset", "xfill",
    YGYOTO_ASTROBJ_GENERIC_KW,
    "get_skypos", "get_txyz", "get_prime", "get_coord", "get_cartesian",
    "startrace",
    0
  };

  YGYOTO_WORKER_INIT(Astrobj, Star, knames, YGYOTO_ASTROBJ_GENERIC_KW_N+25);

  YGYOTO_WORKER_SET_UNIT;
  YGYOTO_WORKER_GETSET_DOUBLE2_UNIT(radius); 
  YGYOTO_WORKER_GETSET_OBJECT2(metric,Metric);

  /* INITCOORD */
  if ((iarg=kiargs[++k])>=0) { //initcoord
    iarg+=*rvset;
    if (yarg_nil(iarg)) {
      if ((*rvset)++) y_error(rmsg);
      state_t v(8);
      (*ao)->getInitialCoord(v);
      long dims[]= {1, 8};
      double *coord = ypush_d(dims);
      for (int i=0; i<8; ++i) coord[i]=v[i];
    } else {
      long ptot=1, vtot=1;
      double *pos=ygeta_d(iarg, &ptot, 0);
      if (ptot<4) y_error("POS should have at least 4 elements");
      int n;
      double coord[8];
      for (n=0; n<4; ++n) coord[n]=pos[n];

      double *v=NULL, tdot0;

      if (yarg_number(piargs[0]+*rvset)) {
	if ((*paUsed)++) y_error(pmsg);
	v=ygeta_d(piargs[0]+*rvset, &vtot, 0);
	if (vtot!=3) y_error("V should have 3 elements");
	if (!(*ao)->metric())
	  y_error("Please set metric before setting initial condition");
	tdot0=(*ao)->metric()->SysPrimeToTdot(pos, v);
	coord[4]=tdot0;
	for (n=0; n<3; ++n) coord[5+n]=v[n]*tdot0;
      } else if (ptot==8) {
	for (n=4; n<8; ++n) coord[n]=pos[n];
      } else if (ptot==7) {
	v=pos+4;
	tdot0=(*ao)->metric()->SysPrimeToTdot(pos, v);
	coord[4]=tdot0;
	for (n=0; n<3; ++n) coord[5+n]=v[n]*tdot0;
      } else y_error("Not enough information to set initial condition");

      (*ao)->setInitialCondition(coord);
    }
  }
 
  YGYOTO_WORKER_GETSET_OBJECT2(spectrum,Spectrum);
  YGYOTO_WORKER_GETSET_OBJECT2(opacity,Spectrum);
  YGYOTO_WORKER_GETSET_DOUBLE2_UNIT(delta);
  YGYOTO_WORKER_GETSET_LONG2(adaptive);
  YGYOTO_WORKER_GETSET_DOUBLE2(deltaMaxOverRadius);
  YGYOTO_WORKER_GETSET_DOUBLE2(deltaMaxOverDistance);
  YGYOTO_WORKER_GETSET_LONG2( maxiter );
  YGYOTO_WORKER_GETSET_STRING2( integrator );
  YGYOTO_WORKER_GETSET_DOUBLE2( deltaMin );
  YGYOTO_WORKER_GETSET_DOUBLE2( Worldline::deltaMax );
  YGYOTO_WORKER_GETSET_DOUBLE2( deltaMaxOverR );
  YGYOTO_WORKER_GETSET_DOUBLE2( absTol );
  YGYOTO_WORKER_GETSET_DOUBLE2( relTol );
  YGYOTO_WORKER_RUN( reset() );
  YGYOTO_WORKER_RUN( xFill(ygets_d(iarg+*rvset), false) );
 
  YGYOTO_WORKER_CALL_GENERIC(Astrobj);

  k += YGYOTO_ASTROBJ_GENERIC_KW_N;

  // SPECIFIC GET KEYWORDS
  if ((iarg=kiargs[++k])>=0) { // get skypos
    if ((*rvset)++) y_error(rmsg);
    if ((*paUsed)++) y_error(pmsg);

    if (!yarg_Screen(iarg)) y_error("Expecting gyoto_Screen argument");
    SmartPointer<Screen> *screen = yget_Screen(iarg);

    int nel =(*ao)->get_nelements();
    long dims[] = {2, nel, 3};
    double * data=ypush_d(dims);

    (*ao)->getSkyPos(*screen, data, data+nel, data+2*nel);
  }

  if ((iarg=kiargs[++k])>=0) { // get_txyz
    if ((*rvset)++) y_error(rmsg);
    int nel =(*ao)->get_nelements();
      
    long dims[] = {2, nel, 4};
    double * data=ypush_d(dims);

    (*ao)->get_t(data);
    (*ao)->get_xyz(data+nel, data+2*nel, data+3*nel);
  }

  if ((iarg=kiargs[++k])>=0) { // get_prime
    if ((*rvset)++) y_error(rmsg);
    int nel =(*ao)->get_nelements();
      
    long dims[] = {2, nel, 3};
    double * data=ypush_d(dims);

    (*ao)->get_prime(data, data+nel, data+2*nel);
  }

  if ((iarg=kiargs[++k])>=0) { // get_coord
    if ((*rvset)++) y_error(rmsg);
    // Two cases : get_coord=1 or get_coord=dates. 
    // We will recognize the first case if the parameter is the
    // integer scalar 1, any other type or value is a date
    // specification.
    int argt = yarg_number(iarg);
    long ntot[1] = {1};
    long dims[Y_DIMSIZE];
    double* dates = NULL;
    double dummy[]= {1.};
    if (yarg_nil(iarg)) {
      // void is same as 1 in this context
      argt=1;
      *ntot=1;
      dims[0]=0;
      dates=dummy;
    } else dates = ygeta_d(iarg, ntot, dims);

    if (debug())
      cerr << "DEBUG: gyoto_Star(get_coord=array("
	   << (argt==1?"long":"double")
	   << ", " << *ntot 
	   << ") (rank=" << dims[0]
	   << ", dates[0]=" << dates[0] << ")" << endl;

    if (dims[0] == 0 && argt == 1 && dates[0] == 1) {
      if(debug())
	cerr << "DEBUG: retrieving all already computed coordinates" << endl;
      int nel =(*ao)->get_nelements();
      long ddims[] = {2, nel, 8};
      double * data = ypush_d(ddims);

      (*ao)->getCoord(data, data+nel, data+2*nel, data+3*nel);
      (*ao)->get_dot(data+4*nel, data+5*nel, data+6*nel, data+7*nel);
    } else {
      if(debug())
	cerr << "DEBUG: retrieving coordinates for specified date(s)"<< endl;
      if (dims[0] > Y_DIMSIZE-2)
	y_errorn("gyoto_Star(get_coord=dates): DATES must have at most %ld "
		 "dimensions", Y_DIMSIZE-2);
      dims[0] += 1;
      dims[dims[0]] = 7;
      double * data = ypush_d(dims);
      size_t nel = *ntot;
      (*ao) -> getCoord(dates, nel,
			data, data+ nel,data+2* nel,data+3* nel,
			data+4* nel, data+5* nel, data+6* nel);
    }
  }
  if ((iarg=kiargs[++k])>=0) { // get_cartesian
    if ((*rvset)++) y_error(rmsg);
    long ntot[1] = {1};
    long dims[Y_DIMSIZE];
    double* dates = ygeta_d(iarg, ntot, dims);
    if(debug())
      cerr << "DEBUG: gyoto_Star(get_cartesian=dates)"<< endl;
    if (dims[0] > Y_DIMSIZE-2)
      y_errorn("gyoto_Star(get_cartesian=dates): DATES must have at most %ld "
	       "dimensions", Y_DIMSIZE-2);
    dims[0] += 1;
    dims[dims[0]] = 6;
    double * data = ypush_d(dims);
    size_t nel = *ntot;
    (*ao) -> getCartesian(dates, nel,
			  data,         data +   nel, data + 2*nel,
			  data + 3*nel, data + 4*nel, data + 5*nel);
    
  }

  /* STARTRACE */
  if ((iarg=kiargs[++k])>=0) {
    if (*rvset)  y_error(rmsg);
    if (*paUsed) y_error(pmsg);
    double tmin=ygets_d(iarg);
    double tmax=ygets_d(piargs[0]);
    ++(*rvset);
    ++(*paUsed);
     *ypush_Astrobj() = new StarTrace(**ao, tmin, tmax);
  }
  
  
}



extern "C" {
  void Y__gyoto_Star_register_as_Astrobj(){
    ygyoto_Astrobj_register("Star",&ygyoto_Star_eval);
  }

  // Generic constructor/accessor
  void
  Y_gyoto_Star(int argc)
  {
    GYOTO_DEBUG << endl;
    YGYOTO_CONSTRUCTOR_INIT2(Astrobj, Astrobj::Generic, Star, astrobj);
    if ((*ao)->kind().compare("Star"))
      y_error("Expecting Astrobj of kind Star");
    ygyoto_Star_eval(ao, argc);
  }

}
