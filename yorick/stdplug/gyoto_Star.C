/*
    Copyright 2011 Thibaut Paumard

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

#include "ygyoto.h"
#include "yapi.h"
#include <GyotoStar.h>
#ifdef GYOTO_USE_XERCES
#include <GyotoFactory.h>
#endif

#include <iostream>
using namespace std;

using namespace Gyoto;
using namespace Gyoto::Astrobj;

// on_eval worker
void ygyoto_Star_eval(Gyoto::SmartPointer<Gyoto::Astrobj::Generic>* ao_, int argc) {
  if (debug()) cerr << "in ygyoto_Star_eval" << endl;
  int rvset[1]={0}, paUsed[1]={0}, constructor=0;

  // If needed, create the object.
  if (!ao_) { // Constructor mode
    constructor=1;
    ao_ = ypush_Astrobj();
  } else *ypush_Astrobj()=*ao_;

  // Parse arguments
  static char const * knames[]={
    "radius", "metric", "initcoord", "spectrum", "opacity", "reset", "xfill",
    YGYOTO_ASTROBJ_GENERIC_KW,
    "get_skypos", "get_txyz", "get_coord", "get_cartesian",
    0
  };
#define nkw 11
  static long kglobs[YGYOTO_ASTROBJ_GENERIC_KW_N+nkw+1];
  int kiargs[YGYOTO_ASTROBJ_GENERIC_KW_N+nkw];
  int piargs[]={-1,-1,-1,-1};
  yarg_kw_init(const_cast<char**>(knames), kglobs, kiargs);
  int iarg=argc, parg=0;
  while (iarg>=1) {
    iarg = yarg_kw(iarg, kglobs, kiargs);
    if (iarg>=1) {
      if (parg<4) piargs[parg++]=iarg--;
      else y_error("gyoto_Astrobj takes at most 4 positional arguments");
    }
  }
  if (debug()) {
    cerr <<   "DEBUG: gyoto_Star parameters:" << endl;
    for (int i=0; i<4; ++i)
      cerr << "DEBUG:         piargs[" << i << "]="
	   << piargs[i] << endl;
    for (int i=0; i<YGYOTO_ASTROBJ_GENERIC_KW_N+nkw; ++i)
      cerr << "DEBUG:         kiargs[" << i << "]="
	   << kiargs[i] 
	   << " (" << knames[i] << ")" << endl;
  }

  // Constructor mode from XML file
  if (constructor) {
    if (yarg_string(piargs[0])) {
#ifdef GYOTO_USE_XERCES
      try { *ao_ = Factory(ygets_q(piargs[0])).getAstrobj(); } 
      YGYOTO_STD_CATCH;
      *paUsed=1;
#else
      y_error("This GYOTO was compiled without XERCES: no xml i/o");
#endif
    } else *ao_ = new Star();
  }
  SmartPointer<Star> *ao = (SmartPointer<Star> *)ao_;

  // Process specific keywords
  int k=-1;
  char const * rmsg="Cannot set return value more than once";
  char const * pmsg="Cannot use positional argument more than once";

  //// MEMBERS ////
  /* RADIUS */
  if ((iarg=kiargs[++k])>=0) {
    iarg+=*rvset;
    if (yarg_nil(iarg)) {
      if ((*rvset)++) y_error(rmsg);
      ypush_double((*ao)->getRadius());
    } else
      (*ao)->setRadius(ygets_d(iarg)) ;
  }

  /* METRIC */
  if ((iarg=kiargs[++k])>=0) {
    iarg+=*rvset;
    if (yarg_nil(iarg)) {
      if ((*rvset)++) y_error(rmsg);
      *ypush_Metric() = (*ao)->getMetric();
    } else
      (*ao)->setMetric(*yget_Metric(iarg)) ;
  }

  /* INITCOORD */
  if ((iarg=kiargs[++k])>=0) { //initcoord
    iarg+=*rvset;
    if (yarg_nil(iarg)) {
      if ((*rvset)++) y_error(rmsg);
      long dims[]= {1, 8};
      double *coord = ypush_d(dims);
      (*ao)->getInitialCoord(coord);
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
	if (!(*ao)->getMetric())
	  y_error("Please set metric before setting initial condition");
	tdot0=(*ao)->getMetric()->SysPrimeToTdot(pos, v);
	coord[4]=tdot0;
	for (n=0; n<3; ++n) coord[5+n]=v[n]*tdot0;
      } else if (ptot==8) {
	for (n=4; n<8; ++n) coord[n]=pos[n];
      } else if (ptot==7) {
	v=pos+4;
	tdot0=(*ao)->getMetric()->SysPrimeToTdot(pos, v);
	coord[4]=tdot0;
	for (n=0; n<3; ++n) coord[5+n]=v[n]*tdot0;
      } else y_error("Not enough information to set initial condition");

      (*ao)->setInitialCondition(coord);
    }
  }
 
  /* SPECTRUM */
  if ((iarg=kiargs[++k])>=0) {
    if (yarg_nil(iarg)) {
      SmartPointer<Spectrum::Generic> * sp = ypush_Spectrum();
      *sp = (*ao) -> getSpectrum();
    } else {
      (*ao) -> setSpectrum ( *yget_Spectrum(iarg) );
    }
  }

  /* OPACITY */
  if ((iarg=kiargs[++k])>=0) {
    if (yarg_nil(iarg)) {
      SmartPointer<Spectrum::Generic> * sp = ypush_Spectrum();
      *sp = (*ao) -> getOpacity();
    } else {
      (*ao) -> setOpacity ( *yget_Spectrum(iarg) );
    }
  }

  //// METHODS ////
  if ((iarg=kiargs[++k])>=0) (*ao)->reset();
  if ((iarg=kiargs[++k])>=0) (*ao)->xFill(ygets_d(iarg+*rvset));
 
  // GENERIC WORKER
    ygyoto_Astrobj_generic_eval(ao_, kiargs+k+1, piargs, rvset, paUsed);

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

  if (yarg_true(kiargs[++k])) { // get_txyz
    if ((*rvset)++) y_error(rmsg);
    int nel =(*ao)->get_nelements();
      
    long dims[] = {2, nel, 4};
    double * data=ypush_d(dims);

    (*ao)->get_t(data);
    (*ao)->get_xyz(data+nel, data+2*nel, data+3*nel);
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
    double* dates = ygeta_d(iarg, ntot, dims);

    if (debug())
      cerr << "DEBUG: gyoto_Star(get_coord=array("
	   << (yarg_number(iarg)==1?"long":"double")
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
  
  
}



extern "C" {
  void Y__gyoto_Star_register_as_Astrobj(){
    ygyoto_Astrobj_register("Star",&ygyoto_Star_eval);
  }

  // Generic constructor/accessor
  void
  Y_gyoto_Star(int argc)
  {
    if (debug()) cerr << "In Y_gyoto_Star" << endl;
    SmartPointer<Astrobj::Generic> *ao = NULL;
    if (yarg_Astrobj(argc-1)) {
      ao = yget_Astrobj(--argc);
      if ((*ao)->getKind().compare("Star"))
	y_error("Expecting Astrobj of kind Star");
    }
    ygyoto_Star_eval(ao, argc);
  }

  // STAR CLASS
  void
  Y_gyoto_Star_get_t(int n)
  {
    if (n!=1) y_error("gyoto_Star_get_t takes exactly 1 argument");
    SmartPointer<Astrobj::Generic> *astrobj=yget_Astrobj(0);
    if (strcmp((*astrobj)->getKind().c_str(),"Star")) y_error("first argument must be a GYOTO Star object ");
    SmartPointer<Star> st (*astrobj);
    int nel ;
    try { nel = st->get_nelements(); }
    YGYOTO_STD_CATCH ;
    long dims[] = {1, nel};
    double * t=ypush_d(dims);
    try { st->get_t(t); }
    YGYOTO_STD_CATCH ;
  }

  void
  Y_gyoto_Star_xFill(int n)
  {
    if (n!=2) y_error("gyoto_Star_xFill takes exactly 2 argument");
    SmartPointer<Astrobj::Generic> *astrobj=yget_Astrobj(n-1);
    if (strcmp((*astrobj)->getKind().c_str(),"Star")) y_error("first argument must be a GYOTO Star object ");
    double tlim=ygets_d(n-2);
    SmartPointer<Star> st (*astrobj);
    try { st->xFill(tlim); }
    YGYOTO_STD_CATCH ;
  }

  void
  Y_gyoto_Star_getSkyPos(int n)
  {
    if (n<2) y_error("gyoto_Star_get_xyz takes at least 3 argument");
    if (n>4) y_error("gyoto_Star_get_xyz takes at most 5 arguments");
    SmartPointer<Astrobj::Generic> *astrobj=yget_Astrobj(n-1);
    if (strcmp((*astrobj)->getKind().c_str(),"Star")) y_error("first argument must be a GYOTO Star object ");

    SmartPointer<Screen> *screen = yget_Screen(n-2);

    long da_ref=yget_ref(n-3);
    long dd_ref=yget_ref(n-4);
    long dD_ref=yget_ref(n-5);

    SmartPointer<Star> st (*astrobj);
    int nel ;
    try { nel = st->get_nelements(); }
    YGYOTO_STD_CATCH ;
    long dims[] = {1, nel};
    double * da=ypush_d(dims);
    double * dd=ypush_d(dims);
    double * dD=ypush_d(dims);

    try { st->getSkyPos(*screen, da, dd, dD); }
    YGYOTO_STD_CATCH ;

    yput_global(dD_ref,0);
    yarg_drop(1);
    yput_global(dd_ref,0);
    yarg_drop(1);
    yput_global(da_ref,0);

  }

  void
  Y_gyoto_Star_get_xyz(int n)
  {
    if (n<2) y_error("gyoto_Star_get_xyz takes at least 2 argument");
    if (n>4) y_error("gyoto_Star_get_xyz takes at most 4 arguments");
    SmartPointer<Astrobj::Generic> *astrobj=yget_Astrobj(n-1);
    if (strcmp((*astrobj)->getKind().c_str(),"Star")) y_error("first argument must be a GYOTO Star object ");
    long x_ref=yget_ref(n-2);
    long y_ref=yget_ref(n-3);
    long z_ref=yget_ref(n-4);

    SmartPointer<Star> st (*astrobj);
    int nel ;
    try { nel = st->get_nelements(); }
    YGYOTO_STD_CATCH ;
    long dims[] = {1, nel};
    double * x=ypush_d(dims);
    double * y=ypush_d(dims);
    double * z=ypush_d(dims);

    try { st->get_xyz(x,y,z); }
    YGYOTO_STD_CATCH ;

    yput_global(z_ref,0);
    yarg_drop(1);
    yput_global(y_ref,0);
    yarg_drop(1);
    yput_global(x_ref,0);

  }

  void
  Y_gyoto_Star_get_coord(int n)
  {
    if (n<2) y_error("gyoto_Star_get_coord takes at least 2 argument");
    if (n>5) y_error("gyoto_Star_get_coord takes at most 5 arguments");
    SmartPointer<Astrobj::Generic> *astrobj=yget_Astrobj(n-1);
    if (strcmp((*astrobj)->getKind().c_str(),"Star")) y_error("first argument must be a GYOTO Star object ");
    long x0_ref=yget_ref(n-2);
    long x1_ref=yget_ref(n-3);
    long x2_ref=yget_ref(n-4);
    long x3_ref=yget_ref(n-5);

    SmartPointer<Star> st (*astrobj);
    int nel ;
    try { nel = st->get_nelements(); }
    YGYOTO_STD_CATCH ;
    long dims[] = {1, nel};
    double * x0=ypush_d(dims);
    double * x1=ypush_d(dims);
    double * x2=ypush_d(dims);
    double * x3=ypush_d(dims);

    try { st->getCoord(x0, x1, x2, x3); }
    YGYOTO_STD_CATCH ;

    yput_global(x3_ref,0);
    yarg_drop(1);
    yput_global(x2_ref,0);
    yarg_drop(1);
    yput_global(x1_ref,0);
    yarg_drop(1);
    yput_global(x0_ref,0);

  }

  void
  Y_gyoto_Star_get_dot(int n)
  {
    if (n<2) y_error("gyoto_Star_get_dot takes at least 2 argument");
    if (n>5) y_error("gyoto_Star_get_dot takes at most 5 arguments");
    SmartPointer<Astrobj::Generic> *astrobj=yget_Astrobj(n-1);
    if (strcmp((*astrobj)->getKind().c_str(),"Star")) y_error("first argument must be a GYOTO Star object ");
    long x0_ref=yget_ref(n-2);
    long x1_ref=yget_ref(n-3);
    long x2_ref=yget_ref(n-4);
    long x3_ref=yget_ref(n-5);

    SmartPointer<Star> st (*astrobj);
    int nel ;
    try { nel = st->get_nelements(); }
    YGYOTO_STD_CATCH ;
    long dims[] = {1, nel};
    double * x0=ypush_d(dims);
    double * x1=ypush_d(dims);
    double * x2=ypush_d(dims);
    double * x3=ypush_d(dims);

    try { st->get_dot(x0, x1, x2, x3); }
    YGYOTO_STD_CATCH ;

    yput_global(x3_ref,0);
    yarg_drop(1);
    yput_global(x2_ref,0);
    yarg_drop(1);
    yput_global(x1_ref,0);
    yarg_drop(1);
    yput_global(x0_ref,0);
  }

  void
  Y_gyoto_Star_get_prime(int n)
  {
    if (n<2) y_error("gyoto_Star_get_prime takes at least 2 argument");
    if (n>4) y_error("gyoto_Star_get_prime takes at most 4 arguments");
    SmartPointer<Astrobj::Generic> *astrobj=yget_Astrobj(n-1);
    if (strcmp((*astrobj)->getKind().c_str(),"Star")) y_error("first argument must be a GYOTO Star object ");
    long x1_ref=yget_ref(n-2);
    long x2_ref=yget_ref(n-3);
    long x3_ref=yget_ref(n-4);

    SmartPointer<Star> st (*astrobj);
    int nel ;
    try { nel = st->get_nelements(); }
    YGYOTO_STD_CATCH ;
    long dims[] = {1, nel};
    double * x1=ypush_d(dims);
    double * x2=ypush_d(dims);
    double * x3=ypush_d(dims);

    try { st->get_prime(x1, x2, x3); }
    YGYOTO_STD_CATCH ;

    yput_global(x3_ref,0);
    yarg_drop(1);
    yput_global(x2_ref,0);
    yarg_drop(1);
    yput_global(x1_ref,0);

  }

  void
  Y_is_gyoto_Star(int argc)
  {
    if (!yarg_Astrobj(0)) {
      ypush_long(0);
      return;
    }
    SmartPointer<Astrobj::Generic> *ao = yget_Astrobj(0);
    ypush_long((*ao)->getKind() == "Star");
  }



}
