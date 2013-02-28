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

#include "GyotoPhoton.h"
#include "GyotoFactory.h"
#include "GyotoDefs.h"
#include "ygyoto.h"
#include "yapi.h"

using namespace Gyoto;
using namespace std;

#define YGYOTO_PHOTON_GENERIC_KW "metric", "initcoord", "astrobj",	\
    "spectro", "tmin",							\
    "xfill", "save_txyz", "xmlwrite", "is_hit",				\
    "get_txyz", "get_coord", "get_cartesian", "clone"
#define YGYOTO_PHOTON_GENERIC_KW_N 13

void ygyoto_Photon_generic_eval(Gyoto::SmartPointer<Gyoto::Photon>* ph,
				 int *kiargs, int *piargs, int *rvset, int *paUsed) {
  if (debug()) cerr << "\nDEBUG: In ygyoto_Photon_generic_eval: ";
  int k=-1, iarg;
  char const * rmsg="Cannot set return value more than once";
  char const * pmsg="Cannot use positional argument more than once";

  //// MEMBERS ////
  /* METRIC */
  if ((iarg=kiargs[++k])>=0) {
    iarg+=*rvset;
    if (debug()) cerr << "metric=";
    if (yarg_nil(iarg)) {
      if ((*rvset)++) y_error(rmsg);
      *ypush_Metric()=(*ph)->getMetric();
    } else
      (*ph)->setMetric(*yget_Metric(iarg)) ;
    if (debug()) cerr << "... ";
  }

  /* INITCOORD */
  if ((iarg=kiargs[++k])>=0) {
    iarg+=*rvset;
    if (debug()) cerr << "initcoord=";
    if (yarg_nil(iarg)) { // Getting initcoord
      if (debug()) cerr << "     get_initcoord=1" << endl;
      if ((*rvset)++) y_error(rmsg);
      long dims[]= {1, 8};
      double *coord = ypush_d(dims);
      (*ph)->getInitialCoord(coord);
    } else {          // Setting initcoord
      long ptot=1;
      double coord[8];

      if (yarg_number(iarg)) { // initcoord=coord or initcoord=pos, vel
	double *pos=ygeta_d(iarg, &ptot, 0);
	if (ptot!=4 && ptot !=8)
	  y_error("POS should have either 4 or 8 elements.");
	for (int ii=0; ii<ptot; ++ii) coord[ii]=pos[ii];
	if (ptot==4) {
	  if ((*paUsed)++) y_error(pmsg);
	  long vtot=1;
	  double *vel=ygeta_d(piargs[0]+*rvset, &vtot, 0);
	  if (vtot==4) for (int ii=0; ii<4; ++ii) coord[ii+4]=vel[ii];
	  else if (vtot==3) {
	    for (int ii=0; ii<3; ++ii) coord[ii+5]=vel[ii];
	    if (!(*ph)->getMetric())
	      y_error("METRIC should have been set already");
	    (*ph)->getMetric()->nullifyCoord(coord);
	  } else y_error("VEL should have 3 or 4 elements");
	}
      } else {
	SmartPointer<Screen> sc = NULL;
	if (yarg_Scenery(iarg)) { // initcoord=scenery,i,j or senery,da,dd
	  SmartPointer<Scenery> *scenery = yget_Scenery(iarg);
	  (*ph) -> setMetric((*scenery)->getMetric());
	  (*ph) -> setAstrobj((*scenery)->getAstrobj());
	  sc = (*scenery)->getScreen();
	} else sc = *yget_Screen(iarg); //initcoord=screen,i,j or screen, da, dd
	if (yarg_number(piargs[0]+*rvset)==2) {
	  double da = ygets_d(piargs[0]+*rvset);
	  double dd = ygets_d(piargs[1]+*rvset);
	  if (debug()) cerr << "screen, dalpha="<<da <<", ddelta="<<dd <<endl;
	  sc -> getRayCoord(da, dd, coord);
	} else {
	  size_t i = ygets_l(piargs[0]+*rvset);
	  size_t j = ygets_l(piargs[1]+*rvset);
	  if (debug()) cerr << "screen, i="<<i <<", j="<<j << endl;
	  sc -> getRayCoord(i, j, coord);
	}
	(*ph) -> setSpectrometer(sc->getSpectrometer());
      }

      if (debug()) {
	cerr << "DEBUG:        initcoord=[" ;
	for (int ii=0; ii<7; ++ii) cerr << coord[ii] << ", ";
	cerr << coord[7] ;
      }
      (*ph)->setInitCoord(coord) ;
      if (debug()) cerr << "]" << endl;
    }
    if (debug()) cerr << "... " << endl;
  }

  /* ASTROBJ */
  if ((iarg=kiargs[++k])>=0) {
    iarg+=*rvset;
    if (debug()) cerr << "astrobj=";
    if (yarg_nil(iarg)) {
      if ((*rvset)++) y_error(rmsg);
      *ypush_Astrobj() = (*ph)->getAstrobj();
    } else
      (*ph) -> setAstrobj( *yget_Astrobj(iarg) ) ;
    if (debug()) cerr << "... " << endl;
  }

  /* SPECTRO */
  if ((iarg=kiargs[++k])>=0) {
    iarg+=*rvset;
    if (yarg_nil(iarg)) {
      if ((*rvset)++) y_error(rmsg);
      *ypush_Spectrometer() = (*ph) -> getSpectrometer();
    } else {
      (*ph) -> setSpectrometer(*yget_Spectrometer(iarg));
    }
  }

  /* TMIN */
  if ((iarg=kiargs[++k])>=0) {
    iarg+=*rvset;
    if (yarg_nil(iarg)) {
      if ((*rvset)++) y_error(rmsg);
      ypush_double((*ph) -> getTmin());
    } else {
      (*ph) -> setTmin(ygets_d(iarg));
    }
  }

  //// METHODS ////
  // SUBROUTINE-LIKE //
  // xfill=tlim: integrate geodesic
  if ((iarg=kiargs[++k])>=0) {
    iarg+=*rvset;
    if (debug()) cerr << "xfill=";
    double tlim = ygets_d(iarg);
    if (debug()) cerr << tlim;
    (*ph)->xFill(tlim);
    if (debug()) cerr << "... " << endl;
  }

  // save_txyz=filename, t1, mass_sun, distance_kpc, unit, screen
  if ((iarg=kiargs[++k])>=0) {
    iarg+=*rvset;
    if (debug()) cerr << "     save_txyz=";
    int myarg=0;
    if (debug()) cerr << "save_txyz=";
    char*  filename          = ygets_q(iarg);
    double t1                = 0.;//ygets_d(piargs[0]+*rvset);
    double mass_sun          = 1.;//ygets_d(piargs[1]+*rvset);
    double distance_kpc      = 1.;//ygets_d(piargs[2]+*rvset);
    string unit              = "geometrical";//ygets_q(piargs[3]+*rvset);
    SmartPointer<Screen>* sc = NULL; //yget_Screen(piargs[4]+*rvset);

    if (yarg_number(piargs[myarg]+*rvset)) {
      t1= ygets_d(piargs[myarg++]+*rvset);
      if (yarg_number(piargs[myarg]+*rvset)) {
	mass_sun= ygets_d(piargs[myarg++]+*rvset);
	if (yarg_number(piargs[myarg]+*rvset)) {
	  distance_kpc= ygets_d(piargs[myarg++]+*rvset);
	}
      }
    }
    if (yarg_string(piargs[myarg]+*rvset))
      unit = ygets_q(piargs[myarg++]+*rvset);
    if (yarg_Screen(piargs[myarg]+*rvset))
      sc = yget_Screen(piargs[myarg++]+*rvset);

    (*ph) -> save_txyz ( filename, t1, mass_sun, distance_kpc, unit,
			 (sc?*sc:SmartPointer<Screen>(NULL)) );

    if (debug()) cerr << filename << endl ;
  }

  // Save to file
  if ((iarg=kiargs[++k])>=0) { // xmlwrite
    iarg+=*rvset;
    if (debug()) cerr << "     xmlwrite=";
#ifdef GYOTO_USE_XERCES
    char *filename=ygets_q(iarg);
    Factory(*ph).write(filename);
    if (debug()) cerr << filename << endl;
#else
    y_error("This GYOTO was compiled without XERCES: no xml i/o");
#endif
  }

  // FUNCTION-LIKE //

  /* is_hit= */
  if ((iarg=kiargs[++k])>=0) { // is_hit=whatever
    if (debug()) cerr << "     is_hit=" << endl;
    if ((*rvset)++) y_error(rmsg);
    Astrobj::Properties junk;
    ypush_int((*ph)->hit(&junk));
  }

  /* GET_TXYZ */
  if ((iarg=kiargs[++k])>=0) { // get_txyz=
    if (debug()) cerr << "     get_txyz=1" << endl;
    if ((*rvset)++) y_error(rmsg);
    int nel =(*ph)->get_nelements();
      
    long dims[] = {2, nel, 4};
    double * data=ypush_d(dims);

    (*ph)->get_t(data);
    (*ph)->get_xyz(data+nel, data+2*nel, data+3*nel);
  }

  if ((iarg=kiargs[++k])>=0) { // get_coord
    if (debug()) cerr << "     get_coord=dates" << endl;
    if ((*rvset)++) y_error(rmsg);
    // Two cases : get_coord=1 or get_coord=dates. 
    // We will recognize the first case if the parameter is the
    // integer scalar 1, any other type or value is a date
    // specification.
    int argt = 1;yarg_number(iarg);
    long ntot[1] = {1};
    long dims[Y_DIMSIZE] = {0};
    double * dates = NULL;

    if (!yarg_nil(iarg)) {
      argt = yarg_number(iarg);
      dates=ygeta_d(iarg, ntot, dims);
    } else {
      dates = ypush_d(dims);
      dates[0]=1.;
    }

    if (debug())
      cerr << "DEBUG: gyoto_Photon(get_coord=array("
	   << (yarg_number(iarg)==1?"long":"double")
	   << ", " << *ntot 
	   << ") (rank=" << dims[0]
	   << ", dates[0]=" << dates[0] << ")" << endl;

    if (dims[0] == 0 && argt == 1 && dates[0] == 1) {
      if(debug())
	cerr << "DEBUG: retrieving all already computed coordinates" << endl;
      int nel =(*ph)->get_nelements();
      long ddims[] = {2, nel, 8};
      double * data = ypush_d(ddims);

      (*ph)->getCoord(data, data+nel, data+2*nel, data+3*nel);
      (*ph)->get_dot(data+4*nel, data+5*nel, data+6*nel, data+7*nel);
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
      (*ph) -> getCoord(dates, nel,
			data, data+ nel,data+2* nel,data+3* nel,
			data+4* nel, data+5* nel, data+6* nel);
    }
  }

  if ((iarg=kiargs[++k])>=0) { // get_cartesian
    if (debug()) cerr << "     get_cartesian=dates" << endl;
    if ((*rvset)++) y_error(rmsg);
    long ntot[1] = {1};
    long dims[Y_DIMSIZE];
    double* dates = ygeta_d(iarg, ntot, dims);
    if (dims[0] > Y_DIMSIZE-2)
      y_errorn("gyoto_Star(get_cartesian=dates): DATES must have at most %ld "
	       "dimensions", Y_DIMSIZE-2);
    dims[0] += 1;
    dims[dims[0]] = 6;
    double * data = ypush_d(dims);
    size_t nel = *ntot;
    (*ph) -> getCartesian(dates, nel,
			  data,         data +   nel, data + 2*nel,
			  data + 3*nel, data + 4*nel, data + 5*nel);
    
  }

  /* CLONE */
  if ((iarg=kiargs[++k])>=0) {
    if ((*rvset)++) y_error(rmsg);
    *ypush_Photon() = (*ph)->clone();
  } 

  if (debug()) cerr << endl;
}


extern "C" {

// PHOTON CLASS

// Photon()
// setInitialCondition(Metric*, Astrobj*, coord[8], sys);
// setDelta(double)
// hit(double)
  typedef struct gyoto_Photon { SmartPointer<Photon> photon; } gyoto_Photon;
  void gyoto_Photon_free(void *obj) {
    if (((gyoto_Photon*)obj)->photon) {
      ((gyoto_Photon*)obj)->photon=NULL;
    } else printf("null pointer\n");
  }
  void gyoto_Photon_print(void *obj) {
#ifdef GYOTO_USE_XERCES
    string rest="", sub="";
    size_t pos=0, len=0;
    try {rest = Factory(((gyoto_Photon*)obj)->photon).format();}
    YGYOTO_STD_CATCH;
    while (len=rest.length())  {
      sub=rest.substr(0, pos=rest.find_first_of("\n",0));
      rest=rest.substr(pos+1, len-1);
      y_print( sub.c_str(),1 );
    }
#else
    y_print("GYOTO photon",0);
#endif
  }
  void gyoto_Photon_eval(void *obj, int argc) {
    // If no parameters, return pointer
    if (argc==1 && yarg_nil(0)) {
      ypush_long((long)((gyoto_Photon*)obj)->photon());
      return;
    }

    SmartPointer<Photon> *ph = &(((gyoto_Photon*)obj)->photon);

    static char const * knames[]={
      YGYOTO_PHOTON_GENERIC_KW, 0
    };
    static long kglobs[YGYOTO_PHOTON_GENERIC_KW_N+1];
    int kiargs[YGYOTO_PHOTON_GENERIC_KW_N];
    int piargs[]={-1,-1,-1,-1,-1};
    // push default return value: the photon itsef
    *ypush_Photon() = *ph;
    yarg_kw_init(const_cast<char**>(knames), kglobs, kiargs);
   
    int iarg=argc, parg=0;
    while (iarg>=1) {
      iarg = yarg_kw(iarg, kglobs, kiargs);
      if (iarg>=1) {
	if (parg<5) piargs[parg++]=iarg--;
	else y_error("gyoto_Photon takes at most 5 positional arguments");
      }
    }

    int rvset[1]={0}, paUsed[1]={0};
    ygyoto_Photon_generic_eval(ph, kiargs, piargs, rvset, paUsed);


  }
  static y_userobj_t gyoto_Photon_obj =
    {const_cast<char*>("gyoto_Photon"), &gyoto_Photon_free, &gyoto_Photon_print, &gyoto_Photon_eval, 0, 0};
  
  // Generic constructor/accessor
  void
  Y_gyoto_Photon(int argc)
  {
    int rvset[1]={0}, paUsed[1]={0}, builder=0;
    SmartPointer<Photon> *ph = NULL;
    if (yarg_Photon(argc-1)) { // photon on the stack
      ph = yget_Photon(--argc);
      *ypush_Photon() = *ph; // push back photon on the stack
    } else { // constructor mode
      ph = ypush_Photon();
      builder=1;
    }

    static char const * knames[]={
      YGYOTO_PHOTON_GENERIC_KW, 0
    };
    static long kglobs[YGYOTO_PHOTON_GENERIC_KW_N+1];
    int kiargs[YGYOTO_PHOTON_GENERIC_KW_N];
    int piargs[]={-1,-1,-1,-1};
    // push default return value: need to drop before pushing another one
    yarg_kw_init(const_cast<char**>(knames), kglobs, kiargs);
   
    int iarg=argc, parg=0;
    while (iarg>=1) {
      iarg = yarg_kw(iarg, kglobs, kiargs);
      if (iarg>=1) {
	if (parg<4) piargs[parg++]=iarg--;
	else y_error("gyoto_Photon takes at most 4 positional arguments");
      }
    }
    
    // if rvset==1, constructor mode:
    if (builder) {
      if (yarg_string(piargs[0])) {
#ifdef GYOTO_USE_XERCES
	// 	try { *ao_ = Factory(ygets_q(piargs[0])).getAstrobj(); } 
	// 	YGYOTO_STD_CATCH;
	*paUsed=1;
#else
	y_error("This GYOTO was compiled without XERCES: no xml i/o");
#endif
      } else *ph = new Photon();
    }

    ygyoto_Photon_generic_eval(ph, kiargs, piargs, rvset, paUsed);

  }

  void
  Y_gyoto_Photon_new(int n)
  {
    gyoto_Photon * obj=(gyoto_Photon*)ypush_obj(&gyoto_Photon_obj, sizeof(gyoto_Photon));
    try { obj->photon = new Photon(); }
    YGYOTO_STD_CATCH ;
  }

  void
  Y_gyoto_Photon_setInitialCondition(int n)
  {  //(Metric*, Astrobj*, coord[8], sys);
    gyoto_Photon  *phobj =(gyoto_Photon*) yget_obj(n-1, &gyoto_Photon_obj);
    SmartPointer<Metric::Generic> *gg = yget_Metric(n-2);
    SmartPointer<Astrobj::Generic> *astrobj = yget_Astrobj(n-3);
    SmartPointer<Screen> *screen = yget_Screen(n-4);

    if (n==5) {
      long ntot=1;
      double * coord = ygeta_d(n-4, &ntot, NULL);
      if (ntot < 4) y_error("coord must have at least 4 elements");
      try {
	(phobj->photon)->setInitialCondition(*gg, *astrobj, coord);
      } YGYOTO_STD_CATCH ;
    } else if (n==6) {
      double d_alpha = ygets_d(n-5);
      double d_delta = ygets_d(n-6);
      try {
	(phobj->photon)->setInitialCondition(*gg, *astrobj, *screen, d_alpha, d_delta);
      } YGYOTO_STD_CATCH ;
    } else y_error("gyoto_Photon_setInitialCondition takes either 4 or 7 arguments");
  }

  void
  Y_gyoto_Photon_setDelta(int n)
  {  //(Metric*, Astrobj*, coord[8], sys);
    gyoto_Photon * phobj =(gyoto_Photon*)yget_obj(n-1, &gyoto_Photon_obj);
    double delta=ygets_d(n-2);
    try {
      (phobj->photon)->setDelta(delta);
    } YGYOTO_STD_CATCH ;
  }

  void
  Y_gyoto_Photon_hit(int n)
  {  //(Metric*, Astrobj*, coord[8], sys);
    gyoto_Photon * phobj =(gyoto_Photon*)yget_obj(n-1, &gyoto_Photon_obj);
    double tlim=ygets_d(n-2);
    ypush_int((phobj->photon)->hit(NULL));
  }

}

SmartPointer<Photon>* yget_Photon(int iarg) {
  return &(((gyoto_Photon*)yget_obj(iarg, &gyoto_Photon_obj))->photon);
}

SmartPointer<Photon>* ypush_Photon() {
  return &(((gyoto_Photon*)ypush_obj(&gyoto_Photon_obj, sizeof(gyoto_Photon)))->photon);
}

int yarg_Photon(int iarg) {
  return yget_obj(iarg,0)==gyoto_Photon_obj.type_name;
}



