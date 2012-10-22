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
#include <yapi.h>
#include <pstdlib.h>
#ifdef GYOTO_USE_XERCES
#include <GyotoFactory.h>
#endif

#include <iostream>
using namespace std;


using namespace Gyoto;

extern "C" {
  typedef struct gyoto_Screen { SmartPointer<Screen> screen; } gyoto_Screen;
  void gyoto_Screen_free(void *obj);
  void gyoto_Screen_print(void *obj);
  void gyoto_Screen_eval(void *obj, int n);
  static y_userobj_t gyoto_Screen_obj =
    {const_cast<char*>("gyoto_Screen"), &gyoto_Screen_free, &gyoto_Screen_print, &gyoto_Screen_eval, 0, 0};

  // SCREEN CLASS
  void gyoto_Screen_free(void *obj) {
    if (((gyoto_Screen*)obj)->screen) {
      ((gyoto_Screen*)obj)->screen=NULL;
    } else printf("null pointer\n");
  }

  void gyoto_Screen_print(void *obj) {
#ifdef GYOTO_USE_XERCES
    string rest="", sub="";
    try { rest = Factory(((gyoto_Screen*)obj)->screen).format(); }
    YGYOTO_STD_CATCH;
    size_t pos=0, len;
    while (len=rest.length())  {
      sub=rest.substr(0, pos=rest.find_first_of("\n",0));
      rest=rest.substr(pos+1, len-1);
      y_print( sub.c_str(),1 );
    }
#else
    y_print("GYOTO screen ",0);
    //    SmartPointer<Screen> gg = ((gyoto_Screen*)obj)->screen;
    //    y_print(gg->getKind().c_str(),0);
#endif
  }

  void gyoto_Screen_eval(void *obj, int argc) {
    int rvset[1]={0}, paUsed[1]={0};
    int k=-1;
    char const * rmsg="Cannot set return value more than once";
    char const * pmsg="Cannot use positional argument more than once";

    SmartPointer<Screen> *screen = NULL; //&((gyoto_Screen*)obj)->screen;

    if (!obj) {
      obj = ypush_obj(&gyoto_Screen_obj, sizeof(gyoto_Screen));
      screen = &(((gyoto_Screen*)obj)->screen);
      *screen = new Screen();
    } else if (argc==1 && yarg_nil(0)) { // If no parameters, return pointer
      ypush_long( (long) ((gyoto_Screen*)obj)->screen() );
      return;
    } else {
      screen = &(((gyoto_Screen*)obj)->screen);
      *ypush_Screen()=*screen;
    }
      
    static char const * knames[]={
      "unit",
      "metric",
      "time","fov","resolution",
      "distance", "dmax", "inclination", "paln", "argument",
      "projection", "observerpos",
      "spectro",
      "skycoord",  "raycoord",
      "xmlwrite", "clone",
      0
    };
#define nkw 17
    static long kglobs[nkw+1];
    int kiargs[nkw];
    int piargs[]={-1,-1,-1,-1};
    // push default return value: need to drop before pushing another one
    yarg_kw_init(const_cast<char**>(knames), kglobs, kiargs);
    char * unit=NULL;
      
    int iarg=argc, parg=0;
    while (iarg>=1) {
      iarg = yarg_kw(iarg, kglobs, kiargs);
      if (iarg>=1) {
	if (parg<4) piargs[parg++]=iarg--;
	else y_error("gyoto_Metric takes at most 4 positional arguments");
      }
    }

    /* UNIT */
    if ((iarg=kiargs[++k])>=0) {
      iarg+=*rvset;
      GYOTO_DEBUG << "get unit" << endl;
      unit = ygets_q(iarg);
    }

    /* METRIC */
    if ((iarg=kiargs[++k])>=0) {
      iarg+=*rvset;
      if (yarg_nil(iarg)) { // get 
	if ((*rvset)++) y_error("Only one return value possible");
	*ypush_Metric() = (*screen)->getMetric();
      } else                // set
	(*screen)->setMetric(*yget_Metric(kiargs[k]));
    }

    /* OBSERVING TIME */
    if ((iarg=kiargs[++k])>=0) {
      iarg+=*rvset;
      if (yarg_nil(iarg)) { // get 
	if ((*rvset)++) y_error("Only one return value possible");
	ypush_double((*screen)->getTime(unit?unit:""));
      } else
	(*screen) -> setTime(ygets_d(iarg), unit?unit:"");
    }

    /* FIELD OF VIEW (fov) */
    if ((iarg=kiargs[++k])>=0) {
      iarg+=*rvset;
      if (yarg_nil(iarg)) { // get 
	if ((*rvset)++) y_error("Only one return value possible");
	ypush_double((*screen)->getFieldOfView());
      } else
	(*screen) -> setFieldOfView(ygets_d(iarg));
    }

    /* RESOLUTION */
    if ((iarg=kiargs[++k])>=0) {
      iarg+=*rvset;
      if (yarg_nil(iarg)) { // get 
	if ((*rvset)++) y_error("Only one return value possible");
	ypush_long((*screen)->getResolution());
      } else
	(*screen) -> setResolution  (ygets_l(iarg));
    }

    /* DISTANCE */
    if ((iarg=kiargs[++k])>=0) {
      iarg+=*rvset;
      if (yarg_nil(iarg)) { // get 
	if ((*rvset)++) y_error(rmsg);
	ypush_double((*screen)->getDistance(unit?unit:""));
      } else
	(*screen) -> setDistance    (ygets_d(iarg), unit?unit:"") ;
    }

    /* DMAX */ 
    if ((iarg=kiargs[++k])>=0) {
      iarg+=*rvset;
      if (yarg_nil(iarg)) { // get 
	if ((*rvset)++) y_error(rmsg);
	ypush_double((*screen)->getDmax());
      } else
	(*screen) -> setDmax (ygets_d(iarg)) ;
    }

    /* INCLINATION */ 
    if ((iarg=kiargs[++k])>=0) {
      iarg+=*rvset;
      if (yarg_nil(iarg)) { // get 
	if ((*rvset)++) y_error(rmsg);
	ypush_double((*screen)->getInclination(unit?unit:""));
      } else
	(*screen) -> setInclination (ygets_d(iarg), unit?unit:"") ;
    }

    /* POSITION ANGLE OF THE LINE OF NODES (paln) */
    if ((iarg=kiargs[++k])>=0) {
      iarg+=*rvset;
      if (yarg_nil(iarg)) { // get 
	if ((*rvset)++) y_error(rmsg);
	ypush_double((*screen)->getPALN(unit?unit:""));
      } else
	(*screen) -> setPALN        (ygets_d(iarg), unit?unit:"") ;
    }

    /* ARGUMENT */
    if ((iarg=kiargs[++k])>=0) {
      iarg+=*rvset;
      if (yarg_nil(iarg)) { // get 
	if ((*rvset)++) y_error(rmsg);
	ypush_double((*screen)->getArgument(unit?unit:""));
      } else
	(*screen) -> setArgument    (ygets_d(iarg), unit?unit:"") ;
    }
      
    /* PROJECTION */
    if ((iarg=kiargs[++k])>=0) { // Set Projection
      iarg+=*rvset;
      long ntot;
      double *proj=ygeta_d(iarg, &ntot, NULL);
      switch (ntot) {
      case 4:
	(*screen)->setProjection(proj[0], proj[1], proj[2], proj[3]);
	break;
      case 3:
	(*screen)->setProjection(proj[0], proj[1], proj[2]);
	break;
      }
    }

    /* OBSERVERPOS */
    if ((iarg=kiargs[++k])>=0) {
      iarg+=*rvset;
      if (yarg_nil(iarg)) { // get_observerpos
	if ((*rvset)++) y_error(rmsg);
	long dims[] = {1,4};
	double * coord=ypush_d(dims);
	(*screen)->getObserverPos(coord);
      } else { // Set ObserverPos
	long ntot;
	double * pos = ygeta_d(iarg, &ntot, NULL);
	if (ntot<4) y_error("POS must have at least 4 elements");
	(*screen) -> setObserverPos(pos);
      }
    }

    /* SPECTRO */
    if ((iarg=kiargs[++k])>=0) {
      iarg+=*rvset;
      if (yarg_nil(iarg)) {
	if ((*rvset)++) y_error(rmsg);
	*ypush_Spectrometer() = (*screen) -> getSpectrometer();
      } else {
	(*screen) -> setSpectrometer(*yget_Spectrometer(iarg));
      }
    }


    ///// METHODS //////
    /* SKYCOORD METHOD */
    if ((iarg=kiargs[++k])>=0) { // skycoord
      if ((*rvset)++) y_error(rmsg);
      long ntot=1;
      double *pos=ygeta_d(iarg, &ntot, NULL);
      if (ntot<4) y_error("POS argument should have at lest 4 elements");
	
      long dims[] = {1, 3};
      double * skypos=ypush_d(dims);
	
      (*screen)->coordToXYZ(pos, skypos);
    }
      
    /* RAYCOORD METHOD */
    if ((iarg=kiargs[++k])>=0) { // raycoord
      if ((*rvset)++) y_error(rmsg);
	
      long ntot=1;
      double *pos=ygeta_d(iarg, &ntot, NULL);
      if (ntot<2) y_error("X_Y argument should have at lest 4 elements");
	
      long dims[] = {1,8};
      yarg_drop(1);
      double * coord=ypush_d(dims);
	
      (*screen)->getRayCoord(pos[0], pos[1], coord);
    }
      
    // Save to file
    if ((iarg=kiargs[++k])>=0) { // xmlwrite
      iarg+=*rvset;
#ifdef GYOTO_USE_XERCES
      char *filename=ygets_q(kiargs[k]);
      Factory((*screen)).write(filename);
#else
      y_error("This GYOTO was compiled without XERCES: no xml i/o");
#endif
    }

    /* CLONE */
    if ((iarg=kiargs[++k])>=0) {
      if ((*rvset)++) y_error(rmsg);
      *ypush_Screen() = (*screen)->clone();
    } 
  }
}

// PUBLIC API

SmartPointer<Screen> *yget_Screen(int iarg) {
  return &((gyoto_Screen*)yget_obj(iarg, &gyoto_Screen_obj))->screen;
}
SmartPointer<Screen> *ypush_Screen() {
  gyoto_Screen* obj = (gyoto_Screen*)ypush_obj(&gyoto_Screen_obj, sizeof(gyoto_Screen));
  return &(obj->screen);
}

int yarg_Screen(int iarg) {
  return yget_obj(iarg,0)==gyoto_Screen_obj.type_name;
}


// YAPI FUNCTIONS

extern "C" {

  void Y_gyoto_Screen(int argc) {
    void *obj = NULL;
    if (yarg_Screen(argc)) obj = yget_obj(--argc, &gyoto_Screen_obj) ;
    gyoto_Screen_eval(obj, argc);
  }


}
